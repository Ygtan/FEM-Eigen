#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: nonlinear cable_solver.exe (fem_3d_nl.hpp) end-to-end smoke.

Distinct from ``compare_cable_arch_3d.py``: this exercises the geometric-
nonlinear co-rotational chain (tension-only cables, unstressed length L0,
sag / Ernst modulus, cable-length iteration, COMPLETION closure).  There is
no linear ANSYS baseline for this chain, so we check internal consistency
(convergence / closure / tension / cable-length closure) and compare the
fresh run against the frozen ``*_cable.txt``.

All runtime output is ASCII.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MSYS_BIN = r"C:\msys64\mingw64\bin"
STEM = "cb3d05_completion"


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def parse_cable_txt(path: Path) -> dict:
    meta: dict[str, float] = {}
    cables: list[dict] = []
    sec = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("["):
            sec = line.strip("[]").upper()
            continue
        if sec:
            if sec == "CABLE" and len(line.split()) >= 7:
                p = line.split()
                try:
                    cables.append({
                        "elem": int(p[0]), "anchor": int(p[1]),
                        "T": float(p[2]), "L0": float(p[3]),
                        "chord": float(p[4]), "Ltarget": float(p[5]),
                        "dL": float(p[6]),
                    })
                except ValueError:
                    pass
            continue
        for tok in line.split():
            if "=" not in tok:
                continue
            k, _, v = tok.partition("=")
            try:
                meta[k] = float(v)
            except ValueError:
                pass
    return {"meta": meta, "cables": cables}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="FEM-Eigen build dir")
    ap.add_argument("--smoke", default=str(HERE / "nonlinear_smoke"))
    a = ap.parse_args()

    exe = Path(a.bin) / "cable_solver.exe"
    smoke = Path(a.smoke)
    dat = smoke / (STEM + ".dat")
    frozen = smoke / (STEM + "_cable.txt")
    print("== cable_solver nonlinear smoke ==")
    print("   exe: %s" % exe)
    print("   dat: %s" % dat)
    if not exe.exists():
        print("[FAIL] cable_solver.exe not found: %s" % exe)
        return 1
    if not dat.exists() or not frozen.exists():
        print("[FAIL] archived .dat / frozen _cable.txt missing")
        return 1

    errs: list[str] = []
    npass = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal npass
        if ok:
            npass += 1
            print("[PASS] %-46s %s" % (name, detail))
        else:
            errs.append(name + " " + detail)
            print("[FAIL] %-46s %s" % (name, detail))

    with tempfile.TemporaryDirectory(prefix="cbsmoke_") as td:
        work = Path(td)
        shutil.copy2(dat, work / dat.name)
        p = subprocess.run([str(exe), str(work / dat.name)], cwd=str(work),
                           env=env_with_msys(), capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=600)
        fresh = work / (STEM + "_cable.txt")
        if not fresh.exists():
            print("[FAIL] no _cable.txt produced (rc=%d) %s" %
                  (p.returncode, (p.stdout or "")[-400:]))
            return 1
        cur = parse_cable_txt(fresh)
        old = parse_cable_txt(frozen)

    m, c = cur["meta"], cur["cables"]
    om, oc = old["meta"], old["cables"]
    check("converged == 1", m.get("converged", 0) == 1,
          "iterations=%g elapsed_ms=%g" % (m.get("iterations", -1),
                                         m.get("elapsed_ms", -1)))
    check("closure_err < 1e-6 m", m.get("closure_err", 1) < 1e-6,
          "closure_err=%.6e" % m.get("closure_err", float("nan")))
    check("tension_err < 1e-6", m.get("tension_err", 1) < 1e-6,
          "tension_err=%.6e" % m.get("tension_err", float("nan")))
    check("cable count matches frozen", len(c) == len(oc) == 48,
          "fresh=%d frozen=%d" % (len(c), len(oc)))
    if c:
        dL = max(abs(x["dL"]) for x in c)
        check("max |dL| < 1e-6 m (length iteration)", dL < 1e-6,
              "max |dL|=%.6e m" % dL)
        check("all cable tensions > 0 (tension only)",
              all(x["T"] > 0 for x in c),
              "min=%.1f N max=%.1f N" % (min(x["T"] for x in c),
                                        max(x["T"] for x in c)))
    sv, av = m.get("support_v_sum_kN", 0), m.get("applied_v_kN", 1)
    check("support vs applied load within 1%",
          abs(sv - av) / max(abs(av), 1.0) < 1e-2,
          "support=%.1f applied=%.1f kN rel=%.3e" %
          (sv, av, abs(sv - av) / max(abs(av), 1.0)))

    # Compare against frozen result: key scalars + every cable force / L0 / dL
    for k in ("closure_err", "tension_err", "tension_sum_kN",
              "tension_min_kN", "tension_max_kN", "support_v_sum_kN"):
        x, y = m.get(k), om.get(k)
        if x is None or y is None:
            check("frozen key %s present" % k, False, "missing")
            continue
        d = abs(x - y)
        check("frozen %s matches" % k, d <= 1e-6 + 1e-3 * abs(y),
              "fresh=%.8g frozen=%.8g" % (x, y))
    by_elem = {x["elem"]: x for x in oc}
    worst = 0.0
    nmatch = 0
    for x in c:
        y = by_elem.get(x["elem"])
        if not y:
            continue
        nmatch += 1
        for fld in ("T", "L0", "dL"):
            d = abs(x[fld] - y[fld])
            scale = max(abs(y[fld]), 1e-9)
            worst = max(worst, d / scale if fld != "dL" else
                        d / max(scale, 1e-9))
    check("48 frozen cable rows reproduced", nmatch == 48,
          "matched=%d worst T/L0 rel=%.3e" % (nmatch, worst))
    check("frozen cable tensions within 1e-3", worst <= 1e-3,
          "worst rel=%.3e" % worst)

    print("== Result: %d PASS, %d FAIL ==" % (npass, len(errs)))
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
