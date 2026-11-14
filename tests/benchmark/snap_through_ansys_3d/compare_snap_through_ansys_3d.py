#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: CB3D17 post-buckling / snap-through (arc-length) vs ANSYS 18.0.

Offline by default (no ANSYS licence needed):

  * re-run ``snap_probe.exe`` on the archived ``*.snap`` and verify it
    reproduces the frozen ``*_fem.csv`` value by value;
  * re-derive all 1686 golden values from the frozen ANSYS CSVs (guards drift);
  * check every ANSYS substep load factor against the independent analytic
    ``P_log(v)`` anchor, and the reaction-derived lambda against ANSYS TIME;
  * check the fresh FEM arc-length path against its analytic anchor
    (engineering / logarithmic strain), the quadratic-fit limit point,
    and the descending branch;
  * compare the FEM arc-length path with the ANSYS ARCLEN path by quadratic
    interpolation at every ANSYS substep v;
  * check archived ANSYS logs are 0 error.

``--live-ansys`` re-runs all 3 LINK180 jobs in a temporary directory and
diffs every CSV cell against the frozen copies.  All output is ASCII.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from build_reference_snap import (  # noqa: E402
    SNAP_FILE, FEM_CSV, check,
)

MSYS_BIN = r"C:\msys64\mingw64\bin"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
TOL_LIVE_CSV = 1.0e-6        # absolute
TOL_LIVE_CSV_REL = 1.0e-10   # relative


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def read_fem_csv(path: Path) -> tuple[dict, dict]:
    pts: dict[str, list[list[float]]] = {}
    summ: dict[str, dict[str, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        q = line.split(",")
        if q[0] == "P" and len(q) >= 9:
            pts.setdefault(q[1], []).append([float(x) for x in q[2:]])
        elif q[0] == "S" and len(q) >= 4:
            summ.setdefault(q[1], {})[q[2]] = float(q[3])
    return pts, summ


def run_fem_probe(exe: Path, snap: Path, work: Path):
    shutil.copy2(snap, work / snap.name)
    out = work / "snap_fresh.csv"
    p = subprocess.run([str(exe), snap.name, out.name], cwd=str(work),
                       env=env_with_msys(), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=900)
    if p.returncode != 0 or not out.exists():
        raise RuntimeError("snap_probe rc=%d: %s %s" %
                           (p.returncode, p.stdout, p.stderr))
    return read_fem_csv(out)


def reproduction(fresh, frozen) -> tuple[int, float, list[str]]:
    n, worst, errs = 0, 0.0, []
    if set(fresh[0]) != set(frozen[0]) or set(fresh[1]) != set(frozen[1]):
        return 0, 0.0, ["FEM case ids differ"]
    for cid in frozen[0]:
        a, b = fresh[0][cid], frozen[0][cid]
        if len(a) != len(b):
            errs.append("{}: point count {} vs {}".format(cid, len(a), len(b)))
            continue
        for ra, rb in zip(a, b):
            for x, y in zip(ra, rb):
                n += 1
                d = abs(x - y)
                worst = max(worst, d)
                if d > 1e-12 + 1e-12 * abs(y):
                    errs.append("{} P: {:.12g} vs {:.12g}".format(cid, x, y))
    for cid in frozen[1]:
        a, b = fresh[1][cid], frozen[1][cid]
        if set(a) != set(b):
            errs.append("{}: summary keys differ".format(cid))
            continue
        for k in b:
            n += 1
            d = abs(a[k] - b[k])
            worst = max(worst, d)
            if d > 1e-12 + 1e-12 * abs(b[k]):
                errs.append("{} S {}: {:.12g} vs {:.12g}".format(
                    cid, k, a[k], b[k]))
    return n, worst, errs[:6]


def interp2(curve, vq: float) -> float:
    pts = sorted(curve, key=lambda x: x[0])
    if len(pts) < 3:
        return 0.0
    i = min(range(len(pts)), key=lambda k: abs(pts[k][0] - vq))
    j = max(0, min(i - 1, len(pts) - 3))
    xs = [pts[j + k][0] for k in range(3)]
    ys = [pts[j + k][1] for k in range(3)]
    if any(abs(xs[a] - xs[b]) < 1e-14
           for a in range(3) for b in range(a + 1, 3)):
        return 0.0
    s = 0.0
    for a in range(3):
        t = ys[a]
        for b in range(3):
            if a != b:
                t *= (vq - xs[b]) / (xs[a] - xs[b])
        s += t
    return s


def check_logs(doc: dict, cases: Path) -> tuple[list[str], list[str]]:
    lines, errs = [], []
    for cm in doc["cases"]:
        path = cases / cm["out"]
        if not path.exists():
            errs.append("missing ANSYS log " + cm["out"])
            continue
        txt = path.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
        ne = int(m.group(1)) if m else txt.count("*** ERROR ***")
        m = re.search(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
        nw = int(m.group(1)) if m else txt.count("*** WARNING ***")
        lines.append("   {:<16s} log: {} error / {} warning".format(
            cm["id"], ne, nw))
        if ne != 0:
            errs.append("{} errors={}".format(cm["id"], ne))
    return lines, errs


def read_raw_csv(path: Path) -> list[list[float]]:
    out = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        p = line.split()
        if len(p) >= 5:
            try:
                out.append([float(x) for x in p[:5]])
            except ValueError:
                pass
    return out


def live_ansys(doc: dict, cases: Path, work: Path,
               ansys: Path) -> tuple[list[str], list[str]]:
    lines = ["-- live ANSYS re-run (temporary dir, frozen CSV untouched) --"]
    errs: list[str] = []
    total, worst = 0, 0.0
    for cm in doc["cases"]:
        job, inp = cm["id"], cm["apdl"]
        src = cases / inp
        if not src.exists():
            errs.append("missing APDL " + inp)
            continue
        shutil.copy2(src, work / inp)
        p = subprocess.run(
            [str(ansys), "-b", "-j", job, "-i", inp, "-o", job + "_live.out"],
            cwd=str(work), env=env_with_msys(), capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=1800)
        got = work / (job + "_curve.csv")
        ref = cases / cm["csv"]
        if p.returncode != 0 or not got.exists() or not ref.exists():
            errs.append("{} rc={} csv={}".format(job, p.returncode, got.exists()))
            continue
        a, b = read_raw_csv(got), read_raw_csv(ref)
        if len(a) != len(b):
            errs.append("{} row count {} vs {}".format(job, len(a), len(b)))
            continue
        jw = 0.0
        for ra, rb in zip(a, b):
            for x, y in zip(ra, rb):
                total += 1
                d = abs(x - y)
                jw = max(jw, d)
                worst = max(worst, d)
                if d > TOL_LIVE_CSV + TOL_LIVE_CSV_REL * abs(y):
                    errs.append("{} CSV {:.12g} vs {:.12g}".format(job, x, y))
        lines.append("   {:<16s} rows={:4d} max|d|={:.3e}".format(
            job, len(a), jw))
    lines.append("   total {} CSV values, max|d|={:.3e}".format(total, worst))
    return lines, errs[:6]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", type=Path, required=True)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys-exe", type=Path, default=Path(DEFAULT_ANSYS))
    a = ap.parse_args()

    doc = json.loads((a.reference / "cb3d17_snap_through_ansys.json").read_text(
        encoding="utf-8"))
    cfg = doc["config"]
    crit = doc["criteria"]
    tol = doc["tolerances"]
    lines: list[str] = []
    errs: list[str] = []

    lines.append("CB3D17 post-buckling / snap-through: FEM arc-length vs ANSYS ARCLEN")
    lines.append("   {} golden checks, {} cases, {} criteria".format(
        len(doc["checks"]), len(doc["cases"]), len(crit)))

    # 1) golden self-consistency
    if not check(doc):
        errs.append("golden reference drifted from frozen ANSYS CSV")

    # 2) FEM reproduction
    exe = a.bin / "snap_probe.exe"
    with tempfile.TemporaryDirectory(prefix="snap_bench_") as td:
        work = Path(td)
        try:
            fresh = run_fem_probe(exe, a.cases / SNAP_FILE, work)
        except Exception as exc:                            # noqa: BLE001
            errs.append("FEM probe: {}".format(exc))
            fresh = ({}, {})
        frozen = read_fem_csv(a.cases / FEM_CSV)
        if fresh[0]:
            n, worst, rep = reproduction(fresh, frozen)
            lines.append("   FEM reproduction: {} values, worst|d|={:.3e}".format(
                n, worst))
            errs.extend(rep)

        # 3) ANSYS frozen values vs independent analytic anchor
        worst_dp = 0.0
        worst_lt = 0.0
        nchk = 0
        for c in doc["checks"]:
            if c["key"] == "LAM" and c["analytic"] is not None:
                worst_dp = max(worst_dp,
                               abs(c["ref"] - c["analytic"]) * cfg["PREF"])
                nchk += 1
            elif c["key"] == "TIME" and c["analytic"] is not None:
                worst_lt = max(worst_lt, abs(c["ref"] - c["analytic"]))
        lines.append("   ANSYS vs analytic P_log(v): {} substeps, "
                     "worst|dP|={:.4e} N (limit {:.1e})".format(
                         nchk, worst_dp, tol["dp_ansys"]))
        if worst_dp > tol["dp_ansys"]:
            errs.append("ANSYS analytic worst dP {:.4e} > {:.1e}".format(
                worst_dp, tol["dp_ansys"]))
        lines.append("   ANSYS lambda_reaction vs TIME: worst|d|={:.3e} "
                     "(limit {:.1e})".format(worst_lt, tol["lambda_time"]))
        if worst_lt > tol["lambda_time"]:
            errs.append("lambda/TIME mismatch {:.3e}".format(worst_lt))

        # 4) fresh FEM path criteria
        if fresh[1]:
            fem_worst = max(s["worst_dP"] for s in fresh[1].values())
            fit_worst = max(abs(s["lambda_max_fit"]
                                - s["analytic_lambda_cr"])
                            for s in fresh[1].values())
            all_conv = all(s["converged"] > 0.5 and s["error_empty"] > 0.5
                           for s in fresh[1].values())
            lines.append("   FEM vs analytic: worst|dP|={:.4e} N, "
                         "worst|lambda*_fit-analytic|={:.3e}, "
                         "converged={}".format(fem_worst, fit_worst, all_conv))
            if fem_worst > tol["dp_fem"]:
                errs.append("FEM analytic worst dP {:.4e}".format(fem_worst))
            if fit_worst > tol["extremum"]:
                errs.append("FEM limit-point fit {:.3e}".format(fit_worst))
            if not all_conv:
                errs.append("some FEM arc-length paths did not converge")

        # 5) FEM arc-length curve vs ANSYS ARCLEN curve (quadratic interp)
        curve_worst = 0.0
        dcases = 0
        if fresh[0]:
            for cm in doc["cases"]:
                if cm["keyopt2"] != 1 or not cm["fem"]:
                    continue
                fcur = [(r[2], r[1]) for r in fresh[0][cm["fem"]]]
                acur = doc["curves"][cm["id"]]
                fvmax = max(v for v, _ in fcur)
                d = max(abs(lam - interp2(fcur, v))
                        for v, lam in acur
                        if 0.0 <= v <= min(2 * cfg["H"], fvmax))
                curve_worst = max(curve_worst, d)
                dcases += 1
        lines.append("   FEM arc-length vs ANSYS ARCLEN: {} curve pairs, "
                     "worst|dlambda|={:.3e} (limit {:.1e})".format(
                         dcases, curve_worst, tol["curve"]))
        if dcases < 2 or curve_worst > tol["curve"]:
            errs.append("FEM/ANSYS curve mismatch {:.3e}".format(curve_worst))

        # 6) descending branch + strain contrast
        main_d = doc["derived"]["sn_log_fine"]
        desc = crit["descending_branch"]
        if not (main_d["lambda_min_fit"] < desc["lambda_min_limit"]
                and main_d["v_max"] >= desc["v_min_limit"]):
            errs.append("descending branch not traversed")
        lines.append("   descending branch: lambda_min_fit={:+.7f} at v={:.6f}, "
                     "v_max={:.6f}".format(
                         main_d["lambda_min_fit"],
                         main_d["v_at_lambda_min_fit"], main_d["v_max"]))
        sc = crit["strain_measure_contrast"]["value"]
        lines.append("   strain-measure contrast KEYOPT(2)=0 vs 1: {:+.3f}% "
                     "(constitutive, not algorithmic)".format(sc))
        if not (crit["strain_measure_contrast"]["lo"] <= sc
                <= crit["strain_measure_contrast"]["hi"]):
            errs.append("strain contrast out of expected band")

        # 7) archived ANSYS logs
        lg, le = check_logs(doc, a.cases)
        lines.extend(lg)
        errs.extend(le)

        # 8) optional live ANSYS
        if a.live_ansys:
            lv, le = live_ansys(doc, a.cases, work, a.ansys_exe)
            lines.extend(lv)
            errs.extend(le)

    lines.append("")
    if errs:
        lines.append("[FAIL] CB3D17 snap-through benchmark")
        for e in errs[:12]:
            lines.append("   " + e)
    else:
        lines.append("[PASS] CB3D17 snap-through benchmark: "
                     "{} checks, {} criteria, {} ANSYS jobs".format(
                         len(doc["checks"]), len(crit), len(doc["cases"])))
    print("\n".join(lines))
    return 1 if errs else 0


if __name__ == "__main__":
    raise SystemExit(main())
