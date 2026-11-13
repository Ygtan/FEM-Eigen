#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: 3D hybrid-element (BEAM + CABLE/TRUSS) ANSYS benchmark.

For every ``reference/*.json``:
  1. verify the archived case (dat / inp / frozen ANSYS CSV) is complete;
  2. re-run the CURRENT build's ``example_3d_report.exe`` on the archived
     .dat, parse the fresh .report.txt, and compare every check item against
     the frozen ANSYS reference values;
  3. optionally (``--live-ansys``) re-run the archived APDL with ANSYS 18.0
     and compare every generated CSV row/column against the frozen CSV.

The companion nonlinear chain (``cable_solver.exe`` / fem_3d_nl.hpp) is
checked by ``compare_cable_solver_smoke.py`` and registered separately.

All runtime output is ASCII so CTest logs / PowerShell / redirects stay clean.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_ANSYS = Path(
    r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
MSYS_BIN = r"C:\msys64\mingw64\bin"

FORCE_KEYS = ("N", "Vy", "Vz", "T", "My", "Mz")


def read_csv_table(path: Path) -> dict[int, list[float]]:
    rows: dict[int, list[float]] = {}
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        toks = raw.replace(",", " ").split()
        if len(toks) < 2:
            continue
        try:
            vals = [float(t) for t in toks]
        except ValueError:
            continue
        rows[int(round(vals[0]))] = vals[1:]
    return rows


def parse_fem_report(path: Path) -> dict:
    """Self-contained parser for example_3d_report.exe output.

    No BridgeStudio import on purpose: FEM-Eigen CTest must not depend on
    an external project on disk.
    """
    disp: dict[int, list[float]] = {}
    force: dict[int, dict[str, float]] = {}
    kinds: dict[int, str] = {}
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "节点位移" in line or "node displacement" in line.lower():
            section = "disp"
            continue
        if "单元内力" in line or "element force" in line.lower():
            section = "force"
            continue
        if "单元应力" in line or "支座反力" in line:
            section = ""
            continue
        toks = line.replace(",", " ").split()
        try:
            key = int(toks[0])
        except (ValueError, IndexError):
            continue
        if section == "disp":
            try:
                vals = [float(v) for v in toks[1:7]]
            except ValueError:
                continue
            if vals:
                disp[key] = (vals + [0.0] * 6)[:6]
        elif section == "force":
            kv: dict[str, float] = {}
            for t in toks[2:]:
                if "=" in t:
                    k, _, v = t.partition("=")
                    try:
                        kv[k] = float(v)
                    except ValueError:
                        pass
            if any(k in kv for k in FORCE_KEYS):
                force[key] = {k: kv.get(k, 0.0) for k in FORCE_KEYS}
                kinds[key] = toks[1] if len(toks) > 1 else ""
    return {"disp": disp, "force": force, "kinds": kinds}


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def run(cmd: list[str], cwd: Path) -> tuple[int, str]:
    p = subprocess.run([str(c) for c in cmd], cwd=str(cwd), env=env_with_msys(),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=600)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def archive_checks(case: dict, cases_dir: Path) -> list[str]:
    errs: list[str] = []
    stem = case["id"]
    must = [stem + ".dat", stem + ".inp", stem + ".report.txt",
           stem + ".model.json", case["csv"]["disp"], case["csv"]["reac"],
           case["csv"]["beam"], case["csv"]["link"]]
    for name in must:
        if not (cases_dir / name).exists():
            errs.append("missing %s" % name)
    if errs:
        return errs
    for kind, want in case["expected_rows"].items():
        got = len(read_csv_table(cases_dir / case["csv"][kind]))
        if got != want:
            errs.append("%s has %d rows, expects %d" %
                        (case["csv"][kind], got, want))
    return errs


def live_ansys(case: dict, cases_dir: Path, ansys_exe: Path,
               work: Path) -> tuple[list[str], list[str]]:
    """Re-run archived APDL and diff every CSV against the frozen copies."""
    errs: list[str] = []
    lines: list[str] = []
    stem = case["id"]
    inp = cases_dir / (stem + ".inp")
    shutil.copy2(inp, work / inp.name)
    rc, out = run([ansys_exe, "-b", "-j", stem, "-i", inp.name,
                  "-o", stem + "_live.out"], work)
    if rc != 0:
        return ["ANSYS batch failed rc=%d: %s" % (rc, out[-400:])], lines
    log = (work / (stem + "_live.out")).read_text(
        encoding="utf-8", errors="replace")
    import re
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nerr = int(m.group(1)) if m else 0
    m = re.search(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nwarn = int(m.group(1)) if m else 0
    if nerr:
        errs.append("ANSYS %s_live.out: %d error(s)" % (stem, nerr))
    lines.append("  LIVE ANSYS %-26s log: %d error / %d warning" %
                 (stem, nerr, nwarn))
    total = 0
    worst = 0.0
    for kind, frozen_name in case["csv"].items():
        fresh = work / (stem + "_%s.csv" % kind)
        if not fresh.exists():
            errs.append("LIVE %s: %s not produced" % (kind, fresh.name))
            continue
        a = read_csv_table(fresh)
        b = read_csv_table(cases_dir / frozen_name)
        if set(a) != set(b):
            errs.append("LIVE %s: row keys differ (%d vs %d)" %
                        (kind, len(a), len(b)))
            continue
        for k in sorted(a):
            if len(a[k]) != len(b[k]):
                errs.append("LIVE %s row %d: column count differs" % (kind, k))
                continue
            for c, (x, y) in enumerate(zip(a[k], b[k])):
                total += 1
                d = abs(x - y)
                scale = max(abs(x), abs(y), 1.0)
                worst = max(worst, d / scale)
                if d > 1e-6 + 1e-3 * abs(y):
                    errs.append("LIVE %s row %d col %d: %g vs %g" %
                                (kind, k, c + 1, x, y))
    lines.append("  LIVE ANSYS %-26s %d values, max rel diff %.3e  %s" %
                 (stem, total, worst, "PASS" if not errs else "FAIL"))
    # force equilibrium of the freshly computed reactions
    reac = read_csv_table(work / (stem + "_reac.csv"))
    fz = sum(reac.get(nd, [0.0] * 6)[2] for nd in case["support_nodes"])
    want = case["support_fz_sum"]
    if abs(fz - want) > 1e-3 * max(abs(want), 1.0):
        errs.append("LIVE %s reaction sum %.6g vs frozen %.6g" % (stem, fz, want))
    else:
        lines.append("  LIVE ANSYS %-26s support sum FZ %.6g N  PASS" % (stem, fz))
    return errs, lines


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="FEM-Eigen build dir")
    ap.add_argument("--cases", default=str(HERE / "cases"))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys-exe", default=str(DEFAULT_ANSYS))
    a = ap.parse_args()

    bin_dir = Path(a.bin)
    fem_exe = bin_dir / "example_3d_report.exe"
    cases_dir = Path(a.cases)
    ref_dir = Path(a.reference)
    ansys_exe = Path(a.ansys_exe)

    print("== cable_arch_3d benchmark ==")
    print("   exe: %s" % fem_exe)
    print("   cases: %s" % cases_dir)
    if not fem_exe.exists():
        print("[FAIL] solver not found: %s" % fem_exe)
        return 1
    refs = sorted(ref_dir.glob("*_ansys.json"))
    if not refs:
        print("[FAIL] no reference JSON in %s" % ref_dir)
        return 1

    npass = nfail = 0
    all_errs: list[str] = []
    for rf in refs:
        case = json.loads(rf.read_text(encoding="utf-8"))
        stem = case["id"]
        print("== %s ==" % stem)
        errs = archive_checks(case, cases_dir)
        if errs:
            for e in errs:
                print("[FAIL] archive: %s" % e)
            nfail += len(errs)
            all_errs.extend(errs)
            continue

        with tempfile.TemporaryDirectory(prefix="carch_") as td:
            work = Path(td)
            shutil.copy2(cases_dir / (stem + ".dat"), work / (stem + ".dat"))
            rc, out = run([fem_exe, stem + ".dat"], work)
            if rc != 0:
                print("[FAIL] solver rc=%d %s" % (rc, out[-300:]))
                nfail += 1
                continue
            report = work / (stem + ".report.txt")
            if not report.exists():
                print("[FAIL] no .report.txt produced")
                nfail += 1
                continue
            fem = parse_fem_report(report)

            for chk in case["checks"]:
                name = "%s:%s" % (stem, chk["name"])
                k, col = chk["key"], chk["col"]
                if chk["kind"] == "disp":
                    row = fem["disp"].get(k)
                    got = row[col] if row else None
                elif chk["kind"] == "beam":
                    f = fem["force"].get(k)
                    got = f[FORCE_KEYS[col]] if f else None
                else:  # link
                    f = fem["force"].get(k)
                    got = f["N"] if f else None
                if got is None:
                    nfail += 1
                    all_errs.append("%s: no FEM value" % name)
                    print("[FAIL] %-40s no FEM value" % name)
                    continue
                if chk.get("abs"):
                    got = abs(got)
                ref = chk["ref"]
                d = abs(got - ref)
                tol = chk["atol"] + chk["rtol"] * abs(ref)
                if d <= tol:
                    npass += 1
                else:
                    nfail += 1
                    all_errs.append("%s got=%g ref=%g d=%g tol=%g" %
                                    (name, got, ref, d, tol))
                    print("[FAIL] %-40s got=%-14.8g ref=%-14.8g d=%-10.3g tol=%g" %
                          (name, got, ref, d, tol))

            if a.live_ansys:
                if not ansys_exe.exists():
                    print("[FAIL] ANSYS exe not found: %s" % ansys_exe)
                    nfail += 1
                else:
                    lerrs, llines = live_ansys(case, cases_dir, ansys_exe, work)
                    for s in llines:
                        print(s)
                    all_errs.extend(lerrs)
                    nfail += len(lerrs)

        print("   %d checks" % len(case["checks"]))

    print("== Result: %d/%d PASS, %d FAIL ==" % (npass, npass + nfail, nfail))
    if all_errs:
        print("   first failures:")
        for e in all_errs[:8]:
            print("     " + e)
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
