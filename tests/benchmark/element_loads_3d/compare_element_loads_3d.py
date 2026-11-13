#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: re-run FEM-Eigen 3D on E-stage .dat and compare with ANSYS.

Offline by default (no ANSYS licence needed): the golden values are frozen
from real ANSYS 18.0 BEAM4 batch CSVs. ``--live-ansys`` additionally re-runs
the archived APDL and compares its CSVs cell-by-cell with the frozen ones.

Kernel-level closed-form validation lives in verify_kernel_loads.py and is
registered as a separate CTest case.
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

from verify_kernel_loads import MSYS_BIN, parse_report  # noqa: E402

FORCE_KEYS = ("N", "Vy", "Vz", "T", "My", "Mz")
ARCHIVE = [
    "{stem}.dat", "{stem}.inp", "{stem}.model.json", "{stem}.report.txt",
    "{stem}_ansys_disp.csv", "{stem}_ansys_reac.csv", "{stem}_ansys_elem.csv",
    "{stem}_ansys.out", "{stem}_compare_report.txt", "{stem}_manifest.txt",
]
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"


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


def archive_checks(doc: dict, cases: Path) -> list[str]:
    errs = []
    stem = doc["id"]
    for pat in ARCHIVE:
        if not (cases / pat.format(stem=stem)).exists():
            errs.append("missing archive file " + pat.format(stem=stem))
    for kind, want in doc.get("expected_rows", {}).items():
        p = cases / doc["csv"][kind]
        if not p.exists():
            continue
        got = len(read_csv_table(p))
        if got != want:
            errs.append("%s has %d rows, expected %d" % (p.name, got, want))
    return errs


def extract(parsed: tuple[dict, dict], chk: dict):
    disp, force = parsed
    if chk["kind"] == "disp":
        row = disp.get(chk["key"])
        if row is None or len(row) <= chk["col"]:
            return None
        return row[chk["col"]]
    if chk["kind"] == "force":
        row = force.get(chk["key"])
        if row is None or chk["col"] >= len(FORCE_KEYS):
            return None
        return row.get(FORCE_KEYS[chk["col"]])
    return None


def run_fem(exe: Path, dat: Path, work: Path) -> tuple[dict, dict]:
    local = work / dat.name
    shutil.copy2(dat, local)
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    r = subprocess.run([str(exe), local.name], cwd=str(work), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=300)
    if r.returncode != 0:
        raise RuntimeError("solver rc=%d: %s" % (r.returncode, r.stderr[-300:]))
    rp = work / (dat.stem + ".report.txt")
    if not rp.exists():
        raise RuntimeError("no report.txt produced")
    return parse_report(rp.read_text(encoding="utf-8", errors="replace"))


def compare_one(doc: dict, exe: Path, cases: Path, work: Path):
    stem = doc["id"]
    lines = ["== %s ==" % stem]
    errs = archive_checks(doc, cases)
    parsed = run_fem(exe, cases / doc["dat"], work)
    npass = nfail = 0
    for chk in doc["checks"]:
        name = "%s:%s" % (stem, chk["name"])
        got = extract(parsed, chk)
        if got is None:
            lines.append("[FAIL] %-40s missing value" % name)
            nfail += 1
            continue
        if chk.get("abs"):
            got = abs(got)
        ref = chk["ref"]
        d = abs(got - ref)
        tol = chk.get("atol", 1e-9) + chk.get("rtol", 1e-3) * abs(ref)
        if d <= tol:
            npass += 1
        else:
            lines.append("[FAIL] %-40s got=%-15.8g ref=%-15.8g d=%g tol=%g" %
                         (name, got, ref, d, tol))
            nfail += 1
    lines.append("    %d/%d checks passed" % (npass, npass + nfail))
    return lines, errs, nfail


def live_ansys(doc: dict, cases: Path, work: Path, ansys: Path):
    """Re-run archived APDL, compare every cell with the frozen CSV."""
    stem = doc["id"]
    lines = ["== %s LIVE ANSYS cross-check ==" % stem]
    errs = []
    shutil.copy2(cases / doc["apdl"], work / doc["apdl"])
    r = subprocess.run([str(ansys), "-b", "-j", stem + "_live",
                        "-i", doc["apdl"], "-o", stem + "_live.out"],
                       cwd=str(work), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600)
    if r.returncode != 0:
        return lines, ["ANSYS rc=%d" % r.returncode]
    log = work / (stem + "_live.out")
    txt = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    for tag in ("WARNING", "ERROR"):
        for ln in txt.splitlines():
            if "NUMBER OF %s MESSAGES ENCOUNTERED=" % tag in ln:
                val = ln.split("=")[-1].strip()
                if val not in ("0", ""):
                    errs.append("ANSYS %s count = %s" % (tag, val))
                break
    for kind in ("disp", "reac", "elem"):
        fresh_name = "%s_%s.csv" % (stem, kind)
        fresh = work / fresh_name
        if not fresh.exists():
            errs.append("LIVE missing " + fresh_name)
            continue
        a = read_csv_table(fresh)
        b = read_csv_table(cases / doc["csv"][kind])
        if set(a) != set(b):
            errs.append("%s row keys differ" % kind)
            continue
        dmax = 0.0
        scale = 0.0
        for k in a:
            if len(a[k]) != len(b[k]):
                errs.append("%s row %d length differs" % (kind, k))
                break
            for x, y in zip(a[k], b[k]):
                dmax = max(dmax, abs(x - y))
                scale = max(scale, abs(y))
        ok = dmax <= 1e-3 * scale + 1e-6
        lines.append("    %-5s maxdiff=%-12.6g scale=%-12.6g %s" %
                     (kind, dmax, scale, "PASS" if ok else "FAIL"))
        if not ok:
            errs.append("LIVE %s out of tolerance" % kind)
    return lines, errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys", default=DEFAULT_ANSYS)
    a = ap.parse_args()
    exe = Path(a.bin) / "example_3d_report.exe"
    cases = Path(a.cases)
    ref = Path(a.reference)
    if not exe.exists():
        print("[FAIL] solver not found: %s" % exe)
        return 1
    docs = [json.loads(p.read_text(encoding="utf-8"))
            for p in sorted(ref.glob("*_ansys.json"))]
    if not docs:
        print("[FAIL] no reference JSON found in %s" % ref)
        return 1
    print("== 3D element-load / selfweight benchmark ==")
    all_errs: list[str] = []
    total_fail = 0
    with tempfile.TemporaryDirectory(prefix="fem3d_eload_") as td:
        work = Path(td)
        for doc in docs:
            try:
                lines, errs, nfail = compare_one(doc, exe, cases, work)
            except Exception as exc:                       # noqa: BLE001
                print("[FAIL] %s exception: %s" % (doc["id"], exc))
                all_errs.append(str(exc))
                total_fail += 1
                continue
            print("\n".join(lines))
            all_errs.extend(errs)
            total_fail += nfail
            if a.live_ansys:
                l2, e2 = live_ansys(doc, cases, work, Path(a.ansys))
                print("\n".join(l2))
                all_errs.extend(e2)
    for e in all_errs:
        print("[FAIL] " + e)
    print("Result: %d case(s), %d failure(s)" % (len(docs), total_fail + len(all_errs)))
    return 0 if (total_fail == 0 and not all_errs) else 1


if __name__ == "__main__":
    sys.exit(main())
