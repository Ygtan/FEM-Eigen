#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: FEM-Eigen nonlinear completion state vs frozen ANSYS NLGEOM.

Offline by default (no ANSYS licence needed):

  * re-run ``cable_solver.exe`` (fem_3d_nl.hpp co-rotational chain) on the
    archived COMPLETION .dat;
  * compare all 300 [U] DOF, 16 [CABLE] forces and 9 [REACTION] components
    against the golden values frozen from a real ANSYS 18.0
    BEAM4 + LINK180 + NLGEOM + INISTATE batch run;
  * verify the fresh FEM result reproduces the archived ``*_cable.txt``;
  * verify every frozen CSV cell still matches the reference JSON.

``--live-ansys`` additionally re-runs the archived APDL and diffs every CSV
cell against the frozen copies.  The one documented ANSYS warning
(LINK180 KEYOPT(3)=1 "undocumented option") is an expected baseline, not a
failure; it is checked explicitly rather than silently ignored.

All runtime output is ASCII.
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
MSYS_BIN = r"C:\msys64\mingw64\bin"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


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


def parse_cable_result(path: Path) -> dict:
    """FEM completion state: [U] / [CABLE] / [REACTION] + header scalars."""
    u: dict[tuple[int, int], float] = {}
    cables: dict[int, dict] = {}
    reaction: dict[tuple[int, int], float] = {}
    meta: dict[str, float] = {}
    sec = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("["):
            sec = line.strip("[]").upper()
            continue
        parts = line.split()
        if sec == "":
            for tok in parts:
                if "=" not in tok:
                    continue
                k, _, v = tok.partition("=")
                try:
                    meta[k] = float(v)
                except ValueError:
                    pass
            continue
        try:
            if sec == "U" and len(parts) >= 3:
                u[(int(parts[0]), int(parts[1]))] = float(parts[2])
            elif sec == "CABLE" and len(parts) >= 7:
                cables[int(parts[0])] = {
                    "anchor": int(parts[1]), "T": float(parts[2]),
                    "L0": float(parts[3]), "chord": float(parts[4]),
                }
            elif sec == "REACTION" and len(parts) >= 3:
                reaction[(int(parts[0]), int(parts[1]))] = float(parts[2])
        except ValueError:
            continue
    return {"u": u, "cables": cables, "reaction": reaction, "meta": meta}


def run_fem(exe: Path, dat: Path, work: Path) -> dict:
    local = work / dat.name
    shutil.copy2(dat, local)
    p = subprocess.run([str(exe), str(local)], cwd=str(work),
                       env=env_with_msys(), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600)
    if p.returncode != 0:
        raise RuntimeError("cable_solver rc=%d: %s" %
                           (p.returncode, (p.stderr or p.stdout or "")[-300:]))
    fresh = work / (dat.stem + "_cable.txt")
    if not fresh.exists():
        raise RuntimeError("no %s produced" % fresh.name)
    return parse_cable_result(fresh)


def extract(res: dict, chk: dict):
    kind = chk["kind"]
    if kind == "disp":
        return res["u"].get((chk["key"], chk["dof"]))
    if kind == "reac":
        return res["reaction"].get((chk["key"], chk["dof"]))
    if kind == "force":
        c = res["cables"].get(chk["key"])
        return None if c is None else c["T"]
    return None


def archive_checks(doc: dict, cases: Path, tables: dict) -> list[str]:
    errs: list[str] = []
    for name in doc.get("archive", []):
        if not (cases / name).exists():
            errs.append("missing archive file " + name)
    for kind, want in doc.get("expected_rows", {}).items():
        got = len(tables.get(kind, {}))
        if got != want:
            errs.append("%s CSV has %d rows, expected %d" % (kind, got, want))
    frozen = parse_cable_result(cases / doc["fem_result"])
    for key, want in doc.get("fem_rows", {}).items():
        got = len(frozen["u"] if key == "u" else
                  frozen["cables"] if key == "cables" else
                  frozen["reaction"])
        if got != want:
            errs.append("frozen %s has %d rows, expected %d" % (key, got, want))
    counts: dict[str, int] = {}
    for chk in doc["checks"]:
        counts[chk["kind"]] = counts.get(chk["kind"], 0) + 1
    for kind, want in doc.get("check_counts", {}).items():
        if counts.get(kind, 0) != want:
            errs.append("checks[%s]=%d, expected %d" %
                        (kind, counts.get(kind, 0), want))
    return errs


def frozen_integrity(doc: dict, tables: dict) -> list[str]:
    """Every ref must equal the archived ANSYS CSV cell (guard against drift)."""
    errs: list[str] = []
    for chk in doc["checks"]:
        kind = chk["kind"]
        src = {"disp": "disp", "reac": "reac", "force": "link"}[kind]
        row = tables[src].get(chk["key"])
        if row is None or len(row) <= chk["col"]:
            errs.append("%s: %s CSV row %d col %d missing" %
                        (chk["name"], src, chk["key"], chk["col"] + 1))
            continue
        val = abs(row[chk["col"]]) if chk.get("abs") else row[chk["col"]]
        if abs(val - chk["ref"]) > 1e-9 + 1e-12 * abs(chk["ref"]):
            errs.append("%s: frozen CSV %g != ref %g" %
                        (chk["name"], val, chk["ref"]))
    return errs


def frozen_reproduction(fresh: dict, frozen: dict) -> tuple[int, list[str], float]:
    """Tight reproducibility check: re-run FEM == archived *_cable.txt."""
    errs: list[str] = []
    n = 0
    worst = 0.0
    for key, label in (("u", "[U]"), ("reaction", "[REACTION]")):
        for k, old in frozen[key].items():
            new = fresh[key].get(k)
            if new is None:
                errs.append("%s missing key %s" % (label, k))
                continue
            n += 1
            d = abs(new - old)
            tol = 1e-9 + 1e-6 * abs(old)
            worst = max(worst, d / max(abs(old), 1e-30))
            if d > tol:
                errs.append("%s key %s: %g vs frozen %g" % (label, k, new, old))
    for eid, old in frozen["cables"].items():
        new = fresh["cables"].get(eid)
        if new is None:
            errs.append("[CABLE] missing elem %d" % eid)
            continue
        n += 1
        d = abs(new["T"] - old["T"])
        tol = 1e-6 + 1e-6 * abs(old["T"])
        worst = max(worst, d / max(abs(old["T"]), 1e-30))
        if d > tol:
            errs.append("[CABLE] e%d: %g vs frozen %g" %
                       (eid, new["T"], old["T"]))
    return n, errs, worst


def compare_one(doc: dict, exe: Path, cases: Path, work: Path):
    stem = doc["id"]
    lines = ["== %s ==" % stem]
    tables = {k: read_csv_table(cases / v) for k, v in doc["csv"].items()}
    errs = archive_checks(doc, cases, tables)
    errs.extend(frozen_integrity(doc, tables))
    fresh = run_fem(exe, cases / doc["dat"], work)
    frozen = parse_cable_result(cases / doc["fem_result"])

    md = fresh["meta"]
    if int(md.get("converged", 0)) != int(doc["converged"]):
        errs.append("converged=%g, expected %d" %
                    (md.get("converged", 0), doc["converged"]))
    ce = float(md.get("closure_err", float("nan")))
    te = float(md.get("tension_err", float("nan")))
    if not (ce < 1e-6 and te < 1e-6):
        errs.append("closure/tension not converged: %.3e / %.3e" % (ce, te))
    lines.append("   FEM converged=%d iters=%g closure=%.3e tension=%.3e" %
                 (int(md.get("converged", 0)), md.get("iterations", -1), ce, te))

    nrep, rerrs, rworst = frozen_reproduction(fresh, frozen)
    errs.extend(rerrs[:5])
    lines.append("   FEM reproduction: %d values, worst rel=%.3e %s" %
                 (nrep, rworst, "PASS" if not rerrs else "FAIL"))

    npass = nfail = 0
    worst_ratio = 0.0
    for chk in doc["checks"]:
        got = extract(fresh, chk)
        if got is None:
            errs.append("%s: missing FEM value" % chk["name"])
            nfail += 1
            continue
        if chk.get("abs"):
            got = abs(got)
        ref = chk["ref"]
        scale = chk.get("scale", abs(ref))
        tol = chk.get("atol", 1e-6) + chk.get("rtol", 1e-3) * scale
        d = abs(got - ref)
        worst_ratio = max(worst_ratio, d / max(tol, 1e-30))
        if d <= tol:
            npass += 1
        else:
            nfail += 1
            if nfail <= 8:
                lines.append("[FAIL] %-34s got=%-15.8g ref=%-15.8g "
                             "d=%g tol=%g" % (chk["name"], got, ref, d, tol))
    lines.append("   ANSYS golden: %d/%d checks passed, worst d/tol=%.3f" %
                 (npass, npass + nfail, worst_ratio))
    lines.append("   group scales: %s" %
                 ", ".join("%s=%.6g" % (k, v)
                           for k, v in sorted(doc["group_scale"].items())))
    return lines, errs, nfail


def live_ansys(doc: dict, cases: Path, work: Path, ansys: Path):
    stem = doc["apdl_stem"]
    lines = ["== %s LIVE ANSYS cross-check ==" % stem]
    errs: list[str] = []
    inp = cases / doc["apdl"]
    shutil.copy2(inp, work / inp.name)
    p = subprocess.run([str(ansys), "-b", "-j", stem, "-i", inp.name,
                        "-o", stem + "_live.out"], cwd=str(work),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=900)
    if p.returncode != 0:
        return lines, ["ANSYS rc=%d: %s" %
                       (p.returncode, (p.stdout or p.stderr or "")[-300:])]
    log = (work / (stem + "_live.out")).read_text(
        encoding="utf-8", errors="replace")
    exp = doc.get("ansys_log_expect", {})
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nerr = int(m.group(1)) if m else 0
    m = re.search(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nwarn = int(m.group(1)) if m else 0
    lines.append("   log: %d error / %d warning (expected %d / %d)" %
                 (nerr, nwarn, exp.get("errors", 0), exp.get("warnings", 0)))
    if nerr != exp.get("errors", 0):
        errs.append("ANSYS error count %d, expected %d" %
                    (nerr, exp.get("errors", 0)))
    if nwarn != exp.get("warnings", 0):
        errs.append("ANSYS warning count %d, expected %d" %
                    (nwarn, exp.get("warnings", 0)))
    sig = exp.get("warning_signature")
    if sig and sig not in log:
        errs.append("expected ANSYS warning signature missing: %s" % sig)

    # APDL writes <apdl_stem>_{disp,reac,link}.csv; the archive inserts
    # "_ansys" before the suffix, hence the name mapping below.
    total = 0
    worst = 0.0
    for kind, frozen_name in doc["csv"].items():
        fresh = work / ("%s_%s.csv" % (stem, kind))
        if not fresh.exists():
            errs.append("LIVE %s: %s not produced" % (kind, fresh.name))
            continue
        a, b = read_csv_table(fresh), read_csv_table(cases / frozen_name)
        if set(a) != set(b):
            errs.append("LIVE %s row keys differ (%d vs %d)" %
                        (kind, len(a), len(b)))
            continue
        for k in sorted(a):
            if len(a[k]) != len(b[k]):
                errs.append("LIVE %s row %d column count differs" % (kind, k))
                continue
            for x, y in zip(a[k], b[k]):
                total += 1
                d = abs(x - y)
                worst = max(worst, d)
                if d > 1e-6 + 1e-3 * abs(y):
                    errs.append("LIVE %s row %d: %g vs %g" % (kind, k, x, y))
    lines.append("   %d CSV values, max abs diff %.3e  %s" %
                 (total, worst, "PASS" if not errs else "FAIL"))
    return lines, errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="FEM-Eigen build dir")
    ap.add_argument("--cases", default=str(HERE / "cases"))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys", default=DEFAULT_ANSYS)
    a = ap.parse_args()
    exe = Path(a.bin) / "cable_solver.exe"
    cases, ref = Path(a.cases), Path(a.reference)
    if not exe.exists():
        print("[FAIL] cable_solver.exe not found: %s" % exe)
        return 1
    docs = [json.loads(p.read_text(encoding="utf-8"))
            for p in sorted(ref.glob("*_ansys.json"))]
    if not docs:
        print("[FAIL] no reference JSON found in %s" % ref)
        return 1
    print("== nonlinear completion vs ANSYS NLGEOM benchmark ==")
    all_errs: list[str] = []
    total_fail = 0
    with tempfile.TemporaryDirectory(prefix="femnl_") as td:
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
    print("Result: %d case(s), %d failure(s)" %
          (len(docs), total_fail + len(all_errs)))
    return 0 if (total_fail == 0 and not all_errs) else 1


if __name__ == "__main__":
    sys.exit(main())
