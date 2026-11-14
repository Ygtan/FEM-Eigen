#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: CB3D13 balanced self-anchored suspension + erection stages.

Offline by default (no ANSYS licence needed):

  * re-run ``cable_solver.exe`` on the archived ``*SUSP_BRIDGE + *ERECTION``
    .dat and verify it reproduces the frozen ``*_cable.txt`` value by value,
    including all four new ``[ERECTION_*]`` sections;
  * compare every [A] full-bridge displacement / cable force / reaction and
    every [C] stage shape / segment force / horizontal force against golden
    values frozen from real ANSYS 18.0 batch CSVs;
  * re-derive every [C] golden value from the frozen CSV (guards drift);
  * re-check the [B] balanced-side-span criterion from the fresh FEM run.

``--live-ansys`` additionally re-runs the full-bridge APDL and all five
stage APDLs, and diffs every CSV cell against the frozen copies.  The two
documented ``KEYOPT(3)=1`` warnings of the full-bridge model are an expected
baseline; the stage chains use KEYOPT(3)=0 and must be warning-free.

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
sys.path.insert(0, str(HERE))

from build_reference_erection import (  # noqa: E402
    parse_fem_result, read_csv_table, stage_ansys, stage_apdl_zinit,
)

MSYS_BIN = r"C:\msys64\mingw64\bin"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def run_fem(exe: Path, dat: Path, work: Path) -> dict:
    shutil.copy2(dat, work / dat.name)
    p = subprocess.run([str(exe), dat.name], cwd=str(work),
                       env=env_with_msys(), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600)
    if p.returncode != 0:
        raise RuntimeError("cable_solver rc=%d: %s" %
                           (p.returncode, (p.stderr or p.stdout or "")[-300:]))
    fresh = work / (dat.stem + "_cable.txt")
    if not fresh.exists():
        raise RuntimeError("no %s produced" % fresh.name)
    return parse_fem_result(fresh)


def extract(res: dict, chk: dict):
    k = chk["kind"]
    if k == "disp":
        return res["u"].get((chk["key"], chk["dof"]))
    if k == "reac":
        return res["reaction"].get((chk["key"], chk["dof"]))
    if k == "force":
        c = res["cables"].get(chk["key"])
        return None if c is None else c["T"]
    st = chk.get("stage")
    if k == "stage_shape":
        row = res["shape"].get(st)
        return None if row is None or chk["key"] >= len(row) else row[chk["key"]]
    if k == "stage_force":
        row = res["force"].get(st)
        return None if row is None or chk["key"] >= len(row) else row[chk["key"]]
    if k == "stage_H":
        row = res["erection"].get(st)
        return None if row is None else row.get("H")
    return None


def stage_tables(doc: dict, cases: Path) -> dict[int, dict]:
    out: dict[int, dict] = {}
    for st in doc["stages"]:
        sf = doc["stage_files"][str(st)]
        zi = stage_apdl_zinit(cases / sf["apdl"])
        d = read_csv_table(cases / sf["disp"])
        lk = read_csv_table(cases / sf["link"])
        out[st] = {"z_init": zi, "disp": d, "link": lk,
                   "ansys": stage_ansys(d, lk, xs=doc["chain_x"], z_init=zi)}
    return out


def archive_checks(doc: dict, cases: Path, tables: dict,
                   stages: dict[int, dict]) -> list[str]:
    errs: list[str] = []
    for name in doc.get("archive", []):
        if not (cases / name).exists():
            errs.append("missing archive file " + name)
    for kind, want in doc.get("expected_rows", {}).items():
        got = len(tables.get(kind, {}))
        if got != want:
            errs.append("%s CSV has %d rows, expected %d" % (kind, got, want))
    frozen = parse_fem_result(cases / doc["fem_result"])
    fr = doc.get("fem_rows", {})
    for key in ("u", "cables", "reaction", "erection", "chain"):
        got = len(frozen[key]) if key != "chain" else len(frozen["chain"])
        if got != fr.get(key):
            errs.append("frozen %s has %d rows, expected %s" %
                        (key, got, fr.get(key)))
    for key in ("shape", "load", "force"):
        for st, want in fr.get(key, {}).items():
            got = len(frozen[key].get(int(st), []))
            if got != want:
                errs.append("frozen %s[%s] has %d cols, expected %d" %
                           (key, st, got, want))
    counts: dict[str, int] = {}
    for chk in doc["checks"]:
        counts[chk["kind"]] = counts.get(chk["kind"], 0) + 1
    for kind, want in doc.get("check_counts", {}).items():
        if counts.get(kind, 0) != want:
            errs.append("checks[%s]=%d, expected %d" %
                        (kind, counts.get(kind, 0), want))
    n = len(doc["chain_x"])
    for st, sd in stages.items():
        if len(sd["z_init"]) != n:
            errs.append("stage %d APDL has %d nodes, expected %d" %
                        (st, len(sd["z_init"]), n))
        if len(sd["disp"]) != n or len(sd["link"]) != n - 1:
            errs.append("stage %d CSV rows %d/%d, expected %d/%d" %
                        (st, len(sd["disp"]), len(sd["link"]), n, n - 1))
    return errs


def frozen_integrity(doc: dict, tables: dict,
                     stages: dict[int, dict]) -> list[str]:
    """Every ref must equal a frozen ANSYS CSV cell (or its stated derivation)."""
    errs: list[str] = []
    srcmap = {"disp": "disp", "reac": "reac", "force": "link"}
    for chk in doc["checks"]:
        k = chk["kind"]
        if k in srcmap:
            row = tables[srcmap[k]].get(chk["key"])
            if row is None or len(row) <= chk["col"]:
                errs.append("%s: CSV row %d col %d missing" %
                            (chk["name"], chk["key"], chk["col"] + 1))
                continue
            val = abs(row[chk["col"]]) if chk.get("abs") else row[chk["col"]]
        elif k == "stage_shape":
            sd = stages[chk["stage"]]
            row = sd["disp"].get(chk["key"] + 1)
            if row is None or len(row) < 3:
                errs.append("%s: stage disp row missing" % chk["name"])
                continue
            val = sd["z_init"][chk["key"]] + row[2]
        elif k == "stage_force":
            sd = stages[chk["stage"]]
            row = sd["link"].get(chk["key"] + 1)
            if row is None or not row:
                errs.append("%s: stage link row missing" % chk["name"])
                continue
            val = abs(row[0])
        elif k == "stage_H":
            val = abs(stages[chk["stage"]]["ansys"]["H"])
        else:
            errs.append("%s: unknown kind %s" % (chk["name"], k))
            continue
        if abs(val - chk["ref"]) > 1e-9 + 1e-12 * abs(chk["ref"]):
            errs.append("%s: frozen CSV %g != ref %g" %
                        (chk["name"], val, chk["ref"]))
    return errs


def frozen_reproduction(fresh: dict, frozen: dict) -> tuple[int, list[str], float]:
    errs: list[str] = []
    n = 0
    worst = 0.0

    def cmp(label: str, a, b, tol: float) -> None:
        nonlocal n, worst
        if a is None or b is None:
            errs.append("%s missing value" % label)
            return
        n += 1
        d = abs(a - b)
        worst = max(worst, d / max(abs(b), 1e-30))
        if d > tol:
            errs.append("%s: %g vs frozen %g" % (label, a, b))

    for key, label in (("u", "[U]"), ("reaction", "[REACTION]")):
        for k, old in frozen[key].items():
            cmp("%s %s" % (label, k), fresh[key].get(k), old,
                1e-9 + 1e-6 * abs(old))
    for eid, old in frozen["cables"].items():
        new = fresh["cables"].get(eid)
        cmp("[CABLE] e%d T" % eid, None if new is None else new["T"],
            old["T"], 1e-6 + 1e-6 * abs(old["T"]))
    for st, row in frozen["erection"].items():
        f = fresh["erection"].get(st, {})
        for fld in ("H", "sag", "midz", "hmin", "hmax"):
            cmp("[ERECTION] s%d %s" % (st, fld), f.get(fld), row.get(fld),
                1e-6 + 1e-6 * abs(row.get(fld, 0.0)))
    for key in ("shape", "load", "force"):
        for st, old in frozen[key].items():
            new = fresh[key].get(st, [])
            if len(new) != len(old):
                errs.append("[%s] stage %d length %d vs %d" %
                           (key.upper(), st, len(new), len(old)))
                continue
            for k, v in enumerate(old):
                cmp("[%s] s%d k%d" % (key.upper(), st, k), new[k], v,
                    1e-9 + 1e-6 * abs(v))
    if len(frozen["chain"]) != len(fresh["chain"]):
        errs.append("[ERECTION_CHAIN] length differs")
    else:
        for k, (oid, ox) in enumerate(frozen["chain"]):
            nid, nx = fresh["chain"][k]
            cmp("[CHAIN] k%d id" % k, float(nid), float(oid), 0)
            cmp("[CHAIN] k%d x" % k, nx, ox, 1e-9 + 1e-6 * abs(ox))
    return n, errs, worst


def balance_check(fresh: dict, doc: dict) -> tuple[str, list[str]]:
    exp = doc["balance_expect"]
    nodes = exp["end_nodes"]
    got = [fresh["reaction"].get((nd, 3), 0.0) for nd in nodes]
    ratio = max(abs(v) for v in got) / (
        exp["q_deck_N_per_m"] * exp["main_span_m"])
    errs: list[str] = []
    for nd, a, b in zip(nodes, got, exp["end_rz_N"]):
        if abs(a - b) > 1e-3 + 1e-6 * abs(b):
            errs.append("[B] end node %d Rz %g != frozen %g" % (nd, a, b))
    if ratio > exp["tol"]:
        errs.append("[B] max|Rz|/(q*L) = %.3e > %.3e" % (ratio, exp["tol"]))
    return ("   [B] balanced side span: max|Rz|/(q*L)=%.3e (<=%.1e) %s" %
            (ratio, exp["tol"], "PASS" if not errs else "FAIL"), errs)


def compare_one(doc: dict, exe: Path, cases: Path, work: Path):
    lines = ["== %s ==" % doc["id"]]
    tables = {k: read_csv_table(cases / v) for k, v in doc["csv"].items()}
    stages = stage_tables(doc, cases)
    errs = archive_checks(doc, cases, tables, stages)
    errs.extend(frozen_integrity(doc, tables, stages))
    fresh = run_fem(exe, cases / doc["dat"], work)
    frozen = parse_fem_result(cases / doc["fem_result"])

    md = fresh["meta"]
    if int(md.get("converged", 0)) != int(doc["converged"]):
        errs.append("converged=%g, expected %d" %
                    (md.get("converged", 0), doc["converged"]))
    ce = float(md.get("closure_err", float("nan")))
    te = float(md.get("tension_err", float("nan")))
    if not (ce < doc["closure_tol"] and te < doc["tension_tol"]):
        errs.append("closure/tension %.3e / %.3e exceed %.3e / %.3e" %
                    (ce, te, doc["closure_tol"], doc["tension_tol"]))
    lines.append("   FEM converged=%d iters=%g closure=%.3e tension=%.3e" %
                 (int(md.get("converged", 0)), md.get("iterations", -1), ce, te))

    nrep, rerrs, rworst = frozen_reproduction(fresh, frozen)
    errs.extend(rerrs[:6])
    lines.append("   FEM reproduction: %d values (incl. ERECTION_*), "
                 "worst rel=%.3e %s" %
                 (nrep, rworst, "PASS" if not rerrs else "FAIL"))

    bline, berrs = balance_check(fresh, doc)
    lines.append(bline)
    errs.extend(berrs)

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
        tol = chk.get("atol", 1e-6) + chk.get("rtol", 1e-3) * chk["scale"]
        d = abs(got - ref)
        worst_ratio = max(worst_ratio, d / max(tol, 1e-30))
        if d <= tol:
            npass += 1
        else:
            nfail += 1
            if nfail <= 8:
                lines.append("[FAIL] %-30s got=%-15.8g ref=%-15.8g "
                             "d=%g tol=%g" %
                             (chk["name"], got, ref, d, tol))
    lines.append("   ANSYS golden: %d/%d checks passed, worst d/tol=%.3f" %
                 (npass, npass + nfail, worst_ratio))
    lines.append("   group scales: %s" %
                 ", ".join("%s=%.6g" % (k, v)
                           for k, v in sorted(doc["group_scale"].items())))
    return lines, errs, nfail


def check_log(log: str, expect: dict, tag: str) -> tuple[str, list[str]]:
    errs: list[str] = []
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nerr = int(m.group(1)) if m else 0
    m = re.search(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nwarn = int(m.group(1)) if m else 0
    if nerr != expect.get("errors", 0):
        errs.append("%s ANSYS error count %d, expected %d" %
                    (tag, nerr, expect.get("errors", 0)))
    if nwarn != expect.get("warnings", 0):
        errs.append("%s ANSYS warning count %d, expected %d" %
                    (tag, nwarn, expect.get("warnings", 0)))
    sig = expect.get("warning_signature")
    if sig and sig not in log:
        errs.append("%s expected warning signature missing" % tag)
    return "   %s log: %d error / %d warning" % (tag, nerr, nwarn), errs


def live_ansys(doc: dict, cases: Path, work: Path, ansys: Path):
    lines = ["== %s LIVE ANSYS cross-check ==" % doc["id"]]
    errs: list[str] = []
    jobs = [(doc["apdl"], doc["apdl_stem"],
             [v for v in doc["csv"].values()], "full")]
    for st in doc["stages"]:
        sf = doc["stage_files"][str(st)]
        jobs.append((sf["apdl"], Path(sf["apdl"]).stem,
                     [sf["disp"], sf["link"]], "stage%d" % st))
    total = 0
    worst = 0.0
    for inp_name, job, csv_names, tag in jobs:
        inp = cases / inp_name
        shutil.copy2(inp, work / inp_name)
        for cn in csv_names:
            p = work / cn
            if p.exists():
                p.unlink()
        p = subprocess.run([str(ansys), "-b", "-j", job, "-i", inp_name,
                            "-o", job + "_live.out"], cwd=str(work),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=900)
        if p.returncode != 0:
            errs.append("LIVE %s ANSYS rc=%d" % (tag, p.returncode))
            continue
        log = (work / (job + "_live.out")).read_text(
            encoding="utf-8", errors="replace")
        expect_key = "full" if tag == "full" else "stage"
        s, e = check_log(log, doc["ansys_log_expect"][expect_key], tag)
        lines.append(s)
        errs.extend(e)
        for cn in csv_names:
            fresh = work / cn
            if not fresh.exists():
                errs.append("LIVE %s: %s not produced" % (tag, cn))
                continue
            a, b = read_csv_table(fresh), read_csv_table(cases / cn)
            if set(a) != set(b):
                errs.append("LIVE %s %s row keys differ" % (tag, cn))
                continue
            for k in sorted(a):
                if len(a[k]) != len(b[k]):
                    errs.append("LIVE %s %s row %d columns differ" % (tag, cn, k))
                    continue
                for x, y in zip(a[k], b[k]):
                    total += 1
                    d = abs(x - y)
                    worst = max(worst, d)
                    if d > 1e-6 + 1e-3 * abs(y):
                        errs.append("LIVE %s %s row %d: %g vs %g" %
                                    (tag, cn, k, x, y))
    lines.append("   %d CSV values across %d jobs, max abs diff %.3e  %s" %
                 (total, len(jobs), worst, "PASS" if not errs else "FAIL"))
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
    print("== balanced suspension + erection stages vs ANSYS NLGEOM benchmark ==")
    all_errs: list[str] = []
    total_fail = 0
    with tempfile.TemporaryDirectory(prefix="femerect_") as td:
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
