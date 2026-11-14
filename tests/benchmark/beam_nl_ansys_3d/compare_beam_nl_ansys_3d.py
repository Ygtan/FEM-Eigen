#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: CB3D16 large-rotation beam-column (CR kernel) vs ANSYS 18.0.

Offline by default (no ANSYS licence needed):

  * re-run ``beam_nl_probe.exe`` on the archived ``*.beam`` and verify it
    reproduces the frozen ``*_fem.csv`` value by value;
  * re-derive all golden values from the frozen ANSYS CSVs (guards drift);
  * compare every ANSYS nodal displacement/rotation with the freshly
    computed FEM co-rotational solution using
    ``|got-ref| <= atol + rtol * per-DOF scale``;
  * check analytic criteria (full-circle closure, tip rotation,
    P-Delta amplification, beam-column exact solution);
  * check archived ANSYS logs are 0 error / 0 warning.

``--live-ansys`` re-runs all 8 BEAM188 jobs in a temporary directory and
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

from build_reference_beamnl import (  # noqa: E402
    derive, exact_values, parse_beam, read_disp_csv, read_fem_csv,
)

MSYS_BIN = r"C:\msys64\mingw64\bin"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def run_fem_probe(exe: Path, beam: Path, work: Path):
    shutil.copy2(beam, work / beam.name)
    out = work / "beamnl_fresh.csv"
    p = subprocess.run([str(exe), beam.name, out.name], cwd=str(work),
                       env=env_with_msys(), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600)
    if p.returncode != 0 or not out.exists():
        raise RuntimeError("beam_nl_probe rc=%d: %s %s" %
                           (p.returncode, p.stdout, p.stderr))
    return read_fem_csv(out)


def reproduction(fresh, frozen) -> tuple[int, float, list[str]]:
    """逐值比对 FEM CLI 重跑结果与冻结 CSV。"""
    n = 0
    worst = 0.0
    errs: list[str] = []
    if set(fresh[0]) != set(frozen[0]) or set(fresh[1]) != set(frozen[1]):
        return 0, 0.0, ["FEM case ids differ"]
    for cid in frozen[0]:
        a, b = fresh[0][cid], frozen[0][cid]
        if set(a) != set(b):
            errs.append("{}: node ids differ".format(cid))
            continue
        for i in b:
            for x, y in zip(a[i], b[i]):
                n += 1
                d = abs(x - y)
                worst = max(worst, d)
                if d > 1e-12 + 1e-12 * abs(y):
                    errs.append("{} node {}: {:.12g} vs {:.12g}".format(cid, i, x, y))
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
                errs.append("{} S {}: {:.12g} vs {:.12g}".format(cid, k, a[k], b[k]))
    return n, worst, errs[:6]


def check_log(path: Path, expect: dict, tag: str) -> tuple[str, list[str]]:
    txt = path.read_text(encoding="utf-8", errors="replace")
    nerr = len(re.findall(r"\*\*\* ERROR \*\*\*", txt))
    nwarn = len(re.findall(r"\*\*\* WARNING \*\*\*", txt))
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
    if m:
        nerr = int(m.group(1))
    m = re.search(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
    if m:
        nwarn = int(m.group(1))
    errs = []
    if nerr != expect["errors"]:
        errs.append("%s errors=%d, expected %d" % (tag, nerr, expect["errors"]))
    if nwarn != expect["warnings"]:
        errs.append("%s warnings=%d, expected %d" %
                    (tag, nwarn, expect["warnings"]))
    return "   %s log: %d error / %d warning" % (tag, nerr, nwarn), errs


def live_ansys(doc: dict, cases: Path, work: Path,
               ansys: Path) -> tuple[list[str], list[str]]:
    lines = ["-- live ANSYS re-run (temporary dir, frozen CSV untouched) --"]
    errs: list[str] = []
    total = 0
    worst = 0.0
    for cm in doc["cases"]:
        cid, job = cm["id"], cm["id"]
        src = cases / cm["apdl"]
        shutil.copy2(src, work / src.name)
        out = work / (job + ".out")
        p = subprocess.run([str(ansys), "-b", "-j", job, "-i", src.name,
                            "-o", out.name], cwd=str(work),
                           env=env_with_msys(), capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=1200)
        fresh = work / cm["csv"]["disp"]
        frozen = cases / cm["csv"]["disp"]
        if p.returncode != 0 or not fresh.exists():
            errs.append("{}: rc={} csv missing".format(cid, p.returncode))
            continue
        a, b = read_disp_csv(fresh), read_disp_csv(frozen)
        if set(a) != set(b):
            errs.append("{}: node ids differ".format(cid))
            continue
        for i in b:
            for x, y in zip(a[i], b[i]):
                total += 1
                d = abs(x - y)
                worst = max(worst, d)
                # ANSYS NLGEOM 用 Newton 力/位移收敛判据，重跑同一输入会有
                # 1e-8 量级的末位波动（实测 max 9.6e-9）；这里按
                # NLGEOM 数值复现口径给容差，而不是逐位一致口径。
                if d > 1e-8 + 1e-5 * abs(y):
                    errs.append("{} node {}: {:.10g} vs {:.10g}".format(cid, i, x, y))
    lines.append("   %d CSV values across %d jobs, max abs diff %.3e  %s" %
                 (total, len(doc["cases"]), worst,
                  "PASS" if not errs else "FAIL"))
    return lines, errs[:8]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="FEM-Eigen build dir")
    ap.add_argument("--cases", default=str(HERE / "cases"))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys", default=DEFAULT_ANSYS)
    a = ap.parse_args()

    exe = Path(a.bin) / "beam_nl_probe.exe"
    cases, ref = Path(a.cases), Path(a.reference)
    print("== CB3D16 large-rotation beam-column vs ANSYS 18.0 BEAM188 benchmark ==")
    if not exe.exists():
        print("[FAIL] required exe not found: %s" % exe)
        return 1
    docs = [json.loads(f.read_text(encoding="utf-8"))
            for f in sorted(ref.glob("*_ansys.json"))]
    if not docs:
        print("[FAIL] no reference JSON found in %s" % ref)
        return 1

    all_errs: list[str] = []
    nfail = 0
    with tempfile.TemporaryDirectory(prefix="fem_beamnl_") as td:
        work = Path(td)
        for doc in docs:
            try:
                fresh = run_fem_probe(exe, cases / doc["beam_file"], work)
                frozen = read_fem_csv(cases / doc["fem_csv"])
            except Exception as exc:                        # noqa: BLE001
                print("[FAIL] FEM probe: %s" % exc)
                return 1
            nrep, rworst, rerrs = reproduction(fresh, frozen)
            print("== %s ==" % doc["id"])
            print("   FEM reproduction: %d values, worst abs diff %.3e  %s" %
                  (nrep, rworst, "PASS" if not rerrs else "FAIL"))
            all_errs.extend(rerrs)

            # 防漂移：从冻结 CSV 重算全部黄金值
            groups, _ = parse_beam(cases / doc["beam_file"])
            rebuilt: dict[tuple, float] = {}
            scales: dict[tuple, float] = {}
            for cm in doc["cases"]:
                c = dict(id=cm["id"], group=cm["group"], kind=cm["kind"],
                         lam=cm["lam"])
                d = derive(c, groups[cm["group"]])
                for ck in d["checks"]:
                    rebuilt[(ck["case"], ck["key"])] = ck["ref"]
                    scales[(ck["case"], ck["key"])] = ck["scale"]

            nodal = fresh[0]
            ncheck = 0
            worst_ratio = 0.0
            worst_key = ""
            for ck in doc["checks"]:
                key = (ck["case"], ck["key"])
                if key not in rebuilt or \
                   abs(rebuilt[key] - ck["ref"]) > 1e-12 + 1e-12 * abs(ck["ref"]):
                    all_errs.append("golden value drift: %s/%s" % key)
                    continue
                node = int(ck["key"].split("_")[1])
                dof = int(ck["key"].split("_")[2]) - 1
                try:
                    got = nodal[ck["case"]][node][dof]
                except (KeyError, IndexError):
                    all_errs.append("missing FEM value: %s/%s" % key)
                    continue
                ncheck += 1
                tol = ck["atol"] + ck["tol"] * ck["scale"]
                ratio = abs(got - ck["ref"]) / tol
                if ratio > worst_ratio:
                    worst_ratio, worst_key = ratio, "%s/%s" % key
                if ratio > 1.0:
                    all_errs.append("%s/%s: FEM %.9g vs ANSYS %.9g, d/tol=%.3f" %
                                    (ck["case"], ck["key"], got, ck["ref"], ratio))
            print("   ANSYS golden vs FEM: %d/%d checks, worst d/tol=%.4f (%s)  %s" %
                  (ncheck, len(doc["checks"]), worst_ratio, worst_key,
                   "PASS" if ncheck == len(doc["checks"]) else "FAIL"))
            if ncheck != len(doc["checks"]):
                nfail += 1

            # 解析标尺判据（从冻结 CSV 重算，防止 criteria 值被改写）
            derived = {}
            for cm in doc["cases"]:
                c = dict(id=cm["id"], group=cm["group"], kind=cm["kind"],
                         lam=cm["lam"])
                derived[cm["id"]] = derive(c, groups[cm["group"]])["derived"]
            cir = [c for c in doc["cases"] if c["kind"] == "CIRCLE"]
            pd = [c for c in doc["cases"] if c["kind"] == "PDELTA"]
            cur = {
                "full_circle_closure": next(
                    derived[c["id"]]["closure"] for c in cir if c["lam"] == 1.0),
                "max_tip_rotation_rel_error": max(
                    derived[c["id"]]["rot_rel"] for c in cir),
                "pdelta_amplification": [derived[c["id"]]["amp"] for c in pd],
                "pdelta_y_vs_exact_max_rel": max(
                    derived[c["id"]]["y_rel"] for c in pd),
            }
            crit = doc["criteria"]
            crit_errs: list[str] = []
            for k in ("full_circle_closure", "max_tip_rotation_rel_error",
                      "pdelta_y_vs_exact_max_rel"):
                if cur[k] > crit[k]["limit"]:
                    crit_errs.append("%s %.4e > %.4e" % (k, cur[k], crit[k]["limit"]))
                if abs(cur[k] - crit[k]["value"]) > 1e-12 + 1e-12 * abs(crit[k]["value"]):
                    crit_errs.append("%s criterion value drift" % k)
            amps = cur["pdelta_amplification"]
            if crit["pdelta_amplification"]["monotonic"] and \
               not all(amps[i + 1] > amps[i] for i in range(len(amps) - 1)):
                crit_errs.append("P-Delta amplification not monotonic")
            if amps[-1] < crit["pdelta_amplification"]["min_final"]:
                crit_errs.append("P-Delta final amplification %.3f < %.3f" %
                                 (amps[-1], crit["pdelta_amplification"]["min_final"]))
            print("   analytic criteria: %d checked  %s" %
                  (len(crit), "PASS" if not crit_errs else "FAIL"))
            print("      full-circle closure %.3e m | max rot error %.3e | "
                  "P-Delta amp 1 -> %.3f | y vs exact %.3e" %
                  (cur["full_circle_closure"], cur["max_tip_rotation_rel_error"],
                   amps[-1], cur["pdelta_y_vs_exact_max_rel"]))
            all_errs.extend(crit_errs)

            # ANSYS 日志
            for cm in doc["cases"]:
                p = cases / cm["out"]
                if not p.exists():
                    all_errs.append("missing ANSYS log %s" % p.name)
                    continue
                line, e = check_log(p, doc["ansys_log_expect"], cm["id"])
                all_errs.extend(e)
            print("   %d archived ANSYS logs checked for 0 error / 0 warning" %
                  len(doc["cases"]))

            if a.live_ansys:
                l2, e2 = live_ansys(doc, cases, work, Path(a.ansys))
                print("\n".join(l2))
                all_errs.extend(e2)

    for e in all_errs:
        print("[FAIL] " + e)
    print("Result: %d case(s), %d failure(s)" % (len(docs), nfail + len(all_errs)))
    return 0 if not all_errs and nfail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
