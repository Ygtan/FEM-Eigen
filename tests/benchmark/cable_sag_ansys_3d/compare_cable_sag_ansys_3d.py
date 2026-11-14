#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: CB3D14 cable sag / Ernst + CB3D15 W>0 empty-cable state.

Offline by default (no ANSYS licence needed):

  * re-run ``cable_sag_probe.exe`` on the archived ``*.sag`` and verify it
    reproduces the frozen ``*_fem.csv`` value by value;
  * re-derive all 573 golden values from the frozen ANSYS CSVs (guards drift);
  * compare every ANSYS H / sag / segment force with the FEM kernel closed
    form on a per-group relative criterion (rtol=2e-3);
  * verify the Ernst tangent ``dT/dL0`` from a central difference of ANSYS
    H over the L0 sweep, and the ``W=0 -> Eeq/E=1`` degeneration;
  * re-run ``cable_solver.exe`` on the CB3D15 ``W>0`` .dat, reproduce the
    frozen ``*_cable.txt`` value by value, and check the empty-cable stage
    is a real catenary (H >> 100 N, sag > 0).

``--live-ansys`` re-runs all 7 LINK180 chain jobs and diffs every CSV cell
against the frozen copies.  All runtime output is ASCII.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from build_reference_sag import (  # noqa: E402
    derive, parse_bridge_stage0, parse_sag, read_csv_table, read_fem_csv,
    z_init,
)

MSYS_BIN = r"C:\msys64\mingw64\bin"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def run_fem_probe(exe: Path, sag: Path, work: Path) -> dict:
    shutil.copy2(sag, work / sag.name)
    out = work / "probe_fresh.csv"
    p = subprocess.run([str(exe), sag.name, out.name], cwd=str(work),
                       env=env_with_msys(), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=300)
    if p.returncode != 0 or not out.exists():
        raise RuntimeError("cable_sag_probe rc=%d: %s %s" %
                           (p.returncode, p.stdout, p.stderr))
    return read_fem_csv(out)


def run_fem_bridge(exe: Path, dat: Path, work: Path) -> Path:
    shutil.copy2(dat, work / dat.name)
    p = subprocess.run([str(exe), dat.name], cwd=str(work),
                       env=env_with_msys(), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600)
    fresh = work / (dat.stem + "_cable.txt")
    if not fresh.exists():
        raise RuntimeError("cable_solver rc=%d: %s %s" %
                           (p.returncode, p.stdout, p.stderr))
    return fresh


def numeric_lines(path: Path) -> dict[int, list[float]]:
    """逐行抽浮点数；跳过耗时行（每次运行必然不同）。"""
    out: dict[int, list[float]] = {}
    for i, raw in enumerate(path.read_text(encoding="utf-8",
                                           errors="replace").splitlines(), 1):
        if "elapsed_ms" in raw or raw.lstrip().startswith("#"):
            continue
        vals = []
        for tok in raw.replace(",", " ").split():
            try:
                vals.append(float(tok))
            except ValueError:
                continue
        if vals:
            out[i] = vals
    return out


def reproduction(fresh: Path, frozen: Path) -> tuple[int, float, list[str]]:
    a, b = numeric_lines(fresh), numeric_lines(frozen)
    errs: list[str] = []
    n = 0
    worst = 0.0
    if set(a) != set(b):
        return 0, 0.0, ["line structure differs: %d vs %d" % (len(a), len(b))]
    for k in sorted(a):
        if len(a[k]) != len(b[k]):
            errs.append("line %d: %d vs %d numbers" % (k, len(a[k]), len(b[k])))
            continue
        for x, y in zip(a[k], b[k]):
            n += 1
            d = abs(x - y)
            worst = max(worst, d)
            if d > 1e-9 + 1e-12 * abs(y):
                errs.append("line %d: %.12g vs %.12g" % (k, x, y))
    return n, worst, errs[:6]


def fem_value(fem: dict, cid: str, key: str):
    row = fem[cid]
    if key == "H":
        return row["T"]
    if key == "f":
        return row["f"]
    if key.startswith("T_"):
        # FEM 用一根索的单一 T；ANSYS 分段索链的张力沿索变化（见 README）。
        return row["T"]
    return None


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


def compare_one(doc: dict, probe: Path, solver: Path, cases: Path,
                 work: Path):
    lines = ["== %s ==" % doc["id"]]
    errs: list[str] = []
    nfail = 0

    # 1) FEM 内核 CLI 逐位复现
    fresh = run_fem_probe(probe, cases / doc["sag_file"], work)
    frozen = read_fem_csv(cases / doc["fem_csv"])
    if set(fresh) != set(frozen):
        errs.append("FEM probe case ids differ")
        nfail += 1
    else:
        nrep, rworst = 0, 0.0
        for cid in frozen:
            if set(fresh[cid]) != set(frozen[cid]):
                errs.append("%s FEM columns differ" % cid)
                nfail += 1
                continue
            for k, v in frozen[cid].items():
                nrep += 1
                d = abs(fresh[cid][k] - v)
                rworst = max(rworst, d / max(abs(v), 1e-12))
                if d > 1e-6 + 1e-12 * abs(v):
                    errs.append("%s %s got=%r ref=%r" %
                                (cid, k, fresh[cid][k], v))
                    nfail += 1
        lines.append("   FEM probe reproduction: %d values, worst rel=%.3e %s" %
                     (nrep, rworst, "PASS" if rworst < 1e-9 else "FAIL"))

    # 2) 冻结 CSV -> 黄金值 逐项重算（防漂移）
    cfg, sag_cases = parse_sag(cases / doc["sag_file"])
    derived_all = {}
    for cs in sag_cases:
        derived_all[cs["id"]] = derive(cs, cfg)
    ninteg = 0
    for chk in doc["checks"]:
        cid, key = chk["case"], chk["key"]
        d = derived_all[cid]
        got = d["H"] if key == "H" else (
            d["f"] if key == "f" else d["T"][int(key[2:]) - 1])
        ninteg += 1
        if abs(got - chk["ref"]) > 1e-9 + 1e-12 * abs(chk["ref"]):
            errs.append("golden drift %s %s" % (cid, key))
            nfail += 1
    lines.append("   frozen ANSYS CSV re-derivation: %d/%d checks %s" %
                 (ninteg, len(doc["checks"]),
                  "PASS" if ninteg == len(doc["checks"]) else "FAIL"))

    # 3) ANSYS 黄金值 vs FEM 内核闭式，分组相对判据
    scales: dict[str, float] = {}
    for chk in doc["checks"]:
        got = fem_value(fresh, chk["case"], chk["key"])
        if got is None:
            errs.append("no FEM value for %s %s" % (chk["case"], chk["key"]))
            nfail += 1
            continue
        scales[chk["group"]] = max(scales.get(chk["group"], 0.0),
                                   max(abs(got), abs(chk["ref"])))
    npass = 0
    worst = 0.0
    for chk in doc["checks"]:
        got = fem_value(fresh, chk["case"], chk["key"])
        ref = chk["ref"]
        scale = scales.get(chk["group"], 1.0)
        tol = 1e-9 + chk["tol"] * scale
        r = abs(got - ref) / tol
        worst = max(worst, r)
        if r > 1.0:
            errs.append("%s %s: got=%.9g ref=%.9g d/tol=%.3f" %
                        (chk["case"], chk["key"], got, ref, r))
            nfail += 1
        else:
            npass += 1
    lines.append("   ANSYS golden: %d/%d checks passed, worst d/tol=%.3f" %
                 (npass, len(doc["checks"]), worst))
    lines.append("   group scales: " + ", ".join(
        "%s=%.6g" % (k, v) for k, v in sorted(scales.items())))

    # 4) Ernst 切线
    t = doc["tangent"]
    ana = fresh[t["center"]]["dTdL0"]
    rel = abs(t["ansys_fd"] - ana) / max(abs(ana), 1.0)
    lines.append("   Ernst tangent: ANSYS_fd=%.6e  FEM=%.6e  rel=%.3e %s" %
                 (t["ansys_fd"], ana, rel, "PASS" if rel <= t["tol"] else "FAIL"))
    if rel > t["tol"]:
        errs.append("Ernst tangent rel=%.3e > %.3e" % (rel, t["tol"]))
        nfail += 1

    # 5) W=0 直杆退化
    dg = doc["degeneration"]
    ratio = fresh[dg["case"]]["Eeq_ratio"]
    ok = abs(ratio - dg["expect"]) <= dg["tol"]
    lines.append("   W=0 degeneration: Eeq/E=%.9f %s" %
                 (ratio, "PASS" if ok else "FAIL"))
    if not ok:
        errs.append("W=0 Eeq/E != 1")
        nfail += 1

    # 6) ANSYS 日志签名
    for cm in doc["cases"]:
        op = cases / (cm["id"] + "_ansys.out")
        if not op.exists():
            errs.append("%s: frozen ANSYS .out missing" % cm["id"])
            nfail += 1
            continue
        s, e = check_log(op, doc["ansys_log_expect"], cm["id"])
        errs.extend(e)
        nfail += len(e)

    # 7) CB3D15 空缆态：重跑 FEM + 判据
    b = doc["bridge"]
    bdat = cases / b["dat"]
    bfrozen = cases / b["cable"]
    try:
        bfresh = run_fem_bridge(solver, bdat, work)
        nrep, rworst, rerrs = reproduction(bfresh, bfrozen)
        errs.extend(rerrs[:6])
        nfail += len(rerrs)
        b0 = parse_bridge_stage0(bfresh)
        crit = b["criteria"]
        bok = (b0["H"] > crit["H_min"] and b0["sag"] > crit["sag_min"])
        lines.append("   CB3D15 W>0 bridge: %d values, worst rel=%.3e; "
                     "stage0 H=%.3f kN sag=%.6f m %s" %
                     (nrep, rworst, b0["H"] / 1e3, b0["sag"],
                      "PASS" if bok else "FAIL"))
        if not bok:
            errs.append("CB3D15 empty-cable stage is degenerate")
            nfail += 1
    except Exception as exc:                                # noqa: BLE001
        errs.append("CB3D15 exception: %s" % exc)
        nfail += 1

    return lines, errs, nfail


def live_ansys(doc: dict, cases: Path, work: Path, ansys: Path):
    lines = ["== %s LIVE ANSYS cross-check ==" % doc["id"]]
    errs: list[str] = []
    total = 0
    worst = 0.0
    for cm in doc["cases"]:
        job = cm["id"]
        shutil.copy2(cases / cm["apdl"], work / cm["apdl"])
        for name in cm["csv"].values():
            p = work / name
            if p.exists():
                p.unlink()
        r = subprocess.run([str(ansys), "-b", "-j", job, "-i", cm["apdl"],
                            "-o", job + "_live.out"], cwd=str(work),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=900)
        if r.returncode != 0:
            errs.append("LIVE %s ANSYS rc=%d" % (job, r.returncode))
            continue
        log = (work / (job + "_live.out")).read_text(
            encoding="utf-8", errors="replace")
        s, e = check_log(work / (job + "_live.out"),
                         doc["ansys_log_expect"], "LIVE-" + job)
        lines.append(s)
        errs.extend(e)
        for name in cm["csv"].values():
            fresh = work / name
            if not fresh.exists():
                errs.append("LIVE %s: %s not produced" % (job, name))
                continue
            a, b = read_csv_table(fresh), read_csv_table(cases / name)
            if set(a) != set(b):
                errs.append("LIVE %s %s row keys differ" % (job, name))
                continue
            for k in sorted(a):
                if len(a[k]) != len(b[k]):
                    errs.append("LIVE %s %s row %d columns differ" % (job, name, k))
                    continue
                for x, y in zip(a[k], b[k]):
                    total += 1
                    d = abs(x - y)
                    worst = max(worst, d)
                    if d > 1e-6 + 1e-3 * abs(y):
                        errs.append("LIVE %s %s row %d: %g vs %g" %
                                    (job, name, k, x, y))
    lines.append("   %d CSV values across %d jobs, max abs diff %.3e  %s" %
                 (total, len(doc["cases"]), worst,
                  "PASS" if not errs else "FAIL"))
    return lines, errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="FEM-Eigen build dir")
    ap.add_argument("--cases", default=str(HERE / "cases"))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys", default=DEFAULT_ANSYS)
    a = ap.parse_args()

    exe = Path(a.bin)
    probe, solver = exe / "cable_sag_probe.exe", exe / "cable_solver.exe"
    cases, ref = Path(a.cases), Path(a.reference)
    print("== cable sag / Ernst + W>0 empty-cable state vs ANSYS 18.0 benchmark ==")
    for p in (probe, solver):
        if not p.exists():
            print("[FAIL] required exe not found: %s" % p)
            return 1
    docs = [json.loads(f.read_text(encoding="utf-8"))
            for f in sorted(ref.glob("*_ansys.json"))]
    if not docs:
        print("[FAIL] no reference JSON found in %s" % ref)
        return 1

    all_errs: list[str] = []
    nfail = 0
    with tempfile.TemporaryDirectory(prefix="fem_sag_") as td:
        work = Path(td)
        for doc in docs:
            try:
                lines, errs, nf = compare_one(doc, probe, solver, cases, work)
            except Exception as exc:                       # noqa: BLE001
                print("[FAIL] %s exception: %s" % (doc["id"], exc))
                return 1
            print("\n".join(lines))
            all_errs.extend(errs)
            nfail += nf
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
