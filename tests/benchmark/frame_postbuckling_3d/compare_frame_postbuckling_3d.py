#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 21 项 CTest：CB3D19 CR 梁刚接框架后屈曲（离线）。

与 CB3D17（两杆 TRUSS 桁架 + 解析 P(v) 闭式）不同，本题必须用梁单元：
顶点刚接带来弯曲刚度，路径是「第一极限点 -> 下降段 -> 反向荷载 ->
再硬化」。独立锚点换成纯 Python 2D 线性刚架刚度矩阵 k0（本 benchmark
自带第三套实现，与 BridgeStudio 侧、frame_probe.exe 内核三方互证）。

检查项：
  1. 黄金值自洽（reference JSON <-> 冻结 cases/）；
  2. FEM CLI 复现：临时目录重跑 frame_probe.exe，逐值与冻结 CSV 一致；
  3. 独立线性刚度锚点 k0 与 FEM 内核精确切线严格对拍；
  4. FEM 路径自洽：收敛 / 残差 / 下降段 / 终止位移 / 第一极限点；
  5. ANSYS 冻结曲线：反力反算 λ vs TIME、初始斜率、下降段、0 error；
  6. FEM（fresh）与 ANSYS（frozen）完整路径三点二次插值逐点对拍；
  7. 负例 zero_reference_load 不静默产出路径；
  8. 可选 --live-ansys：临时目录重跑 2 个 ANSYS job，逐 CSV 单元 diff，
     不覆写任何冻结文件。

用法：
    python compare_frame_postbuckling_3d.py --bin <build_dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_reference_frame import (  # noqa: E402
    ANSYS_JOBS, CASE_ID, FEM_CASES, FEM_CSV, NEG_CASE, PARAMS_FILE, P_KEYS,
    build, cfg_of, check, curve_of, first_limit_point, interp_lambda,
    linear_frame_stiffness, parse_ansys, parse_fem,
)

REF_FILE = "cb3d19_frame_postbuckling.json"
EXE = "frame_probe.exe"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
TOL_LIVE_CSV = 1.0e-6
TOL_LIVE_CSV_REL = 1.0e-10

FAILURES: list[str] = []


def ok(cond: bool, msg: str) -> None:
    if not cond:
        FAILURES.append(msg)


def reproduction(frozen: dict, fresh: dict) -> tuple[float, int]:
    """逐值比较两个 parse_fem() 结果，返回 (worst|d|, nvalues)。"""
    worst, n = 0.0, 0
    fa, fb = frozen[0], fresh[0]
    fs, gs = frozen[1], fresh[1]
    ok(set(fa) == set(fb), "FEM P cases differ after CLI reproduction")
    for cid in set(fa) & set(fb):
        ok(len(fa[cid]) == len(fb[cid]),
           "%s point count %d vs %d" % (cid, len(fa[cid]), len(fb[cid])))
        for ra, rb in zip(fa[cid], fb[cid]):
            for x, y in zip(ra, rb):
                worst = max(worst, abs(x - y))
                n += 1
    ok(set(fs) == set(gs), "FEM S cases differ after CLI reproduction")
    for cid in set(fs) & set(gs):
        ok(set(fs[cid]) == set(gs[cid]), "%s S keys differ" % cid)
        for k in set(fs[cid]) & set(gs[cid]):
            worst = max(worst, abs(fs[cid][k] - gs[cid][k]))
            n += 1
    ok(frozen[2] == fresh[2], "FEM E (error) rows differ after CLI reproduction")
    return worst, n


def live_ansys(cases: Path, cfg: dict, ansys_exe: str, tmp: Path) -> None:
    """临时目录重跑 ANSYS，逐 CSV 单元对比；绝不写回 cases/。"""
    exe = Path(ansys_exe)
    if not exe.is_file():
        FAILURES.append("live-ansys: ANSYS180.exe not found: %s" % exe)
        return
    for job in ANSYS_JOBS:
        shutil.copy2(cases / (job + ".inp"), tmp / (job + ".inp"))
        p = subprocess.run(
            [str(exe), "-b", "-j", job, "-i", job + ".inp",
             "-o", job + "_live.out"],
            cwd=str(tmp), timeout=1800, capture_output=True,
            text=True, errors="replace")
        ok(p.returncode == 0, "live-ansys %s exit=%d" % (job, p.returncode))
        got_path = tmp / (job + "_curve.csv")
        if not got_path.is_file():
            FAILURES.append("live-ansys %s produced no curve CSV" % job)
            continue
        got = parse_ansys(got_path, cfg)
        ref = parse_ansys(cases / (job + "_curve.csv"), cfg)
        ok(len(got) == len(ref), "live-ansys %s rows %d vs %d" % (
            job, len(got), len(ref)))
        worst = 0.0
        nval = 0
        for a, b in zip(got, ref):
            for x, y in zip(a, b):
                ok(abs(x - y) <= TOL_LIVE_CSV + TOL_LIVE_CSV_REL * abs(y),
                   "live-ansys %s cell drift %.3e" % (job, x - y))
                worst = max(worst, abs(x - y))
                nval += 1
        print("[PASS] live-ansys %s: %d values, max|d|=%.3e" % (
            job, nval, worst))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="dir containing frame_probe.exe")
    ap.add_argument("--cases", required=True, help="frozen cases directory")
    ap.add_argument("--reference", required=True, help="reference JSON directory")
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys-exe", default=DEFAULT_ANSYS)
    a = ap.parse_args()

    bin_dir, cases = Path(a.bin), Path(a.cases)
    ref_file = Path(a.reference) / REF_FILE
    if not ref_file.is_file():
        print("[FAIL] missing reference: %s" % ref_file)
        return 1
    doc = json.loads(ref_file.read_text(encoding="utf-8"))
    cfg = doc["config"]
    tol = doc["tolerances"]

    # ---- 1. 黄金值自洽 ----
    good, msg = check(doc, cases_dir=cases)
    ok(good, "golden self-consistency: " + msg)
    if good:
        print("[PASS] golden: %s" % msg)

    # ---- 2. FEM CLI 复现 ----
    exe = bin_dir / EXE
    fresh = None
    if not exe.is_file():
        FAILURES.append("missing FEM CLI: %s" % exe)
    else:
        with tempfile.TemporaryDirectory(prefix="cb3d19_fem_") as td:
            tmp = Path(td)
            shutil.copy2(cases / PARAMS_FILE, tmp / PARAMS_FILE)
            p = subprocess.run([str(exe), PARAMS_FILE, FEM_CSV],
                               cwd=str(tmp), timeout=900,
                               capture_output=True, text=True, errors="replace")
            ok(p.returncode == 0, "frame_probe exit=%d" % p.returncode)
            out_csv = tmp / FEM_CSV
            if not out_csv.is_file():
                FAILURES.append("frame_probe produced no CSV")
            else:
                fresh = parse_fem(out_csv)
                frozen = parse_fem(cases / FEM_CSV)
                worst, nval = reproduction(frozen, fresh)
                ok(worst <= tol["reproduce"],
                   "FEM reproduction worst|d|=%.3e > %.1e" % (
                       worst, tol["reproduce"]))
                print("[PASS] FEM CLI reproduction: %d values, worst|d|=%.3e" % (
                    nval, worst))
    if fresh is None:
        fresh = parse_fem(cases / FEM_CSV)
    fpts, fsum, ferr = fresh

    # ---- 3. 独立线性刚度锚点（第三套实现）----
    k0 = linear_frame_stiffness(cfg)
    ok(abs(k0 - doc["k0_analytic"]) <= 1.0e-6,
       "independent k0 %.9e vs golden %.9e" % (k0, doc["k0_analytic"]))
    for cid in FEM_CASES:
        kf = fsum.get(cid, {}).get("k0_fem", 0.0)
        rel = abs(kf - k0) / k0
        ok(rel <= tol["k0_rel"],
           "%s k0_fem rel=%.3e > %.1e" % (cid, rel, tol["k0_rel"]))
    print("[PASS] independent linear-frame anchor: k0=%.9e, "
          "k0/PREF=%.9f, worst rel=%.3e" % (
              k0, k0 / cfg["PREF"],
              max(abs(fsum.get(c, {}).get("k0_fem", 0.0) - k0) / k0
                  for c in FEM_CASES)))

    # ---- 4. FEM 路径自洽 ----
    for cid in FEM_CASES:
        s = fsum.get(cid, {})
        pts_rows = fpts.get(cid, [])
        ok(len(pts_rows) >= 400, "%s too few points: %d" % (cid, len(pts_rows)))
        ok(s.get("converged", 0) > 0.5 and s.get("error_empty", 0) > 0.5,
           "%s not converged" % cid)
        ok(s.get("max_residual_rel", 9.9) <= tol["res_rel"],
           "%s max_res_rel=%.3e" % (cid, s.get("max_residual_rel", 9.9)))
        ok(s.get("lambda_min", 0.0) < 0.0,
           "%s no descending branch" % cid)
        ok(s.get("v_final", 0.0) >= 0.99 * cfg["VMAX"],
           "%s v_final=%.6f" % (cid, s.get("v_final", 0.0)))
        ok(max(r[3] for r in pts_rows) <= tol["res_rel"] * cfg["PREF"],
           "%s pointwise residual too large" % cid)
        print("[PASS] FEM %-8s %4d points: lambda_cr_fit=%.9f, "
              "lambda_min=%.6f, v_final=%.6f, maxRES=%.3e" % (
                  cid, len(pts_rows), s.get("lambda_max_fit", 0.0),
                  s.get("lambda_min", 0.0), s.get("v_final", 0.0),
                  max(r[3] for r in pts_rows)))

    # ---- 7. 负例 ----
    ok(ferr.get(NEG_CASE) == "zero_reference_load",
       "negative case error=%r" % ferr.get(NEG_CASE))
    ok(not fpts.get(NEG_CASE), "negative case produced equilibrium points")
    ok(doc["strings"].get("bad_noload/error") == "zero_reference_load",
       "golden negative error string drifted")
    ok(doc["strings"].get("report/result") == "ALL CHECKS PASSED",
       "frozen compare report is not ALL CHECKS PASSED")
    print("[PASS] negative case: ERROR=zero_reference_load, no points")

    # ---- 5/6. ANSYS 冻结曲线 + FEM(fresh) vs ANSYS 完整路径 ----
    for job in ANSYS_JOBS:
        rows = parse_ansys(cases / (job + "_curve.csv"), cfg)
        ok(len(rows) >= 40, "%s too few ANSYS substeps: %d" % (job, len(rows)))
        if not rows:
            continue
        cur = [(v, lr) for v, lr, lt, iset in rows]
        worst_lt = max(abs(lr - lt) for v, lr, lt, iset in rows)
        ok(worst_lt <= tol["lambda_time"],
           "%s worst|lambda_r-lambda_t|=%.3e" % (job, worst_lt))
        small = [x for x in rows if 0.0 < x[0] <= 0.002 * cfg["H"]]
        ok(bool(small), "%s no small-displacement substeps" % job)
        if small:
            slope = small[-1][1] / small[-1][0]
            srel = abs(slope - k0 / cfg["PREF"]) / (k0 / cfg["PREF"])
            ok(srel <= tol["slope_rel"],
               "%s initial slope rel=%.3e" % (job, srel))
        vf, lf, _, _ = first_limit_point(cur)
        lmin = min(y for x, y in cur)
        ok(lmin < 0.0, "%s no descending branch" % job)
        ok(doc["strings"].get(job + "/log", "").startswith("errors=0"),
           "%s ANSYS log: %s" % (job, doc["strings"].get(job + "/log")))

        fem_cid = "fine" if job == "fr_fine" else "coarse"
        fs = fsum.get(fem_cid, {})
        ok(abs(lf - fs.get("lambda_max_fit", -9.9)) <= tol["extremum"]
           and abs(vf - fs.get("v_at_lambda_max_fit", -9.9)) <= tol["extremum"],
           "%s first limit point ANSYS=(%.7f,%.7f) FEM=(%.7f,%.7f)" % (
               job, vf, lf, fs.get("v_at_lambda_max_fit", -9.9),
               fs.get("lambda_max_fit", -9.9)))
        fc = curve_of(fpts, fem_cid)
        worst = 0.0
        wv = 0.0
        for v, lr, lt, iset in rows:
            if 0.0 <= v <= cfg["VMAX"]:
                d = abs(lr - interp_lambda(fc, v))
                if d > worst:
                    worst, wv = d, v
        ok(worst <= tol["curve"],
           "%s FEM-vs-ANSYS worst|dlambda|=%.3e @v=%.5f" % (job, worst, wv))
        print("[PASS] ANSYS %-9s %3d substeps: lambda_cr=%.9f, "
              "lambda_min=%.6f, worst|dlambda|=%.3e @v=%.5f, "
              "worst|lam_r-time|=%.3e" % (
                  job, len(rows), lf, lmin, worst, wv, worst_lt))

    # ---- 8. criteria 与冻结值一致（check 已比对，这里再显式核对关键项）----
    for key in ("k0_rel_fem", "worst_curve_fr_fine", "worst_curve_fr_coarse",
                "lambda_cr_ansys_fr_fine", "lambda_min_fine"):
        ok(key in doc["criteria"], "missing criterion %s" % key)
    ok(doc["criteria"]["worst_curve_fr_fine"] <= tol["curve"]
       and doc["criteria"]["worst_curve_fr_coarse"] <= tol["curve"],
       "golden curve criteria exceed tolerance")
    print("[PASS] criteria: %d entries, %d checks, %d strings" % (
        len(doc["criteria"]), len(doc["checks"]), len(doc["strings"])))

    if a.live_ansys:
        with tempfile.TemporaryDirectory(prefix="cb3d19_live_") as td:
            live_ansys(cases, cfg, a.ansys_exe, Path(td))

    if FAILURES:
        print("\n[FAIL] %d check(s) failed:" % len(FAILURES))
        for f in FAILURES:
            print("  - " + f)
        return 1
    print("\n[PASS] frame_postbuckling_3d benchmark: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
