#!/usr/bin/env python3
"""第 20 项 CTest：B3 `[STABILITY]` 主链路 opt-in 对标（离线）。

与 CB3D17 的关系：那里用 ANSYS `ARCLEN` 验证了**独立探针**
`snap_probe.exe`；这里验证把同一内核以 `STABILITY=1` 接入
`cable_solver.exe` 主程序后，路径逐值同源、默认全关零污染、
解析锚点成立、ERROR 双函数不静默。本项**不再跑 ANSYS**。

检查项：
  1. 黄金值自洽（reference JSON ↔ 冻结 cases/*_stability.txt）；
  2. 默认零污染：新 exe 跑既有 CB3D10 算例，输出除 elapsed_ms 外
     与冻结 _cable.txt 逐字节一致，且不含 [STABILITY]、仍含
     [BUCKLING] ERROR=nonlinear_state_not_supported；
  3. FEM CLI 重跑 5 个 .dat 到临时目录，与冻结输出逐字节一致；
  4. 独立闭式解析锚点 P(v) 与 λ·PREF 逐点对比、极值点二次插值、
     下降段 / P(H)=0 / P(2H)=0 / 逐点残差判据；
  5. 3 类 ERROR 负例码正确、不出平衡点、[BUCKLING] 段仍在；
  6. 与 CB3D17 冻结 FEM 路径（log_fine / eng_fine）逐值对拍。

用法：
    python compare_stable_path_3d.py --bin <build_dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_reference_stable import (  # noqa: E402
    A, B, CASE_SPEC, E, H, L0, PREF, analytic_critical, build,
    p_of_v, parse_stability,
)

REF_FILE = "cb3d18_stable_path.json"
EXE = "cable_solver.exe"
TOL_DP = 1.0e-4             # N，|λ·PREF − P(v)| 上限
TOL_LAMBDA_DISCRETE = 1.0e-3   # 离散 λ_max 与解析 λ_cr
TOL_LAMBDA_FIT = 5.0e-4       # 三点二次插值后的极值误差
TOL_RES = 2.002e-3          # N = STABTOL(1e-9) × PREF(2e6)，留 0.1% 余量
TOL_V_FINAL = 5.0e-3        # |v_final − 1.0|（VMAX=1.0）
TOL_CURVE = 0.0             # 与 CB3D17 冻结 FEM 路径逐位一致

FAILURES: list[str] = []


def ok(cond: bool, msg: str) -> None:
    if not cond:
        FAILURES.append(msg)


def run_exe(exe: Path, dat: Path, cwd: Path) -> tuple[int, str]:
    p = subprocess.run([str(exe), dat.name], cwd=str(cwd),
                       capture_output=True, text=True, errors="replace")
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def quad_extremum(curve: list[tuple[float, float]], i: int) -> tuple[float, float]:
    if i <= 0 or i + 1 >= len(curve):
        return curve[i]
    x0, y0 = curve[i - 1]
    x1, y1 = curve[i]
    x2, y2 = curve[i + 1]
    s1 = (y1 - y0) / (x1 - x0)
    s2 = (y2 - y1) / (x2 - x1)
    d2 = (s2 - s1) / (x2 - x0)
    if abs(d2) < 1.0e-14:
        return x1, y1
    xs = 0.5 * (x0 + x1) - s1 / (2.0 * d2)
    ys = y0 + s1 * (xs - x0) + d2 * (xs - x0) * (xs - x1)
    return xs, ys


def normalized_cable_txt(txt: str) -> list[str]:
    """去掉与物理无关的运行耗时行，其余逐行保留用于零污染比对。"""
    return [re.sub(r"elapsed_ms=[0-9.eE+-]+", "elapsed_ms=*", ln)
            for ln in txt.splitlines()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="含 cable_solver.exe 的构建目录")
    ap.add_argument("--cases", required=True, help="冻结 cases/ 目录")
    ap.add_argument("--reference", required=True, help="冻结 reference/ 目录")
    a = ap.parse_args()

    exe = Path(a.bin) / EXE
    cases = Path(a.cases)
    refdir = Path(a.reference)
    if not exe.is_file():
        print(f"[FATAL] 找不到 {exe}")
        return 1

    # ---- 1. 黄金值自洽 -------------------------------------------------
    doc = json.loads((refdir / REF_FILE).read_text(encoding="utf-8"))
    fresh = build()
    ok(doc["checks"] == fresh["checks"] and doc["strings"] == fresh["strings"],
       "黄金值与冻结输出不一致")
    ok(all(abs(doc["criteria"][k] - v) <= 1e-12 * max(1.0, abs(v))
           for k, v in fresh["criteria"].items()),
       "解析 criteria 与闭式解不一致")
    n_checks = len(doc["checks"])
    print(f"[PASS] {n_checks} checks + {len(doc['strings'])} strings + "
          f"{len(doc['criteria'])} criteria consistent with frozen output")

    with tempfile.TemporaryDirectory(prefix="stable_path_") as td:
        tmp = Path(td)

        # ---- 2. 默认零污染（行为级，用既有 CB3D10 算例）---------------
        bench = cases.parent.parent
        base_dat = bench / "nonlinear_ansys_3d/cases/cb3d10_completion.dat"
        base_txt = bench / "nonlinear_ansys_3d/cases/cb3d10_completion_cable.txt"
        if base_dat.is_file() and base_txt.is_file():
            off = tmp / "cb3d10_off.dat"
            shutil.copy2(base_dat, off)
            rc, _ = run_exe(exe, off, tmp)
            got = tmp / "cb3d10_off_cable.txt"
            if not got.is_file():
                FAILURES.append("默认路径未产出 _cable.txt")
            else:
                txt = got.read_text(encoding="utf-8", errors="replace")
                ok(rc == 0, "既有算例默认路径退出码非 0")
                ok("[STABILITY]" not in txt, "默认路径泄漏了 [STABILITY] 段")
                ok("ERROR=nonlinear_state_not_supported" in txt,
                   "默认路径丢失 [BUCKLING] ERROR 段")
                frozen = base_txt.read_text(encoding="utf-8", errors="replace")
                ok(normalized_cable_txt(txt) == normalized_cable_txt(frozen),
                   "默认路径输出与冻结 _cable.txt 不一致（除 elapsed_ms）")
                print("[PASS] default STABILITY=0: no [STABILITY] section, "
                      "[BUCKLING] ERROR kept, cable.txt byte-identical "
                      "(excluding elapsed_ms)")
        else:
            FAILURES.append("缺少既有 CB3D10 零污染对照算例")

        # ---- 3. FEM CLI 重跑 5 个 .dat --------------------------------
        worst_repro = 0
        for stem, _, _ in CASE_SPEC:
            shutil.copy2(cases / f"{stem}.dat", tmp / f"{stem}.dat")
            rc, _ = run_exe(exe, tmp / f"{stem}.dat", tmp)
            got = tmp / f"{stem}_stability.txt"
            frozen = (cases / f"{stem}_stability.txt").read_bytes()
            if not got.is_file():
                FAILURES.append(f"{stem} 未产出 _stability.txt")
                continue
            cur = got.read_bytes()
            if cur != frozen:
                FAILURES.append(f"{stem} 重跑输出与冻结归档不一致")
                worst_repro = max(worst_repro, 1)
            expected_rc = 0 if stem in ("von_mises_log", "von_mises_eng") else 3
            ok(rc == expected_rc, f"{stem} exit={rc} != {expected_rc}")
        print(f"[PASS] FEM CLI reproduction: 5 cases, "
              f"byte-identical frozen output (worst diff flag={worst_repro})")

        # ---- 4. 解析锚点 / 极值点 / 下降段 / 残差 ----------------------
        summary: list[str] = []
        for stem, tag, strain in (("von_mises_log", "log", 1),
                                  ("von_mises_eng", "eng", 0)):
            r = parse_stability(tmp / f"{stem}_stability.txt")
            pts = r["points"]
            ok(len(pts) >= 200, f"{stem} 平衡点过少：{len(pts)}")
            curve = [(p["v"], p["lambda"]) for p in pts]
            worst_dp = max(abs(p["lambda"] * PREF - p_of_v(p["v"], strain))
                           for p in pts if 0.0 <= p["v"] <= 2.0 * H)
            ok(worst_dp <= TOL_DP,
               f"{stem} |λ·PREF−P(v)|={worst_dp:.3e} > {TOL_DP}")
            lams = [p["lambda"] for p in pts]
            lam_max, lam_min = max(lams), min(lams)
            i_max = lams.index(lam_max)
            _, lam_fit = quad_extremum(curve, i_max)
            v_cr, p_cr = analytic_critical(strain)
            lam_cr = p_cr / PREF
            ok(abs(lam_max - lam_cr) <= TOL_LAMBDA_DISCRETE,
               f"{stem} 离散 λ_max 误差 {abs(lam_max-lam_cr):.3e}")
            ok(abs(lam_fit - lam_cr) <= TOL_LAMBDA_FIT,
               f"{stem} 插值 λ_max 误差 {abs(lam_fit-lam_cr):.3e}")
            ok(lam_min < -0.5, f"{stem} 无下降段：λ_min={lam_min}")
            ok(abs(pts[-1]["v"] - 1.0) <= TOL_V_FINAL,
               f"{stem} v_final={pts[-1]['v']} 偏离 VMAX=1.0")
            ok(max(p["res"] for p in pts) <= TOL_RES,
               f"{stem} 逐点残差超 STABTOL·|Fref|={TOL_RES}")
            ok(abs(p_of_v(H, strain)) <= 1e-9
               and abs(p_of_v(2.0 * H, strain)) <= 1e-9,
               f"{stem} 解析 P(H)/P(2H) 不为 0")
            summary.append(
                f"   {stem}: {len(pts)} points, worst|dP|={worst_dp:.3e} N, "
                f"lambda_cr={lam_cr:.9f}, lambda_max={lam_max:.9f}, "
                f"lambda_max_fit={lam_fit:.9f} (err {abs(lam_fit-lam_cr):.2e}), "
                f"lambda_min={lam_min:.6f}, v_final={pts[-1]['v']:.6f}, "
                f"maxRES={max(p['res'] for p in pts):.3e} N")
        print("[PASS] analytic oracle: P(v) / critical point / descending branch")
        for s in summary:
            print(s)

        # ---- 5. ERROR 负例 ---------------------------------------------
        expect_err = {"von_mises_bad_ctrl": "invalid_control_dof",
                      "von_mises_bad_method": "unsupported_path_method",
                      "von_mises_no_load": "zero_reference_load"}
        for stem, code in expect_err.items():
            r = parse_stability(tmp / f"{stem}_stability.txt")
            ok(r["error"] == code, f"{stem} ERROR={r['error']} != {code}")
            ok(r["points"] == [] and r["npoint"] == 0,
               f"{stem} ERROR 段却产出平衡点")
            ok(r["has_buckling_error"] == 1,
               f"{stem} 丢失 [BUCKLING] ERROR 段")
            print(f"[PASS] {stem}: ERROR={code}, no points, BUCKLING kept")

        # ---- 6. 与 CB3D17 冻结 FEM 路径逐值对拍 ------------------------
        fem_csv = bench / ("snap_through_ansys_3d/cases/"
                           "cb3d17_snap_through_arc_length_fem.csv")
        if fem_csv.is_file():
            worst = 0.0
            for stem, case in (("von_mises_log", "log_fine"),
                               ("von_mises_eng", "eng_fine")):
                ref = [l.strip().split(",") for l in
                       fem_csv.read_text(encoding="utf-8").splitlines()
                       if l.startswith(f"P,{case},")]
                cur = parse_stability(tmp / f"{stem}_stability.txt")["points"]
                ok(len(ref) == len(cur),
                   f"{stem} 与 L 阶段点数不一致 {len(ref)} vs {len(cur)}")
                for arow, brow in zip(ref, cur):
                    for x, y in ((float(arow[3]), brow["lambda"]),
                                 (float(arow[4]), brow["v"]),
                                 (float(arow[5]), brow["res"]),
                                 (float(arow[6]), float(brow["iter"])),
                                 (float(arow[7]), brow["ds"]),
                                 (float(arow[8]), brow["pivot"])):
                        worst = max(worst, abs(x - y))
            ok(worst <= TOL_CURVE,
               f"与 CB3D17 冻结路径 worst|d|={worst:.3e} > {TOL_CURVE}")
            print(f"[PASS] vs CB3D17 snap_probe frozen path: "
                  f"2 curves, worst|d|={worst:.3e}")
        else:
            FAILURES.append("缺少 CB3D17 冻结 FEM 路径 CSV")

    if FAILURES:
        print("\n[FAIL] " + f"{len(FAILURES)} 项未通过：")
        for f in FAILURES:
            print("  - " + f)
        return 1
    print("\n[PASS] stable_path_3d benchmark: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
