#!/usr/bin/env python3
"""第 22 项 CTest：成桥态起点全过程稳定（STABASIS=COMPLETION）对标（离线）。

与第 20 项 CTest 的关系：那里验证 `STABILITY=1` 在**参考构形**上接入
`cable_solver.exe`；这里验证第十九次 ① 把路径起点换成**成桥态收敛位移
u0 + 恒定荷载 F0**（`R = F0 + λ·Finc − Fint(u)`），且两种 BASIS 分列、
默认零污染、独立闭式锚点成立。本项**不跑 ANSYS**（真实桥梁两段对标
属阶段 ②）。

检查项：
  1. 黄金值自洽（reference JSON ↔ 冻结 cases/*）；
  2. FEM CLI 重跑 6 个 .dat 到临时目录，输出逐字节一致、退出码符合预期；
  3. 独立闭式 `λ(v)=[P(v)−F0]/PINC` 逐点对拍 `MODE=STATIC` / `TENSION`
     两条成桥态路径；初始斜率对拍独立共转切线 `K_yy`；
  4. 起点自洽：`STATE CTRL0 ≡ [U](2,2)`、`F0_NORM ≡ |F0|`、
     `FINC_NORM ≡ PINC`、`RES0/|F0| < 1e-6`；
  5. 两种 BASIS 分列：正常段 `completion_state`、未知基准 ERROR 段
     `BASIS=unknown` 不伪装成参考构形、默认输出无 `[STABILITY]`；
  6. 先走完成桥态闭环：`[NEWTON]` / `[LENGTH]` / `[BUCKLING]` 段仍在，
     TENSION 的 `[LENGTH] KIND=cable_length STATUS=ok` 在位；
  7. 默认零污染：既有 CB3D10 输出逐位一致、M 阶段 von Mises
     `_stability.txt` 逐字节一致、`vshape_default` 成桥位移仍是 −V0；
  8. 4 类 ERROR 负例（含 `no_stability_section`）不出平衡点、不静默降级；
  9. 独立导数锚点 `K_yy(V0)` 与中心差分 `dP/dv` 相对差 < 1e-8。

用法：
    python compare_completion_stability_3d.py --bin <build_dir> \
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
from build_reference_completion import (  # noqa: E402
    B, CASE_SPEC, E, F0, H, L0, PINC, TSTAR, V0,
    build, check, k_yy_of_v, lambda_of_v, numeric_dP_dv,
    N_of_v, P_of_v, parse_result,
)

REF_FILE = "cb3d20_completion_stability.json"
EXE = "cable_solver.exe"

TOL_LAMBDA_STATIC = 1.0e-9      # STATIC：闭式逐点应到机器精度
TOL_LAMBDA_TENSION = 1.0e-6    # TENSION：索长迭代命中 T* 的微差
TOL_DP = 0.1                   # N，|λ·PINC − (P(v)−F0)| 上限（覆盖 TENSION）
TOL_SLOPE_REL = 1.0e-2         # 第一点割线斜率 vs 解析 dλ/dv
TOL_CTRL0 = 1.0e-12            # CTRL0 ≡ [U](2,2)
TOL_F0_REL = 1.0e-6            # F0_NORM ≡ |F0|
TOL_RES_REL = 1.0e-6           # RES0/|F0|：起点闭合残差
TOL_KYY_REL = 1.0e-8           # 独立共转切线 vs 中心差分
TOL_LIVE_CSV = 0.0             # 冻结输出逐字节

FAILURES: list[str] = []


def ok(cond: bool, msg: str) -> None:
    if not cond:
        FAILURES.append(msg)


def run_exe(exe: Path, dat: Path, cwd: Path) -> tuple[int, str]:
    p = subprocess.run([str(exe), dat.name], cwd=str(cwd),
                       capture_output=True, text=True, errors="replace")
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def normalized(txt: str) -> list[str]:
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
    good, msg = check(doc, cases)
    ok(good, "黄金值与冻结输出不一致")
    print(f"[PASS] {msg} consistent with frozen output")
    crit = doc["criteria"]
    ok(abs(crit["analytic_P_at_V0"] - crit["analytic_F0"])
       <= 1e-9 * max(1.0, abs(crit["analytic_F0"])),
       "闭式 P(V0) ≠ F0（成桥态平衡不成立）")
    dk = abs(crit["analytic_k_yy_at_V0"] - crit["analytic_dP_dv_at_V0"])
    ok(dk <= TOL_KYY_REL * abs(crit["analytic_dP_dv_at_V0"]),
       f"独立共转切线 K_yy 与中心差分不符：{dk:.3e}")
    print("[PASS] independent oracle: P(V0)=F0, "
          f"|K_yy - dP/dv|={dk:.3e} N/m "
          f"(rel {dk/abs(crit['analytic_dP_dv_at_V0']):.3e})")

    with tempfile.TemporaryDirectory(prefix="completion_stab_") as td:
        tmp = Path(td)

        # ---- 2. FEM CLI 重跑 6 个 .dat --------------------------------
        expected_rc = {
            "vshape_static": 0, "vshape_tension": 0, "vshape_no_load": 0,
            "vshape_bad_basis": 3, "vshape_unconverged": 3,
            "vshape_default": 0,
        }
        repro = {}
        for stem, suffix, _, _, _ in CASE_SPEC:
            shutil.copy2(cases / f"{stem}.dat", tmp / f"{stem}.dat")
            rc, _ = run_exe(exe, tmp / f"{stem}.dat", tmp)
            got = tmp / f"{stem}{suffix}"
            frozen = cases / f"{stem}{suffix}"
            if not got.is_file():
                FAILURES.append(f"{stem} 未产出 {stem}{suffix}")
                continue
            cur = got.read_bytes()
            frozen_b = frozen.read_bytes()
            if suffix == "_cable.txt":
                # _cable.txt 头部有 elapsed_ms 运行耗时，按项目既有约定
                # 归一化后逐行比对；其余字段（含 [STABILITY] 全段）逐字节。
                repro[stem] = normalized(cur.decode("utf-8", "replace")) == \
                              normalized(frozen_b.decode("utf-8", "replace"))
            else:
                # _stability.txt 不含计时字段，直接逐字节。
                repro[stem] = (cur == frozen_b)
            ok(repro[stem], f"{stem} 重跑输出与冻结归档不一致")
            ok(rc == expected_rc[stem],
               f"{stem} exit={rc} != {expected_rc[stem]}")
        print(f"[PASS] FEM CLI reproduction: {len(repro)} cases, "
              f"byte-identical frozen output "
              f"({sum(repro.values())}/{len(repro)})")

        # ---- 3. 独立闭式 λ(v) 逐点对拍 ---------------------------------
        summary: list[str] = []
        for stem, suffix, _, _, _ in CASE_SPEC:
            if stem not in ("vshape_static", "vshape_tension"):
                continue
            r = parse_result(tmp / f"{stem}{suffix}")
            pts = r["points"]
            ok(len(pts) >= 200, f"{stem} 平衡点过少：{len(pts)}")
            worst_lam = 0.0
            worst_dp = 0.0
            for q in pts:
                v = V0 - q["ctrl"]           # CTRL 是相对成桥态的 y 增量
                worst_lam = max(worst_lam, abs(q["lambda"] - lambda_of_v(v)))
                worst_dp = max(worst_dp,
                               abs(q["lambda"] * PINC - (P_of_v(v) - F0)))
            tol = (TOL_LAMBDA_STATIC if stem == "vshape_static"
                   else TOL_LAMBDA_TENSION)
            ok(worst_lam <= tol,
               f"{stem} worst|λ−λ_closed|={worst_lam:.3e} > {tol}")
            ok(worst_dp <= TOL_DP,
               f"{stem} worst|λ·PINC−(P(v)−F0)|={worst_dp:.3e} N > {TOL_DP}")
            # 初始割线斜率 vs 独立共转切线 dλ/dv
            p1 = pts[1]
            # CTRL 为负（低点继续下移），故增量位移 = -CTRL。
            slope = p1["lambda"] / (-p1["ctrl"]) if p1["ctrl"] else 0.0
            dslope = k_yy_of_v(V0) / PINC
            ok(abs(slope - dslope) <= TOL_SLOPE_REL * dslope,
               f"{stem} 初始割线斜率 {slope:.6f} vs 解析 {dslope:.6f}")
            # ---- 4. 起点自洽 ------------------------------------------
            ok(abs(r["ctrl0"] - r["u_ctrl"]) <= TOL_CTRL0,
               f"{stem} CTRL0={r['ctrl0']} ≠ [U](2,2)={r['u_ctrl']}")
            ok(abs(r["f0_norm"] - F0) <= TOL_F0_REL * F0,
               f"{stem} F0_NORM={r['f0_norm']} ≠ |F0|={F0}")
            ok(abs(r["finc_norm"] - PINC) <= 1e-9 * PINC,
               f"{stem} FINC_NORM={r['finc_norm']} ≠ PINC")
            ok(r["f0_norm"] > 0.0 and r["res0"] / r["f0_norm"] < TOL_RES_REL,
               f"{stem} RES0/|F0|={r['res0']/r['f0_norm']:.3e} 过大")
            ok(max(q["res"] for q in pts) <= 1e-3 * F0,
               f"{stem} 逐点残差 {max(q['res'] for q in pts):.3e} N 过大")
            # ---- 6. 先走完成桥态闭环 ----------------------------------
            ok(r["has_newton"] == 1 and r["has_buckling_error"] == 1,
               f"{stem} 丢失 [NEWTON] 或 [BUCKLING] ERROR 段")
            ok(r["basis"] == "completion_state",
               f"{stem} BASIS={r['basis']} != completion_state")
            summary.append(
                f"   {stem}: {len(pts)} points, worst|dlam|={worst_lam:.3e}, "
                f"worst|dP|={worst_dp:.3e} N, lambda_final={r['lambda_final']:.9f}, "
                f"ctrl0={r['ctrl0']:.9f}, res0={r['res0']:.3e}, "
                f"slope0={slope:.6f} vs {dslope:.6f}")
        print("[PASS] closed-form lambda(v) / initial tangent / "
              "start-state self-consistency")
        for s in summary:
            print(s)
        ok(parse_result(tmp / "vshape_tension_cable.txt")["has_length"] == 1,
           "MODE=TENSION 的 [LENGTH] 索长迭代过程段丢失")

        # ---- 5. 两种 BASIS 分列 ---------------------------------------
        bad = parse_result(tmp / "vshape_bad_basis_stability.txt")
        ok(bad["basis"] == "unknown" and not bad["points"],
           "未知基准被静默降级成参考构形 / 却产出平衡点")
        dflt = (tmp / "vshape_default_cable.txt").read_text(
            encoding="utf-8", errors="replace")
        ok("[STABILITY]" not in dflt, "默认路径泄漏了 [STABILITY] 段")
        ok("ERROR=nonlinear_state_not_supported" in dflt,
           "默认路径丢失 [BUCKLING] ERROR 段")
        print("[PASS] BASIS separation: completion_state / unknown / no section")

        # ---- 7. 默认零污染（跨算例） ----------------------------------
        bench = cases.parent.parent
        base_dat = bench / "nonlinear_ansys_3d/cases/cb3d10_completion.dat"
        base_txt = bench / "nonlinear_ansys_3d/cases/cb3d10_completion_cable.txt"
        if base_dat.is_file() and base_txt.is_file():
            shutil.copy2(base_dat, tmp / "cb3d10_off.dat")
            rc, _ = run_exe(exe, tmp / "cb3d10_off.dat", tmp)
            got = tmp / "cb3d10_off_cable.txt"
            if got.is_file():
                txt = got.read_text(encoding="utf-8", errors="replace")
                ok(rc == 0, "既有 CB3D10 默认路径退出码非 0")
                ok("[STABILITY]" not in txt,
                   "既有 CB3D10 默认路径泄漏了 [STABILITY] 段")
                ok(normalized(txt) == normalized(base_txt.read_text(
                    encoding="utf-8", errors="replace")),
                   "既有 CB3D10 输出与冻结 _cable.txt 不一致（除 elapsed_ms）")
        else:
            FAILURES.append("缺少既有 CB3D10 零污染对照算例")

        # M 阶段参考构形基准逐字节回归（证明新键不污染第十七次口径）
        vm_dat = bench / "stable_path_3d/cases/von_mises_log.dat"
        vm_txt = bench / "stable_path_3d/cases/von_mises_log_stability.txt"
        if vm_dat.is_file() and vm_txt.is_file():
            shutil.copy2(vm_dat, tmp / "vm_log.dat")
            rc, _ = run_exe(exe, tmp / "vm_log.dat", tmp)
            got = tmp / "vm_log_stability.txt"
            ok(got.is_file() and got.read_bytes() == vm_txt.read_bytes(),
               "M 阶段 von Mises REFERENCE 基准输出漂移")
            ok(rc == 0, "M 阶段 von Mises 退出码非 0")
        else:
            FAILURES.append("缺少 M 阶段 von Mises 冻结参考输出")

        # *STABILITY_LOAD 不进 m.loads：默认态位移必须仍是 -V0
        ud = parse_result(tmp / "vshape_default_cable.txt")["u_ctrl"]
        ok(ud is not None and abs(ud + V0) <= 1e-9,
           f"默认态 U(2,2)={ud} ≠ -V0（增量荷载混进了 m.loads）")
        print("[PASS] zero pollution: CB3D10 / von Mises REFERENCE / "
              "*STABILITY_LOAD not in m.loads")

        # ---- 8. ERROR 负例 --------------------------------------------
        expect = {
            "vshape_no_load": ("zero_incremental_load", "_cable.txt"),
            "vshape_bad_basis": ("unsupported_stability_basis",
                                 "_stability.txt"),
            "vshape_unconverged": ("completion_not_converged", "_cable.txt"),
            "vshape_default": ("no_stability_section", "_cable.txt"),
        }
        for stem, (code, suffix) in expect.items():
            r = parse_result(tmp / f"{stem}{suffix}")
            ok(r["error"] == code, f"{stem} ERROR={r['error']} != {code}")
            ok(r["points"] == [] and r["npoint"] == 0,
               f"{stem} ERROR 段却产出平衡点")
            ok(r["has_buckling_error"] == 1,
               f"{stem} 丢失 [BUCKLING] ERROR 段")
            print(f"[PASS] {stem}: ERROR={code}, no points, BUCKLING kept")

    if FAILURES:
        print(f"\n[FAIL] {len(FAILURES)} 项未通过：")
        for f in FAILURES:
            print("  - " + f)
        return 1
    print("\n[PASS] completion_stability_3d benchmark: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
