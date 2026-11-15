#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 23 项 CTest：真实斜拉桥「成桥态起点全过程稳定」两段对标（离线）。

与第 22 项的关系：那里只在 FEM 侧用 V 形双索单自由度闭式解验证
`STABASIS=COMPLETION` 的力学自洽；这里把同一个主链路 opt-in 放到**真实
双塔三跨斜拉桥**（50 节点 / 63 单元 / 16 索），并与 ANSYS 18.0
「step1 成桥态恒载 + step2 增量荷载」两段求解逐子步对拍。

检查项：
  1. 黄金值自洽（reference JSON ↔ 冻结 cases/*）；
  2. FEM CLI 重跑 2 个 .dat，输出除 elapsed_ms 外逐行一致、退出码 0；
  3. 起点自洽：`STATE CTRL0 ≡ [U](1016,3)`、
     `RES0/|F0| < 1e-6`、`RES0` 与 `[NEWTON] RESLAST` 同量级；
  4. FEM 路径单调硬化、λ_final 覆盖 ANSYS 全部子步；
  5. FEM↔ANSYS 逐子步：`λ=(TIME−1)*LAMMAX`、
     `ΔUZ=UZ−UZ_step1` 对 FEM `CTRL` 插值，worst ≤ 1e-4 m；
     初始线性段斜率相对差 ≤ 2e-2；
  6. ARCLEN job 与荷载控制主 job 终点 UZ 一致（TIME 不作 λ 换算）；
  7. 同源导出独立锚点：16 根索 `L0*(1+EPS0)` ≡ 设计线形弦长、
     16 条 `INISTATE ... EPEL,X,eps` ≡ EPS0、
     APDL `F_end` ≡ `F_orig − LAMMAX*PINC`；
  8. BASIS 分列：主 job 荷载控制、arc job ARCLEN，两个 job 都恰有
     2 个 SOLVE（step1 + step2），不得拿参考构形 ARCLEN 冒充成桥态；
  9. 默认零污染：cs_default 无 `[STABILITY]`、既有 CB3D10 与 M 阶段
     von Mises REFERENCE 输出逐位不变；
 10. 冻结对比报告 `ALL CHECKS PASSED`、两个 ANSYS job 0 error；
 11. 可选 --live-ansys：临时目录重跑 2 个 ANSYS job，逐 CSV 单元 diff，
     绝不写回冻结 cases/。

用法：
    python compare_completion_bridge_3d.py --bin <build_dir> \
        --cases <cases_dir> --reference <reference_dir> [--live-ansys]
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_reference_bridge import (  # noqa: E402
    ANSYS_JOBS, DEFAULT_STEM, REF, STEM, build, check, interp_ctrl,
    parse_apdl, parse_curve, parse_log, parse_stability, split_curve,
)

EXE = "cable_solver.exe"
DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"

TOL_RES0_REL = 1.0e-6
TOL_RES0_NEWTON = 10.0
TOL_CTRL0 = 1.0e-10
TOL_CURVE = 1.0e-4
TOL_SLOPE_REL = 2.0e-2
TOL_ENDPOINT = 1.0e-4
TOL_START_UZ = 1.0e-4
TOL_SAME_SOURCE = 1.0e-9       # 几何 / eps 同源精度
TOL_F_END = 1.0e-5             # APDL 只写 15 位有效数字，N
TOL_LIVE_CSV = 1.0e-6
TOL_LIVE_CSV_REL = 1.0e-10

FAILURES: list[str] = []


def ok(cond: bool, msg: str) -> None:
    if not cond:
        FAILURES.append(msg)


def normalized(txt: str) -> list[str]:
    return [re.sub(r"elapsed_ms=[0-9.eE+-]+", "elapsed_ms=*", ln)
            for ln in txt.splitlines()]


def run_exe(exe: Path, dat: Path, cwd: Path) -> tuple[int, str]:
    p = subprocess.run([str(exe), dat.name], cwd=str(cwd),
                       capture_output=True, text=True, errors="replace",
                       timeout=900)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def live_ansys(cases: Path, ansys_exe: str, tmp: Path) -> None:
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
        got = parse_curve(got_path)
        ref = parse_curve(cases / (job + "_curve.csv"))
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
    ap.add_argument("--bin", required=True, help="含 cable_solver.exe 的构建目录")
    ap.add_argument("--cases", required=True, help="冻结 cases/ 目录")
    ap.add_argument("--reference", required=True, help="冻结 reference/ 目录")
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys-exe", default=DEFAULT_ANSYS)
    a = ap.parse_args()

    exe = Path(a.bin) / EXE
    cases = Path(a.cases)
    ref = Path(a.reference) / REF.name
    if not exe.is_file():
        print(f"[FATAL] 找不到 {exe}")
        return 1
    if not ref.is_file():
        print(f"[FATAL] 找不到 {ref}")
        return 1

    # ---- 1. 黄金值自洽 -------------------------------------------------
    doc = json.loads(ref.read_text(encoding="utf-8"))
    good, msg = check(doc, cases)
    ok(good, "黄金值与冻结证据不一致")
    print(f"[PASS] {msg} consistent with frozen evidence")
    crit = doc["criteria"]
    strings = doc["strings"]
    EPS0, PINC, LAMMAX = (float(crit["eps0"]), float(crit["pinc"]),
                          float(crit["lammax"]))
    CN, CD = int(crit["ctrl_node"]), int(crit["ctrl_dof"])

    # ---- 2. 独立同源锚点（纯 Python）-----------------------------------
    ok(crit["n_cables"] == 16 and crit["n_cable_l0"] == 16
       and crit["n_inistate"] == 16,
       "索 / L0 / INISTATE 数量不齐：%s" % [
           crit["n_cables"], crit["n_cable_l0"], crit["n_inistate"]])
    ok(crit["cable_chord_dev_max"] <= TOL_SAME_SOURCE,
       "L0*(1+EPS0) 与设计线形弦长不同源：%.3e" %
       crit["cable_chord_dev_max"])
    ok(crit["inistate_eps_dev_max"] <= TOL_SAME_SOURCE,
       "INISTATE eps 与 EPS0 不同源：%.3e" % crit["inistate_eps_dev_max"])
    ok(crit["f_end_apdl_dev"] <= TOL_F_END,
       "APDL F_end 与 F_orig-LAMMAX*PINC 不符：%.3e N" %
       crit["f_end_apdl_dev"])
    print("[PASS] independent same-source anchors: "
          "16 cables L0*(1+EPS0)=chord (dev %.3e), "
          "INISTATE eps dev %.3e, F_end dev %.3e N" % (
              crit["cable_chord_dev_max"], crit["inistate_eps_dev_max"],
              crit["f_end_apdl_dev"]))

    # ---- 3. 起点自洽 & FEM 路径 ----------------------------------------
    ok(strings["fem/basis"] == "completion_state",
       "BASIS=%s != completion_state" % strings["fem/basis"])
    ok(strings["fem/error"] == "none" and strings["fem/status"] == "ok",
       "FEM 稳定段 ERROR/STATUS 异常")
    ok(crit["res0_over_f0"] <= TOL_RES0_REL,
       "RES0/|F0|=%.3e 过大" % crit["res0_over_f0"])
    ok(crit["res0_over_reslast"] <= TOL_RES0_NEWTON,
       "RES0/[NEWTON] RESLAST=%.3f 不同量级" % crit["res0_over_reslast"])
    ok(abs(crit["fem_ctrl0"] - crit["ansys_uz_step1"]) <= TOL_START_UZ,
       "FEM CTRL0 与 ANSYS step1 UZ 差 %.3e m" % crit["start_uz_gap"])
    ok(crit["fem_lambda_final"] >= LAMMAX - 1e-8,
       "FEM 路径未覆盖 LAMMAX：%.9g < %.9g" % (
           crit["fem_lambda_final"], LAMMAX))
    print("[PASS] start-state self-consistency: "
          "RES0=%.4g N (F0=%.6g N, RESLAST=%.4g N), "
          "CTRL0=%.9g m, ANSYS step1=%.9g m (gap %.3e), "
          "lambda_final=%.9g" % (
              crit["fem_res0"], crit["fem_f0_norm"], crit["fem_reslast"],
              crit["fem_ctrl0"], crit["ansys_uz_step1"],
              crit["start_uz_gap"], crit["fem_lambda_final"]))

    # ---- 4. 冻结 ANSYS 曲线 vs FEM 路径 --------------------------------
    frozen_fem = parse_stability(cases / (STEM + "_cable.txt"), CN, CD)
    fem_path = [(q["lambda"], q["ctrl"]) for q in frozen_fem["points"]]
    ok(len(fem_path) == int(crit["fem_npoint"]) and len(fem_path) >= 200,
       "FEM 平衡点过少：%d" % len(fem_path))
    lams = [p[0] for p in fem_path]
    ok(all(lams[i] >= lams[i - 1] for i in range(1, len(lams))),
       "FEM 路径非单调硬化（本算例口径）")

    rows_load = parse_curve(cases / "cs_load_curve.csv")
    uz0, s1, s2 = split_curve(rows_load)
    worst, wlam, cnt = 0.0, 0.0, 0
    for iset, tnow, uz, _ in s2:
        lam = (tnow - 1.0) * LAMMAX
        if lam < -1e-8 or lam > lams[-1] + 1e-8:
            continue
        d = abs((uz - uz0) - interp_ctrl(fem_path, lam))
        cnt += 1
        if d > worst:
            worst, wlam = d, lam
    ok(cnt >= 20 and worst <= TOL_CURVE,
       "FEM↔ANSYS 逐子步 worst|dUZ|=%.3e m (n=%d) > %.1e" % (
           worst, cnt, TOL_CURVE))
    ok(abs(worst - crit["worst_duz"]) <= 1e-12,
       "重算 worst|dUZ| 与黄金值不符")
    ok(crit["initial_slope_rel"] <= TOL_SLOPE_REL,
       "初始斜率相对差 %.3e > %.1e" % (
           crit["initial_slope_rel"], TOL_SLOPE_REL))
    rows_arc = parse_curve(cases / "cs_arc_curve.csv")
    gap = abs(rows_arc[-1][2] - rows_load[-1][2])
    ok(gap <= TOL_ENDPOINT,
       "ARCLEN / 荷载控制终点 UZ 差 %.3e m > %.1e" % (gap, TOL_ENDPOINT))
    print("[PASS] FEM vs ANSYS two-step: %d substeps, "
          "worst|dUZ|=%.3e m @ lambda=%.6g, "
          "initial slope rel=%.3e, endpoint gap=%.3e m" % (
              cnt, worst, wlam, crit["initial_slope_rel"], gap))

    # ---- 5. 两段 APDL 语义 / 日志 / 报告 --------------------------------
    for job in ANSYS_JOBS:
        ad = parse_apdl(cases / (job + ".inp"))
        er, wa = parse_log(cases / (job + "_ansys.out"))
        ok(ad["n_solve"] == 2, "%s SOLVE 数=%d（应为 step1+step2 两段）" % (
            job, ad["n_solve"]))
        ok(ad["n_loads"] == 2, "%s F,1016,FZ 数=%d（应为 2）" % (
            job, ad["n_loads"]))
        ok(ad["n_inistate"] == 16, "%s INISTATE 数=%d" % (
            job, ad["n_inistate"]))
        ok(er == 0, "%s ANSYS 日志 %d error" % (job, er))
        ok(strings["log/%s/summary" % job] == "errors=%d warnings=%d" % (er, wa),
           "%s 日志摘要与黄金值不符" % job)
    ok(not parse_apdl(cases / "cs_load.inp")["arclen"],
       "主 job 不应启用 ARCLEN（λ 用 TIME 线性换算）")
    ok(parse_apdl(cases / "cs_arc.inp")["arclen"],
       "交叉验证 job 缺少 ARCLEN,ON")
    ok(strings["apdl/cs_load/step2"] == "load_control"
       and strings["apdl/cs_arc/step2"] == "ARCLEN",
       "两段求解类型登记不符")
    ok(strings["report/result"] == "ALL CHECKS PASSED",
       "冻结对比报告结论不是 ALL CHECKS PASSED")
    ok(strings["cfg/mode"] == "STATIC",
       "MODE 不是 STATIC（COMPLETION 会改写 L0 导致起点不同源）")
    print("[PASS] two-step APDL semantics: 2 SOLVE / 2 F / 16 INISTATE "
          "per job, load-control vs ARCLEN, 0 error, report ALL CHECKS PASSED")

    with tempfile.TemporaryDirectory(prefix="completion_bridge_") as td:
        tmp = Path(td)

        # ---- 6. FEM CLI 离线重跑 ---------------------------------------
        for stem in (STEM, DEFAULT_STEM):
            shutil.copy2(cases / (stem + ".dat"), tmp / (stem + ".dat"))
            rc, _ = run_exe(exe, tmp / (stem + ".dat"), tmp)
            got = tmp / (stem + "_cable.txt")
            frozen = cases / (stem + "_cable.txt")
            if not got.is_file():
                FAILURES.append("%s 未产出 %s_cable.txt" % (stem, stem))
                continue
            ok(rc == 0, "%s exit=%d != 0" % (stem, rc))
            ok(normalized(got.read_text(encoding="utf-8", errors="replace"))
               == normalized(frozen.read_text(encoding="utf-8",
                                              errors="replace")),
               "%s 重跑输出与冻结归档不一致" % stem)
        print("[PASS] FEM CLI reproduction: 2 cases, "
              "byte-identical frozen _cable.txt (except elapsed_ms)")

        # ---- 7. 默认零污染（含跨算例）---------------------------------
        fresh_fem = parse_stability(tmp / (STEM + "_cable.txt"), CN, CD)
        dflt = parse_stability(tmp / (DEFAULT_STEM + "_cable.txt"), CN, CD)
        ok(fresh_fem["basis"] == "completion_state"
           and len(fresh_fem["points"]) == len(fem_path),
           "重跑主算例 BASIS / 点数漂移")
        ok(not dflt["has_section"] and not dflt["points"],
           "cs_default 泄漏了 [STABILITY] 段")
        ok(dflt["has_buckling_error"] == 1 and dflt["has_newton"] == 1,
           "cs_default 丢失 [BUCKLING] ERROR / [NEWTON]")
        ok(strings["default/error"] == "no_stability_section",
           "默认负例错误码登记不符")
        bench = cases.parent.parent
        base_dat = bench / "nonlinear_ansys_3d/cases/cb3d10_completion.dat"
        base_txt = bench / "nonlinear_ansys_3d/cases/cb3d10_completion_cable.txt"
        if base_dat.is_file() and base_txt.is_file():
            shutil.copy2(base_dat, tmp / "cb3d10_off.dat")
            rc, _ = run_exe(exe, tmp / "cb3d10_off.dat", tmp)
            got = tmp / "cb3d10_off_cable.txt"
            ok(got.is_file() and rc == 0, "既有 CB3D10 默认路径重跑失败")
            if got.is_file():
                txt = got.read_text(encoding="utf-8", errors="replace")
                ok("[STABILITY]" not in txt,
                   "既有 CB3D10 默认路径泄漏了 [STABILITY] 段")
                ok(normalized(txt) == normalized(base_txt.read_text(
                    encoding="utf-8", errors="replace")),
                   "既有 CB3D10 输出与冻结归档不一致")
        else:
            FAILURES.append("缺少既有 CB3D10 零污染对照算例")
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
        print("[PASS] zero pollution: default / CB3D10 / "
              "von Mises REFERENCE unchanged")

        # ---- 8. 可选 live-ansys ----------------------------------------
        if a.live_ansys:
            live_ansys(cases, a.ansys_exe, tmp)

    if FAILURES:
        print(f"\n[FAIL] {len(FAILURES)} 项未通过：")
        for f in FAILURES:
            print("  - " + f)
        return 1
    print("\n[PASS] completion_bridge_3d benchmark: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
