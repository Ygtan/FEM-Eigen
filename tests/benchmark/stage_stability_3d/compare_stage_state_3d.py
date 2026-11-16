#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 26 项 CTest：CB3D23 施工阶段中间态全过程稳定，纯离线比对。

    python compare_stage_state_3d.py --bin <build dir> \
        --cases <frozen cases dir> --reference <golden dir>

判据分五层：
  1. 15 份冻结证据的 sha256 / 字节数（不放宽）；
  2. 黄金值自洽：用同一套解析函数从冻结证据重新派生，逐字段核对；
  3. 物理判据：BASIS=stage_state、STAGE=中间态、起点闭合 RES0/|F0|、
     CTRL0 != 0、逐阶段递进、路径单调收敛；
  4. 同源锚点 + ANSYS 两段求解证据（APDL N/E/F/D 与 .dat 逐项一致，
     两个 job 0 error、λ=TIME−1）与 FEM↔ANSYS 逐点对拍；
  5. 当前二进制离线重放 `<stem>_stage_stability.txt`，逐字节一致
     （冻结证据由 mingw64 工具链产出，见 tests/benchmark/README.md）。

不重跑 ANSYS，不写回冻结证据，不修改任何源文件。
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

from build_reference_stage_state import (  # noqa: E402
    ANSYS_CASES, FILES, LAMMAX, PINC, REF_NAME, STAGE_ID, STEM,
    CTRLTARGET, BAD_STEM, ansys_errors, ansys_step2_curve, apdl_anchor,
    build, compare_curve, fem_metrics, normalized, read_curve_csv,
    read_text, sha256_file,
)

TOL_RES0_REL = 1.0e-9
TOL_CTRL0_ABS = 1.0e-6
TOL_LAMBDA_PT = 6.0e-3
TOL_LAMBDA_MEAN = 2.5e-3
EPS = 1.0e-12

FAILED: list[str] = []
CHECKS = [0]


def ok(cond: bool, msg: str) -> None:
    CHECKS[0] += 1
    if not cond:
        FAILED.append(msg)


def close(a: float, b: float, tol: float = 1e-12) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", default=str(HERE / "cases"))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    args = ap.parse_args()
    cases, refdir, bindir = Path(args.cases), Path(args.reference), Path(args.bin)

    refpath = refdir / REF_NAME
    ok(refpath.is_file(), f"黄金值缺失: {refpath}")
    if not refpath.is_file():
        return 1
    ref = json.loads(refpath.read_text(encoding="utf-8"))

    # ---- 1. 冻结证据 ----------------------------------------------------
    ok(ref.get("case") == "CB3D23", "黄金值 case 标记不符")
    ok(ref.get("basis") == "stage_state", "黄金值 basis 不是 stage_state")
    ok(ref.get("n_files") == len(FILES) == 15, "黄金值文件数不是 15")
    for name in FILES:
        p = cases / name
        ok(p.is_file(), f"缺冻结文件 {name}")
        if not p.is_file():
            continue
        meta = ref["files"].get(name, {})
        ok(meta.get("sha256") == sha256_file(p), f"{name} sha256 与黄金值不符")
        ok(meta.get("bytes") == p.stat().st_size, f"{name} 字节数与黄金值不符")

    # ---- 2. 黄金值自洽：重新派生 ---------------------------------------
    fresh = build(cases)
    ok(len(fresh["files"]) == len(ref["files"]), "冻结文件清单长度不符")
    for name in FILES:
        ok(fresh["files"][name]["sha256"] == ref["files"][name]["sha256"],
           f"{name} 重新派生 sha256 不一致")
    for key in ("lammax", "pinc", "ctrltarget", "stage_id"):
        ok(fresh[key] == ref[key], f"黄金值 {key} 与冻结证据不符")

    f_ref, f_new = ref["fem"], fresh["fem"]
    for key in ("nstep", "npoint", "n_descending", "head", "state", "stage"):
        ok(f_new[key] == f_ref[key], f"FEM {key} 与黄金值不符")
    for key in ("res0_rel", "lambda_end", "ctrl_end", "lambda_max"):
        ok(close(f_new[key], f_ref[key], 1e-11), f"FEM {key} 与黄金值不符")
    ok(len(f_new["stages"]) == len(f_ref["stages"]), "[STAGES] 阶段数不符")
    for a, b in zip(f_new["stages"], f_ref["stages"]):
        ok(a == b, f"[STAGES] {a.get('ID')} 与黄金值不符")
    ok(f_new["points"] == f_ref["points"], "FEM 逐点路径与黄金值不符")
    ok(len(f_ref["points"]) == f_ref["npoint"] >= 100,
       "黄金值逐点路径长度自洽性失败")

    for job in ANSYS_CASES:
        a_ref, a_new = ref["ansys"][job], fresh["ansys"][job]
        for key in ("n_rows", "n_step2", "errors", "n_solve_inp", "has_arclen"):
            ok(a_new[key] == a_ref[key], f"{job} {key} 与黄金值不符")
        ok(close(a_new["u0"], a_ref["u0"], 1e-12), f"{job} u0 与黄金值不符")
        ok(close(a_new["lam_end"], a_ref["lam_end"], 1e-12),
           f"{job} lam_end 与黄金值不符")
        ok(a_new["curve"] == a_ref["curve"], f"{job} step2 曲线与黄金值不符")
        c_ref, c_new = ref["comparisons"][job], fresh["comparisons"][job]
        for key in ("n", "worst", "mean", "lo", "hi", "coverage", "fem_span"):
            ok(close(c_new[key], c_ref[key], 1e-10),
               f"{job} 对拍 {key} 与黄金值不符")
        an_ref, an_new = ref["anchors"][job], fresh["anchors"][job]
        for key in ("n_nodes_dat", "n_nodes_inp", "n_elems_inp", "coord_max_dev",
                    "n_fixes", "n_solve", "loads_step1", "loads_step2",
                    "elems_inp", "fixes", "nodes_dat"):
            ok(an_new[key] == an_ref[key], f"{job} 锚点 {key} 与黄金值不符")
    for key in ("negative", "config", "report"):
        ok(fresh[key] == ref[key], f"黄金值 {key} 与冻结证据不符")

    # ---- 3. FEM 物理判据 -------------------------------------------------
    head, state, stage = f_ref["head"], f_ref["state"], f_ref["stage"]
    stages = f_ref["stages"]
    ok(head.get("BASIS") == "stage_state", "BASIS 不是 stage_state")
    ok(head.get("STATUS") == "ok" and head.get("CONVERGED") == "1",
       "FEM 弧长路径未收敛")
    ok(head.get("PATHMETHOD") == "SPHERE", "路径方法不是球面弧长法")
    ok(close(stage.get("STAGE", -1), STAGE_ID), "STAGE= 不是指定中间态 2")
    ok(stage.get("NELEM", 0) == len(ref["config"]["stage_active_elements"]),
       "STAGE NELEM 与阶段 2 激活单元数不符")
    ok(stage.get("NDOF", 0) == 12, "STAGE NDOF 不是 12")
    ok([int(s["ID"]) for s in stages] == [1, 2, 3], "[STAGES] 阶段序列不符")
    ok(len(stages) > STAGE_ID, "选中阶段不是中间态（后面没有阶段）")
    by_id = {int(s["ID"]): s for s in stages}
    ok(by_id[STAGE_ID]["CONVERGED"] == 1, "阶段 2 共转列式平衡未收敛")
    ok(by_id[1]["NELEM"] < by_id[STAGE_ID]["NELEM"], "阶段 1→2 激活单元无递进")
    ok(by_id[STAGE_ID]["NDLOAD"] == len(ref["config"]["stage_constant_loads"]),
       "阶段 2 累计阶段荷载条数不符")
    ok(by_id[3]["ID"] == 3, "阶段 3 缺失")
    ok(f_ref["res0_rel"] < TOL_RES0_REL,
       "起点闭合残差 RES0/|F0| 超容差: %.3g" % f_ref["res0_rel"])
    ctrl0 = state.get("CTRL0", 0.0)
    ok(abs(ctrl0) > 1e-4, "起点 CTRL0 近似 0，疑似退化为参考构形")
    ok(state.get("F0_NORM", 0) > 0 and state.get("FINC_NORM", 0) > 0,
       "F0 / Finc 范数缺失")
    ok(f_ref["npoint"] >= 100 and f_ref["n_descending"] == 0,
       "FEM 路径点数不足或 λ 非单调")
    ok(close(abs(f_ref["ctrl_end"]), CTRLTARGET, 0.1),
       "FEM 终点 |ΔUY| 与 CTRLTARGET 不符")
    ok(f_ref["lambda_end"] > 1.0, "FEM 终点 λ 过小")

    # 负例：未知 STAGEID 必须显式报错，不得静默降级。
    neg = ref["negative"]
    ok(neg["error"] == "stage_not_found", "负例未写 stage_not_found")
    ok(neg["basis"] == "stage_state", "负例 BASIS 不是 stage_state")
    ok(neg["has_points"] is False, "负例不应产出平衡点")

    # ---- 4. 同源锚点 / ANSYS 两段 / 逐点对拍 -----------------------------
    cfg = ref["config"]
    const_loads = {int(k): v for k, v in cfg["stage_constant_loads"].items()}
    active = cfg["stage_active_elements"]
    ok(cfg["stage_id"] == STAGE_ID and cfg["pinc"] == PINC, "cfg.json 不同源")
    ok(active == [1, 2, 3] and len(active) == 3, "阶段 2 激活单元清单不符")
    ok(cfg["fixed_x"] == [2, 3], "侧向导向约束与配置不符")
    for job in ANSYS_CASES:
        job_ref = ref["ansys"][job]
        an = ref["anchors"][job]
        ok(an["n_nodes_dat"] == 4 and an["n_nodes_inp"] == 4,
           f"{job} 节点数不是 4")
        ok(an["coord_max_dev"] < 1e-12, f"{job} APDL 节点坐标与 .dat 不同源")
        ok(an["n_elems_inp"] == 3, f"{job} APDL 单元数不是 3")
        ok(an["n_fixes"] == 10, f"{job} 约束条数不是 10")
        ok(an["n_solve"] == 2 and job_ref["n_solve_inp"] == 2,
           f"{job} 不是两段求解")
        ok(job_ref["errors"] == 0, f"{job} ANSYS 有 error")
        ok(job_ref["n_step2"] >= 20, f"{job} step2 点数过少")
        ok(close(job_ref["lam_end"], LAMMAX, 1e-6),
           f"{job} 终点 λ 不是 {LAMMAX}")
        ok(abs(job_ref["u0"] - ctrl0) < TOL_CTRL0_ABS,
           f"{job} step1 UY 与 FEM CTRL0 不符")
        ok(job_ref["has_arclen"] == (job == "st_arc"),
           f"{job} ARCLEN 开关不符")
        s1 = {int(k): v for k, v in an["loads_step1"].items()}
        s2 = {int(k): v for k, v in an["loads_step2"].items()}
        ok(s1 == const_loads, f"{job} step1 恒定荷载与阶段 2 累计荷载不同源")
        ok(all(close(s2[k], const_loads[k] + LAMMAX * PINC, 1e-9)
               for k in const_loads),
           f"{job} step2 终点荷载不满足 F_orig+LAMMAX·PINC")
        inp_txt = read_text(cases / f"{job}.inp")
        ok("NLGEOM,ON" in inp_txt and "TIME,1" in inp_txt,
           f"{job} 缺 NLGEOM / step1 TIME=1")
        ok(f"TIME,{1 + LAMMAX:.12g}" in inp_txt, f"{job} 缺 step2 TIME")
        for nid in cfg["fixed_x"]:
            ok(f"D,{nid},UX,0.0" in inp_txt,
               f"{job} 缺 D,{nid},UX 侧向导向约束")
        out_txt = read_text(cases / f"{job}_ansys.out")
        ok(ansys_errors(out_txt) == 0, f"{job} .out 中 error 计数非 0")

        c = ref["comparisons"][job]
        ok(c["n"] >= 100, f"{job} 共同区间逐点样本过少")
        ok(c["coverage"] > 0.99, f"{job} 共同区间覆盖不足")
        ok(c["worst"] < TOL_LAMBDA_PT,
           f"{job} worst|Δλ|={c['worst']:.4g} 超容差")
        ok(c["mean"] < TOL_LAMBDA_MEAN,
           f"{job} mean|Δλ|={c['mean']:.4g} 超容差")
    ok(ref["report"]["result"] == "ALL CHECKS PASSED"
       and ref["report"]["failed"] == 0
       and ref["report"]["checks"] >= 70,
       "BS 侧 compare_report 不是全绿")

    # ---- 5. 当前二进制离线重放 -------------------------------------------
    exe = bindir / "cable_solver.exe"
    ok(exe.is_file(), f"缺 cable_solver.exe: {exe}")
    if exe.is_file():
        for stem in (STEM, BAD_STEM):
            with tempfile.TemporaryDirectory(prefix="cb3d23_replay_") as td:
                tmp = Path(td)
                shutil.copy2(cases / f"{stem}.dat", tmp / f"{stem}.dat")
                try:
                    r = subprocess.run([str(exe), f"{stem}.dat"], cwd=str(tmp),
                                       timeout=600, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True,
                                       errors="replace")
                    ok(r.returncode in (0, 3),
                       f"{stem} 重放退出码异常 {r.returncode}: {r.stdout[-200:]}")
                except Exception as e:  # noqa: BLE001
                    ok(False, f"{stem} 重放异常: {e}")
                    continue
                got = tmp / f"{stem}_stage_stability.txt"
                ok(got.is_file(), f"{stem} 重放未产出结果文件")
                if not got.is_file():
                    continue
                frozen = read_text(cases / f"{stem}_stage_stability.txt")
                fresh_txt = normalized(got.read_text(encoding="utf-8",
                                                     errors="replace"))
                if fresh_txt != frozen:
                    # 给出首个差异行，便于判断是工具链浮点漂移还是真回归。
                    a, b = frozen.splitlines(), fresh_txt.splitlines()
                    diff = next((i for i in range(min(len(a), len(b)))
                                 if a[i] != b[i]), min(len(a), len(b)))
                    ok(False, f"{stem} 离线重放与冻结证据不一致 "
                              f"(行 {diff + 1}, 行数 {len(a)}/{len(b)})；"
                              f"确认使用 mingw64 规范工具链")
                else:
                    ok(True, "")

    print("=" * 78)
    print("第 26 项 CTest：CB3D23 施工阶段中间态全过程稳定（BASIS=stage_state）")
    print("=" * 78)
    print("[PASS] 冻结证据 15 份 sha256 / 字节数与黄金值一致")
    print("[PASS] 黄金值自洽：FEM %d 点 / [STAGES] %d 阶段 / 2 条 ANSYS 曲线"
          % (f_ref["npoint"], len(stages)))
    print("[PASS] 物理判据：STAGE=2 NELEM=3 NDOF=12 RES0/|F0|=%.3g "
          "CTRL0=%.9g（非零初态）" % (f_ref["res0_rel"], ctrl0))
    for job in ANSYS_CASES:
        c = ref["comparisons"][job]
        print("[PASS] %s: %d 点逐点对拍 worst|Δλ|=%.4g mean=%.4g 覆盖 %.1f%%"
              % (job, c["n"], c["worst"], c["mean"], 100 * c["coverage"]))
    print("[PASS] 当前二进制离线重放 2 个 case 逐字节一致")
    print("[%s] stage_stability_3d benchmark: %d checks, %d failed"
          % ("PASS" if not FAILED else "FAIL", CHECKS[0], len(FAILED)))
    for f in FAILED:
        print("  [FAIL] " + f)
    return 0 if not FAILED else 1


if __name__ == "__main__":
    raise SystemExit(main())
