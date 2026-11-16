#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry: 第 28 项 —— 第二十一次力学①位移控制 / 线搜索（opt-in）。

离线（无需 ANSYS / 求解器许可）：
  * 在临时目录用当前 snap_probe.exe 重跑冻结 .snap，逐值复现 _fem.csv；
  * 三条正常路径逐点对齐**独立闭式** P(v)（不读求解器）；
  * line search 与无 line search 走同一平衡路径（逐点 |Δλ|≈0），
    且 line search 确被触发（cuts>0）、无 line search 不触发（cuts==0）；
  * 细步 λ_max 比粗步更接近闭式 λ_cr；
  * 负例（控制 DOF 指向固定端）必须报错且不收敛。

用法：
    python compare_displacement_control_3d.py --bin <dir_with_snap_probe.exe> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from build_reference import (  # noqa: E402
    CASES_IDS, NORMAL, analytic, p_of_v, sha256,
)
from build_reference import PREF, H  # noqa: E402

MSYS_BIN = r"C:\msys64\mingw64\bin"
TOL_REPRO = 0.0            # 冻结 CSV 逐位复现（同工具链）
SNAP = "cb3d25_displacement_control.snap"
FROZEN = "cb3d25_displacement_control_fem.csv"
REF_JSON = "cb3d25_displacement_control.json"

RESULTS: list[tuple[bool, str, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((bool(ok), name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def read_csv(path: Path) -> tuple[dict, dict]:
    pts, summ = {}, {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        q = line.split(",")
        if q[0] == "P" and len(q) >= 9:
            pts.setdefault(q[1], []).append([float(x) for x in q[2:]])
        elif q[0] == "S" and len(q) >= 4:
            summ.setdefault(q[1], {})[q[2]] = float(q[3])
    return pts, summ


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()

    bin_dir = Path(args.bin).resolve()
    cases = Path(args.cases).resolve()
    ref_dir = Path(args.reference).resolve()
    exe = bin_dir / "snap_probe.exe"
    snap = cases / SNAP
    frozen = cases / FROZEN
    ref = json.loads((ref_dir / REF_JSON).read_text(encoding="utf-8"))

    print("=" * 78)
    print("第 28 项 CTest：位移控制 / 线搜索（B3 第二十一次力学①）")
    print("=" * 78)
    check("snap_probe.exe 与冻结输入在档", exe.is_file() and snap.is_file(),
          f"exe={exe.is_file()} snap={snap.is_file()}")

    # 冻结 CSV sha256 与 reference 一致（ guards drift）
    check("冻结 FEM CSV sha256 = reference",
          sha256(frozen) == ref["frozen_sha256"])

    # 用当前二进制重跑
    with tempfile.TemporaryDirectory(prefix="dispctl_") as td:
        work = Path(td)
        shutil.copy2(snap, work / SNAP)
        fresh_csv = work / "fresh.csv"
        p = subprocess.run([str(exe), SNAP, fresh_csv.name], cwd=str(work),
                           env=env_with_msys(), capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=900)
        check("snap_probe 退出码 0", p.returncode == 0, f"rc={p.returncode}")
        if not fresh_csv.is_file():
            check("重跑 CSV 生成", False, "missing")
            return summary()

        fresh_pts, fresh_summ = read_csv(fresh_csv)
        fresh_bytes = fresh_csv.read_bytes()

    frozen_pts, frozen_summ = read_csv(frozen)

    # 1) 逐位复现冻结 CSV（字节级）
    same_bytes = fresh_bytes == frozen.read_bytes()
    check("当前二进制重跑 = 冻结 CSV（逐字节）", same_bytes)

    # 2) 冻结值逐值一致（结构哈希已覆盖，这里逐点再核一次数值）
    n_vals, worst_repro = 0, 0.0
    for cid in CASES_IDS:
        fp, zp = fresh_pts.get(cid, []), frozen_pts.get(cid, [])
        if len(fp) != len(zp):
            check(f"{cid} 点数 {len(fp)} vs {len(zp)}", False)
            continue
        for a, b in zip(fp, zp):
            for x, y in zip(a, b):
                n_vals += 1
                worst_repro = max(worst_repro, abs(x - y))
    check(f"重跑 vs 冻结逐点数值一致（{n_vals} 值）",
          worst_repro <= TOL_REPRO, f"worst={worst_repro:.3e}")
    for cid in CASES_IDS:
        for k, v in ref["per_case"][cid].items():
            if k == "points":
                continue
            got = fresh_summ.get(cid, {}).get(k)
            if isinstance(v, float):
                ok = got is not None and abs(got - v) <= 1e-12 * max(1.0, abs(v))
            else:
                ok = got == v
            if not ok:
                check(f"{cid} summary {k} vs reference", False, f"{got} vs {v}")

    # 3) 三条正常路径：收敛 / 闭合 / 闭式 P(v) 逐点锚点
    th = ref["thresholds"]
    ana = analytic()
    for cid in NORMAL:
        pc = ref["per_case"][cid]
        s = fresh_summ.get(cid, {})
        check(f"{cid} 收敛且无错误", int(s.get("converged", 0)) == 1
              and int(s.get("error_empty", 0)) == 1)
        check(f"{cid} 走到目标 v=1", abs(s.get("v_final", 0) - 1.0) < 1e-9,
              f"v_final={s.get('v_final')}")
        worst_dp = 0.0
        for row in fresh_pts.get(cid, []):
            lam, v = row[1], row[2]   # P 行: step, lambda, v, residual, iters, ds, pivot
            if 0.0 <= v <= 2.0 * H:
                worst_dp = max(worst_dp, abs(lam * PREF - p_of_v(v)))
        check(f"{cid} 闭式 P(v) 逐点锚点",
              worst_dp <= th["tol_dP_abs"], f"worst_dP={worst_dp:.3e} N")
        check(f"{cid} worst_dP_rel 达标",
              s.get("worst_dP_rel", 1.0) <= th["tol_dP_rel"],
              f"rel={s.get('worst_dP_rel')}")

    # 4) line search 与无 line search 同路径 + 触发断言
    pl = fresh_pts.get("dc_plain", [])
    ll = fresh_pts.get("dc_ls", [])
    if len(pl) == len(ll) and pl:
        worst_eq = max(abs(a[1] - b[1]) for a, b in zip(pl, ll))
        check("line search 与无 line search 逐点同路径",
              worst_eq <= th["tol_ls_eq"], f"worst|Δλ|={worst_eq:.3e}")
    else:
        check("line search 路径点数匹配", False)
    check("无 line search：0 次回溯 / alpha=1",
          fresh_summ["dc_plain"]["ls_total_cuts"] == 0
          and abs(fresh_summ["dc_plain"]["ls_min_alpha"] - 1.0) < 1e-12)
    check("有 line search：确有回溯被触发",
          fresh_summ["dc_ls"]["ls_total_cuts"] > 0
          and fresh_summ["dc_ls"]["ls_steps_cut"] > 0,
          f"cuts={fresh_summ['dc_ls']['ls_total_cuts']}")

    # 5) 细步 λ_max 更接近闭式 λ_cr
    coarse = fresh_summ["dc_plain"]["lambda_max"]
    fine = fresh_summ["dc_ls_fine"]["lambda_max"]
    lcr = ana["lambda_cr"]
    check("粗步 λ_max 在采样误差内接近 λ_cr",
          abs(coarse - lcr) <= th["tol_lcr_coarse"],
          f"|Δ|={abs(coarse-lcr):.3e}")
    check("细步 λ_max 在采样误差内接近 λ_cr",
          abs(fine - lcr) <= th["tol_lcr_fine"],
          f"|Δ|={abs(fine-lcr):.3e}")
    check("细步比粗步更接近 λ_cr", abs(fine - lcr) < abs(coarse - lcr))

    # 6) 负例：控制 DOF 指向固定端 → 报错、不收敛、无路径推进
    bad = fresh_summ.get("dc_badctrl", {})
    check("负例 singular_control_dof：error 置位 / 不收敛",
          int(bad.get("error_empty", 1)) == 0
          and int(bad.get("converged", 1)) == 0,
          f"error_empty={bad.get('error_empty')} conv={bad.get('converged')}")
    check("负例无路径推进（npoints=1）",
          int(bad.get("npoints", -1)) == 1)

    return summary()


def summary() -> int:
    n = sum(1 for ok, _, _ in RESULTS if ok)
    print("\n" + "=" * 78)
    print(f"  {n}/{len(RESULTS)} checks passed")
    print("  [PASS] CB3D25 displacement control + line search benchmark"
          if n == len(RESULTS) else
          "  [FAIL] CB3D25 displacement control benchmark")
    print("=" * 78)
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
