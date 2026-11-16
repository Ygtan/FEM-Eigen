#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 28 项 CTest 黄金值：第二十一次力学①位移控制 / 线搜索。

数据源只有两类：
  1. cases/cb3d25_displacement_control_fem.csv —— 用 mingw64 工具链冻结的
     snap_probe.exe（opt-in CTRLMODE/LINESEARCH）输出；
  2. von Mises 两杆浅桁架的**闭式解析** P(v)（独立 oracle，不读求解器）。

闭式：L(v)=sqrt(B^2+(H-v)^2), N=EA(L-L0)/L0, P(v)=-2N(H-v)/L；
极限点由黄金分割求 P 最负处，lambda_cr = P_cr/PREF。

用法：
    python build_reference.py            # 生成 reference JSON
    python build_reference.py --check    # 校验冻结文件在档且 JSON 一致
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
FEM_CSV = CASES / "cb3d25_displacement_control_fem.csv"
REF = HERE / "reference" / "cb3d25_displacement_control.json"

# von Mises 两杆（与 snap_probe.cpp / CB3D17 同源）
E, A, B, H, PREF = 2.1e11, 0.001, 2.0, 0.5, 2.0e6
L0 = math.hypot(B, H)

CASES_IDS = ["dc_plain", "dc_ls", "dc_ls_fine", "dc_badctrl"]
NORMAL = ["dc_plain", "dc_ls", "dc_ls_fine"]

# 判据阈值
TOL_DP_ABS = 1.0e-6        # 闭式 P(v) 逐点绝对误差上限（N）
TOL_DP_REL = 1.0e-12        # worst_dP_rel 上限
TOL_LS_EQ = 1.0e-10         # line search 与无 line search 的逐点 |Δλ|
TOL_LCR_COARSE = 5.0e-4     # 粗步 λ_max 与闭式 λ_cr 差（采样误差）
TOL_LCR_FINE = 1.0e-4       # 细步 λ_max 与闭式 λ_cr 差


def p_of_v(v: float) -> float:
    L = math.hypot(B, H - v)
    N = E * A * (L - L0) / L0
    return -2.0 * N * (H - v) / L


def analytic() -> dict:
    lo, hi = 0.0, H
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.3819660112501051
        m2 = lo + (hi - lo) * 0.6180339887498949
        if p_of_v(m1) < p_of_v(m2):
            lo = m1
        else:
            hi = m2
    vcr = 0.5 * (lo + hi)
    pcr = p_of_v(vcr)
    return {"v_cr": vcr, "P_cr": pcr, "lambda_cr": pcr / PREF}


def read_csv(path: Path) -> tuple[dict, dict]:
    pts, summ = {}, {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        q = line.split(",")
        if q[0] == "P" and len(q) >= 9:
            pts.setdefault(q[1], []).append([float(x) for x in q[2:]])
        elif q[0] == "S" and len(q) >= 4:
            summ.setdefault(q[1], {})[q[2]] = float(q[3])
    return pts, summ


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build() -> dict:
    pts, summ = read_csv(FEM_CSV)
    ana = analytic()
    out = {
        "schema": "displacement_control_3d/v1",
        "mechanic": "B3 fifth-wave (mechanic 1): displacement control + line search",
        "frozen_sha256": sha256(FEM_CSV),
        "analytic": ana,
        "thresholds": {
            "tol_dP_abs": TOL_DP_ABS, "tol_dP_rel": TOL_DP_REL,
            "tol_ls_eq": TOL_LS_EQ,
            "tol_lcr_coarse": TOL_LCR_COARSE, "tol_lcr_fine": TOL_LCR_FINE,
        },
        "per_case": {},
    }
    for cid in CASES_IDS:
        s = summ.get(cid, {})
        cpoints = pts.get(cid, [])
        out["per_case"][cid] = {
            "npoints": int(s.get("npoints", 0)),
            "steps": int(s.get("steps", 0)),
            "converged": int(s.get("converged", 0)),
            "error_empty": int(s.get("error_empty", 0)),
            "v_final": s.get("v_final", 0.0),
            "lambda_max": s.get("lambda_max", 0.0),
            "worst_dP": s.get("worst_dP", 0.0),
            "ctrlmode": int(s.get("ctrlmode", 0)),
            "linesearch": int(s.get("linesearch", 0)),
            "ctrlstep": s.get("ctrlstep", 0.0),
            "ls_total_cuts": int(s.get("ls_total_cuts", 0)),
            "ls_steps_cut": int(s.get("ls_steps_cut", 0)),
            "ls_min_alpha": s.get("ls_min_alpha", 1.0),
            "points": [[float(r[0]), float(r[1])] for r in cpoints],  # (lambda, v)
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    ref = build()
    if args.check:
        if not FEM_CSV.is_file() or not REF.is_file():
            print("[FAIL] frozen CSV / reference JSON missing")
            return 1
        saved = json.loads(REF.read_text(encoding="utf-8"))
        ok = saved.get("frozen_sha256") == ref["frozen_sha256"]
        print(("[PASS]" if ok else "[FAIL]") + " displacement control frozen CSV sha256")
        return 0 if ok else 1
    REF.parent.mkdir(parents=True, exist_ok=True)
    REF.write_text(json.dumps(ref, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[OK] wrote {REF.relative_to(HERE)}")
    print(f"     analytic lambda_cr = {ref['analytic']['lambda_cr']:.10f}")
    for cid in CASES_IDS:
        pc = ref["per_case"][cid]
        print(f"     {cid:12s} conv={pc['converged']} errEmpty={pc['error_empty']} "
              f"pts={pc['npoints']} ls_cuts={pc['ls_total_cuts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
