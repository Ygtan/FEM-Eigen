#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 29 项 CTest 黄金值：第二十一次力学②初始缺陷 / 缺陷敏感性。

数据源只有两类：
  1. cases/cb3d26_initial_imperfection_fem.csv —— 冻结的 snap_probe.exe
     （opt-in IMPERF/IMPMODE/IMPAMP）输出；
  2. von Mises 两杆浅桁架的**闭式解析** P(v) —— compare 脚本独立实现，
     不读求解器、不读本 JSON 的解析值（本文件只存冻结摘要与 sha256）。

闭式（缺陷后几何，无应力 L0' = hypot(B, H−a)）：
  h = H − a,  L(v) = sqrt(B² + (h−v)²),  N = EA(L − L0')/L0'
  P(v) = −2 N (h−v)/L,  极值 dP/dv = 0 由黄金分割求, λ_cr = P(v_cr)/PREF。

用法：
    python build_reference.py            # 生成 reference JSON
    python build_reference.py --check    # 校验冻结文件在档且 sha256 一致
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
SNAP = CASES / "cb3d26_initial_imperfection.snap"
FEM_CSV = CASES / "cb3d26_initial_imperfection_fem.csv"
REF = HERE / "reference" / "cb3d26_initial_imperfection.json"

E, A, B, H, PREF = 2.1e11, 0.001, 2.0, 0.5, 2.0e6
CASE_IDS = ["imp_off", "imp_a0", "imp_a1", "imp_a2", "imp_a3",
            "imp_anti", "imp_badmode"]
SWEEP = ["imp_a0", "imp_a1", "imp_a2", "imp_a3"]

TOL_LCR_FIT = 2.0e-4      # λ_max_fit 与闭式 λ_cr 的差（细步二次插值）
TOL_DP_ABS = 1.0e-5       # 闭式 P(v) 逐点绝对误差上限（N）
TOL_DP_REL = 1.0e-11
TOL_SYM_DX = 1.0e-12      # 对称模态顶点 X 偏移上限
TOL_ANTI_DY = 1.0e-12     # 反对称模态顶点 Y 偏移上限


def closed_lcr(amp: float) -> dict:
    """独立闭式：缺陷后几何的极值荷载因子。"""
    h = H - amp
    l0 = math.hypot(B, h)

    def p_of_v(v: float) -> float:
        L = math.hypot(B, h - v)
        return -2.0 * (E * A * (L - l0) / l0) * (h - v) / L

    lo, hi = 0.0, h
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.3819660112501051
        m2 = lo + (hi - lo) * 0.6180339887498949
        if p_of_v(m1) < p_of_v(m2):
            lo = m1
        else:
            hi = m2
    vcr = 0.5 * (lo + hi)
    return {"amp": amp, "h_eff": h, "L0_eff": l0,
            "v_cr": vcr, "lambda_cr": p_of_v(vcr) / PREF}


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
    out = {
        "schema": "initial_imperfection_3d/v1",
        "mechanic": "B3 fifth-wave (mechanic 2): initial imperfection + sensitivity",
        "frozen_sha256": sha256(FEM_CSV),
        "snap_sha256": sha256(SNAP),
        "model": {"E": E, "A": A, "B": B, "H": H, "PREF": PREF},
        "thresholds": {
            "tol_lcr_fit": TOL_LCR_FIT, "tol_dP_abs": TOL_DP_ABS,
            "tol_dP_rel": TOL_DP_REL, "tol_sym_dx": TOL_SYM_DX,
            "tol_anti_dy": TOL_ANTI_DY,
        },
        "per_case": {},
    }
    for cid in CASE_IDS:
        s = summ.get(cid, {})
        out["per_case"][cid] = {
            "npoints": int(s.get("npoints", 0)),
            "steps": int(s.get("steps", 0)),
            "converged": int(s.get("converged", 0)),
            "error_empty": int(s.get("error_empty", 0)),
            "lambda_max": s.get("lambda_max", 0.0),
            "lambda_max_fit": s.get("lambda_max_fit", 0.0),
            "imp_ok": int(s.get("imp_ok", 0)),
            "imp_lambda": s.get("imp_lambda", 0.0),
            "imp_dx2": s.get("imp_dx2", 0.0),
            "imp_dy2": s.get("imp_dy2", 0.0),
            "points": [[float(r[0]), float(r[1])] for r in pts.get(cid, [])],
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    ref = build()
    if args.check:
        ok = (FEM_CSV.is_file() and SNAP.is_file() and REF.is_file() and
              json.loads(REF.read_text(encoding="utf-8")).get("frozen_sha256")
              == ref["frozen_sha256"])
        print(("[PASS]" if ok else "[FAIL]") + " initial imperfection frozen files")
        return 0 if ok else 1
    REF.parent.mkdir(parents=True, exist_ok=True)
    REF.write_text(json.dumps(ref, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[OK] wrote {REF.relative_to(HERE)}")
    for cid in CASE_IDS:
        pc = ref["per_case"][cid]
        print(f"  {cid:12s} conv={pc['converged']} imp_ok={pc['imp_ok']} "
              f"lam_fit={pc['lambda_max_fit']:.9f} dx2={pc['imp_dx2']:+.4f} "
              f"dy2={pc['imp_dy2']:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
