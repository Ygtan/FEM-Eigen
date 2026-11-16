#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 29 项 CTest：第二十一次力学②初始缺陷形态 / 幅值 / 敏感性。

检查：
  1. 当前 snap_probe.exe 重跑冻结 .snap → 与冻结 CSV 逐字节 / 逐值一致；
  2. 对称模态幅值扫描 a0→a3 的 λ_max_fit 与**独立闭式** λ_cr 对齐（细步二次插值）；
  3. 闭式 P(v) 逐点锚点（缺陷后几何，不读求解器）；
  4. λ_max 随幅值单调下降；
  5. 反对称模态的对称破缺（顶点 X 偏移、Y 不偏移）；
  6. IMPMODE 越界负例（imp_ok=0，不静默）；
  7. IMPERF=0 默认零污染：关闭路径与零幅值路径逐位一致、且不输出诊断键。

闭式独立实现：H'=H−a, L0'=hypot(B,H'), P(v)=−2EA(L−L0')/L0'·(H'−v)/L，
黄金分割求极值。模型参数从冻结 .snap 的 G 行读取，不读 reference JSON。
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

SWEEP = ["imp_a0", "imp_a1", "imp_a2", "imp_a3"]
ALL_IDS = ["imp_off", "imp_a0", "imp_a1", "imp_a2", "imp_a3",
           "imp_anti", "imp_badmode"]


def parse_snap(path: Path) -> tuple[dict, dict]:
    groups, cases = {}, {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line[0] in "#!":
            continue
        t = line.split()
        kv = {}
        for x in t[2:]:
            if "=" in x:
                k, v = x.split("=", 1)
                try:
                    kv[k] = float(v)
                except ValueError:
                    kv[k] = v
        if t[0] == "G":
            groups[t[1]] = kv
        elif t[0] == "C":
            c = dict(kv)
            c["gid"] = t[2]          # C <cid> <gid> KEY=...
            cases[t[1]] = c
    return groups, cases


def read_csv(path: Path) -> tuple[dict, dict]:
    pts, summ = {}, {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        q = line.split(",")
        if q[0] == "P" and len(q) >= 9:
            pts.setdefault(q[1], []).append([float(x) for x in q[2:]])
        elif q[0] == "S" and len(q) >= 4:
            summ.setdefault(q[1], {})[q[2]] = float(q[3])
    return pts, summ


def closed_lcr(g: dict, amp: float, strain: int = 0) -> tuple[float, float]:
    E, A, B, H, PREF = g["E"], g["A"], g["B"], g["H"], g["PREF"]
    h = H - amp
    l0 = math.hypot(B, h)

    def p(v: float) -> float:
        L = math.hypot(B, h - v)
        N = E * A * (L - l0) / l0
        return -2.0 * N * (h - v) / L

    lo, hi = 0.0, h
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.3819660112501051
        m2 = lo + (hi - lo) * 0.6180339887498949
        if p(m1) < p(m2):
            lo = m1
        else:
            hi = m2
    vcr = 0.5 * (lo + hi)
    return vcr, p(vcr) / PREF


def p_of_v(g: dict, amp: float, v: float) -> float:
    E, A, B, H = g["E"], g["A"], g["B"], g["H"]
    h = H - amp
    l0 = math.hypot(B, h)
    L = math.hypot(B, h - v)
    return -2.0 * (E * A * (L - l0) / l0) * (h - v) / L


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()

    bin_dir = Path(args.bin)
    cases_dir = Path(args.cases)
    ref_dir = Path(args.reference)
    exe = bin_dir / "snap_probe.exe"
    snap = cases_dir / "cb3d26_initial_imperfection.snap"
    frozen = cases_dir / "cb3d26_initial_imperfection_fem.csv"
    ref = json.loads((ref_dir / "cb3d26_initial_imperfection.json").read_text(encoding="utf-8"))

    ok_all = True
    passed = 0
    total = 0

    def ck(cond: bool, name: str, detail: str = "") -> None:
        nonlocal ok_all, passed, total
        total += 1
        if cond:
            passed += 1
            print(f"[PASS] {name}")
        else:
            ok_all = False
            print(f"[FAIL] {name}: {detail}")

    groups, cases = parse_snap(snap)
    g = groups[cases[ALL_IDS[0]]["gid"]]
    th = ref["thresholds"]

    # 当前 exe 重跑冻结输入（在上下文内读完所有内容）
    with tempfile.TemporaryDirectory() as td:
        fresh = Path(td) / "fresh.csv"
        proc = subprocess.run([str(exe), str(snap), str(fresh)],
                              capture_output=True, text=True)
        ck(proc.returncode == 0, "snap_probe 重跑冻结输入", proc.stderr[:200])
        fresh_bytes = fresh.read_bytes() if fresh.exists() else b""
    frozen_bytes = frozen.read_bytes()
    ck(fresh_bytes == frozen_bytes, "冻结 CSV 逐字节复现",
       f"fresh={len(fresh_bytes)} frozen={len(frozen_bytes)}")

    fresh_pts, fresh_summ = read_csv(frozen)  # 逐字节一致 → 同一内容

    # 逐 case 基础健全性
    for cid in ALL_IDS:
        s = fresh_summ.get(cid, {})
        ck(s.get("converged", 0) == 1.0 and s.get("error_empty", 0) == 1.0,
           f"{cid} 收敛且无错误", f"conv={s.get('converged')} err={s.get('error_empty')}")
        ck(len(fresh_pts.get(cid, [])) > 1, f"{cid} 路径非空",
           f"n={len(fresh_pts.get(cid, []))}")

    # 闭式锚点：对称模态幅值扫描
    amps = []
    lcrs = []
    for cid in SWEEP:
        s = fresh_summ[cid]
        amp = -s["imp_dy2"]           # 对称模态顶点下移量 = 幅值
        amps.append(amp)
        ck(abs(s["imp_dx2"]) <= th["tol_sym_dx"], f"{cid} 对称模态无横向偏移",
           f"dx2={s['imp_dx2']}")
        vcr, lcr = closed_lcr(g, amp)
        lcrs.append(lcr)
        lfit = s["lambda_max_fit"]
        ck(abs(lfit - lcr) <= th["tol_lcr_fit"],
           f"{cid} λ_max_fit 对齐闭式 λ_cr",
           f"fit={lfit:.9f} closed={lcr:.9f} diff={abs(lfit-lcr):.2e}")
        # 闭式 P(v) 逐点锚点
        worst = 0.0
        pcr = p_of_v(g, amp, vcr)          # 向下极值荷载为正
        for row in fresh_pts[cid]:
            lam, v = row[1], row[2]        # P 行: step, lambda, v, ...
            if 0.0 <= v <= 2.0 * (g["H"] - amp):
                worst = max(worst, abs(lam * g["PREF"] - p_of_v(g, amp, v)))
        rel = worst / pcr if pcr > 0 else 0.0
        ck(worst <= th["tol_dP_abs"] and rel <= th["tol_dP_rel"],
           f"{cid} 闭式 P(v) 逐点锚点", f"worst={worst:.3e} rel={rel:.3e}")

    # 单调性：幅值越大 λ_max 越低
    mono = all(lcrs[i] > lcrs[i + 1] for i in range(len(lcrs) - 1))
    ck(mono, "幅值→λ_max 严格单调下降",
       " > ".join(f"{x:.6f}" for x in lcrs) + f" (amps={amps})")
    lfits = [fresh_summ[c]["lambda_max_fit"] for c in SWEEP]
    ck(all(lfits[i] > lfits[i + 1] for i in range(len(lfits) - 1)),
       "FEM λ_max_fit 同步单调下降", " > ".join(f"{x:.6f}" for x in lfits))

    # 对称破缺：反对称模态顶点 X 偏移、Y 不偏移
    anti = fresh_summ["imp_anti"]
    anti_amp = cases["imp_anti"]["IMPAMP"]
    ck(abs(anti["imp_dx2"] + anti_amp) <= 1.0e-12,
       "反对称模态顶点 X 偏移=幅值", f"dx2={anti['imp_dx2']} amp={anti_amp}")
    ck(abs(anti["imp_dy2"]) <= th["tol_anti_dy"],
       "反对称模态顶点 Y 不偏移", f"dy2={anti['imp_dy2']}")

    # 负例：IMPMODE 越界 → imp_ok=0，不静默
    ck(fresh_summ["imp_badmode"]["imp_ok"] == 0.0,
       "IMPMODE 越界负例 imp_ok=0", f"imp_ok={fresh_summ['imp_badmode']['imp_ok']}")

    # 默认零污染：IMPERF=0 不输出诊断键；与零幅值路径逐位一致
    ck("imperf" not in fresh_summ["imp_off"],
       "IMPERF=0 不输出缺陷诊断键")
    off_pts = fresh_pts["imp_off"]
    a0_pts = fresh_pts["imp_a0"]
    ck(off_pts == a0_pts, "IMPERF=0 与零幅值路径逐位一致",
       f"off={len(off_pts)} a0={len(a0_pts)}")

    print(f"{'ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED'} "
          f"({passed}/{total})")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
