#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成力学⑤动力 snap-through benchmark 的冻结参考 JSON。

参考值全部由**独立闭式 / 独立积分**产生，不读取 FEM 的 T 行位移/速度：
  - von Mises 浅拱静力极限 P(v)、v_cr：黄金分割求极值；
  - 等效质量：一致质量 m=2·rhoA·L0/3，集中质量 m=rhoA·L0（闭式）；
  - 时程锚点：单自由度 ODE  m·v'' + F_spring(v) = P(t)  的 RK4 积分，
    步长 1e-5 s（FEM 的 1/5），每 20 个 FEM 步取一个样本冻结；
  - FEM CSV / 输入文件 sha256 冻结，供 compare 做逐字节重放。
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

B = 2.0
H = 0.5
E = 2.1e11
A = 1.0e-3
DENS = 7850.0
PREF_SNAP = 1010157.0
PREF_LOW = 713052.0
DT_FEM = 5.0e-5
TEND = 0.02
SAMPLE = 20          # 每 20 个 FEM 步冻结一个 RK4 样本
RK4_DT = DT_FEM / 5.0


def p_of_v(v):
    L0 = math.hypot(B, H)
    L = math.hypot(B, H - v)
    return 2.0 * E * A * (L0 - L) * (H - v) / (L0 * L)


def limit_v():
    lo, hi = 0.0, H
    r1 = 0.3819660112501051
    r2 = 0.6180339887498949
    for _ in range(200):
        m1 = lo + (hi - lo) * r1
        m2 = lo + (hi - lo) * r2
        if p_of_v(m1) < p_of_v(m2):
            lo = m1
        else:
            hi = m2
    return 0.5 * (lo + hi)


def restoring(v):
    L0 = math.hypot(B, H)
    L = math.hypot(B, H - v)
    return 2.0 * E * A * (L0 - L) * (H - v) / (L0 * L)


def rk4_step_load(m, P, nsteps):
    """阶跃荷载 SDOF：m·v'' + F(v) = P, v(0)=v'(0)=0。返回 {step:(v,vel)}。"""
    out = {0: (0.0, 0.0)}
    v, vel = 0.0, 0.0
    dt = RK4_DT
    sub = 5  # 每个 FEM 步内 5 个 RK4 子步

    def acc(vv):
        return (P - restoring(vv)) / m

    for i in range(1, nsteps + 1):
        for _ in range(sub):
            k1v = vel
            k1a = acc(v)
            k2v = vel + 0.5 * dt * k1a
            k2a = acc(v + 0.5 * dt * k1v)
            k3v = vel + 0.5 * dt * k2a
            k3a = acc(v + 0.5 * dt * k2v)
            k4v = vel + dt * k3a
            k4a = acc(v + dt * k3v)
            v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
            vel += dt / 6.0 * (k1a + 2.0 * k2a + 2.0 * k3a + k4a)
        out[i] = (v, vel)
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()
    cases = Path(args.cases)
    refdir = Path(args.reference)
    refdir.mkdir(parents=True, exist_ok=True)

    inp = cases / "cb3d29_dynamic_snapthrough.dyn"
    fem = cases / "cb3d29_dynamic_snapthrough_fem.csv"

    L0 = math.hypot(B, H)
    vcr = limit_v()
    pcr = p_of_v(vcr)
    rhoA = DENS * A
    m_cons = 2.0 * rhoA * L0 / 3.0
    m_lump = rhoA * L0
    nsteps = int(round(TEND / DT_FEM))

    n_snap = rk4_step_load(m_cons, PREF_SNAP, nsteps)
    n_lump = rk4_step_load(m_lump, PREF_SNAP, nsteps)
    n_low = rk4_step_load(m_cons, PREF_LOW, nsteps)

    def sample(series):
        return [
            [round(i * DT_FEM, 10), series[i][0], -series[i][1]]
            for i in range(0, nsteps + 1, SAMPLE)
        ]

    ref = {
        "meta": {
            "case": "cb3d29_dynamic_snapthrough",
            "mechanic": "B3-5 dynamic snap-through (Newmark-beta, opt-in)",
            "integrator": "independent SDOF RK4, substep=5, dt=1e-5 s",
            "note": "v 向下为正；vel 为 FEM 的 Y 向速度（向下为负）",
        },
        "sha256": {
            "input": sha256(inp),
            "fem_csv": sha256(fem),
        },
        "analytic": {
            "B": B, "H": H, "E": E, "A": A, "dens": DENS,
            "L0": L0, "rhoA": rhoA,
            "v_cr": vcr, "P_cr": pcr,
            "PREF_snap": PREF_SNAP, "PREF_low": PREF_LOW,
            "ratio_snap": PREF_SNAP / pcr,
            "ratio_low": PREF_LOW / pcr,
            "m_consistent": m_cons, "m_lumped": m_lump,
            "dt": DT_FEM, "t_end": TEND, "steps": nsteps,
            "static_snap_v_2nd_root": 2.0 * H,
        },
        "rk4": {
            "dyn_snap": sample(n_snap),
            "dyn_lump": sample(n_lump),
            "dyn_low": sample(n_low),
        },
    }
    out = refdir / "cb3d29_dynamic_snapthrough.json"
    out.write_text(json.dumps(ref, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print(f"v_cr={vcr:.10f}  P_cr={pcr:.6f}  P_snap/P_cr={PREF_SNAP/pcr:.6f}"
          f"  P_low/P_cr={PREF_LOW/pcr:.6f}")
    print(f"m_consistent={m_cons:.6f}  m_lumped={m_lump:.6f}")
    print(f"RK4 samples per case: {len(ref['rk4']['dyn_snap'])}")


if __name__ == "__main__":
    main()
