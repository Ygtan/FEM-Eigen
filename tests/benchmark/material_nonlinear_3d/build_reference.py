#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 30 项 CTest 黄金值：第二十一次力学③材料非线性本构接入 3D CR 梁。

数据源只有两类：
  1. cases/cb3d27_material_nl_fem.csv —— 冻结的 mat_nl_probe.exe 输出
     （section.hpp materialStress/sectionResponse + fem_3d_nl.hpp
      crFiberLocalForce/crFiberLocalTangent，Gauss 纤维离散）；
  2. 独立 Python 实现（本文件 + compare 脚本共用）的纤维截面积分 ——
     不调用任何 C++ 代码，用 numpy.polynomial.legendre.leggauss 生成
     与 C++ Newton 求根**独立**的 Gauss 点/权重，逐条重算 SECT/FBEAM
     期望值。

闭式/独立锚点：
  * 单轴应力-应变：ELASTIC N=E·A·ε；STEEL |Eε|≤fy 弹性否则 N=±fy·A；
    CONCRETE 受拉 Eε>ft 退出 N=0，受压线弹性。
  * 矩形截面塑性铰：高斯纤维离散逐点 materialStress 后积分 N/My/Mz
    与 3×3 截面切线 D（对 y²,z² 二阶精确，弹性段退化为 E·Iz/E·Iy）。
  * CR 局部力：q[0/6]=∓N，q[5/11]=∓/+Mz(端)，q[4/10]=∓/+My(端)，
    端剪力由无分布荷载平衡 q[1]=(q5+q11)/L、q[2]=−(q4+q10)/L，
    扭转 q[3]=GJ/L·(θx1−θx2) 保持线弹性。
  * 弹性逐位对齐：ELASTIC 组纤维力/切线与解析 Kl12·d 逐元素差 <1e-6 相对。

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

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
MAT = CASES / "cb3d27_material_nl.mat"
FEM_CSV = CASES / "cb3d27_material_nl_fem.csv"
REF = HERE / "reference" / "cb3d27_material_nl.json"

# 与 .mat 输入严格一致（build_reference 不读输入文本，靠这些常量独立复算）
E_ST = 2.1e11
E_CO = 3.0e10
A = 0.02
IY = 1.66666666666667e-5
IZ = 6.66666666666667e-5
G_ST = 8.07692307692308e10
G_CO = 1.15384615384615e10
J = 2.66666666666667e-4
FY = 2.35e8
FT = 2.0e6
NFY, NFZ = 8, 4
H_RECT = math.sqrt(12.0 * IZ / A)   # 0.2
B_RECT = math.sqrt(12.0 * IY / A)   # 0.1
LREF = 1.0

GROUPS = {
    "gel": dict(kind="ELASTIC", E=E_ST, G=G_ST),
    "gst": dict(kind="STEEL", E=E_ST, G=G_ST, fy=FY),
    "gco": dict(kind="CONCRETE", E=E_ST, G=G_CO, ft=FT),
}
# gco 在 .mat 里 E=3.0e10，单轴段用它；纤维积分统一看 kind。
GROUPS["gco"]["E"] = E_CO

# S 行：(sid, gid, eps, kz, ky)
SECT_CASES = [
    ("el_t1", "gel", 1e-4, 0.0, 0.0),
    ("el_c1", "gel", -5e-5, 0.0, 0.0),
    ("st_e1", "gst", 5e-4, 0.0, 0.0),
    ("st_e2", "gst", 1.5e-3, 0.0, 0.0),
    ("st_c2", "gst", -1.5e-3, 0.0, 0.0),
    ("co_t1", "gco", 3e-5, 0.0, 0.0),
    ("co_t2", "gco", 1e-4, 0.0, 0.0),
    ("co_c1", "gco", -1e-4, 0.0, 0.0),
    ("el_k1", "gel", 0.0, 1e-3, 0.0),
    ("st_k1", "gst", 0.0, 5e-3, 0.0),
    ("st_k2", "gst", 0.0, 2e-2, 0.0),
    ("st_k3", "gst", 0.0, 5e-2, 0.0),
    ("el_ky", "gel", 0.0, 0.0, 2e-3),
    ("st_ky", "gst", 0.0, 0.0, 2e-2),
]

# F 行：(fid, gid, 12 维 deff)
FIBER_CASES = [
    ("f_el1", "gel", [0, 0, 0, 5e-4, 1e-3, 1e-3, 1e-3, 0, 0, -5e-4, -1e-3, 2e-3]),
    ("f_el2", "gel", [0, 0, 0, 0, 0, 2e-3, 0, 0, 0, 0, 0, -2e-3]),
    ("f_st1", "gst", [0, 0, 0, 0, 2e-3, 2e-3, 0, 0, 0, 0, -2e-3, -2e-3]),
    ("f_st2", "gst", [0, 0, 0, 0, 1e-2, 1e-2, 0, 0, 0, 0, -1e-2, -1e-2]),
    ("f_co1", "gco", [0, 0, 0, 0, 0, 0, 1e-3, 0, 0, 0, 0, 0]),
]

TOL_REL = 1.0e-6        # 截面/纤维值相对误差（Gauss 离散与 C++ 同源但独立实现）
TOL_ABS = 1.0e-6        # 近零值绝对误差（力 N）
TOL_EQ = 1.0e-6         # 弹性段纤维 vs Kl12 相对误差


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gauss(n: int):
    """独立生成 [-1,1] 上 n 点 Gauss 点/权重（numpy，与 C++ Newton 求根无关）。"""
    x, w = leggauss(n)
    return x, w


def material_stress(g: dict, eps: float):
    """与 section.hpp materialStress 同语义的独立实现：返回 (应力, 切线)。"""
    E = g["E"]
    if g["kind"] == "STEEL":
        s = E * eps
        if abs(s) <= g["fy"]:
            return s, E
        return math.copysign(g["fy"], eps), 0.0
    if g["kind"] == "CONCRETE":
        if eps >= 0.0 and g["ft"] > 0.0:
            s = E * eps
            if s <= g["ft"]:
                return s, E
            return 0.0, 0.0
        return E * eps, E
    return E * eps, E


def fiber_grid():
    gy, wy = gauss(NFY)
    gz, wz = gauss(NFZ)
    halfA = 0.25 * H_RECT * B_RECT
    fibers = []
    for i, (gx, wx) in enumerate(zip(gy, wy)):
        y = 0.5 * H_RECT * gx
        for j, (gzj, wzj) in enumerate(zip(gz, wz)):
            z = 0.5 * B_RECT * gzj
            fibers.append((y, z, halfA * wx * wzj))
    return fibers


def section_response(g: dict, eps0: float, kz: float, ky: float) -> dict:
    """独立截面积分：逐纤维 materialStress 后积分 N/My/Mz 与 3×3 切线。"""
    N = My = Mz = 0.0
    D11 = D12 = D13 = D22 = D23 = D33 = 0.0
    for y, z, area in fiber_grid():
        eps = eps0 - kz * y + ky * z
        sig, Et = material_stress(g, eps)
        N += sig * area
        Mz += -sig * y * area
        My += sig * z * area
        D11 += Et * area
        D12 += -Et * y * area
        D13 += Et * z * area
        D22 += Et * y * y * area
        D23 += -Et * y * z * area
        D33 += Et * z * z * area
    return dict(N=N, My=My, Mz=Mz, D11=D11, D12=D12, D13=D13,
                D22=D22, D23=D23, D33=D33)


def fiber_local_force(g: dict, d: list) -> list:
    """独立 CR 局部力（与 crFiberLocalForce 同语义）。"""
    L = LREF
    GJ = g["G"] * J
    eps0 = d[6] / L
    kz0 = -(4.0 * d[5] + 2.0 * d[11]) / L
    kz1 = (2.0 * d[5] + 4.0 * d[11]) / L
    ky0 = -(4.0 * d[4] + 2.0 * d[10]) / L
    ky1 = (2.0 * d[4] + 4.0 * d[10]) / L
    r0 = section_response(g, eps0, kz0, ky0)
    r1 = section_response(g, eps0, kz1, ky1)
    q = [0.0] * 12
    q[0] = -r0["N"]; q[6] = r0["N"]
    q[5] = -r0["Mz"]; q[11] = r1["Mz"]
    q[4] = -r0["My"]; q[10] = r1["My"]
    q[1] = (q[5] + q[11]) / L; q[7] = -q[1]
    q[2] = -(q[4] + q[10]) / L; q[8] = -q[2]
    q[3] = GJ / L * (d[3] - d[9]); q[9] = -q[3]
    return q


def kl12_local(g: dict) -> np.ndarray:
    """独立解析 12×12 梁局部刚度（Lref=1，与 beam3DLocalStiffness 同式）。"""
    E = g["E"]
    K = np.zeros((12, 12))
    EAL = E * A / LREF
    kz = E * IZ / LREF ** 3
    ky = E * IY / LREF ** 3
    GJL = g["G"] * J / LREF
    K[0, 0] = EAL; K[0, 6] = -EAL; K[6, 0] = -EAL; K[6, 6] = EAL
    K[1, 1] = 12*kz; K[1, 5] = 6*kz; K[1, 7] = -12*kz; K[1, 11] = 6*kz
    K[5, 1] = 6*kz; K[5, 5] = 4*kz; K[5, 7] = -6*kz; K[5, 11] = 2*kz
    K[7, 1] = -12*kz; K[7, 5] = -6*kz; K[7, 7] = 12*kz; K[7, 11] = -6*kz
    K[11, 1] = 6*kz; K[11, 5] = 2*kz; K[11, 7] = -6*kz; K[11, 11] = 4*kz
    K[2, 2] = 12*ky; K[2, 4] = -6*ky; K[2, 8] = -12*ky; K[2, 10] = -6*ky
    K[4, 2] = -6*ky; K[4, 4] = 4*ky; K[4, 8] = 6*ky; K[4, 10] = 2*ky
    K[8, 2] = -12*ky; K[8, 4] = 6*ky; K[8, 8] = 12*ky; K[8, 10] = 6*ky
    K[10, 2] = -6*ky; K[10, 4] = 2*ky; K[10, 8] = 6*ky; K[10, 10] = 4*ky
    K[3, 3] = GJL; K[3, 9] = -GJL; K[9, 3] = -GJL; K[9, 9] = GJL
    return K


IDX = [3, 4, 5, 6, 9, 10, 11]


def tang_sub(g: dict, d: list) -> list:
    """独立中心差分切线（材料敏感子块，行/列 idx）。"""
    q0 = fiber_local_force(g, d)
    K = np.zeros((12, 12))
    for j in range(12):
        h = 1e-7 * (abs(d[j]) + 1e-3)
        dp = list(d); dp[j] += h
        dm = list(d); dm[j] -= h
        qp = fiber_local_force(g, dp)
        qm = fiber_local_force(g, dm)
        K[:, j] = [(qp[i] - qm[i]) / (2*h) for i in range(12)]
    K = 0.5 * (K + K.T)
    return [[K[i, j] for j in IDX] for i in IDX]


def expected() -> dict:
    out = {"sect": {}, "fbeam": {}, "ftang": {}, "feq": {}}
    for sid, gid, eps, kz, ky in SECT_CASES:
        out["sect"][sid] = section_response(GROUPS[gid], eps, kz, ky)
    for fid, gid, d in FIBER_CASES:
        out["fbeam"][fid] = fiber_local_force(GROUPS[gid], d)
        out["ftang"][fid] = tang_sub(GROUPS[gid], d)
        if GROUPS[gid]["kind"] == "ELASTIC":
            q = np.array(out["fbeam"][fid])
            qe = kl12_local(GROUPS[gid]) @ np.array(d)
            Kt = np.array(out["ftang"][fid])
            Kr = kl12_local(GROUPS[gid])[np.ix_(IDX, IDX)]
            out["feq"][fid] = {
                "maxdq": float(np.max(np.abs(q - qe))),
                "maxdK": float(np.max(np.abs(Kt - Kr))),
                "qscale": float(np.max(np.abs(qe))),
                "Kscale": float(np.max(np.abs(Kr))),
            }
    return out


def build() -> dict:
    exp = expected()
    return {
        "schema": "material_nonlinear_3d/v1",
        "mechanic": "B3 fifth-wave (mechanic 3): material nonlinearity in 3D CR beam",
        "frozen_sha256": sha256(FEM_CSV),
        "input_sha256": sha256(MAT),
        "model": {"E_STEEL": E_ST, "E_CONCRETE": E_CO, "A": A,
                  "IY": IY, "IZ": IZ, "J": J, "FY": FY, "FT": FT,
                  "h_rect": H_RECT, "b_rect": B_RECT, "nfy": NFY,
                  "nfz": NFZ, "Lref": LREF},
        "thresholds": {"tol_rel": TOL_REL, "tol_abs": TOL_ABS,
                       "tol_eq": TOL_EQ},
        "expected": exp,
    }


def check() -> None:
    if not FEM_CSV.exists():
        raise SystemExit(f"缺冻结文件: {FEM_CSV}")
    if not MAT.exists():
        raise SystemExit(f"缺输入文件: {MAT}")
    ref = build()
    got_csv = sha256(FEM_CSV)
    if got_csv != ref["frozen_sha256"]:
        raise SystemExit(f"冻结 CSV sha256 不符: {got_csv}")
    print(f"[OK] 冻结文件在档，sha256={got_csv[:16]}…")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    if args.check:
        check()
    else:
        REF.parent.mkdir(parents=True, exist_ok=True)
        REF.write_text(json.dumps(build(), ensure_ascii=False, indent=2),
                       encoding="utf-8")
        print(f"已生成 {REF}")
