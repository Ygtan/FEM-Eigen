#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 37 项 benchmark 的独立闭式 oracle 与冻结脚本（B5 第二十六次）。

力学口径（3D 欧拉-伯努利梁，局部 x 轴沿单元）：
  * 无应力曲率 κ0 的固端力（局部）：f_eq(θ₁) = −EI·κ0，f_eq(θ₂) = +EI·κ0；
    经 Tᵀ 转全局后并入 elementEquivalentLoad3D()，故
        K·u = f_ext + f_eq,     端内力 = K·u − f_eq = EI(κ − κ0)。
  * 悬臂（节点 1 固结、节点 2 自由、零外荷载）：K·u0 = f_eq 精确成立
        u0(v2) = κ0z·L²/2,  u0(θz2) = κ0z·L   （绕 z）
        u0(w2) = −κ0y·L²/2, u0(θy2) = +κ0y·L （绕 y，因 θy = −dw/dx）
    且端内力恒零（自平衡）。
  * 两端固结：U = 0，截面弯矩 M = −EI·κ0（常数），扭矩 T = 0。
  * 简支 + 跨中集中力 P：内力与 κ0 无关（M_sup = 0、M_mid = PL/4），
    位移叠加 v_mid = −P·L³/(48EI) − κ0z·L²/8，支座转角 ±(PL²/16EI + κ0L/2)。

本脚本把每个算例的 .report.txt / .curvature.txt 冻结进 reference/，并把解析值
与 sha256 一起写进 reference/unstressed_curvature_3d.json。compare 脚本只读
冻结证据、只跑当前二进制，不依赖本脚本。

用法：
    python build_reference.py --bin <build_dir> [--cases ...] [--reference ...]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REFERENCE = HERE / "reference"

# ---- 算例材料常数（与 cases/*.dat 严格一致）----
E = 3.45e10
A = 1.0
IY = 1.0
IZ = 1.0
J = 1.0
G = 1.0e10
L = 10.0
K0 = 0.001
P = 1000.0

# 正常算例（需要 --curvature-detail 的与不需要的分开，便于零污染对照）
NORMAL = [
    ("uc01_cantilever_z", True),
    ("uc02_cantilever_y", True),
    ("uc03_fixed_z", True),
    ("uc04_fixed_y", True),
    ("uc05_span_k0", True),
    ("uc06_span_base", False),   # 默认路径零污染基线：不带开关、无 κ0
    ("uc09_torsion_j", True),
]
HARD_NEGATIVE = {
    "uc07_bad_truss": "只支持 BEAM 单元",
    "uc08_bad_nan": "不是有限数",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def run(bin_dir: Path, dat: Path, curvature: bool, outdir: Path) -> subprocess.CompletedProcess:
    """在 outdir 里跑一份拷贝，避免污染 cases/。"""
    work = outdir / dat.stem
    work.mkdir(parents=True, exist_ok=True)
    local = work / dat.name
    shutil.copy2(dat, local)
    cmd = [str(bin_dir / "example_3d_report.exe"), local.name]
    if curvature:
        cmd.append("--curvature-detail")
    return subprocess.run(cmd, cwd=work, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")


def parse_report(txt: str) -> dict:
    """解析 .report.txt：节点位移 6 分量 + 单元内力 key=value。"""
    disp: dict[str, list[float]] = {}
    force: dict[str, dict[str, float]] = {}
    section = ""
    for raw in txt.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "节点位移" in line:
            section = "disp"
            continue
        if "单元内力" in line:
            section = "force"
            continue
        toks = line.split()
        if section == "disp":
            try:
                disp[toks[0]] = [float(v) for v in toks[1:7]]
            except (ValueError, IndexError):
                continue
        elif section == "force":
            kv = {}
            for t in toks[2:]:
                if "=" in t:
                    k, _, v = t.partition("=")
                    try:
                        kv[k] = float(v)
                    except ValueError:
                        pass
            force[toks[0]] = kv
    return {"disp": disp, "force": force}


def parse_curvature(txt: str) -> dict:
    """解析 .curvature.txt（[UNSTRESSED_CURVATURE] 段）。"""
    out: dict = {"basis": None, "error": None, "nelem": None, "rows": {}}
    for raw in txt.splitlines():
        line = raw.strip()
        if not line or line.startswith("["):
            continue
        if line.startswith("BASIS="):
            out["basis"] = line.split("=", 1)[1]
            continue
        if line.startswith("ERROR="):
            out["error"] = line.split("=", 1)[1]
            continue
        if line.startswith("NELEM="):
            out["nelem"] = int(line.split("=", 1)[1])
            continue
        toks = line.split()
        kv = {}
        for t in toks[1:]:
            if "=" in t:
                k, _, v = t.partition("=")
                kv[k] = float(v)
        out["rows"][toks[0]] = kv
    return out


def oracle() -> dict:
    """独立闭式解。单位 m / N / rad，与 .dat 的 SI 口径一致。"""
    EIz = E * IZ
    EIy = E * IY
    return {
        "constants": {"E": E, "A": A, "IY": IY, "IZ": IZ, "J": J, "G": G,
                      "L": L, "K0": K0, "P": P},
        "uc01_cantilever_z": {
            "disp": {"2": [0.0, K0 * L * L / 2.0, 0.0, 0.0, 0.0, K0 * L]},
            "force": {"1": {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                            "My": 0.0, "Mz": 0.0}},
        },
        "uc02_cantilever_y": {
            "disp": {"2": [0.0, 0.0, -K0 * L * L / 2.0, 0.0, K0 * L, 0.0]},
            "force": {"1": {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                            "My": 0.0, "Mz": 0.0}},
        },
        "uc03_fixed_z": {
            "disp": {"1": [0.0] * 6, "2": [0.0] * 6},
            "force": {"1": {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                            "My": 0.0, "Mz": -EIz * K0}},
        },
        "uc04_fixed_y": {
            "disp": {"1": [0.0] * 6, "2": [0.0] * 6},
            "force": {"1": {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                            "My": -EIy * K0, "Mz": 0.0}},
        },
        # 简支 + 跨中集中力：内力与 κ0 无关；位移叠加
        "uc05_span_k0": {
            "disp": {
                "1": [0.0, 0.0, 0.0, 0.0, 0.0,
                      -(P * L * L / (16.0 * EIz) + K0 * L / 2.0)],
                "2": [0.0, -(P * L ** 3 / (48.0 * EIz) + K0 * L * L / 8.0),
                      0.0, 0.0, 0.0, 0.0],
                "3": [0.0, 0.0, 0.0, 0.0, 0.0,
                      +(P * L * L / (16.0 * EIz) + K0 * L / 2.0)],
            },
            "force": {
                "1": {"N": 0.0, "Vy": P / 2.0, "Vz": 0.0, "T": 0.0,
                      "My": 0.0, "Mz": 0.0},
                "2": {"N": 0.0, "Vy": -P / 2.0, "Vz": 0.0, "T": 0.0,
                      "My": 0.0, "Mz": P * L / 4.0},
            },
        },
        "uc06_span_base": {
            "disp": {
                "1": [0.0, 0.0, 0.0, 0.0, 0.0, -(P * L * L / (16.0 * EIz))],
                "2": [0.0, -(P * L ** 3 / (48.0 * EIz)), 0.0, 0.0, 0.0, 0.0],
                "3": [0.0, 0.0, 0.0, 0.0, 0.0, +(P * L * L / (16.0 * EIz))],
            },
            "force": {
                "1": {"N": 0.0, "Vy": P / 2.0, "Vz": 0.0, "T": 0.0,
                      "My": 0.0, "Mz": 0.0},
                "2": {"N": 0.0, "Vy": -P / 2.0, "Vz": 0.0, "T": 0.0,
                      "My": 0.0, "Mz": P * L / 4.0},
            },
        },
        "uc09_torsion_j": {
            "disp": {"1": [0.0] * 6, "2": [0.0] * 6},
            "force": {"1": {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                            "My": 0.0, "Mz": -EIz * K0}},
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", default=str(CASES))
    ap.add_argument("--reference", default=str(REFERENCE))
    args = ap.parse_args()
    bin_dir = Path(args.bin)
    cases = Path(args.cases)
    ref = Path(args.reference)
    ref.mkdir(parents=True, exist_ok=True)

    import tempfile
    frozen: dict[str, str] = {}
    actual: dict[str, dict] = {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for stem, curvature in NORMAL:
            cp = run(bin_dir, cases / f"{stem}.dat", curvature, tmp)
            if cp.returncode != 0:
                print(f"[FAIL] {stem}: {cp.stderr.strip()}")
                return 1
            rep = tmp / stem / f"{stem}.report.txt"
            if not rep.exists():
                print(f"[FAIL] {stem}: 未生成 report")
                return 1
            shutil.copy2(rep, ref / f"{stem}.report.txt")
            frozen[f"{stem}.report.txt"] = sha256(ref / f"{stem}.report.txt")
            actual[stem] = parse_report(rep.read_text(encoding="utf-8"))
            cur = tmp / stem / f"{stem}.curvature.txt"
            if curvature:
                if not cur.exists():
                    print(f"[FAIL] {stem}: 未生成 curvature")
                    return 1
                shutil.copy2(cur, ref / f"{stem}.curvature.txt")
                frozen[f"{stem}.curvature.txt"] = sha256(
                    ref / f"{stem}.curvature.txt")
                actual[stem]["curvature"] = parse_curvature(
                    cur.read_text(encoding="utf-8"))
            elif cur.exists():
                print(f"[FAIL] {stem}: 默认路径不应生成 curvature 文件")
                return 1

    # 冻结前先用 oracle 自检，避免把错值冻进去
    orc = oracle()
    worst = 0.0
    for stem, _ in NORMAL:
        for node, vals in orc[stem]["disp"].items():
            got = actual[stem]["disp"].get(node, [])
            for a, b in zip(vals, got):
                worst = max(worst, abs(a - b) / max(1.0, abs(b)))
        for eid, kv in orc[stem]["force"].items():
            got = actual[stem]["force"].get(eid, {})
            for k, v in kv.items():
                worst = max(worst, abs(v - got.get(k, 0.0)) / max(1.0, abs(v)))
    print(f"oracle 自检：最差相对偏差 = {worst:.3e}")
    if worst > 1e-9:
        print("[FAIL] oracle 与实测偏差超限，拒绝冻结")
        return 1

    doc = {
        "benchmark": "unstressed_curvature_3d",
        "item": 37,
        "basis": "unstressed_curvature_fem_only",
        "semantics": "κ0 是制造/设计给定的无应力曲率，程序不反算；"
                     "FEM 自洽口径，不含真实制造参数、不做 ANSYS 对拍",
        "normal": [s for s, _ in NORMAL],
        "curvature_detail": {s: c for s, c in NORMAL},
        "hard_negative": HARD_NEGATIVE,
        "frozen_sha256": frozen,
        "oracle": orc,
    }
    (ref / "unstressed_curvature_3d.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"已冻结 {len(frozen)} 份证据 → {ref}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
