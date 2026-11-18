#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 38 项 CTest：双主缆悬索桥 benchmark（双索面 3D opt-in）。

离线（无需 ANSYS / 网络 / 许可）：
  * 两份冻结 .dat 的 sha256 与 reference 对齐；
  * 当前 cable_solver.exe 重放：双主缆 plane_1 / plane_2 独立找形、
    H 相等、全部索登记（2×(主缆段+吊索)）；
  * 单索面对照输出**不含** susp_plane_count，默认行为逐字节语义不变；
  * 源码契约：附加平面段解析在位、第二平面索不并入 cfg.mainCable 单链、
    BridgeStudio 侧 *.dat 追加 *_2 段与镜像模块在位。

用法：
    python compare_dual_plane_suspension_3d.py --bin <dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS: list[tuple[bool, str, str]] = []
FE_ROOT = HERE.parents[2]
BS_ROOT = FE_ROOT.parent / "BridgeStudio"
KV = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)=([^\s]+)")


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((bool(ok), name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def close(a: float, b: float, rel: float = 1e-6) -> bool:
    return abs(a - b) <= max(1e-9, rel * max(abs(a), abs(b)))


def parse_kv(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith("["):
            break
        for k, v in KV.findall(line):
            out.setdefault(k, v)
    return out


def num(d: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(d.get(key, default))
    except (TypeError, ValueError):
        return default


def replay(bin_dir: Path, case: Path, tmp: Path) -> tuple[dict[str, str], str]:
    work = tmp / case.stem
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case, work / case.name)
    exe = bin_dir / "cable_solver.exe"
    if not exe.exists():
        exe = bin_dir / "build" / "cable_solver.exe"
    cp = subprocess.run([str(exe.resolve()), case.name], cwd=work,
                        capture_output=True, text=True,
                        encoding="utf-8", errors="replace")
    out = work / f"{case.stem}_cable.txt"
    if cp.returncode != 0 or not out.exists():
        raise RuntimeError(f"{case.name} 求解失败 rc={cp.returncode}: "
                           + (cp.stdout or cp.stderr)[-400:])
    return parse_kv(out.read_text(encoding="utf-8")), out.read_text(encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()
    bin_dir, cases, ref = (Path(args.bin), Path(args.cases), Path(args.reference))
    doc = json.loads((ref / "dual_plane_suspension_3d.json").read_text(encoding="utf-8"))

    print("=" * 68)
    print("  第 38 项：双主缆悬索桥 benchmark（双索面 opt-in）")
    print("=" * 68)

    # ---------- 1. 冻结输入 ----------
    print("\n[1] 冻结输入逐字节")
    for name, want in doc["frozen_sha256"].items():
        p = cases / name
        check(f"冻结 .dat 存在：{name}", p.exists())
        if p.exists():
            got = hashlib.sha256(p.read_bytes()).hexdigest()
            check(f"冻结哈希未漂移：{name}", got == want, f"{got[:16]}")

    # ---------- 2. 重放对拍 ----------
    print("\n[2] 当前二进制重放")
    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            dual, dual_raw = replay(bin_dir, cases / "dual_susp_case.dat", tmp)
            single, single_raw = replay(bin_dir, cases / "single_susp_case.dat", tmp)
    except Exception as e:  # noqa: BLE001
        check("双/单索面重放", False, str(e).splitlines()[-1][:250])
        return 1

    dref, sref = doc["dual"], doc["single"]
    check("[双] 收敛", num(dual, "converged") == dref["converged"])
    check("[双] susp_plane_count=2",
          int(num(dual, "susp_plane_count")) == dref["susp_plane_count"])
    h1, h2 = num(dual, "plane_1_H_kN"), num(dual, "plane_2_H_kN")
    check("[双] 两平面 H 相等且命中冻结值",
          h1 > 0 and close(h1, h2) and close(h1, dref["plane_1_H_kN"]),
          f"H1={h1:.1f} H2={h2:.1f} kN")
    check("[双] 索数 = 2×(主缆段 + 吊索)",
          int(num(dual, "cables")) == dref["cables"]
          == 2 * (dref["main_cable_segments"] + dref["hanger_count"]),
          f"cables={dual.get('cables')}")
    check("[双] 矢高命中目标",
          close(num(dual, "sag_m"), dref["sag_m"]),
          f"sag={dual.get('sag_m')} m")

    check("[单] 收敛", num(single, "converged") == sref["converged"])
    check("[单] 输出无 susp_plane_count（默认零污染）",
          "susp_plane_count" not in single,
          f"键存在={single.get('susp_plane_count')}")
    check("[单] 索数 = 主缆段 + 吊索",
          int(num(single, "cables")) == sref["cables"],
          f"cables={single.get('cables')}")
    check("[单] H 命中冻结值",
          close(num(single, "H_kN"), sref["H_kN"]),
          f"H={single.get('H_kN')} kN")

    # ---------- 3. 源码契约 ----------
    print("\n[3] FE / BS 源码契约")
    cpp = (FE_ROOT / "cable_solver.cpp").read_text(encoding="utf-8")
    check("内核解析 *SUSP_BRIDGE_PLANE_n / *MAINCABLE_n / *HANGER_n",
          all(k in cpp for k in ("SUSP_BRIDGE_PLANE_", "MAINCABLE_", "HANGER_",
                                 "struct SuspPlane")))
    check("第二平面索不并入 cfg.mainCable 单链",
          "cfg.planes[pi].mainCable" in cpp
          and "cfg.mainCable.insert(cfg.mainCable.end()" not in cpp)
    check("多平面报告段仅多平面时输出",
          "susp_plane_count=%zu" in cpp and "plane_%zu_H_kN=" in cpp
          and "fit.planeH.size() > 1" in cpp)
    spec = (BS_ROOT / "model_spec_3d.py").read_text(encoding="utf-8")
    mirror = BS_ROOT / "dual_plane_3d.py"
    check("BridgeStudio 追加平面段 + 镜像模块在位",
          "SUSP_BRIDGE_PLANE_{_pi}" in spec
          and mirror.exists()
          and "def mirror_to_dual_plane" in mirror.read_text(encoding="utf-8"))

    n = sum(1 for ok, _, _ in RESULTS if ok)
    print("\n" + "=" * 68)
    print(f"  {n}/{len(RESULTS)} 通过")
    if n < len(RESULTS):
        for ok, name, detail in RESULTS:
            if not ok:
                print(f"    · {name}  {detail}")
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
