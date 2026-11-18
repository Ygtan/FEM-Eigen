#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 40 项 CTest：协作体系（斜拉-悬索）benchmark，阶段一内核侧。

离线（无需 ANSYS / 网络 / 许可）：
  * 冻结协作 .dat 的 sha256 与 reference 对齐；
  * 当前 cable_solver.exe 重放：SYSTEM=COOP → cooperative=true，
    斜拉索目标力保留（找形不清空），四类残差 r_main/r_stay/r_girder/r_tower
    各自命中，[COOPERATIVE] 段字段可读且 coop_ok=1；
  * 负例（动态生成）：协作斜拉索缺省目标力、单元编号冲突 均被显式拒绝；
  * 默认零污染：去掉 SYSTEM=COOP 的同一悬索模型输出不含 [COOPERATIVE] 段，
    既有主缆/吊索索力口径不变；
  * 源码契约：cooperative 独立字段、suspElems 保留集合、拒绝分支在位。

用法：
    python compare_combo_system_3d.py --bin <dir> \
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


def parse_head(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith("["):
            break
        for k, v in KV.findall(line):
            out.setdefault(k, v)
    return out


def parse_section(text: str, name: str) -> dict[str, str]:
    m = re.search(r"\[" + re.escape(name) + r"\]\n(.*?)(?:\n\[|\Z)", text, re.S)
    if not m:
        return {}
    out: dict[str, str] = {}
    for k, v in KV.findall(m.group(1)):
        out.setdefault(k, v)
    return out


def num(d: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(d.get(key, default))
    except (TypeError, ValueError):
        return default


def run_exe(bin_dir: Path, work: Path, case: Path) -> tuple[int, str, str]:
    exe = bin_dir / "cable_solver.exe"
    if not exe.exists():
        exe = bin_dir / "build" / "cable_solver.exe"
    cp = subprocess.run([str(exe.resolve()), case.name], cwd=work,
                        capture_output=True, text=True,
                        encoding="utf-8", errors="replace")
    out = work / f"{case.stem}_cable.txt"
    raw = out.read_text(encoding="utf-8") if out.exists() else ""
    return cp.returncode, cp.stderr, raw


def replay(bin_dir: Path, case: Path, tmp: Path) -> tuple[dict[str, str], dict[str, str], str]:
    work = tmp / case.stem
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case, work / case.name)
    rc, err, raw = run_exe(bin_dir, work, case)
    if rc != 0 or not raw:
        raise RuntimeError(f"{case.name} 求解失败 rc={rc}: " + (err or "")[-400:])
    return parse_head(raw), parse_section(raw, "COOPERATIVE"), raw


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()
    bin_dir, cases, ref = (Path(args.bin), Path(args.cases), Path(args.reference))
    doc = json.loads((ref / "combo_system_3d.json").read_text(encoding="utf-8"))

    print("=" * 68)
    print("  第 40 项：协作体系（斜拉-悬索）benchmark")
    print("=" * 68)

    # ---------- 1. 冻结输入 ----------
    print("\n[1] 冻结输入逐字节")
    case_file = cases / "combo_case.dat"
    check("冻结 .dat 存在", case_file.exists())
    if case_file.exists():
        got = hashlib.sha256(case_file.read_bytes()).hexdigest()
        want = doc["frozen_sha256"]["combo_case.dat"]
        check("冻结哈希未漂移", got == want, got[:16])

    # ---------- 2. 协作重放 ----------
    print("\n[2] 当前二进制重放（协作体系）")
    c = doc["combo"]
    coop_ref = doc["cooperative"]
    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            head, coop, raw = replay(bin_dir, case_file, tmp)
    except Exception as e:  # noqa: BLE001
        check("协作算例重放", False, str(e).splitlines()[-1][:250])
        return 1

    check("收敛", num(head, "converged") == c["converged"])
    check("主缆 H 命中冻结值", close(num(head, "H_kN"), c["H_kN"]),
          f"H={head.get('H_kN')} kN")
    check("矢高命中目标", close(num(head, "sag_m"), c["sag_m"]),
          f"sag={head.get('sag_m')} m")
    check("索数 = 主缆段+吊索+斜拉索",
          int(num(head, "cables")) == c["cables"]
          == c["main_cable_segments"] + c["hanger_count"] + 4,
          f"cables={head.get('cables')}")

    print("\n[3] [COOPERATIVE] 段与四类残差")
    check("[COOPERATIVE] 段存在", bool(coop), str(sorted(coop)[:4]))
    check("coop_system=COOP", coop.get("coop_system") == "COOP")
    check("susp_segments=205 / stay_cables=4 / susp_planes=1",
          int(num(coop, "susp_segments")) == coop_ref["susp_segments"]
          and int(num(coop, "stay_cables")) == coop_ref["stay_cables"]
          and int(num(coop, "susp_planes")) == coop_ref["susp_planes"],
          f"susp={coop.get('susp_segments')} stay={coop.get('stay_cables')}")
    for key, bound in (("r_main", "r_main_le"), ("r_stay", "r_stay_le"),
                       ("r_girder", "r_girder_le"), ("r_tower", "r_tower_le")):
        v = num(coop, key, 9e9)
        check(f"{key} 命中各自容差（≤{coop_ref[bound]:g}）", v <= coop_ref[bound],
              f"{key}={v:.4g}")
    check("coop_ok=1（四类同时满足）",
          int(num(coop, "coop_ok")) == coop_ref["coop_ok"])
    thm, ths = num(coop, "coop_tower_H_main_kN"), num(coop, "coop_tower_H_side_kN")
    check("塔两侧水平力基本平衡（差≤5%）",
          thm > 0 and ths > 0 and abs(thm - ths) / max(thm, ths) <= 0.05,
          f"main={thm:.1f} side={ths:.1f} kN")

    # 斜拉索目标力保留：找形后 T 应接近给定目标（非悬索找形值）
    print("\n[4] 斜拉索目标力保留（找形未覆盖）")
    stay_rows: dict[str, list[str]] = {}
    for line in raw.splitlines():
        mm = re.match(r"^(60[1-4])\s+(\d+)\s+([\d.eE+-]+)", line)
        if mm:
            stay_rows[mm.group(1)] = line.split()
    check("4 根斜拉索全部登记", len(stay_rows) == 4, str(sorted(stay_rows)))
    for eid, target in c["stay_targets_N"].items():
        row = stay_rows.get(eid)
        ok = row is not None and close(float(row[2]), target, 1e-2)
        check(f"索 {eid} 索力命中目标 {target/1e6:g}MN（相对 1%）", ok,
              f"T={row[2] if row else 'NA'}")

    # ---------- 3. 默认零污染 ----------
    print("\n[5] 默认零污染（去 SYSTEM=COOP 走既有悬索路径）")
    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            plain = tmp / "plain_case.dat"
            s = case_file.read_text(encoding="utf-8")
            s = s.replace("SYSTEM=COOP\n", "")
            # 去掉斜拉索负例条目，恢复单悬索模型
            s = re.sub(r"\*TYPE,3\n\*MAT,3\n\*REAL,4\n(\*EN,60[1-4],[^\n]+\n)+", "", s)
            s = re.sub(r"\*CABLE\n(?:#[^\n]*\n)?(?:\*EN[^\n]*\n)?"
                       r"(?:60[1-4] [^\n]+\n)+", "*CABLE\n", s)
            plain.write_text(s, encoding="utf-8", newline="\n")
            rc, err, plain_raw = run_exe(bin_dir, tmp, plain)
            check("既有路径仍收敛", rc == 0 and "converged=1" in plain_raw,
                  f"rc={rc}")
            check("既有路径输出无 [COOPERATIVE] 段", "[COOPERATIVE]" not in plain_raw)
            ph = parse_head(plain_raw)
            check("既有路径索数 = 主缆段 + 吊索（无斜拉索）",
                  int(num(ph, "cables")) == c["main_cable_segments"] + c["hanger_count"],
                  f"cables={ph.get('cables')}")
    except Exception as e:  # noqa: BLE001
        check("默认零污染重放", False, str(e).splitlines()[-1][:250])

    # ---------- 4. 拒绝路径（动态负例）----------
    print("\n[6] 显式拒绝（缺省目标力 / 编号冲突）")
    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            base = case_file.read_text(encoding="utf-8")
            # 负例 A：斜拉索目标力缺省
            negA = tmp / "neg_no_target.dat"
            negA.write_text(base.replace("601 30 2000000.0", "601 30 0"),
                            encoding="utf-8", newline="\n")
            rcA, errA, _ = run_exe(bin_dir, tmp, negA)
            check("缺省目标力被拒绝", rcA != 0 and "必须给目标索力" in errA,
                  f"rc={rcA}")
            # 负例 B：编号冲突——把 *CABLE 里一根斜拉索的单元号改成主缆单元号，
            # 使同一单元既在 *CABLE 又在 *MAINCABLE。注意不能把 601 插进
            # *MAINCABLE 链首：那会先触发主缆链序校验（节点不连续），
            # 走不到协作体系冲突分支，错误文本对不上。
            negB = tmp / "neg_conflict.dat"
            b = base.replace("601 30 2000000.0", "201 30 2000000.0", 1)
            negB.write_text(b, encoding="utf-8", newline="\n")
            rcB, errB, _ = run_exe(bin_dir, tmp, negB)
            check("编号冲突被拒绝", rcB != 0 and "编号冲突" in errB,
                  f"rc={rcB}")
    except Exception as e:  # noqa: BLE001
        check("拒绝路径重放", False, str(e).splitlines()[-1][:250])

    # ---------- 5. 源码契约 ----------
    print("\n[7] FE / BS 源码契约")
    cpp = (FE_ROOT / "cable_solver.cpp").read_text(encoding="utf-8")
    check("Cfg 有独立 cooperative 字段（默认 false）",
          "bool cooperative = false;" in cpp)
    check("*CABLE_BRIDGE 解析 SYSTEM=COOP", 'key == "SYSTEM"' in cpp
          and 'upper(val) == "COOP"' in cpp)
    check("suspElems 保留集合 + 只清悬索索",
          "std::set<int> suspElems" in cpp and "suspElems.count(r.elem)" in cpp)
    check("非协作路径保留旧 clear() 全量重填",
          "} else {\n                // 目标索力全部由找形与荷载推出" in cpp)
    check("[COOPERATIVE] 段仅协作时写", "[COOPERATIVE]" in cpp
          and "if (cfg.cooperative)" in cpp)
    fe = (BS_ROOT / "fem_eigen_3d.py").read_text(encoding="utf-8")
    check("fem_eigen_3d.sections 登记 [COOPERATIVE]", '"cooperative"' in fe)
    spec = (BS_ROOT / "组合体系接入规格.md").read_text(encoding="utf-8")
    check("规格 §3.7 用户口径在位",
          "3.7 用户确认的协作体系力学口径" in spec and "r_tower" in spec)

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
