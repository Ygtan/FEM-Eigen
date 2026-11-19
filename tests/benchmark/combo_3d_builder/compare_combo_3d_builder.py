"""compare_combo_3d_builder.py
=============================
第 43 项 CTest：阶段二 combo_3d.py 建造器 → 内核端到端离线比对。

与第 40 项（手工冻结 .dat，只验内核机制）的区别：本项冻结的是**建造器
产出的 .dat**，验证 §3.7 口径参数落盘、协作残差命中、默认零污染，以及
Python 侧路由（kind=suspension + combo_kind=STAY_SUSP → 第 5 求解键）。

纯离线：只读冻结文件、重放 cable_solver.exe、读源码文本；不写 .dat、
不依赖 pydantic / streamlit。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FE_ROOT = ROOT.parents[2]                       # FEM-Eigen
BS_ROOT = FE_ROOT.parent / "BridgeStudio"
RESULTS: list[tuple[bool, str, str]] = []
KV = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^\s]+)")


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((ok, name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def parse_head(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for k, v in KV.findall(text.split("[", 1)[0]):
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


def run_case(bin_dir: Path, work: Path, case: Path) -> tuple[int, str, str]:
    exe = bin_dir / "cable_solver.exe"
    if not exe.exists():
        exe = bin_dir / "build" / "cable_solver.exe"
    cp = subprocess.run([str(exe.resolve()), case.name], cwd=work,
                        capture_output=True, text=True,
                        encoding="utf-8", errors="replace")
    out = work / f"{case.stem}_cable.txt"
    raw = out.read_text(encoding="utf-8") if out.exists() else ""
    return cp.returncode, cp.stderr, raw


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()
    bin_dir, cases, ref = Path(args.bin), Path(args.cases), Path(args.reference)
    doc = json.loads((ref / "combo_3d_builder.json").read_text(encoding="utf-8"))
    c, lim = doc["combo"], doc["cooperative"]

    print("=" * 68)
    print("  第 43 项：组合体系建造器（阶段二）benchmark")
    print("=" * 68)

    # ---------- 1. 冻结输入 ----------
    print("\n[1] 建造器产出 .dat 逐字节")
    case = cases / "combo_3d_builder.dat"
    check("冻结 .dat 存在", case.exists(), case.name)
    if not case.exists():
        return 1
    got = hashlib.sha256(case.read_bytes()).hexdigest()
    want = doc["frozen_sha256"]["combo_3d_builder.dat"]
    check("冻结哈希未漂移", got == want, got[:16])
    dat = case.read_text(encoding="utf-8")
    check("建造器写入 SYSTEM=COOP", "SYSTEM=COOP" in dat)
    b = doc["builder"]
    for eid in b["stay_elems"]:
        pattern = rf"^{eid}\s+\d+\s+[0-9.]+\s*$"
        check(f"斜拉索 {eid} 目标索力落盘",
              bool(re.search(pattern, dat, re.M)))

    # ---------- 2. 协作算例重放 ----------
    print("\n[2] 冻结算例重放与协作残差")
    tmp = Path(__file__).parent / "_tmp_combo_builder"
    work = tmp / case.stem
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    shutil.copy2(case, work / case.name)
    rc, err, raw = run_case(bin_dir, work, case)
    check("cable_solver 重放成功", rc == 0 and bool(raw), (err or "")[-160:])
    if rc != 0 or not raw:
        return 1
    head, coop = parse_head(raw), parse_section(raw, "COOPERATIVE")
    c = doc["combo"]
    check("收敛", num(head, "converged") == c["converged"])
    H = num(head, "H_kN")
    check("主缆 H 命中冻结值",
          abs(H - c["H_kN"]) <= c["H_kN"] * c["H_rel_tol"], f"{H:.1f} kN")
    check("矢高命中目标", abs(num(head, "sag_m") - c["sag_m"]) < 1e-6)
    check("[COOPERATIVE] 段存在", bool(coop), str(sorted(coop)[:4]))
    check("coop_system=COOP", coop.get("coop_system") == lim["coop_system"])
    for key, want_v in (("susp_segments", c["susp_segments"]),
                        ("stay_cables", c["stay_cables"]),
                        ("susp_planes", c["susp_planes"])):
        check(f"{key}={want_v}", int(num(coop, key)) == want_v,
              f"got={num(coop, key):g}")
    for key in ("r_main", "r_stay", "r_girder", "r_tower"):
        v = num(coop, key)
        check(f"{key} ≤ {lim[key + '_le']:g}", v <= lim[key + "_le"],
              f"{v:.3e}")
    check("coop_ok=1（四类同时满足）",
          int(num(coop, "coop_ok")) == lim["coop_ok"])
    # 塔两侧水平力基本平衡（§3.7「边跨略短、基本能够平衡」）
    hm = num(coop, "coop_tower_H_main_kN")
    hs = num(coop, "coop_tower_H_side_kN")
    check("塔两侧水平力基本平衡（差≤5%）",
          hm > 0 and hs > 0 and abs(hm - hs) / max(hm, hs) <= 0.05,
          f"main={hm:.1f} side={hs:.1f} kN")

    # ---------- 3. 默认零污染 ----------
    print("\n[3] 去掉 SYSTEM=COOP 走既有悬索路径")
    plain = work / "plain.dat"
    plain.write_text(dat.replace("SYSTEM=COOP\n", ""), encoding="utf-8")
    rc2, err2, raw2 = run_case(bin_dir, work, plain)
    check("既有路径仍收敛", rc2 == 0 and "converged=1" in raw2,
          (err2 or "")[-160:])
    check("既有路径输出无 [COOPERATIVE] 段", "[COOPERATIVE]" not in raw2)
    check("既有路径索数 = 主缆+吊索（无斜拉）",
          num(parse_head(raw2), "cables") == c["susp_segments"],
          f"got={num(parse_head(raw2), 'cables'):g}")

    # ---------- 4. BS 源码契约 ----------
    print("\n[4] BridgeStudio 阶段二路由与 §3.7 参数")
    combo_py = (BS_ROOT / "combo_3d.py").read_text(encoding="utf-8")
    check("combo_3d.py 有 build_combo_3d 签名",
          "def build_combo_3d(" in combo_py
          and "-> tuple[ModelSpec3D, dict]" in combo_py)
    check("建造器置 cooperative=True", "cb.cooperative = True" in combo_py)
    for key in ("susp_fraction", "side_stay_ratio", "tower_height_ratio",
                "plane_count", "coop_boundary_node"):
        check(f"§3.7 参数 {key} 已登记", key in combo_py)
    md_py = (BS_ROOT / "model_doc.py").read_text(encoding="utf-8")
    check("model_doc 路由 combo_kind=STAY_SUSP → combo_3d",
          "combo_kind" in md_py and "_build_combo" in md_py
          and "from combo_3d import build_combo_3d" in md_py)
    check("12 格 BUILDERS / 6 kind 未被破坏",
          '("suspension", "3d"): "suspension"' in md_py
          and "BUILDERS: dict[tuple[str, str], str]" in md_py)
    solver_py = (BS_ROOT / "solver_3d.py").read_text(encoding="utf-8")
    check("solver_3d 第 5 键 combo 在位",
          '"combo": SolverSpec' in solver_py and 'return "combo"' in solver_py)
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
