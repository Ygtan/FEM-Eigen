"""compare_beam_arch.py
=====================
第 44 项 CTest：梁拱 / 拱梁组合建造器（阶段二·线性共同受力）端到端离线比对。

冻结的是 ``beam_arch_3d.py`` 产出的两份 ``.dat``（TIED 系杆拱自锚 /
EXTERNAL 拱梁外部推力），重放线性 3D 内核 ``example_3d_grillage.exe``，
验证"自然语言 → 建造器 → 内核"闭环的力学分岔：

- TIED：拱脚与主梁端共用节点，水平推力沿系梁自锚 —— 系梁恒为拉力且
  接近含拱肋自重的三铰拱闭式 ``(q_deck+q2_deck+w_arch)·L²/(8f)``；
- EXTERNAL：拱脚独立固结、推力入基础 —— 系梁轴力为 0，拱肋仍受压。

纯离线：只读冻结文件、重放内核、读源码文本；不写 ``.dat``、不依赖
pydantic / streamlit，也不需要 matplotlib。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FE_ROOT = ROOT.parents[2]                       # FEM-Eigen
BS_ROOT = FE_ROOT.parent / "BridgeStudio"
RESULTS: list[tuple[bool, str, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((ok, name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def run_case(bin_dir: Path, work: Path, case: Path) -> tuple[int, str, str]:
    exe = bin_dir / "example_3d_grillage.exe"
    if not exe.exists():
        exe = bin_dir / "build" / "example_3d_grillage.exe"
    cp = subprocess.run([str(exe.resolve()), str((work / case.name).resolve())],
                        cwd=work, capture_output=True, text=True,
                        encoding="utf-8", errors="replace")
    sm = work / f"{case.stem}_summary.txt"
    raw = sm.read_text(encoding="utf-8") if sm.exists() else ""
    return cp.returncode, cp.stderr, raw


def parse_beam_N(raw: str) -> dict[int, float]:
    """从 ``final stage beam force`` 段读每单元轴力（N）。i/j 取同一值。"""
    m = re.search(r"final stage beam force\n(.*?)(?:\nfinal stage displacement|\Z)",
                  raw, re.S)
    if not m:
        return {}
    return {int(e): float(n)
            for e, n in re.findall(r"^E(\d+) [ij] N=([-+0-9.eE]+)",
                                   m.group(1), re.M)}


def span(a: int, b: int) -> list[int]:
    return list(range(a, b + 1))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()
    bin_dir, cases, ref = Path(args.bin), Path(args.cases), Path(args.reference)
    doc = json.loads((ref / "beam_arch.json").read_text(encoding="utf-8"))
    b, tol = doc["builder"], doc["tol"]
    tied_ref, ext_ref = doc["tied"], doc["external"]
    arch_ids = span(*b["arch_elems"])
    deck_ids = span(*b["deck_elems"])
    hang_ids = span(*b["hanger_elems"])

    print("=" * 68)
    print("  第 44 项：梁拱 / 拱梁组合建造器（线性共同受力）benchmark")
    print("=" * 68)

    # ---------- 1. 冻结输入 ----------
    print("\n[1] 建造器产出 .dat 逐字节")
    files = ("beam_arch_tied.dat", "beam_arch_external.dat")
    texts: dict[str, str] = {}
    for fn in files:
        case = cases / fn
        check(f"冻结 {fn} 存在", case.exists(), case.name)
        if not case.exists():
            return 1
        got = hashlib.sha256(case.read_bytes()).hexdigest()
        check(f"{fn} 哈希未漂移",
              got == doc["frozen_sha256"][fn], got[:16])
        texts[fn] = case.read_text(encoding="utf-8")
    for fn in files:
        check(f"{fn} 不写 SYSTEM=COOP（不套协作残差）",
              "SYSTEM=COOP" not in texts[fn])
        check(f"{fn} 含 BEAM 与 TRUSS 单元定义",
              "BEAM" in texts[fn] and "TRUSS" in texts[fn])
    check("TIED / EXTERNAL 两种型式确实落盘不同",
          texts[files[0]] != texts[files[1]])

    # ---------- 2. 冻结算例重放 ----------
    print("\n[2] 线性 3D 内核重放")
    tmp = Path(__file__).parent / "_tmp_beam_arch"
    work = tmp / "replay"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    N: dict[str, dict[int, float]] = {}
    for fn in files:
        case = cases / fn
        shutil.copy2(case, work / case.name)
        rc, err, raw = run_case(bin_dir, work, case)
        check(f"{fn} 内核重放成功", rc == 0 and bool(raw), (err or "")[-160:])
        N[fn] = parse_beam_N(raw)
        check(f"{fn} final beam force 段有内力", bool(N[fn]),
              f"entries={len(N[fn])}")
        # 规范要求所有 BEAM 单元都有结果（拱肋 + 主梁）
        need = set(arch_ids) | set(deck_ids)
        check(f"{fn} 拱肋+主梁内力齐全",
              need.issubset(N[fn].keys()), f"need={len(need)} got={len(N[fn])}")

    # ---------- 3. 力学分岔 ----------
    print("\n[3] TIED 系杆拱：推力沿系梁自锚")
    nd = N[files[0]]
    deck = [nd[e] / 1000.0 for e in deck_ids]
    arch = [nd[e] / 1000.0 for e in arch_ids]
    check("系梁全跨受拉（N>0）", all(v > 0 for v in deck),
          f"min={min(deck):.1f} kN")
    check("系梁轴力沿跨恒定（≤1 kN 离散）",
          max(deck) - min(deck) <= tol["deck_spread_abs_kN"],
          f"spread={max(deck) - min(deck):.3f} kN")
    H = statistics.mean(deck)
    H_th = b["H_theory_kN"]
    check(f"系梁拉力 ≈ 含拱肋自重闭式 {H_th:.0f} kN（≤2%）",
          abs(H - H_th) / H_th <= tol["H_rel_tol"], f"got={H:.1f} kN")
    check("拱肋全跨受压（N<0）", all(v < -tol["N_floor_kN"] for v in arch),
          f"max={max(arch):.1f} kN")
    check("拱肋左右拱脚轴力对称",
          abs(arch[0] - arch[-1]) / max(abs(arch[0]), abs(arch[-1]))
          <= tol["symmetry_rel_tol"], f"{arch[0]:.0f} / {arch[-1]:.0f} kN")
    check("拱脚轴力命中冻结值",
          abs(arch[0] - tied_ref["arch_foot_N_kN"])
          / abs(tied_ref["arch_foot_N_kN"]) <= tol["H_rel_tol"],
          f"{arch[0]:.0f} kN")

    print("\n[4] EXTERNAL 拱梁组合：推力入基础，系梁不受轴力")
    ne = N[files[1]]
    deck_e = [ne[e] / 1000.0 for e in deck_ids]
    arch_e = [ne[e] / 1000.0 for e in arch_ids]
    check("系梁轴力 ≈ 0（≤1 kN）",
          all(abs(v) <= tol["external_deck_abs_kN"] for v in deck_e),
          f"max|N|={max(abs(v) for v in deck_e):.3f} kN")
    check("拱肋全跨受压（N<0）",
          all(v < -tol["N_floor_kN"] for v in arch_e),
          f"max={max(arch_e):.1f} kN")
    check("拱肋左右拱脚对称",
          abs(arch_e[0] - arch_e[-1]) / max(abs(arch_e[0]), abs(arch_e[-1]))
          <= tol["symmetry_rel_tol"],
          f"{arch_e[0]:.0f} / {arch_e[-1]:.0f} kN")
    check("两型式拱脚压力量级一致（≤5%）",
          abs(arch[0] - arch_e[0]) / max(abs(arch[0]), abs(arch_e[0]))
          <= tol["arch_foot_cross_rel_tol"],
          f"tied={arch[0]:.0f} external={arch_e[0]:.0f} kN")

    # ---------- 5. 默认零污染与源码契约 ----------
    print("\n[5] 默认零污染与 X20 源码契约")
    md_py = (BS_ROOT / "model_doc.py").read_text(encoding="utf-8")
    check("model_doc 路由 BEAM_ARCH → _build_beam_arch",
          "combo_kind" in md_py and "_build_beam_arch" in md_py
          and "from beam_arch_3d import build_beam_arch_3d" in md_py)
    m = re.search(r"BUILDERS: dict\[tuple\[str, str\], str\] = \{(.*?)\n\}",
                  md_py, re.S)
    entries = re.findall(r'^\s*\("([a-z_]+)", "[23]d"\):', m.group(1), re.M) if m else []
    check("BUILDERS 仍 12 格 / 6 kind",
          len(entries) == 12 and len(set(entries)) == 6,
          f"cells={len(entries)} kinds={len(set(entries))}")
    solver_py = (BS_ROOT / "solver_3d.py").read_text(encoding="utf-8")
    check("未新增梁拱专属求解键",
          '"beam_arch"' not in solver_py and "beam_arch" not in solver_py.split("_registry")[0][:2000])
    ba_py = (BS_ROOT / "beam_arch_3d.py").read_text(encoding="utf-8")
    check("beam_arch_3d.py 签名与两种型式常量在位",
          "def build_beam_arch_3d(" in ba_py
          and "TIED_ARCH" in ba_py and "EXTERNAL_ARCH" in ba_py
          and "H_theory" in ba_py)
    spec = (BS_ROOT / "组合体系接入规格.md").read_text(encoding="utf-8")
    check("规格 §4.7 梁拱口径在位",
          "4.7 梁拱 / 拱梁组合" in spec and "含拱肋自重" in spec)

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
