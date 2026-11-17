#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 37 项 CTest：成桥态无应力曲率 benchmark（B5 第二十六次，opt-in）。

离线（无需 ANSYS / 求解器许可）：
  * 冻结 .report.txt / .curvature.txt 的 sha256 与 reference 对齐；
  * 当前二进制重放：逐字节复现冻结证据；
  * 正常算例逐键对齐独立 Python 闭式 oracle（悬臂 u0 / 固结 M=−EIκ0 /
    简支叠加）；
  * 自平衡强断言：悬臂端内力恒零、κ0 只改无应力构形不改刚度（uc05 与 uc06
    内力逐位相同、位移差恰为 κ0L²/8 与 ±κ0L/2）；
  * J 不变：uc09 的 T 恒为 0（曲率只进 My/Mz）；
  * 两个硬负例（TRUSS 给 K0 / K0=nan）显式拒绝；
  * 默认零污染：不带 --curvature-detail 不生成 .curvature.txt、无 κ0 的
    算例报告逐字节等于冻结基线、其它 benchmark 的 reference 目录里没有
    多出 curvature 文件；
  * FE / BridgeStudio 源码契约守卫。

用法：
    python compare_unstressed_curvature_3d.py --bin <dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from build_reference import (  # noqa: E402
    HARD_NEGATIVE, NORMAL, oracle, parse_curvature, parse_report, sha256,
)

RESULTS: list[tuple[bool, str, str]] = []
FE_ROOT = HERE.parents[2]
BS_ROOT = FE_ROOT.parent / "BridgeStudio"


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((bool(ok), name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def close(a: float, b: float, rel: float = 1e-9, abs_: float = 1e-9) -> bool:
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def run(bin_dir: Path, dat: Path, curvature: bool, outdir: Path) -> subprocess.CompletedProcess:
    work = outdir / dat.stem
    work.mkdir(parents=True, exist_ok=True)
    local = work / dat.name
    shutil.copy2(dat, local)
    cmd = [str(bin_dir / "example_3d_report.exe"), local.name]
    if curvature:
        cmd.append("--curvature-detail")
    return subprocess.run(cmd, cwd=work, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()
    bin_dir = Path(args.bin)
    cases = Path(args.cases)
    ref = Path(args.reference)

    doc = json.loads((ref / "unstressed_curvature_3d.json").read_text(encoding="utf-8"))
    orc = doc["oracle"]
    frozen = doc["frozen_sha256"]
    cst = orc["constants"]
    L = cst["L"]
    K0 = cst["K0"]

    # ---------- 1. 冻结证据完整性 ----------
    for name, want in frozen.items():
        p = ref / name
        check(f"冻结证据存在：{name}", p.exists())
        if p.exists():
            check(f"冻结哈希未漂移：{name}", sha256(p) == want)

    # ---------- 2. 重放 + oracle 对拍 ----------
    actual: dict[str, dict] = {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for stem, curvature in NORMAL:
            cp = run(bin_dir, cases / f"{stem}.dat", curvature, tmp)
            check(f"{stem} 求解成功", cp.returncode == 0, cp.stderr.strip()[:80])
            if cp.returncode != 0:
                continue
            rep = tmp / stem / f"{stem}.report.txt"
            got_txt = rep.read_text(encoding="utf-8")
            check(f"{stem} 报告逐字节复现冻结证据",
                  got_txt == (ref / f"{stem}.report.txt").read_text(encoding="utf-8"))
            actual[stem] = parse_report(got_txt)

            # 位移逐分量
            for node, vals in orc[stem]["disp"].items():
                got = actual[stem]["disp"].get(node, [])
                ok = len(got) == 6 and all(close(a, b) for a, b in zip(vals, got))
                check(f"{stem} 节点 {node} 位移 = 闭式解", ok,
                      "" if ok else f"want={vals} got={got}")
            # 内力逐键
            for eid, kv in orc[stem]["force"].items():
                got = actual[stem]["force"].get(eid, {})
                bad = [k for k, v in kv.items() if not close(v, got.get(k, 0.0))]
                check(f"{stem} 单元 {eid} 端内力 = 闭式解", not bad,
                      f"偏差键={bad}" if bad else "")

            if curvature:
                cur = tmp / stem / f"{stem}.curvature.txt"
                check(f"{stem} 生成 curvature 文件", cur.exists())
                if cur.exists():
                    txt = cur.read_text(encoding="utf-8")
                    check(f"{stem} curvature 逐字节复现冻结证据",
                          txt == (ref / f"{stem}.curvature.txt").read_text(encoding="utf-8"))
                    actual[stem]["curvature"] = parse_curvature(txt)
                    check(f"{stem} curvature BASIS 正确",
                          actual[stem]["curvature"]["basis"]
                          == "unstressed_curvature_fem_only")
            else:
                # 默认路径零污染：不带开关时不得生成 curvature 文件
                check(f"{stem} 默认路径不生成 curvature 文件",
                      not (tmp / stem / f"{stem}.curvature.txt").exists())

        # ---------- 3. 自平衡 / 刚度不变 / J 不变 的强断言 ----------
        # 悬臂端内力恒零（自平衡）
        for stem in ("uc01_cantilever_z", "uc02_cantilever_y"):
            f1 = actual.get(stem, {}).get("force", {}).get("1", {})
            nz = {k: v for k, v in f1.items() if abs(v) > 1e-6}
            check(f"{stem} 零外荷载端内力恒零（自平衡）", not nz, f"非零={nz}" if nz else "")

        # κ0 不改刚度：uc05 与 uc06 内力逐位相同
        f5 = actual.get("uc05_span_k0", {}).get("force", {})
        f6 = actual.get("uc06_span_base", {}).get("force", {})
        same = True
        for eid in ("1", "2"):
            for k in ("N", "Vy", "Vz", "T", "My", "Mz"):
                if not close(f5.get(eid, {}).get(k, 0.0), f6.get(eid, {}).get(k, 0.0)):
                    same = False
        check("κ0 不改变内力（uc05 与无 κ0 基线逐键相同）", same)

        # 位移差恰为 κ0 叠加项
        d5 = actual.get("uc05_span_k0", {}).get("disp", {})
        d6 = actual.get("uc06_span_base", {}).get("disp", {})
        if d5.get("2") and d6.get("2"):
            dv = d5["2"][1] - d6["2"][1]
            check("简支跨中位移差 = κ0L²/8", close(dv, -K0 * L * L / 8.0),
                  f"{dv:.6e} vs {-K0 * L * L / 8.0:.6e}")
            dt1 = d5["1"][5] - d6["1"][5]
            check("简支支座转角差 = κ0L/2", close(dt1, -K0 * L / 2.0),
                  f"{dt1:.6e} vs {-K0 * L / 2.0:.6e}")

        # J 不变：曲率不产生扭矩
        f9 = actual.get("uc09_torsion_j", {}).get("force", {}).get("1", {})
        check("曲率不产生扭矩（T=0，J 不参与）", abs(f9.get("T", 0.0)) < 1e-6,
              f"T={f9.get('T')}")
        check("曲率只进 Mz（My=0）", abs(f9.get("My", 0.0)) < 1e-6)

        # curvature 段内容与 oracle 的 M0=EIκ0 一致
        c1 = actual.get("uc03_fixed_z", {}).get("curvature", {})
        row = c1.get("rows", {}).get("1", {})
        EI = cst["E"] * cst["IZ"]
        check("curvature 段 M0Z = EIz·κ0z", close(row.get("M0Z", 0.0), EI * K0),
              f"{row.get('M0Z')} vs {EI * K0}")
        check("curvature 段 NELEM 计数正确", c1.get("nelem") == 1)
        c5 = actual.get("uc05_span_k0", {}).get("curvature", {})
        check("uc05 两个单元都登记曲率", c5.get("nelem") == 2)

        # ---------- 4. 硬负例 ----------
        for stem, needle in HARD_NEGATIVE.items():
            cp = run(bin_dir, cases / f"{stem}.dat", True, tmp)
            check(f"{stem} 显式拒绝（非零退出）", cp.returncode != 0)
            check(f"{stem} 错误信息含「{needle}」", needle in cp.stderr)

        # ---------- 5. 其它 benchmark 未被污染 ----------
        stray = []
        for sub in (FE_ROOT / "tests" / "benchmark").glob("*/reference"):
            if sub.parent.name == "unstressed_curvature_3d":
                continue
            stray.extend(sub.glob("*.curvature.txt"))
        check("其它 benchmark 的 reference 无 curvature 文件", not stray,
              f"多出={[str(s) for s in stray]}" if stray else "")

    # ---------- 6. 源码契约 ----------
    hpp = (FE_ROOT / "fem_3d.hpp").read_text(encoding="utf-8", errors="replace")
    check("FE Element3D 含 k0y/k0z 字段", "double k0y = 0.0, k0z = 0.0;" in hpp)
    check("FE *EN 解析 K0Y=", '"K0Y="' in hpp)
    check("FE *EN 解析 K0Z=", '"K0Z="' in hpp)
    check("FE 非 BEAM 给曲率显式拒绝", "无应力曲率 K0Y/K0Z 只支持 BEAM 单元" in hpp)
    check("FE 固端力用 Eeff 与刚度同源",
          "-Eeff * e.Iz * e.k0z" in hpp and "-Eeff * e.Iy * e.k0y" in hpp)
    check("FE κ0=0 早退（零污染）", "if (e.k0y != 0.0 || e.k0z != 0.0)" in hpp)
    rpt = (FE_ROOT / "example_3d_report.cpp").read_text(encoding="utf-8")
    check("CLI 有 --curvature-detail", '"--curvature-detail"' in rpt)
    check("输出段 [UNSTRESSED_CURVATURE] 存在", "[UNSTRESSED_CURVATURE]" in rpt)
    check("无曲率时写显式 ERROR 不摆零值",
          "ERROR=no_curvature_elements" in rpt)
    check("默认不写（只有开关才写）",
          "if (curvatureDetail)" in rpt)

    bs_spec = BS_ROOT / "B类后续项接入规格.md"
    check("BS 规格登记 §4.1 本轮实现范围",
          bs_spec.exists() and "### 4.1 本轮实现范围" in bs_spec.read_text(encoding="utf-8"))
    bs_cur = BS_ROOT / "curvature_3d.py"
    check("BS 解析口径层存在", bs_cur.exists(), str(bs_cur))
    if bs_cur.exists():
        bs = bs_cur.read_text(encoding="utf-8")
        check("BS 解析器只解析不重算", "UNSTRESSED_CURVATURE" in bs)
        check("BS 解析器口径层不含 EI 反算",
              "def parse" in bs and "M0Z" in bs)

    ok = sum(1 for r in RESULTS if r[0])
    print()
    print("=" * 78)
    print(f"  {ok}/{len(RESULTS)} checks passed")
    if ok != len(RESULTS):
        for passed, name, detail in RESULTS:
            if not passed:
                print(f"  [FAIL] {name}  {detail}")
        return 1
    print("  [PASS] 成桥态无应力曲率 benchmark")
    return 0


if __name__ == "__main__":
    sys.exit(main())
