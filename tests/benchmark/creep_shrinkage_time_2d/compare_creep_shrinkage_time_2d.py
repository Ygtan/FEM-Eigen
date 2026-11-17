#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 35 项 CTest：结构时变徐变/收缩 JTG 3362 规范算例书锚点（B4 第二十四次）。

离线（无需 ANSYS / 求解器许可）：
  * 冻结 .dat / .csv 的 sha256 与 reference 对齐；
  * 当前二进制重放：逐字节复现冻结证据；
  * 每条正常算例逐键对齐独立 Python 闭式 oracle：
    φ0=φRH·βF·βT0、βH 连续 RH、εcs0=εS·βRH（附录 C.2.1 / C.1）；
  * 逐龄期 φ(t)=φ0·[(t−t0)/(βH+t−t0)]^0.3、εcs(t)=εcs0(βs−βs0)、
    βs、AEMM 折减 E_eff=E/(1+χφ) 逐键核对；
  * 与 σl6 口径隔离的硬断言：本链 βH 用实际 RH 连续值，
    **不是** 55/80 代表值（同一 RH=75% 给出 803.25 而非 961.76）；
  * RH 55% / t0 28d 的敏感性：φ0、βH、βT0 单调变化；
  * 四个负例（缺 RH/H/FCM、旧 CURVE 口径）必须显式拒绝；
  * 探针只读内核（不自带力学，不调用拟合函数之外的求解逻辑）。

用法：
    python compare_creep_shrinkage_time_2d.py --bin <dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from build_reference import (  # noqa: E402
    DAT_CASES, NEGATIVE, NORMAL, parse_probe_csv, run, sha256,
)

RESULTS: list[tuple[bool, str, str]] = []
M_KEYS = ("PHI0", "BETAH", "ECS0")
S_KEYS = ("PHI", "ECS", "BETAS", "EEFF")


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((bool(ok), name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def close(a: float, b: float, rel: float, abs_: float) -> bool:
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    a = ap.parse_args()

    bdir = Path(a.bin).resolve()
    cases = Path(a.cases).resolve()
    rdir = Path(a.reference).resolve()
    exe = bdir / "creep_probe_2d.exe"
    ref = json.loads((rdir / "creep_shrinkage_time_2d.json").read_text(encoding="utf-8"))
    th = ref["thresholds"]
    print("=" * 78)
    print("第 35 项 CTest：结构时变徐变 JTG 3362 规范算例书锚点（B4 第二十四次）")
    print("=" * 78)
    check("creep_probe_2d.exe 在档", exe.is_file(), str(exe))

    fresh: dict[str, Path] = {}
    work = Path(tempfile.mkdtemp(prefix="creepcmp_"))
    try:
        # 1) 冻结哈希 + 当前二进制重放 + 负例拒绝
        for cid in DAT_CASES:
            dat = cases / (cid + ".dat")
            check(f"{cid} .dat sha256 = reference",
                  dat.is_file() and sha256(dat) == ref["cases"][cid]["dat_sha256"])
            wd = work / cid
            wd.mkdir()
            shutil.copy2(dat, wd / dat.name)
            out = wd / "fresh.csv"
            p = run(exe, [dat.name, "fresh.csv"], wd)
            if cid in NEGATIVE:
                kw = ref["cases"][cid]["negative_keyword"]
                check(f"{cid} 显式拒绝（退出码非 0）", p.returncode != 0,
                      f"rc={p.returncode}")
                err = p.stderr + p.stdout
                check(f"{cid} 错误信息含关键语义「{kw}」", kw in err, err[-160:])
                check(f"{cid} 负例不产出数值 CSV", not out.is_file() or out.stat().st_size == 0)
                continue
            check(f"{cid} 重跑退出码 0", p.returncode == 0,
                  f"rc={p.returncode} {p.stderr[-120:]}")
            if not out.is_file():
                check(f"{cid} 重跑产出 CSV", False, "missing")
                continue
            frozen = cases / ref["cases"][cid]["frozen"]
            check(f"{cid} 冻结证据 sha256 = reference",
                  frozen.is_file() and sha256(frozen) == ref["frozen_sha256"][cid])
            check(f"{cid} 当前二进制重跑 = 冻结证据（逐字节）",
                  out.read_bytes() == frozen.read_bytes())
            fresh[cid] = out

        # 2) 逐键对齐独立闭式 oracle
        kvmap: dict[str, dict] = {}
        for cid in NORMAL:
            kv = parse_probe_csv(fresh[cid])
            kvmap[cid] = kv
            want = ref["cases"][cid]["oracle"]
            for key in M_KEYS:
                got = kv["M"].get(key)
                ok = got is not None and close(got, want["M"][key], th["rel"], th["abs"])
                check(f"{cid} M/{key} 同独立闭式 oracle", ok,
                      f"{got} / {want['M'][key]}" if got is not None else "missing")
            for age, wv in want["S"].items():
                for key in S_KEYS:
                    got = kv["S"].get(age, {}).get(key)
                    ok = got is not None and close(got, wv[key], th["rel"], th["abs"])
                    check(f"{cid} S/{age}/{key} 同独立闭式 oracle", ok,
                          f"{got} / {wv[key]}" if got is not None else "missing")

        # 3) 口径隔离硬断言：本链 βH 用实际 RH=75% 连续值（803.25），
        #    不是 σl6 代表值 80% 的 961.76；εcs0 用连续 βRH，不是表值 0.310e-3。
        m_anchor = kvmap["creep_anchor"]["M"]
        check("结构时变链 βH 用实际 RH=75% 连续值（≈803.25，非 σl6 代表值 961.76）",
              close(m_anchor["BETAH"], 803.253, 1e-4, 1e-3)
              and abs(m_anchor["BETAH"] - 961.763) > 100.0,
              f"BETAH={m_anchor['BETAH']:.4f}")
        check("结构时变链 εcs0 用连续 βRH（≈3.674e-4，非 σl6 表值 0.310e-3）",
              close(m_anchor["ECS0"], 3.673984375e-4, 1e-6, 1e-9)
              and abs(m_anchor["ECS0"] - 0.310e-3) > 1e-5,
              f"ECS0={m_anchor['ECS0']:.6e}")

        # 4) RH / t0 敏感性：φ0、βH、βT0 随参数单调变化
        m_rh55 = kvmap["creep_rh55"]["M"]
        check("RH 55% 的 φ0 高于 RH 75%（干燥越大徐变越大）",
              m_rh55["PHI0"] > m_anchor["PHI0"],
              f"{m_rh55['PHI0']:.4f} > {m_anchor['PHI0']:.4f}")
        check("RH 55% 的 βH 小于 RH 75%",
              m_rh55["BETAH"] < m_anchor["BETAH"],
              f"{m_rh55['BETAH']:.2f} < {m_anchor['BETAH']:.2f}")
        m_late = kvmap["creep_t0late"]["M"]
        oracle_late = ref["cases"]["creep_t0late"]["oracle"]["M"]
        check("晚加载 t0=28d 的 βT0 小于 t0=7d（加载越晚徐变越小）",
              oracle_late["betaT0"] < ref["cases"]["creep_anchor"]["oracle"]["M"]["betaT0"],
              f"{oracle_late['betaT0']:.4f} < "
              f"{ref['cases']['creep_anchor']['oracle']['M']['betaT0']:.4f}")
        check("晚加载 t0=28d 的 φ0 低于 t0=7d",
              m_late["PHI0"] < m_anchor["PHI0"],
              f"{m_late['PHI0']:.4f} < {m_anchor['PHI0']:.4f}")

        # 5) AEMM 折减单调性：龄期越大 φ 越大、E_eff 越低
        s28 = kvmap["creep_anchor"]["S"]["28"]
        s97 = kvmap["creep_anchor"]["S"]["97"]
        check("φ(97) > φ(28)（徐变随龄期发展）",
              s97["PHI"] > s28["PHI"], f"{s97['PHI']:.4f} > {s28['PHI']:.4f}")
        check("E_eff(97) < E_eff(28)（AEMM 折减随龄期加深）",
              s97["EEFF"] < s28["EEFF"],
              f"{s97['EEFF']:.3e} < {s28['EEFF']:.3e}")
        check("E_eff = E/(1+0.8φ) 交叉（χ=0.8）",
              close(s97["EEFF"], 200.0e9 / (1.0 + 0.8 * s97["PHI"]), th["rel"], th["abs"]),
              f"{s97['EEFF']:.6e}")

        # 6) 探针只读内核（无自带力学、不实现求解器）
        src = (bdir.parent / "creep_probe_2d.cpp").read_text(encoding="utf-8", errors="replace")
        check("creep_probe_2d 只调 readModel2D + code_loads 拟合函数",
              "readModel2D" in src and "creepCoefficientAt" in src
              and "shrinkageStrainAt" in src and "assemble" not in src)
        check("探针无硬编码代表值 55/80（口径隔离）",
              "55.0" not in src and "80.0" not in src)

        # 6b) h 来源语义：*CREEP_MODEL H 是理论厚度 h=2A/u（非梁高），
        #     给定 H 时内核 m.creepH 覆盖单元 hTheory，探针与内核同源。
        #     防止未来有人改动覆盖逻辑后，探针(全桥 H)与内核(单元 hTheory)分叉。
        src_core = (bdir.parent / "fem_2d.cpp").read_text(encoding="utf-8", errors="replace")
        check("内核 creepCoefficientAt/shrinkageStrainAt 用 m.creepH 覆盖传入单元 h",
              src_core.count("m.creepH > 0.0 ? m.creepH : hMm") >= 2,
              "给定 H 时探针与 AEMM 同源")
        check("探针 φ0/βH 与 βs 采样均取 m.creepH（理论厚度，非梁高）",
              "m.creepRh, m.creepH, m.creepFcm" in src
              and "jtg3362ShrinkageBeta(age, ts, m.creepH)" in src)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    passed = sum(1 for ok, _, _ in RESULTS if ok)
    total = len(RESULTS)
    print("=" * 78)
    print(f"  {passed}/{total} checks passed")
    if passed == total:
        print("  [PASS] 结构时变徐变 JTG 3362 规范算例书锚点 benchmark")
        return 0
    print("  [FAIL] 结构时变徐变 benchmark")
    return 1


if __name__ == "__main__":
    sys.exit(main())
