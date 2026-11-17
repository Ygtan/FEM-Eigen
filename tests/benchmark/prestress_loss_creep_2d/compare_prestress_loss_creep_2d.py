#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 34 项 CTest：σl6 收缩徐变附录 C 显式表值口径（B4 第二十三次，opt-in）。

离线（无需 ANSYS / 求解器许可）：
  * 冻结 .dat / .loss.txt / 3D csv 的 sha256 与 reference 对齐；
  * 当前二进制重放：逐字节复现冻结证据；
  * 每条正常算例逐键对齐独立 Python 闭式 oracle（§6.2 + 附录 C）；
  * εcs(t,t0)/φ(t,t0) 发展公式、βH 区间代表值（55/80）、βs/βc 闭式交叉；
  * C50 折减由用户表值承担（程序不拟合、不查表）的语义；
  * 缺 EPS0 / RH 越界 / C2D+C_SKIP 冲突 / 龄期倒挂四个负例必须显式拒绝；
  * 默认路径零污染（ECS/PHI 直给与仓库既有验证算例逐字节一致，
    且不出现 C_SOURCE=/ECS_DEV=/PHI_DEV= 等新字段）；
  * C_SKIP=1 显式不参与：σl6=0、C_SOURCE=2，与「缺值报错」严格区分；
  * 同一组 *PRESTRESS_LOSS 参数 2D 明细段与 3D 探针同值（含附录 C 中间量）。

用法：
    python compare_prestress_loss_creep_2d.py --bin <dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from build_reference import (  # noqa: E402
    DAT_CASES, NEGATIVE, NORMAL, oracle, parse_3d_csv, parse_loss_txt,
    run, sha256,
)

RESULTS: list[tuple[bool, str, str]] = []
COMMON_KEYS = ("SIGMA_CON", "SL1", "SL2", "SL3", "SL4", "SL5", "SL6",
               "SPE", "RATIO")
CREEP_KEYS = ("C_SOURCE", "ECS_DEV", "PHI_DEV", "BETAS", "BETAS0",
              "BETAC", "BETAH", "RHREP")


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((bool(ok), name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def close(a: float, b: float, rel: float, abs_: float) -> bool:
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def data_lines(path: Path) -> list[str]:
    return [ln for ln in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if ln.strip().startswith("P") and " SIGMA_CON=" in ln]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    a = ap.parse_args()

    bdir = Path(a.bin).resolve()
    cases = Path(a.cases).resolve()
    rdir = Path(a.reference).resolve()
    exe2 = bdir / "example_mixed.exe"
    exe3 = bdir / "loss_probe_3d.exe"
    ref = json.loads((rdir / "prestress_loss_creep_2d.json").read_text(encoding="utf-8"))
    th = ref["thresholds"]
    print("=" * 78)
    print("第 34 项 CTest：σl6 收缩徐变附录 C 显式表值（B4 第二十三次）")
    print("=" * 78)
    check("example_mixed.exe / loss_probe_3d.exe 在档",
          exe2.is_file() and exe3.is_file(),
          f"2d={exe2.is_file()} 3d={exe3.is_file()}")

    fresh: dict[str, Path] = {}
    work = Path(tempfile.mkdtemp(prefix="creepcmp_"))
    try:
        # 1) 冻结哈希 + 当前二进制重放
        for cid in DAT_CASES:
            dat = cases / (cid + ".dat")
            check(f"{cid} .dat sha256 = reference",
                  dat.is_file() and sha256(dat) == ref["cases"][cid]["dat_sha256"])
            if cid in NEGATIVE:
                kw = ref["cases"][cid]["negative_keywords"]
                wd = work / (cid + "_neg")
                wd.mkdir()
                shutil.copy2(dat, wd / dat.name)
                out = wd / (dat.stem + ".loss.txt")
                p = run(exe2, [dat.name, "--loss-detail-only"], wd)
                check(f"{cid} 显式拒绝（退出码非 0）", p.returncode != 0,
                      f"rc={p.returncode}")
                # 错误信息在 stderr；内联到 stdout 兼容不同捕获顺序
                err = p.stderr + p.stdout
                check(f"{cid} 错误信息含关键语义 {kw}",
                      all(k in err for k in kw), err[-200:])
                check(f"{cid} 负例不产出明细文件", not out.is_file())
                continue
            check(f"{cid} 重跑退出码 0", True)  # 下面实际重放会验证
            wd = work / cid
            wd.mkdir()
            shutil.copy2(dat, wd / dat.name)
            if cid == "creep_3d":
                p = run(exe3, [dat.name, "fresh.csv"], wd)
                out = wd / "fresh.csv"
            else:
                p = run(exe2, [dat.name, "--loss-detail-only"], wd)
                out = wd / (dat.stem + ".loss.txt")
            check(f"{cid} 重跑退出码 0", p.returncode == 0, f"rc={p.returncode} {p.stderr[-120:]}")
            if not out.is_file():
                check(f"{cid} 重跑产出明细", False, "missing")
                continue
            frozen = cases / ref["cases"][cid]["frozen"]
            check(f"{cid} 冻结证据 sha256 = reference",
                  frozen.is_file() and sha256(frozen) == ref["frozen_sha256"][cid])
            check(f"{cid} 当前二进制重跑 = 冻结证据（逐字节）",
                  out.read_bytes() == frozen.read_bytes())
            fresh[cid] = out

        # 2) 默认零污染：与仓库既有 verify_prestress_loss.dat 重放逐字节一致
        legacy = bdir.parent / "verify_prestress_loss.dat"
        if not legacy.is_file():
            legacy = Path(__file__).resolve().parents[2] / "verify_prestress_loss.dat"
        if legacy.is_file() and "creep_default" in fresh:
            sub = work / "legacy"
            sub.mkdir()
            shutil.copy2(legacy, sub / legacy.name)
            p = run(exe2, [legacy.name, "--loss-detail-only"], sub)
            got = sub / (legacy.stem + ".loss.txt")
            check("默认路径：既有 verify_prestress_loss.dat 重跑成功",
                  p.returncode == 0 and got.is_file(), f"rc={p.returncode}")
            if got.is_file():
                check("默认路径：明细段逐字节一致（零污染）",
                      got.read_bytes() == fresh["creep_default"].read_bytes())
        else:
            check("默认路径：既有验证算例在档", False, str(legacy))

        # 3) 解析 / oracle 逐键对拍
        kvmap: dict[str, dict] = {}
        for cid in NORMAL:
            if cid == "creep_3d":
                kv = parse_3d_csv(fresh[cid])
            else:
                kv, note = parse_loss_txt(fresh[cid])
                check(f"{cid} NOTE 写明单位 Pa（不误标 MPa）", "Pa" in note)
                check(f"{cid} NOTE 如实写明 σl6 表值不拟合兜底",
                      "不调用规范拟合函数兜底" in note)
            kvmap[cid] = kv
            want = ref["cases"][cid]["oracle"]
            keys = COMMON_KEYS + tuple(k for k in CREEP_KEYS if k in want)
            for key in keys:
                if key not in want:
                    continue
                got = kv.get(key)
                ok = got is not None and close(got, want[key], th["rel"], th["abs"])
                check(f"{cid} {key} 同独立闭式 oracle", ok,
                      f"{got} / {want[key]}" if got is not None else "missing")

        # 4) 默认路径不出现新字段（数据行）
        kv_def, _ = parse_loss_txt(fresh["creep_default"])
        check("默认路径：数据行不出现附录 C 新字段",
              not any(k in kv_def for k in CREEP_KEYS))
        check("默认路径：C_SOURCE=0（用户直给）",
              kv_def.get("C_SOURCE", 0) == 0.0 or "C_SOURCE" not in kv_def)

        # 5) C_SKIP 与「缺值报错」严格区分
        kv_skip, _ = parse_loss_txt(fresh["creep_skip"])
        check("C_SKIP：σl6=0 且 C_SOURCE=2（显式不参与）",
              abs(kv_skip.get("SL6", -1)) < th["abs"] and kv_skip.get("C_SOURCE") == 2.0,
              f"SL6={kv_skip.get('SL6')} C_SOURCE={kv_skip.get('C_SOURCE')}")
        check("C_SKIP：不出现附录 C 发展中间量",
              not any(k in kv_skip for k in ("ECS_DEV", "PHI_DEV", "BETAH")))

        # 6) 附录 C 发展公式交叉 + βH 区间代表值切换
        kv_auto = kvmap["creep_auto"]
        check("creep_auto：βH 代表值 80%（RH=75% 属 70~99%）",
              abs(kv_auto.get("RHREP", 0) - 80.0) < th["abs"])
        kv_rh55 = kvmap["creep_auto_rh55"]
        check("creep_auto_rh55：βH 代表值 55%（RH=55% 属 40~70%）",
              abs(kv_rh55.get("RHREP", 0) - 55.0) < th["abs"])
        check("两区间 βH 不同（代表值切换生效）",
              abs(kv_auto["BETAH"] - kv_rh55["BETAH"]) > 1.0,
              f"{kv_auto['BETAH']:.2f} vs {kv_rh55['BETAH']:.2f}")
        # φ = φ0·βc 交叉（用 cfg 表值）
        cfg_auto = ref["cases"]["creep_auto"]["cfg"]
        phi_cross = cfg_auto["PHI0"] * kv_auto["BETAC"]
        check("creep_auto：φ_dev = φ0·βc 交叉",
              close(kv_auto["PHI_DEV"], phi_cross, th["rel"], th["abs"]))
        eps_cross = cfg_auto["EPS0"] * (kv_auto["BETAS"] - kv_auto["BETAS0"])
        check("creep_auto：εcs_dev = εcs0(βs−βs0) 交叉",
              close(kv_auto["ECS_DEV"], eps_cross, th["rel"], th["abs"]))
        # 独立量纲锚点：规范算例书 σl6≈45MPa
        check("creep_auto：σl6 量级与规范算例书 ~45MPa 一致",
              4.0e7 < kv_auto["SL6"] < 5.0e7, f"{kv_auto['SL6']:.4e}")

        # 7) 2D/3D 同参同值（含附录 C 全部中间量）
        for key in ("SL6", "C_SOURCE", "ECS_DEV", "PHI_DEV", "BETAS",
                    "BETAS0", "BETAC", "BETAH", "RHREP"):
            check(f"2D/3D 同参同值：{key}",
                  close(kv_auto.get(key, float("nan")),
                        kvmap["creep_3d"].get(key, float("nan")),
                        th["rel"], th["abs"]),
                  f"{kv_auto.get(key)} / {kvmap['creep_3d'].get(key)}")

        # 8) 探针只读内核（无自带力学）
        src = (bdir.parent / "loss_probe_3d.cpp").read_text(encoding="utf-8", errors="replace")
        check("loss_probe_3d 只调用 readModel3D + resolvePrestressLosses3D",
              "readModel3D" in src and "resolvePrestressLosses3D" in src
              and "jtg3362NominalCreep" not in src)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    passed = sum(1 for ok, _, _ in RESULTS if ok)
    total = len(RESULTS)
    print("=" * 78)
    print(f"  {passed}/{total} checks passed")
    if passed == total:
        print("  [PASS] σl6 附录 C 收缩徐变 benchmark")
        return 0
    print("  [FAIL] σl6 附录 C 收缩徐变 benchmark")
    return 1


if __name__ == "__main__":
    sys.exit(main())
