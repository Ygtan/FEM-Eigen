#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 33 项 CTest：σl2 附录 G 反摩擦（B4 第二十三次，opt-in）。

离线（无需 ANSYS / 求解器许可）：
  * 冻结 .dat / .loss.txt / 3D csv 的 sha256 与 reference 对齐；
  * 当前二进制重放：逐字节复现冻结证据；
  * 每条正常算例逐键对齐独立 Python 闭式 oracle（§6.2 + 附录 G）；
  * 支判别、两支闭式交叉、lf=l 连续性、单位量纲、NOTE 如实化；
  * 缺 ΣΔl / 非 x=0 代表值两个负例必须显式拒绝；
  * 默认路径零污染（无 L2ANTI 时数值行与仓库既有验证算例逐字节一致，
    且不出现 LF=/DSL2D=/DSL2PEAK=/BRANCH=）；
  * 同一组 *PRESTRESS_LOSS 参数 2D 明细段与 3D 探针 12 键同值。

用法：
    python compare_prestress_loss_anti_2d.py --bin <dir> \
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
ANTI_KEYS = ("LF", "DSL2D", "DSL2PEAK", "BRANCH")


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
    exe2 = bdir / "example_mixed.exe"
    exe3 = bdir / "loss_probe_3d.exe"
    ref = json.loads((rdir / "prestress_loss_anti_2d.json").read_text(encoding="utf-8"))
    th = ref["thresholds"]
    print("=" * 78)
    print("第 33 项 CTest：σl2 附录 G 反摩擦（B4 第二十三次）")
    print("=" * 78)
    check("example_mixed.exe / loss_probe_3d.exe 在档",
          exe2.is_file() and exe3.is_file(),
          f"2d={exe2.is_file()} 3d={exe3.is_file()}")

    for cid in DAT_CASES:
        dat = cases / (cid + ".dat")
        check(f"{cid} .dat sha256 = reference",
              dat.is_file() and sha256(dat) == ref["cases"][cid]["dat_sha256"])

    # 冻结证据哈希守卫
    for cid in NORMAL:
        frozen = cases / ref["cases"][cid]["frozen"]
        check(f"{cid} 冻结证据 sha256 = reference",
              frozen.is_file() and sha256(frozen) == ref["frozen_sha256"][cid])

    fresh: dict[str, Path] = {}
    with tempfile.TemporaryDirectory(prefix="lossanti_cmp_") as td:
        work = Path(td)
        for cid in DAT_CASES:
            dat = cases / (cid + ".dat")
            shutil.copy2(dat, work / dat.name)
            if cid == "anti_3d":
                p = run(exe3, [dat.name, "fresh.csv"], work)
                out = work / "fresh.csv"
            else:
                p = run(exe2, [dat.name, "--loss-detail-only"], work)
                out = work / (dat.stem + ".loss.txt")
            if cid in NEGATIVE:
                kw = ref["cases"][cid]["negative_keywords"]
                err = (p.stderr or "") + (p.stdout or "")
                check(f"{cid} 显式拒绝（退出码非 0）", p.returncode != 0,
                      f"rc={p.returncode}")
                check(f"{cid} 错误信息含 {kw[0]} / {kw[1]}",
                      all(k in err for k in kw), err[-180:])
                check(f"{cid} 负例不产出明细文件", not out.is_file())
                continue
            check(f"{cid} 重跑退出码 0", p.returncode == 0, f"rc={p.returncode}")
            if not out.is_file():
                check(f"{cid} 重跑产出明细", False, "missing")
                continue
            frozen = cases / ref["cases"][cid]["frozen"]
            check(f"{cid} 当前二进制重跑 = 冻结证据（逐字节）",
                  out.read_bytes() == frozen.read_bytes())
            fresh[cid] = out

        # 默认零污染：与仓库既有 verify_prestress_loss.dat 现场生成的段逐字节一致
        legacy = bdir.parent / "verify_prestress_loss.dat"
        if not legacy.is_file():
            legacy = Path(__file__).resolve().parents[2] / "verify_prestress_loss.dat"
        default_out = work / (dat.stem + ".loss.txt")  # 占位，下面用独立子目录
        if legacy.is_file():
            sub = work / "legacy"
            sub.mkdir()
            shutil.copy2(legacy, sub / legacy.name)
            p = run(exe2, [legacy.name, "--loss-detail-only"], sub)
            got = sub / (legacy.stem + ".loss.txt")
            check("默认路径：既有 verify_prestress_loss.dat 重跑成功",
                  p.returncode == 0 and got.is_file(), f"rc={p.returncode}")
            if got.is_file() and "anti_default" in fresh:
                a_bytes = fresh["anti_default"].read_bytes()
                b_bytes = got.read_bytes()
                check("默认路径：明细段逐字节一致（零污染）", a_bytes == b_bytes)
        else:
            check("默认路径：既有验证算例在档", False, str(legacy))

        # 解析 / oracle 逐键对拍
        kvmap: dict[str, dict] = {}
        for cid in NORMAL:
            if cid == "anti_3d":
                kv = parse_3d_csv(fresh[cid])
            else:
                kv, note = parse_loss_txt(fresh[cid])
                check(f"{cid} NOTE 写明单位 Pa（不误标 MPa）", "Pa" in note)
            kvmap[cid] = kv
            want = ref["cases"][cid]["oracle"]
            for key in COMMON_KEYS + (ANTI_KEYS if cid != "anti_default"
                                      and cid != "anti_nofiric" else ()):
                if key not in want:
                    continue
                got = kv.get(key)
                ok = got is not None and close(got, want[key], th["rel"], th["abs"])
                check(f"{cid} {key} 同独立闭式 oracle", ok,
                      f"{got} / {want[key]}" if got is not None else "missing")
            # 不启用 / 无摩擦退化时不得输出反摩擦字段
            if cid in ("anti_default", "anti_nofiric"):
                # 只看 P<id> 数据行：NOTE 里会如实提到 BRANCH= 的含义，
                # 那是文案不是数据字段，不能拿整文件文本误判为字段出现。
                rows = [ln for ln in fresh[cid]
                        .read_text(encoding="utf-8", errors="replace").splitlines()
                        if ln.strip().startswith("P") and " SIGMA_CON=" in ln]
                check(f"{cid} 未启用反摩擦：数据行不追加 LF/DSL2D/DSL2PEAK/BRANCH",
                      bool(rows) and all(not any((k + "=") in ln for k in ANTI_KEYS)
                                         for ln in rows))

        # 支判别与量纲
        cfg_le = ref["cases"]["anti_le"]["cfg"]
        kv_le = kvmap["anti_le"]
        kv_gt = kvmap["anti_gt"]
        check("anti_le 走 BRANCH=1 且 lf<l",
              int(kv_le["BRANCH"]) == 1 and kv_le["LF"] < cfg_le["L"],
              f"lf={kv_le['LF']:.6g}")
        check("anti_gt 走 BRANCH=2 且 lf>l",
              int(kv_gt["BRANCH"]) == 2 and kv_gt["LF"] > cfg_le["L"],
              f"lf={kv_gt['LF']:.6g}")
        for cid, kv in (("anti_le", kv_le), ("anti_gt", kv_gt)):
            check(f"{cid} 单位量纲：σcon/σl2 为 Pa、Δσd 为 Pa/m、lf 为 m",
                  1e8 < kv["SIGMA_CON"] < 1e10
                  and 1e6 < kv["SL2"] < 1e9
                  and 1e5 < kv["DSL2D"] < 1e8
                  and 0.1 < kv["LF"] < 100.0)
            peak = (2.0 * kv["DSL2D"] * kv["LF"] if int(kv["BRANCH"]) == 1
                    else kv["DSL2D"] * (kv["LF"] ** 2 / cfg_le["L"] + cfg_le["L"]))
            check(f"{cid} 分支闭式交叉：DSL2PEAK 与 LF/DSL2D 自洽",
                  close(kv["DSL2PEAK"], peak, th["rel"], th["abs"]),
                  f"{kv['DSL2PEAK']:.6g} / {peak:.6g}")
            check(f"{cid} SL2 = DSL2PEAK（张拉端代表值）",
                  close(kv["SL2"], kv["DSL2PEAK"], th["rel"], th["abs"]))

        # 独立闭式在临界 ΣΔl 上两支连续（完全不读求解器）
        cr = oracle(dict(cfg_le, DELTA_L=0.004502096850693447))
        lf_cr = cr["LF"]
        peak1 = 2.0 * cr["DSL2D"] * lf_cr
        peak2 = cr["DSL2D"] * (lf_cr ** 2 / cfg_le["L"] + cfg_le["L"])
        check("lf=l 两支闭式连续（独立 oracle）",
              close(lf_cr, cfg_le["L"], 1e-9, 1e-9)
              and close(peak1, peak2, 1e-9, 1e-9),
              f"lf={lf_cr:.9f} peak={peak1:.6g}/{peak2:.6g}")

        # 2D / 3D 同参同值
        kv2 = kvmap["anti_le"]
        kv3 = kvmap["anti_3d"]
        for key in COMMON_KEYS + ANTI_KEYS:
            ok = key in kv3 and close(kv2[key], kv3[key], th["rel"], th["abs"])
            check(f"2D/3D 同参同值：{key}", ok,
                  f"{kv2.get(key)} / {kv3.get(key)}")

    # 源码契约：3D 探针不得自带力学，只复用内核入口
    probe = (bdir.parent / "loss_probe_3d.cpp")
    if not probe.is_file():
        probe = Path(__file__).resolve().parents[2] / "loss_probe_3d.cpp"
    src = probe.read_text(encoding="utf-8", errors="replace") if probe.is_file() else ""
    check("loss_probe_3d 只调用 readModel3D + resolvePrestressLosses3D",
          "readModel3D(" in src and "resolvePrestressLosses3D(" in src
          and "computePrestressLosses(" not in src)

    n = sum(1 for ok, _, _ in RESULTS if ok)
    print("\n" + "=" * 78)
    print(f"  {n}/{len(RESULTS)} checks passed")
    print("  [PASS] σl2 附录 G 反摩擦 benchmark"
          if n == len(RESULTS) else "  [FAIL] σl2 附录 G 反摩擦 benchmark")
    print("=" * 78)
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
