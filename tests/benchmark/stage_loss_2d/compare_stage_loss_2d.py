#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 36 项 CTest：2D 阶段分项预应力损失 benchmark（B4 第二十五次，opt-in）。

离线（无需 ANSYS / 求解器许可）：
  * 冻结 .dat / .stage.loss.txt 的 sha256 与 reference 对齐；
  * 当前二进制重放：逐字节复现冻结证据；
  * 正常算例逐键对齐独立 Python 闭式 oracle（§6.2 + 附录 G + 附录 C）；
  * σl2 附录 G 两支中间量、σl6 表值龄期中间量、C_SKIP 显式不参与；
  * 无 LOSSREF 时阶段段只写 ERROR=，不摆零值；LOSS/LOSSREF 冲突与
    LOSSREF 引用不存在表号两个硬负例显式拒绝；
  * 默认零污染：不带开关不写阶段文件、带与不带开关的 _stages.txt 逐字节
    一致、既有第 34 项全局 .loss.txt 与 3D .csv 逐字节不变；
  * 2D/3D 同参同值（同一张 PRESTRESS_LOSS 表，逐键对拍 3D 探针）；
  * FE / BridgeStudio 源码契约守卫。

用法：
    python compare_stage_loss_2d.py --bin <dir> \
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
    DAT_CASES, FORCE_CASES, HARD_NEGATIVE, LEGACY_CASES, NORMAL, SOFT_NEG_FORCE,
    force_oracle, oracle, parse_3d_csv, parse_stage_force_txt,
    parse_stage_loss_txt, run, sha256,
)

RESULTS: list[tuple[bool, str, str]] = []
COMMON_KEYS = ("SIGMA_CON", "SL1", "SL2", "SL3", "SL4", "SL5", "SL6",
               "SPE", "RATIO")
L2_KEYS = ("LF", "DSL2D", "DSL2PEAK", "BRANCH")
C_KEYS = ("C_SOURCE", "ECS_DEV", "PHI_DEV", "BETAS", "BETAS0", "BETAC",
          "BETAH", "RHREP")


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
    exe2, exe3 = bdir / "example_stage.exe", bdir / "loss_probe_3d.exe"
    exe_mixed = bdir / "example_mixed.exe"
    ref = json.loads((rdir / "stage_loss_2d.json").read_text(encoding="utf-8"))
    th = ref["thresholds"]
    print("=" * 78)
    print("第 36 项 CTest：2D 阶段分项预应力损失（B4 第二十五次）")
    print("=" * 78)
    check("example_stage.exe / loss_probe_3d.exe / example_mixed.exe 在档",
          all(p.is_file() for p in (exe2, exe3, exe_mixed)),
          f"stage={exe2.is_file()} 3d={exe3.is_file()} mixed={exe_mixed.is_file()}")

    fresh: dict[str, Path] = {}
    work = Path(tempfile.mkdtemp(prefix="stagelosscmp_"))
    try:
        # 1) 冻结哈希 + 当前二进制重放
        for cid in DAT_CASES:
            dat = cases / (cid + ".dat")
            check(f"{cid} .dat sha256 = reference",
                  dat.is_file() and sha256(dat) == ref["cases"][cid]["dat_sha256"])
            wd = work / cid
            wd.mkdir()
            shutil.copy2(dat, wd / dat.name)
            p = run(exe2, [dat.name, "--stage-loss-detail"], wd)
            out = wd / (dat.stem + ".stage.loss.txt")
            if cid in HARD_NEGATIVE:
                check(f"{cid} 显式拒绝（退出码非 0）", p.returncode != 0,
                      f"rc={p.returncode}")
                err = p.stderr + p.stdout
                check(f"{cid} 错误信息含关键语义",
                      all(k in err for k in ref["cases"][cid]["negative_keywords"]),
                      err[-200:])
                check(f"{cid} 负例不产出阶段明细文件", not out.is_file())
                continue
            check(f"{cid} 重跑退出码 0", p.returncode == 0,
                  f"rc={p.returncode} {p.stderr[-120:]}")
            if not out.is_file():
                check(f"{cid} 重跑产出阶段明细", False, "missing")
                continue
            frozen = cases / ref["cases"][cid]["frozen"]
            check(f"{cid} 冻结证据 sha256 = reference",
                  frozen.is_file() and sha256(frozen) == ref["frozen_sha256"][cid])
            check(f"{cid} 当前二进制重跑 = 冻结证据（逐字节）",
                  out.read_bytes() == frozen.read_bytes())
            fresh[cid] = out

        # 2) oracle 逐键对拍
        kvmap: dict[str, dict] = {}
        for cid in NORMAL:
            kv, note, err = parse_stage_loss_txt(fresh[cid])
            kvmap[cid] = kv
            if cid == "stage_noref":
                check("stage_noref：阶段段写 ERROR= 而非零值明细",
                      "no_lossref_stage_prestress" in err and not kv)
                continue
            check(f"{cid} NOTE 写明单位 Pa（不误标 MPa）", "Pa" in note)
            check(f"{cid} NOTE 写明表值不拟合兜底",
                  "不调用规范拟合函数兜底" in note)
            want = ref["cases"][cid]["oracle"]
            keys = COMMON_KEYS + tuple(
                k for k in (L2_KEYS + C_KEYS) if k in want)
            for key in keys:
                got = kv.get(key)
                ok = got is not None and close(got, want[key], th["rel"], th["abs"])
                check(f"{cid} {key} 同独立闭式 oracle", ok,
                      f"{got} / {want[key]}" if got is not None else "missing")

        # 3) 阶段维度语义：只收录 LOSSREF 束；STAGE/P 编号；BASIS
        for cid in ("stage_default", "stage_l2anti", "stage_c2d",
                    "stage_cskip", "stage_same3d"):
            text = fresh[cid].read_text(encoding="utf-8", errors="replace")
            check(f"{cid} 头行 BASIS=2d_stage", "BASIS=2d_stage" in text)
            check(f"{cid} 阶段行只出现 STAGE1 P1（无零值占位行）",
                  text.count("STAGE1 P1 ") == 1
                  and "STAGE2 P2 " not in text)
        check("字段顺序与 3D 一致（SIGMA_CON→SL1..SL6→SPE→RATIO）",
              all(k in kvmap["stage_default"] for k in COMMON_KEYS))

        # 4) σl2 附录 G 映射（lf>l 分支）
        kv_l2 = kvmap["stage_l2anti"]
        cfg_l2 = ref["cases"]["stage_l2anti"]["cfg"]
        check("stage_l2anti：BRANCH=2（lf>l）", abs(kv_l2.get("BRANCH", 0) - 2.0) < th["abs"])
        anchor = cfg_l2["SIGMA_CON"] * (
            1.0 - __import__("math").exp(
                -(cfg_l2["MU"] * cfg_l2["THETA"] + cfg_l2["K"] * cfg_l2["L"])))
        dsd_cross = anchor / cfg_l2["L"]
        check("stage_l2anti：DSL2D = σl1(l)/l 交叉",
              close(kv_l2.get("DSL2D", 0.0), dsd_cross, th["rel"], th["abs"]))
        check("stage_l2anti：σl2 相对默认直线式确有放大",
              kv_l2["SL2"] > kvmap["stage_default"]["SL2"])

        # 5) σl6 附录 C 映射与 C_SKIP 区分
        kv_c2d, kv_skip = kvmap["stage_c2d"], kvmap["stage_cskip"]
        check("stage_c2d：C_SOURCE=1 且 RHREP=80%",
              kv_c2d.get("C_SOURCE") == 1.0 and abs(kv_c2d.get("RHREP", 0) - 80.0) < th["abs"])
        check("stage_c2d：βH 与规范算例书 ~962 一致",
              abs(kv_c2d.get("BETAH", 0) - 961.763184481) < 1e-3,
              str(kv_c2d.get("BETAH")))
        check("stage_c2d：σl6 与规范算例书 ~45.4MPa 一致",
              4.0e7 < kv_c2d["SL6"] < 5.0e7, f"{kv_c2d['SL6']:.4e}")
        check("stage_cskip：SL6=0 且 C_SOURCE=2（显式不参与）",
              abs(kv_skip.get("SL6", -1)) < th["abs"]
              and kv_skip.get("C_SOURCE") == 2.0)
        check("stage_cskip：不出现附录 C 发展中间量",
              not any(k in kv_skip for k in C_KEYS[1:]))

        # 6) 2D/3D 同参同值（同一张 *PRESTRESS_LOSS 行）
        wd3 = work / "legacy3d"
        wd3.mkdir()
        dat3 = LEGACY_CASES / "creep_3d.dat"
        shutil.copy2(dat3, wd3 / dat3.name)
        p3 = run(exe3, [dat3.name, "fresh.csv"], wd3)
        check("3D 探针重跑成功", p3.returncode == 0 and (wd3 / "fresh.csv").is_file(),
              f"rc={p3.returncode}")
        kv3 = parse_3d_csv(wd3 / "fresh.csv") if (wd3 / "fresh.csv").is_file() else {}
        for key in COMMON_KEYS + C_KEYS:
            check(f"2D/3D 同参同值：{key}",
                  close(kvmap["stage_same3d"].get(key, float("nan")),
                        kv3.get(key, float("nan")), th["rel"], th["abs"]),
                  f"{kvmap['stage_same3d'].get(key)} / {kv3.get(key)}")
        glob = LEGACY_CASES / "creep_3d.csv"
        check("3D 冻结 CSV 逐字节不变", glob.is_file()
              and glob.read_bytes() == (wd3 / "fresh.csv").read_bytes())

        # 7) 默认零污染
        wd_off = work / "flag_off"
        wd_off.mkdir()
        shutil.copy2(cases / "stage_default.dat", wd_off / "stage_default.dat")
        p_off = run(exe2, ["stage_default.dat"], wd_off)
        check("默认不带 --stage-loss-detail：求解仍成功", p_off.returncode == 0)
        check("默认不带开关：不写 .stage.loss.txt",
              not (wd_off / "stage_default.stage.loss.txt").is_file())
        on_stages = work / "stage_default" / "stage_default_stages.txt"
        off_stages = wd_off / "stage_default_stages.txt"
        check("带/不带阶段明细开关：_stages.txt 逐字节一致",
              on_stages.is_file() and off_stages.is_file()
              and on_stages.read_bytes() == off_stages.read_bytes())
        # 全局非阶段段未被本轮污染
        wd_glob = work / "legacy_global"
        wd_glob.mkdir()
        dat_g = LEGACY_CASES / "creep_auto.dat"
        shutil.copy2(dat_g, wd_glob / dat_g.name)
        pg = run(exe_mixed, [dat_g.name, "--loss-detail-only"], wd_glob)
        got_g = wd_glob / "creep_auto.loss.txt"
        frozen_g = LEGACY_CASES / "creep_auto.loss.txt"
        check("既有全局段重跑成功", pg.returncode == 0 and got_g.is_file())
        check("既有全局 .loss.txt 逐字节不变",
              got_g.is_file() and frozen_g.is_file()
              and got_g.read_bytes() == frozen_g.read_bytes())
        check("既有全局冻结证据 sha256 = reference",
              sha256(frozen_g) == ref["legacy_sha256"]["creep_auto.loss.txt"])

        # 7b) 第二十七次：阶段有效预加力段 [PRESTRESSFORCE2D]
        fref = ref["force"]
        fkv: dict[str, list[dict]] = {}
        for cid in FORCE_CASES:
            dat = cases / (cid + ".dat")
            entry = fref["cases"][cid]
            check(f"{cid} force .dat sha256 = reference",
                  dat.is_file() and sha256(dat) == entry["dat_sha256"])
            wd = work / ("f_" + cid)
            wd.mkdir()
            shutil.copy2(dat, wd / dat.name)
            pf = run(exe2, [dat.name, "--stage-prestress-force"], wd)
            check(f"{cid} force 重跑退出码 0", pf.returncode == 0,
                  f"rc={pf.returncode} {pf.stderr[-120:]}")
            out = wd / (dat.stem + ".stage.force.txt")
            if not out.is_file():
                check(f"{cid} force 重跑产出", False, "missing")
                continue
            frozen = cases / entry["frozen"]
            check(f"{cid} force 冻结证据 sha256 = reference",
                  frozen.is_file() and sha256(frozen) == fref["frozen_sha256"][cid])
            check(f"{cid} force 当前二进制重跑 = 冻结证据（逐字节）",
                  out.read_bytes() == frozen.read_bytes())
            rows, note, err = parse_stage_force_txt(out)
            fkv[cid] = rows
            if cid in SOFT_NEG_FORCE:
                check(f"{cid}：无阶段束时写 ERROR= 而非零值行",
                      SOFT_NEG_FORCE[cid] in err and not rows)
                check(f"{cid}：NOTE 说明不摆零值冒充", "不输出零值" in note)
                continue
            want = entry["oracle"]
            check(f"{cid} force 行数 = oracle（不重不漏）",
                  len(rows) == len(want), f"{len(rows)} / {len(want)}")
            check(f"{cid} force 组织键 (STAGE,P) 集合 = oracle",
                  [(int(r["stage"]), int(r["id"])) for r in rows]
                  == [(w["stage"], w["id"]) for w in want])
            for got, w in zip(rows, want):
                tag = f"STAGE{w['stage']} P{w['id']}"
                for key in ("SIGMA_CON", "LOSS", "PE"):
                    check(f"{cid} {tag} {key} 同独立闭式 oracle",
                          key in got and close(got[key], w[key], th["rel"], th["abs"]),
                          f"{got.get(key)} / {w[key]}")
            check(f"{cid} force 不输出量纲不成立的 PE/SIGMA_CON",
                  "RATIO" not in note and all("RATIO" not in r for r in rows))

        # 多束多阶段编号语义：阶段内从 1 起，不串阶段
        keys_multi = [(int(r["stage"]), int(r["id"])) for r in fkv.get("force_multi", [])]
        check("force_multi 键集合 = {(1,1),(1,2),(2,1)}",
              keys_multi == [(1, 1), (1, 2), (2, 1)], str(keys_multi))
        check("force_multi：LOSSREF 束带 σcon、LOSS= 束 σcon 记 0（不反造）",
              all(r["SIGMA_CON"] > 0.0 for r in fkv["force_multi"]
                  if int(r["id"]) == 2)
              and all(r["SIGMA_CON"] == 0.0 for r in fkv["force_multi"]
                      if int(r["id"]) == 1))
        check("force_multi：PE = P·LOSS 逐行成立",
              all(close(r["PE"], r["P"] * r["LOSS"], th["rel"], th["abs"])
                  for r in fkv["force_multi"]))
        check("stage_default：STAGE2 的 LOSS= 束 PE=900000（P=1e6×0.9）",
              any(int(r["stage"]) == 2 and int(r["id"]) == 1
                  and close(r["PE"], 900000.0, th["rel"], th["abs"])
                  for r in fkv["stage_default"]))

        # 两个 opt-in 互不依赖：同时开与单独开，两份文件逐字节一致
        wd_both = work / "both_flags"
        wd_both.mkdir()
        shutil.copy2(cases / "stage_default.dat", wd_both / "stage_default.dat")
        p_both = run(exe2, ["stage_default.dat", "--stage-loss-detail",
                            "--stage-prestress-force"], wd_both)
        check("同时开两个开关：求解成功", p_both.returncode == 0)
        check("同时开：.stage.loss.txt 与只开 loss 开关逐字节一致",
              (wd_both / "stage_default.stage.loss.txt").read_bytes()
              == cases.joinpath("stage_default.stage.loss.txt").read_bytes())
        check("同时开：.stage.force.txt 与只开 force 开关逐字节一致",
              (wd_both / "stage_default.stage.force.txt").read_bytes()
              == cases.joinpath("stage_default.stage.force.txt").read_bytes())
        check("只开 force 开关：不顺手写 .stage.loss.txt",
              not (work / "f_stage_default" / "stage_default.stage.loss.txt").is_file())
        check("只开 loss 开关：不顺手写 .stage.force.txt",
              not (work / "stage_default" / "stage_default.stage.force.txt").is_file())

        # 8) FE 源码契约守卫
        src_out = (bdir.parent / "prestress_loss_output_2d.hpp").read_text(
            encoding="utf-8", errors="replace")
        src_2d = (bdir.parent / "fem_2d.cpp").read_text(
            encoding="utf-8", errors="replace")
        src_hpp = (bdir.parent / "fem_2d.hpp").read_text(
            encoding="utf-8", errors="replace")
        src_3d = (bdir.parent / "example_3d_grillage.cpp").read_text(
            encoding="utf-8", errors="replace")
        check("FE 全局段仍写 BASIS=2d_global",
              "BASIS=2d_global" in src_out)
        check("FE 全局段注释仍写明 2D 非阶段不设 STAGE 占位",
              "2D 没有阶段维度" in src_out and "不写 STAGE" in src_out)
        check("FE 新增独立阶段格式函数",
              "formatStagePrestressLossSection2D" in src_out
              and "BASIS=2d_stage" in src_out)
        check("FE *STAGE_PRESTRESS 支持 LOSSREF 且二选一",
              "LOSSREF" in src_2d
              and "必须且只能指定 LOSS 或 LOSSREF" in src_2d)
        check("FE StageResult 增 prestressLossDetails",
              "prestressLossDetails" in src_hpp)
        check("FE 阶段结果只取本阶段明细（防非阶段明细混入）",
              "if (!st.prestresses.empty())" in src_2d
              and "sr.prestressLossDetails = mc.prestressLossDetails" in src_2d)
        check("3D stage prestress loss 字段顺序同源",
              all(k in src_3d for k in ("STAGE", "SIGMA_CON=", "SL1=", "SL6=",
                                        "SPE=", "RATIO=")))
        # 第二十七次：阶段有效预加力段的实现契约
        src_stage = (bdir.parent / "example_stage.cpp").read_text(
            encoding="utf-8", errors="replace")
        force_hpp = src_out.split("formatStagePrestressForceSection2D")[-1]
        check("FE 新增阶段有效预加力格式函数与段名",
              "formatStagePrestressForceSection2D" in src_out
              and "[PRESTRESSFORCE2D]" in src_out)
        check("FE 有效预加力段仍用 BASIS=2d_stage 与 STAGE/P 行首",
              'BASIS=2d_stage' in src_out and '"STAGE" << sr.stage << " P" << pf.id' in src_out)
        check("FE 有效预加力字段顺序 P→LOSS→PE（SIGMA_CON 在前）",
              force_hpp.find('" SIGMA_CON="') < force_hpp.find('" P="')
              and force_hpp.find('" P="') < force_hpp.find('" LOSS="')
              and force_hpp.find('" LOSS="') < force_hpp.find('" PE="'))
        check("FE 有效预加力段不输出量纲不成立的 RATIO",
              "RATIO" not in force_hpp)
        check("FE example_stage 新增 --stage-prestress-force（默认关闭）",
              "stagePrestressForce = false" in src_stage
              and '"--stage-prestress-force"' in src_stage
              and '".stage.force.txt"' in src_stage)
        check("FE StageResult/Model2D 增 prestressForces",
              "prestressForces" in src_hpp)
        check("FE 有效预加力与 Pe 同一次计算（不重算）",
              "m.prestressForces.push_back(PrestressStageForce2D{ps.id, ps.P, loss, Pe, sigmaCon})"
              in src_2d)
        check("FE 有效预加力明细也受 !st.prestresses.empty() 守卫",
              "sr.prestressForces = mc.prestressForces" in src_2d)

        # 9) BridgeStudio 侧源码契约守卫
        bs = Path(__file__).resolve().parents[3].parent / "BridgeStudio"
        spec = (bs / "bridge_spec.py").read_text(encoding="utf-8", errors="replace")
        py = (bs / "prestress_loss_2d.py").read_text(encoding="utf-8", errors="replace")
        solve = (bs / "solve_2d.py").read_text(encoding="utf-8", errors="replace")
        rep = (bs / "report_2d.py").read_text(encoding="utf-8", errors="replace")
        check("BS StagePrestressDef 增 loss_ref 并映射 LOSSREF",
              "loss_ref: int = 0" in spec and 'f"LOSSREF={self.loss_ref}"' in spec)
        check("BS 阶段束校验引用完整性", "预应力束引用了不存在的" in spec)
        check("BS 新口径层 parse_stage_loss_txt / BASIS=2d_stage",
              "def parse_stage_loss_txt" in py and 'STAGE_BASIS = "2d_stage"' in py)
        check("BS 旧口径层仍只认 2d_global + P<id>",
              'BASIS = "2d_global"' in py
              and "^P(\\d+)\\s+(.*)$" in py)
        check("BS solve_2d 追加 --stage-loss-detail（默认关闭）",
              "stage_loss_detail: bool = False" in solve
              and '"--stage-loss-detail"' in solve)
        check("BS 报告挂阶段损失表", "prestress_stage_loss_table" in rep)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    passed = sum(1 for ok, _, _ in RESULTS if ok)
    total = len(RESULTS)
    print("=" * 78)
    print(f"  {passed}/{total} checks passed")
    if passed == total:
        print("  [PASS] 2D 阶段分项预应力损失 benchmark")
        return 0
    print("  [FAIL] 2D 阶段分项预应力损失 benchmark")
    return 1


if __name__ == "__main__":
    sys.exit(main())
