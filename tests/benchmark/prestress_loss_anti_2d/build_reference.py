#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二十三次 B4 · σl2 附录 G 反摩擦 benchmark 的独立闭式 oracle 与冻结器。

设计纪律（与前面 32 项 benchmark 一致）：

* oracle 用纯 Python 写第二份 JTG 3362-2018 §6.2 + 附录 G 闭式，
  **不调用任何求解器、不读求解器输出**，只从 cases/*.dat 的
  ``*PRESTRESS_LOSS`` 行读参数（参数来源唯一，避免脚本里手写两份）；
* build_reference.py 只负责用**当前二进制**生成冻结证据 + reference JSON；
  compare 脚本再用同一二进制重放并逐字节对拍，防止冻结后内核被悄悄改动；
* 负例不产生数值行，reference 只记录期望的退出码语义与错误关键词。

用法：
    python build_reference.py --bin <dir_with_example_mixed_and_loss_probe_3d> \
        [--cases <cases_dir>] [--reference <reference_dir>]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MSYS_BIN = r"C:\msys64\mingw64\bin"

# 正常算例（有冻结 .loss.txt）与负例（只验拒绝语义）
NORMAL = ("anti_default", "anti_le", "anti_gt", "anti_nofiric", "anti_3d")
DAT_CASES = ("anti_default", "anti_le", "anti_gt", "anti_nodl",
             "anti_badrepr", "anti_nofiric", "anti_3d")
NEGATIVE = {"anti_nodl": ("附录 G", "ΣΔl"),
            "anti_badrepr": ("l2ReprX", "代表截面")}
FROZEN_EXT = {"anti_3d": ".csv"}
EXE_2D = "example_mixed.exe"
EXE_3D = "loss_probe_3d.exe"

# 容差：同工具链同二进制按 12 位有效数字写，允许相对 1e-11 的 JSON 往返差。
REL_TOL = 1.0e-11
ABS_TOL = 1.0e-9


# ---------------------------------------------------------------- 参数与闭式

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_cfg(dat: Path) -> dict:
    """只读 cases/*.dat 的 *PRESTRESS_LOSS 行，不读任何求解器输出。"""
    text = dat.read_text(encoding="utf-8", errors="replace")
    line = ""
    for ln in text.splitlines():
        if ln.strip().startswith("*PRESTRESS_LOSS"):
            continue
        if re.match(r"^\s*\d+\s+SIGMA_CON=", ln):
            line = ln
            break
    if not line:
        raise RuntimeError("算例缺少 *PRESTRESS_LOSS 参数行: " + str(dat))
    cfg: dict[str, float] = {}
    for k, v in re.findall(r"\b([A-Z_][A-Z0-9_]*)=([-+0-9.eE]+)", line):
        cfg[k] = float(v)
    return cfg


def oracle(cfg: dict) -> dict:
    """独立闭式：JTG 3362-2018 §6.2 + 附录 G 张拉端 x=0 代表值。

    SI 单位：Pa / m / rad。σl2 两支与 `lf=l` 处连续；无摩擦时反摩擦无定义，
    退化为直线式（BRANCH=0）。σl5 用 6.2.6-1，低 σpe 时可能为负，如实保留。
    """
    out: dict[str, float] = {}
    con = cfg["SIGMA_CON"]
    ep = cfg["EP"]
    ec = cfg["EC"]
    alpha = ep / ec if ec > 0.0 else 0.0

    def sl1_at(x: float) -> float:
        return con * (1.0 - math.exp(-(cfg["MU"] * cfg["THETA"] + cfg["K"] * x)))

    out["SIGMA_CON"] = con
    out["SL1"] = sl1_at(cfg["X"])
    out["SL3"] = 2.0 * (cfg["T2"] - cfg["T1"]) * 1.0e6
    out["SL4"] = alpha * cfg["SIGMA_PC4"]

    dsl2d = sl1_at(cfg["L"]) / cfg["L"]
    if int(cfg.get("L2ANTI", 0)) == 0:
        out["SL2"] = cfg["DELTA_L"] / cfg["L"] * ep
        out["LF"], out["DSL2D"], out["DSL2PEAK"], out["BRANCH"] = 0.0, 0.0, 0.0, 0.0
    else:
        if int(cfg.get("L2REPX", 0)) != 0:
            raise ValueError("独立 oracle 只覆盖 x=0 代表截面")
        if cfg["DELTA_L"] <= 0.0 or ep <= 0.0:
            raise ValueError("反摩擦需要 ΣΔl>0 且 Ep>0")
        out["DSL2D"] = dsl2d
        if not (dsl2d > 0.0):
            out["SL2"] = cfg["DELTA_L"] / cfg["L"] * ep
            out["LF"], out["DSL2PEAK"], out["BRANCH"] = 0.0, 0.0, 0.0
        else:
            lf = math.sqrt(cfg["DELTA_L"] * ep / dsl2d)
            out["LF"] = lf
            if lf <= cfg["L"]:
                out["DSL2PEAK"] = 2.0 * dsl2d * lf
                out["BRANCH"] = 1.0
            else:
                out["DSL2PEAK"] = dsl2d * (lf * lf / cfg["L"] + cfg["L"])
                out["BRANCH"] = 2.0
            out["SL2"] = out["DSL2PEAK"]

    out["SPE"] = con - out["SL1"] - out["SL2"] - out["SL4"]
    out["RATIO"] = out["SPE"] / con
    out["SL5"] = cfg["PSI"] * cfg["ZETA"] * (0.52 * out["SPE"] / cfg["FPK"] - 0.26) * out["SPE"]
    out["SL6"] = 0.9 * (ep * cfg["ECS"] + alpha * cfg["SIGMA_PC6"] * cfg["PHI"]) \
        / (1.0 + 15.0 * cfg["RHO"] * cfg["RHOPS"])
    return out


# ---------------------------------------------------------------- 输出解析

def parse_loss_txt(path: Path) -> tuple[dict, str]:
    """解析 [PRESTRESSLOSS2D] 数据行；返回 (kv, note)。不重算任何力学量。"""
    text = path.read_text(encoding="utf-8", errors="replace")
    note = ""
    row = ""
    for ln in text.splitlines():
        if ln.startswith("NOTE="):
            note = ln[5:]
        if re.match(r"^P\d+\s+SIGMA_CON=", ln.strip()):
            row = ln.strip()
            break
    if not row:
        raise RuntimeError("明细段没有数据行: " + str(path))
    kv: dict[str, float] = {}
    for k, v in re.findall(r"\b([A-Z_][A-Z0-9_]*)=([-+0-9.eE]+)", row):
        kv[k] = float(v)
    return kv, note


def parse_3d_csv(path: Path) -> dict:
    kv: dict[str, float] = {}
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        p = ln.split(",")
        if len(p) == 5 and p[0] == "S":
            kv[p[3]] = float(p[4])
    return kv


# ---------------------------------------------------------------- 冻结

def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def run(exe: Path, args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run([str(exe)] + args, cwd=str(cwd), env=env_with_msys(),
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=180)


def build(bin_dir: Path, cases: Path, ref_dir: Path) -> int:
    exe2 = bin_dir / EXE_2D
    exe3 = bin_dir / EXE_3D
    if not (exe2.is_file() and exe3.is_file()):
        print("错误: 需要已编译的 example_mixed.exe 与 loss_probe_3d.exe")
        return 2
    ref: dict = {"version": 1, "cases": {}, "frozen_sha256": {},
                 "thresholds": {"rel": REL_TOL, "abs": ABS_TOL}}
    for cid in DAT_CASES:
        dat = cases / (cid + ".dat")
        cfg = load_cfg(dat)
        entry: dict = {"dat_sha256": sha256(dat), "cfg": cfg}
        if cid in NEGATIVE:
            entry["negative_keywords"] = list(NEGATIVE[cid])
            ref["cases"][cid] = entry
            continue
        entry["oracle"] = oracle(cfg)
        with tempfile.TemporaryDirectory(prefix="lossanti_") as td:
            work = Path(td)
            shutil.copy2(dat, work / dat.name)
            if cid == "anti_3d":
                p = run(exe3, [dat.name, "fresh.csv"], work)
                if p.returncode != 0:
                    print("3D 探针失败:", cid, p.returncode, p.stderr[-400:])
                    return 3
                fresh = work / "fresh.csv"
            else:
                p = run(exe2, [dat.name, "--loss-detail-only"], work)
                if p.returncode != 0:
                    print("2D 明细失败:", cid, p.returncode, p.stderr[-400:])
                    return 3
                fresh = work / (dat.stem + ".loss.txt")
            ext = FROZEN_EXT.get(cid, ".loss.txt")
            frozen = cases / (cid + ext)
            shutil.copy2(fresh, frozen)
            entry["frozen"] = cid + ext
            ref["frozen_sha256"][cid] = sha256(frozen)
        ref["cases"][cid] = entry
        print("  冻结", cid, ext, "sha256=" + ref["frozen_sha256"][cid][:16])

    ref_dir.mkdir(parents=True, exist_ok=True)
    out = ref_dir / "prestress_loss_anti_2d.json"
    out.write_text(json.dumps(ref, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print("  写入", out)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", default=str(HERE / "cases"))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    a = ap.parse_args()
    return build(Path(a.bin).resolve(), Path(a.cases).resolve(),
                 Path(a.reference).resolve())


if __name__ == "__main__":
    sys.exit(main())
