#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二十三次 B4 · σl6 收缩徐变附录 C benchmark 的独立闭式 oracle 与冻结器。

设计纪律（与前 33 项 benchmark 一致）：

* oracle 用纯 Python 写第二份 JTG 3362-2018 §6.2 + 附录 C 闭式，
  **不调用任何求解器、不读求解器输出**，只从 cases/*.dat 的
  ``*PRESTRESS_LOSS`` 行读参数（参数来源唯一，避免脚本里手写两份）；
  表值 εcs0/φ0 由用户显式录入，程序只按解析发展公式组装，绝不查表拟合。
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

# 正常算例（有冻结 .loss.txt / .csv）与负例（只验拒绝语义）
NORMAL = ("creep_default", "creep_skip", "creep_auto",
          "creep_auto_rh55", "creep_3d")
DAT_CASES = ("creep_default", "creep_skip", "creep_auto", "creep_auto_rh55",
             "creep_3d", "creep_noc2d", "creep_badrh", "creep_conflict",
             "creep_badage")
NEGATIVE = {
    "creep_noc2d": ("EPS0", "缺 EPS0"),
    "creep_badrh": ("RH", "40"),
    "creep_conflict": ("C_SKIP", "C2D"),
    "creep_badage": ("龄期", "t > t0"),
}
FROZEN_EXT = {"creep_3d": ".csv"}
EXE_2D = "example_mixed.exe"
EXE_3D = "loss_probe_3d.exe"

REL_TOL = 1.0e-11
ABS_TOL = 1.0e-9


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
    """独立闭式：JTG 3362-2018 §6.2 + 附录 C。

    c2d=0：σl6 用 ECS/PHI 用户直给（默认路径）。
    cSkip=1：σl6=0，cSource=2。
    c2d=1：表值 εcs0/φ0 + 解析发展系数：
      βs(t,ts)=[(t-ts)/(350(h/100)²+(t-ts))]^0.5
      εcs(t,t0)=εcs0[βs(t,ts)-βs(t0,ts)]
      βH=min(150[1+(1.2RHrep/100)^18](h/100)+250, 1500), RHrep=55 若 RH<70 否则 80
      βc=[(t-t0)/(βH+(t-t0))]^0.3
      φ(t,t0)=φ0·βc
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
    # σl2 本 benchmark 不涉及反摩擦，用默认直线式（DELTA_L/L*Ep）。
    out["SL2"] = cfg["DELTA_L"] / cfg["L"] * ep
    out["SPE"] = con - out["SL1"] - out["SL2"] - out["SL4"]
    out["RATIO"] = out["SPE"] / con
    out["SL5"] = cfg["PSI"] * cfg["ZETA"] * (0.52 * out["SPE"] / cfg["FPK"] - 0.26) * out["SPE"]

    cskip = int(cfg.get("C_SKIP", 0))
    c2d = int(cfg.get("C2D", 0))
    den = 1.0 + 15.0 * cfg["RHO"] * cfg["RHOPS"]
    if c2d == 0 and cskip == 0:
        # 默认路径内核不输出 C_SOURCE 字段（cSource=0 不追加），oracle 此处也不给，
        # 保持与内核输出字段一致，避免「missing」误判。
        out["SL6"] = 0.9 * (ep * cfg["ECS"] + alpha * cfg["SIGMA_PC6"] * cfg["PHI"]) / den
    elif cskip == 1:
        out["SL6"] = 0.0
        out["C_SOURCE"] = 2.0
    else:
        eps0, phi0 = cfg["EPS0"], cfg["PHI0"]
        rh, hmm = cfg["C_RH"], cfg["C_H"]
        t0, ts, t = cfg["C_T0"], cfg["C_TS"], cfg["C_T"]
        hr = hmm / 100.0
        den_s = 350.0 * hr * hr
        beta_s = math.sqrt((t - ts) / (den_s + (t - ts)))
        beta_s0 = math.sqrt((t0 - ts) / (den_s + (t0 - ts)))
        eps_dev = eps0 * (beta_s - beta_s0)
        rh_rep = 55.0 if rh < 70.0 else 80.0
        beta_h = min(150.0 * (1.0 + (1.2 * rh_rep / 100.0) ** 18) * hr + 250.0, 1500.0)
        beta_c = ((t - t0) / (beta_h + (t - t0))) ** 0.3
        phi_dev = phi0 * beta_c
        out["ECS_DEV"] = eps_dev
        out["PHI_DEV"] = phi_dev
        out["BETAS"] = beta_s
        out["BETAS0"] = beta_s0
        out["BETAC"] = beta_c
        out["BETAH"] = beta_h
        out["RHREP"] = rh_rep
        out["SL6"] = 0.9 * (ep * eps_dev + alpha * cfg["SIGMA_PC6"] * phi_dev) / den
        out["C_SOURCE"] = 1.0
    return out


def parse_loss_txt(path: Path) -> tuple[dict, str]:
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
        with tempfile.TemporaryDirectory(prefix="losscreep_") as td:
            work = Path(td)
            shutil.copy2(dat, work / dat.name)
            if cid == "creep_3d":
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
        print("  冻结", cid, ext, "sha256=" + ref["frozen_sha256"][cid][:16],
              "SL6=" + format(entry["oracle"].get("SL6", 0.0), ".6e"))

    ref_dir.mkdir(parents=True, exist_ok=True)
    out = ref_dir / "prestress_loss_creep_2d.json"
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
