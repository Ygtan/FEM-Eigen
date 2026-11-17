#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 36 项 CTest：2D 阶段分项损失 benchmark 的独立闭式 oracle 与冻结器。

纪律与前 35 项一致：

* oracle 是**第二份**纯 Python 实现（JTG 3362-2018 §6.2 + 附录 G 反摩擦 +
  附录 C 显式表值龄期），不调用求解器、不读求解器输出，只从
  cases/*.dat 的 ``*PRESTRESS_LOSS`` 行读参数（唯一参数来源）；
* build_reference.py 用**当前二进制**生成冻结证据；compare 脚本再用同一
  二进制重放并逐字节对拍，防止冻结后内核被悄悄改动；
* 硬负例只记录期望的退出码语义与错误关键词，不产生数值行；
* 软负例（无 LOSSREF）正常求解，但阶段段只写 ERROR=，也冻结下来对拍。

用法：
    python build_reference.py --bin <dir_with_example_stage_and_loss_probe_3d> \
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

NORMAL = ("stage_default", "stage_l2anti", "stage_c2d",
          "stage_cskip", "stage_same3d", "stage_noref")
DAT_CASES = NORMAL + ("stage_badref", "stage_conflict")
HARD_NEGATIVE = {
    "stage_badref": ("不存在",),       # LOSSREF=99
    "stage_conflict": ("只能",),       # LOSS 与 LOSSREF 二选一
}
# 第二十七次：阶段有效预加力段（[PRESTRESSFORCE2D]）的算例集合。
# NORMAL 全部复算一遍（同一算例同时冻结损失段与有效预加力段），
# 另加多束多阶段锚点与「无阶段束」软负例。
FORCE_CASES = NORMAL + ("force_multi", "force_none")
SOFT_NEG_FORCE = {"force_none": "no_stage_prestress"}
EXE_2D = "example_stage.exe"
EXE_3D = "loss_probe_3d.exe"

# 既有第 34 项 benchmark 的冻结证据：用于证明全局段 / 3D 探针未被污染。
LEGACY_CASES = Path(__file__).resolve().parents[1] / "prestress_loss_creep_2d" / "cases"

REL_TOL = 1.0e-11
ABS_TOL = 1.0e-9


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_cfg(dat: Path) -> dict:
    """只读 cases/*.dat 的 *PRESTRESS_LOSS 行，不读任何求解器输出。"""
    line = ""
    for ln in dat.read_text(encoding="utf-8", errors="replace").splitlines():
        if re.match(r"^\s*\d+\s+SIGMA_CON=", ln):
            line = ln
            break
    if not line:
        raise RuntimeError("算例缺少 *PRESTRESS_LOSS 参数行: " + str(dat))
    return {k: float(v) for k, v in
            re.findall(r"\b([A-Z_][A-Z0-9_]*)=([-+0-9.eE]+)", line)}


def oracle(cfg: dict) -> dict:
    """独立闭式：§6.2 六项 + 附录 G 反摩擦 + 附录 C 表值龄期。"""
    out: dict[str, float] = {}
    con, ep, ec = cfg["SIGMA_CON"], cfg["EP"], cfg["EC"]
    alpha = ep / ec if ec > 0.0 else 0.0
    mu, kk, th = cfg["MU"], cfg["K"], cfg["THETA"]

    sl1 = con * (1.0 - math.exp(-(mu * th + kk * cfg["X"])))
    if int(cfg.get("L2ANTI", 0)):
        sl1_anchor = con * (1.0 - math.exp(-(mu * th + kk * cfg["L"])))
        dsd = sl1_anchor / cfg["L"]
        if dsd > 0.0:
            lf = math.sqrt(cfg["DELTA_L"] * ep / dsd)
            if lf <= cfg["L"]:
                sl2 = 2.0 * dsd * lf
                branch = 1.0
            else:
                sl2 = dsd * (lf * lf / cfg["L"] + cfg["L"])
                branch = 2.0
            out.update({"LF": lf, "DSL2D": dsd, "DSL2PEAK": sl2,
                        "BRANCH": branch})
        else:
            sl2 = cfg["DELTA_L"] / cfg["L"] * ep
    else:
        sl2 = cfg["DELTA_L"] / cfg["L"] * ep

    out["SIGMA_CON"] = con
    out["SL1"] = sl1
    out["SL2"] = sl2
    out["SL3"] = 2.0 * (cfg["T2"] - cfg["T1"]) * 1.0e6
    out["SL4"] = alpha * cfg["SIGMA_PC4"]
    spe = con - sl1 - sl2 - out["SL4"]
    out["SPE"] = spe
    out["RATIO"] = spe / con
    out["SL5"] = cfg["PSI"] * cfg["ZETA"] * (0.52 * spe / cfg["FPK"] - 0.26) * spe

    den = 1.0 + 15.0 * cfg["RHO"] * cfg["RHOPS"]
    cskip, c2d = int(cfg.get("C_SKIP", 0)), int(cfg.get("C2D", 0))
    if cskip:
        out["SL6"] = 0.0
        out["C_SOURCE"] = 2.0
    elif c2d:
        hr = cfg["C_H"] / 100.0
        ds = 350.0 * hr * hr
        bs = math.sqrt((cfg["C_T"] - cfg["C_TS"]) /
                       (ds + cfg["C_T"] - cfg["C_TS"]))
        bs0 = math.sqrt((cfg["C_T0"] - cfg["C_TS"]) /
                        (ds + cfg["C_T0"] - cfg["C_TS"]))
        eps_dev = cfg["EPS0"] * (bs - bs0)
        rep = 55.0 if cfg["C_RH"] < 70.0 else 80.0
        bh = min(150.0 * (1.0 + (1.2 * rep / 100.0) ** 18) * hr + 250.0, 1500.0)
        bc = ((cfg["C_T"] - cfg["C_T0"]) /
              (bh + cfg["C_T"] - cfg["C_T0"])) ** 0.3
        out.update({"ECS_DEV": eps_dev, "PHI_DEV": cfg["PHI0"] * bc,
                    "BETAS": bs, "BETAS0": bs0, "BETAC": bc, "BETAH": bh,
                    "RHREP": rep, "C_SOURCE": 1.0})
        out["SL6"] = 0.9 * (ep * eps_dev + alpha * cfg["SIGMA_PC6"] * out["PHI_DEV"]) / den
    else:
        out["SL6"] = 0.9 * (ep * cfg["ECS"] + alpha * cfg["SIGMA_PC6"] * cfg["PHI"]) / den
    return out


def parse_stage_loss_txt(path: Path) -> tuple[dict, str, str]:
    """返回 (首行数据键值, NOTE, 状态行 ERROR= 或 '')。"""
    text = path.read_text(encoding="utf-8", errors="replace")
    note, row, err = "", "", ""
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("NOTE="):
            note = s[5:]
        if s.startswith("ERROR="):
            err = s
        if re.match(r"^STAGE\d+\s+P\d+\s+SIGMA_CON=", s):
            row = s
            break
    kv = {k: float(v) for k, v in
          re.findall(r"\b([A-Z_][A-Z0-9_]*)=([-+0-9.eE]+)", row)} if row else {}
    return kv, note, err


def load_stage_prestress(dat: Path) -> list[dict]:
    """只读 *STAGE_PRESTRESS 行（不读求解器输出）。

    返回按输入顺序编号的束参数；束号是**阶段内**从 1 起的序号，与内核
    ``ps.id = 本阶段已有束数 + 1`` 同一规则。只取 P / LOSS / LOSSREF 三个键：
    ELEM=1-4 这类带 '-' 的域不参与，避免把区间当数值。
    """
    rows: list[dict] = []
    in_block = False
    per_stage: dict[int, int] = {}
    for ln in dat.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s:
            continue
        if s.startswith("*"):
            in_block = s.upper().startswith("*STAGE_PRESTRESS")
            continue
        if not in_block:
            continue
        m = re.match(r"^(\d+)\s+(.*)$", s)
        if not m:
            continue
        stage = int(m.group(1))
        kv = dict(re.findall(r"\b(P|LOSS|LOSSREF)=([-+0-9.eE]+)", m.group(2)))
        per_stage[stage] = per_stage.get(stage, 0) + 1
        rows.append({"stage": stage, "id": per_stage[stage],
                     "P": float(kv.get("P", 0.0)),
                     "LOSS": float(kv.get("LOSS", 1.0)),
                     "LOSSREF": int(float(kv.get("LOSSREF", 0)))})
    return rows


def force_oracle(dat: Path, cfg: dict, ratio: float) -> list[dict]:
    """独立闭式：Pe = P·loss。LOSSREF 束的 loss 取 oracle 的 σpe/σcon。"""
    out: list[dict] = []
    for r in load_stage_prestress(dat):
        if r["LOSSREF"] > 0:
            loss, sc = ratio, cfg["SIGMA_CON"]
        else:
            loss, sc = r["LOSS"], 0.0
        out.append({"stage": r["stage"], "id": r["id"], "SIGMA_CON": sc,
                    "LOSS": loss, "PE": r["P"] * loss})
    return out


def parse_stage_force_txt(path: Path) -> tuple[list[dict], str, str]:
    """返回 (按行序的有效预加力行, NOTE, ERROR 状态行)。"""
    text = path.read_text(encoding="utf-8", errors="replace")
    note, err = "", ""
    rows: list[dict] = []
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("NOTE="):
            note = s[5:]
        elif s.startswith("ERROR="):
            err = s
        elif re.match(r"^STAGE\d+\s+P\d+\s+SIGMA_CON=", s):
            m = re.match(r"^STAGE(\d+)\s+P(\d+)\s+(.*)$", s)
            kv = {k: float(v) for k, v in
                  re.findall(r"\b([A-Z_][A-Z0-9_]*)=([-+0-9.eE]+)", m.group(3))}
            kv["stage"] = float(m.group(1))
            kv["id"] = float(m.group(2))
            rows.append(kv)
    return rows, note, err


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
    exe2, exe3 = bin_dir / EXE_2D, bin_dir / EXE_3D
    if not (exe2.is_file() and exe3.is_file()):
        print("错误: 需要已编译的 example_stage.exe 与 loss_probe_3d.exe")
        return 2
    ref: dict = {"version": 1, "cases": {}, "frozen_sha256": {},
                 "thresholds": {"rel": REL_TOL, "abs": ABS_TOL}}
    for cid in DAT_CASES:
        dat = cases / (cid + ".dat")
        cfg = load_cfg(dat)
        entry: dict = {"dat_sha256": sha256(dat), "cfg": cfg}
        if cid in HARD_NEGATIVE:
            entry["negative_keywords"] = list(HARD_NEGATIVE[cid])
            ref["cases"][cid] = entry
            continue
        entry["oracle"] = oracle(cfg) if cid != "stage_noref" else {}
        with tempfile.TemporaryDirectory(prefix="stageloss_") as td:
            wd = Path(td)
            shutil.copy2(dat, wd / dat.name)
            p = run(exe2, [dat.name, "--stage-loss-detail"], wd)
            if p.returncode != 0:
                print("2D 阶段明细失败:", cid, p.returncode, p.stderr[-400:])
                return 3
            fresh = wd / (dat.stem + ".stage.loss.txt")
            if not fresh.is_file():
                print("2D 阶段明细未产出:", cid)
                return 3
            frozen = cases / (cid + ".stage.loss.txt")
            shutil.copy2(fresh, frozen)
            entry["frozen"] = cid + ".stage.loss.txt"
            ref["frozen_sha256"][cid] = sha256(frozen)
        ref["cases"][cid] = entry
        print("  冻结", cid, "sha256=" + ref["frozen_sha256"][cid][:16])

    # 第二十七次：阶段有效预加力段（[PRESTRESSFORCE2D]）独立冻结。
    # 与损失段分开跑一次（不同开关、不同文件），证明两个 opt-in 互不依赖。
    ref["force"] = {"cases": {}, "frozen_sha256": {}}
    for cid in FORCE_CASES:
        dat = cases / (cid + ".dat")
        cfg = load_cfg(dat)
        entry = {"dat_sha256": sha256(dat)}
        with tempfile.TemporaryDirectory(prefix="stageforce_") as td:
            wd = Path(td)
            shutil.copy2(dat, wd / dat.name)
            p = run(exe2, [dat.name, "--stage-prestress-force"], wd)
            if p.returncode != 0:
                print("2D 阶段有效预加力失败:", cid, p.returncode, p.stderr[-400:])
                return 3
            fresh = wd / (dat.stem + ".stage.force.txt")
            if not fresh.is_file():
                print("2D 阶段有效预加力未产出:", cid)
                return 3
            frozen = cases / (cid + ".stage.force.txt")
            shutil.copy2(fresh, frozen)
            entry["frozen"] = cid + ".stage.force.txt"
            ref["force"]["frozen_sha256"][cid] = sha256(frozen)
        if cid in SOFT_NEG_FORCE:
            entry["soft_negative"] = SOFT_NEG_FORCE[cid]
        else:
            spell = oracle(cfg)
            ratio = spell["SPE"] / spell["SIGMA_CON"]
            entry["oracle"] = force_oracle(dat, cfg, ratio)
        ref["force"]["cases"][cid] = entry
        print("  冻结", cid, "force sha256="
              + ref["force"]["frozen_sha256"][cid][:16])

    # 既有第 34 项冻结证据的哈希：证明全局段与 3D 探针未被本轮污染。
    legacy_global = LEGACY_CASES / "creep_auto.loss.txt"
    legacy_3d = LEGACY_CASES / "creep_3d.csv"
    ref["legacy_sha256"] = {
        "creep_auto.loss.txt": sha256(legacy_global),
        "creep_3d.csv": sha256(legacy_3d),
    }

    ref_dir.mkdir(parents=True, exist_ok=True)
    out = ref_dir / "stage_loss_2d.json"
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
