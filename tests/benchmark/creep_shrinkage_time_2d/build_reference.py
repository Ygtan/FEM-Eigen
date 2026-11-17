#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二十四次 B4 · 结构时变徐变/收缩规范算例书锚点 benchmark 的独立闭式 oracle 与冻结器。

设计纪律（与 σl6 benchmark 一致，但口径刻意不同）：

* oracle 用**纯 Python** 写第二份 JTG 3362-2018 附录 C.2.1 / C.1 连续拟合闭式，
  **不调用任何求解器、不读求解器输出**，只从 cases/*.dat 的 ``*CREEP_MODEL`` 行
  读参数（参数来源唯一）。这是**结构整体时变**口径：φ0=φRH·βF·βT0，
  βH 用实际 RH 连续值，εcs0=εS·βRH——刻意区别于 σl6 的表值 + RH 代表值口径。
* oracle 同时按 AEMM 写出 E_eff=E/(1+χ·φ(t)) 与逐龄期 φ(t)/εcs(t)/βs，
  与 creep_probe_2d.exe 暴露的内核中间量逐键核对。
* build_reference.py 只用**当前二进制**生成冻结证据 + reference JSON；
  compare 脚本再用同一二进制重放并逐字节对拍，防止冻结后内核被悄悄改动。
* 负例（缺 RH/H/FCM、旧 CURVE 口径）只验拒绝语义，不产生数值行。

用法：
    python build_reference.py --bin <dir_with_creep_probe_2d> \
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

NORMAL = ("creep_anchor", "creep_rh55", "creep_t0late")
NEGATIVE = {
    "creep_norh": "RH",
    "creep_nofcm": "FCM",
    "creep_noh": "H",
    "creep_curve": "未给全 RH/H/FCM",
}
DAT_CASES = NORMAL + tuple(NEGATIVE)
EXE = "creep_probe_2d.exe"

REL_TOL = 1.0e-11
ABS_TOL = 1.0e-9


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_cfg(dat: Path) -> dict:
    text = dat.read_text(encoding="utf-8", errors="replace")
    cfg: dict[str, float] = {}
    for ln in text.splitlines():
        if ln.strip().startswith("*CREEP_MODEL"):
            for k, v in re.findall(r"\b([A-Z_][A-Z0-9_]*)=([-+0-9.eE]+)", ln):
                cfg[k] = float(v)
            break
    return cfg


def read_ages(dat: Path) -> list[float]:
    """从 *STAGE 行读阶段末龄期（探针采样点）。"""
    ages: list[float] = []
    for ln in dat.read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"\s*\d+\s+.*\bTIME=([-+0-9.eE]+)", ln)
        if m:
            ages.append(float(m.group(1)))
    return ages


def nominal_creep(rh: float, h: float, fcm: float, t0: float) -> dict:
    """独立复现 code_loads.hpp::jtg3362NominalCreep（附录 C.2.1 + C.2.2）。"""
    rh = min(99.0, max(40.0, rh))
    h = max(h, 1.0)
    fcm = max(fcm, 1.0)
    t0 = max(t0, 0.1)
    phi_rh = 1.0 + (1.0 - rh / 100.0) / (0.46 * (h / 100.0) ** (1.0 / 3.0))
    beta_f = 5.3 / math.sqrt(fcm / 10.0)
    beta_t0 = 1.0 / (0.1 + t0 ** 0.2)
    beta_h = 150.0 * (1.0 + (1.2 * rh / 100.0) ** 18) * (h / 100.0) + 250.0
    beta_h = min(beta_h, 1500.0)
    return {
        "phiRH": phi_rh,
        "betaF": beta_f,
        "betaT0": beta_t0,
        "PHI0": phi_rh * beta_f * beta_t0,
        "BETAH": beta_h,
    }


def shrinkage_nominal(rh: float, fcm: float) -> dict:
    """独立复现 code_loads.hpp::jtg3362ShrinkageNominal（附录 C.1）。"""
    rh = min(99.0, max(40.0, rh))
    fcm = max(fcm, 1.0)
    eps_s = (160.0 + 50.0 * (9.0 - fcm / 10.0)) * 1e-6
    beta_rh = 1.55 * (1.0 - (rh / 100.0) ** 3)
    return {"epsS": eps_s, "betaRH": beta_rh, "ECS0": eps_s * beta_rh}


def beta_s(t: float, ts: float, h: float) -> float:
    dt = max(0.0, t - ts)
    h = max(h, 1.0)
    den = 350.0 * (h / 100.0) ** 2 + dt
    return 0.0 if den <= 0.0 else math.sqrt(dt / den)


def oracle(cfg: dict, ages: list[float]) -> dict:
    rh, h, fcm = cfg["RH"], cfg["H"], cfg["FCM"]
    t0 = cfg["T0"]
    ts = cfg.get("TS", t0)
    power = cfg.get("POWER", 0.3)
    e0 = 200.0e9
    cre = nominal_creep(rh, h, fcm, t0)
    shr = shrinkage_nominal(rh, fcm)
    out: dict = {"M": {**cre, **shr}, "S": {}}
    for age in ages:
        phi = cre["PHI0"] * ((max(0.0, age - t0)) / (cre["BETAH"] + max(0.0, age - t0))) ** power
        ecs = shr["ECS0"] * (beta_s(age, ts, h) - beta_s(t0, ts, h))
        # age key 用 %g 规范化（28.0→"28"），与探针 CSV 输出格式一致。
        out["S"]["%g" % age] = {
            "PHI": phi,
            "ECS": ecs,
            "BETAS": beta_s(age, ts, h),
            "EEFF": e0 / (1.0 + 0.8 * phi),
        }
    return out


def parse_probe_csv(path: Path) -> dict:
    """返回 {'M': {KEY: v}, 'S': {age: {KEY: v}}}。"""
    out: dict = {"M": {}, "S": {}}
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        p = ln.split(",")
        if len(p) == 3 and p[0] == "M":
            out["M"][p[1]] = float(p[2])
        elif len(p) == 4 and p[0] == "S":
            out["S"].setdefault(p[1], {})[p[2]] = float(p[3])
    return out


def env_with_msys() -> dict:
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    return env


def run(exe: Path, args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run([str(exe)] + args, cwd=str(cwd), env=env_with_msys(),
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=180)


def build(bin_dir: Path, cases: Path, ref_dir: Path) -> int:
    exe = bin_dir / EXE
    if not exe.is_file():
        print("错误: 需要已编译的 " + EXE)
        return 2
    ref: dict = {"version": 1, "cases": {}, "frozen_sha256": {},
                 "thresholds": {"rel": REL_TOL, "abs": ABS_TOL}}
    for cid in DAT_CASES:
        dat = cases / (cid + ".dat")
        cfg = load_cfg(dat)
        entry: dict = {"dat_sha256": sha256(dat), "cfg": cfg}
        if cid in NEGATIVE:
            entry["negative_keyword"] = NEGATIVE[cid]
            ref["cases"][cid] = entry
            wd = Path(tempfile.mkdtemp(prefix="creprefneg_"))
            try:
                shutil.copy2(dat, wd / dat.name)
                p = run(exe, [dat.name, "fresh.csv"], wd)
                print("  负例", cid, "rc=", p.returncode,
                      "stderr=", p.stderr.strip()[-70:])
            finally:
                shutil.rmtree(wd, ignore_errors=True)
            continue
        entry["ages"] = read_ages(dat)
        entry["oracle"] = oracle(cfg, entry["ages"])
        with tempfile.TemporaryDirectory(prefix="crepref_") as td:
            wd = Path(td)
            shutil.copy2(dat, wd / dat.name)
            p = run(exe, [dat.name, "fresh.csv"], wd)
            if p.returncode != 0:
                print("  探针失败:", cid, p.returncode, p.stderr[-400:])
                return 3
            fresh = wd / "fresh.csv"
            frozen = cases / (cid + ".csv")
            shutil.copy2(fresh, frozen)
            entry["frozen"] = cid + ".csv"
            ref["frozen_sha256"][cid] = sha256(frozen)
        ref["cases"][cid] = entry
        print("  冻结", cid, "sha256=" + ref["frozen_sha256"][cid][:16],
              "PHI0=" + format(entry["oracle"]["M"]["PHI0"], ".6f"),
              "BETAH=" + format(entry["oracle"]["M"]["BETAH"], ".4f"))
    ref_dir.mkdir(parents=True, exist_ok=True)
    out = ref_dir / "creep_shrinkage_time_2d.json"
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
