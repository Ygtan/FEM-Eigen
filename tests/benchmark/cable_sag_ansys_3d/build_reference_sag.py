#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS CSV 生成 CB3D14 / CB3D15 黄金参考值。

CB3D14 = 单根 CABLE 的索垂度 / Ernst 等效弹模对标。FEM 侧由
``cable_sag_probe.exe`` 调用内核闭式；ANSYS 侧是 80 段 LINK180 索链 +
水平均布节点自重 + NLGEOM，从正弦扰动几何独立收敛。7 个工况：
5 点 L0 扫描、1 个 W=0 直杆基线、1 个 W/2 低重度。

CB3D15 = 与 CB3D13 同参数的自锚式悬索桥，但 cable_dens=7850 (W>0)。
验证空缆 stage 0 不再命中内核 H 二分下界 100 N，矢高 > 0。

每一项 ``ref`` 都来自真实 ANSYS 18.0 批处理 CSV；H / f / 逐段索力是
从 CSV 派生的量，派生口径写在 ``derived`` 字段，并由
``compare_cable_sag_ansys_3d.py`` 从冻结 CSV 重算校验，防黄金值漂移。

用法：
    python build_reference_sag.py            # 生成 reference/*.json
    python build_reference_sag.py --check    # 只校验一致性
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference"

RTOL = 2e-3
CASE_ID = "cb3d14_cable_sag"
BRIDGE_ID = "cb3d15_cable_selfweight"
SAG_FILE = CASE_ID + ".sag"
FEM_CSV = CASE_ID + "_fem.csv"

TITLE = "CB3D14 cable sag / Ernst modulus + CB3D15 W>0 empty-cable state vs ANSYS 18.0"
DESCRIPTION = (
    "Single CABLE element: E=195 GPa, A=0.006 m2, LH=100 m, W=462.051 N/m. "
    "FEM closed form fem_3d.hpp::cableTensionFromChord + "
    "fem_3d_nl.hpp::cableErnstModulus; ANSYS 80-segment LINK180 chain with "
    "horizontal-uniform nodal self weight, NLGEOM,ON, INISTATE,EPEL and a "
    "0.20 m sine-perturbed straight initial geometry. 7 jobs (5-point L0 "
    "sweep, W=0 rod baseline, W/2 low-weight). Ernst tangent is verified by "
    "finite difference of ANSYS H over the L0 sweep. CB3D15 confirms that "
    "W>0 removes the degenerate empty-cable stage (H >> 100 N, sag > 0)."
)

ANSYS_LOG_EXPECT = {"errors": 0, "warnings": 0}


def read_csv_table(path: Path) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        p = line.split()
        if not p:
            continue
        try:
            out[int(float(p[0]))] = [float(v) for v in p[1:]]
        except ValueError:
            continue
    return out


def read_fem_csv(path: Path) -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    keys: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split(",")
        if p[0] == "case":
            keys = p[1:]
            continue
        rows[p[0]] = {k: float(v) for k, v in zip(keys, p[1:])}
    return rows


def parse_sag(path: Path) -> tuple[dict[str, float], list[dict]]:
    """解析 *.sag：全局参数 + 工况列表。"""
    cfg: dict[str, float] = {}
    cases: list[dict] = []
    for line in path.read_text(encoding="ascii").splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        if "=" in t and not t.startswith("C "):
            k, v = t.split("=", 1)
            cfg[k.strip()] = float(v)
        elif t.startswith("C "):
            p = t.split()
            cases.append(dict(id=p[1], W=float(p[2]), L0=float(p[3]),
                              amp=float(p[4])))
    return cfg, cases


def z_init(amp: float, N: int) -> list[float]:
    """与 ansys_bridge_sag.sag_apdl() 完全相同的初始几何扰动。"""
    return [-amp * math.sin(math.pi * k / N) for k in range(N + 1)]


def derive(cs: dict, cfg: dict) -> dict:
    """从冻结 CSV 派生 H / f / 逐段索力 T，口径与 sag_from_ansys 一致。"""
    N = int(cfg["NSEG"])
    dx = cfg["LH"] / N
    z0 = z_init(cs["amp"], N)
    disp = read_csv_table(CASES / (cs["id"] + "_disp.csv"))
    link = read_csv_table(CASES / (cs["id"] + "_link.csv"))
    if any(k not in disp for k in range(1, N + 2)) or \
       any(k not in link for k in range(1, N + 1)):
        raise ValueError("{}: incomplete ANSYS CSV".format(cs["id"]))
    x = [(k - 1) * dx + disp[k][0] for k in range(1, N + 2)]
    z = [z0[k - 1] + disp[k][2] for k in range(1, N + 2)]
    T = [link[k][0] for k in range(1, N + 1)]
    Hs = []
    for k in range(N):
        ddx, ddz = x[k + 1] - x[k], z[k + 1] - z[k]
        Lk = math.hypot(ddx, ddz)
        Hs.append(T[k] * abs(ddx) / Lk if Lk > 0 else 0.0)
    return {"H": sum(Hs) / N, "f": z[0] - min(z), "T": T,
            "H_spread": max(Hs) - min(Hs)}


def parse_bridge_stage0(path: Path) -> dict:
    """从冻结 *_cable.txt 的 [ERECTION] 首行取空缆态 H / sag / midz。"""
    stage = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if s == "[ERECTION]":
            stage = True
            continue
        if stage and s.startswith("["):
            break
        if stage and s:
            p = s.split()
            # 格式: "<stage> <name> H=<N> sag=<m> midz=<m> hmin=.. hmax=.."
            kv = {}
            for tok in p[2:]:
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    kv[k] = float(v)
            return {"stage": int(p[0]), "name": p[1],
                    "H": kv["H"], "sag": kv["sag"], "midz": kv["midz"]}
    raise ValueError("{}: no [ERECTION] stage 0".format(path))


def build() -> dict:
    cfg, cases = parse_sag(CASES / SAG_FILE)
    fem = read_fem_csv(CASES / FEM_CSV)
    checks: list[dict] = []
    case_meta: list[dict] = []
    derived: dict[str, dict] = {}

    for cs in cases:
        cid = cs["id"]
        d = derive(cs, cfg)
        derived[cid] = {"H": d["H"], "f": d["f"],
                        "T": [d["T"][0], d["T"][len(d["T"]) // 2], d["T"][-1]]}
        case_meta.append(dict(
            id=cid, W=cs["W"], L0=cs["L0"], amp=cs["amp"],
            group="rod" if cs["W"] == 0.0 else "sag",
            apdl=cid + ".inp",
            csv={"disp": cid + "_disp.csv", "link": cid + "_link.csv"},
        ))
        checks.append(dict(case=cid, key="H", ref=d["H"],
                           group="force", tol=RTOL, source="ansys_csv"))
        if cs["W"] > 0.0:
            checks.append(dict(case=cid, key="f", ref=d["f"],
                               group="shape", tol=RTOL, source="ansys_csv"))
        for k, t in enumerate(d["T"], 1):
            checks.append(dict(case=cid, key="T_%02d" % k, ref=t,
                               group="segforce", tol=RTOL, source="ansys_csv"))

    # Ernst 切线：ANSYS 的 H 在 L0 扫描上做中心差分（sag_02 / sag_04）
    sweep = [c for c in cases if c["W"] == cases[0]["W"] and
             abs(c["L0"] - cases[0]["L0"]) < 1.0 and c["id"].startswith("sag_0")]
    if len(sweep) < 3:
        raise ValueError("need >= 3 L0 sweep points for the tangent")
    # 用中间三点做中心差分（步长 0.002 m），与 BridgeStudio 报告口径一致；
    # 用最远两点会引入 O(h^2) 曲率误差，把切线相对差从 1.05e-4 放大到 1.32e-4。
    lo, mid, hi = sweep[1], sweep[2], sweep[3]
    fd = (derived[hi["id"]]["H"] - derived[lo["id"]]["H"]) / (hi["L0"] - lo["L0"])
    fem_tangent = fem[mid["id"]]["dTdL0"]
    tangent = dict(pair=[lo["id"], hi["id"]], center=mid["id"],
                   ansys_fd=fd, fem_ernst=fem_tangent,
                   rel=abs(fd - fem_tangent) / abs(fem_tangent),
                   tol=RTOL,
                   description="dT/dL0 = -(1+T/EA)*Eeq*A/L0 (kernel analytic tangent)")

    rod = next(c for c in cases if c["W"] == 0.0)
    degeneration = dict(case=rod["id"], Eeq_ratio=fem[rod["id"]]["Eeq_ratio"],
                        expect=1.0, tol=1e-12)

    bpath = CASES / (BRIDGE_ID + "_cable.txt")
    b0 = parse_bridge_stage0(bpath)
    bridge = dict(id=BRIDGE_ID,
                  dat=BRIDGE_ID + ".dat", cable=bpath.name,
                  stage0=b0,
                  criteria=dict(H_min=1.0e3, sag_min=1.0e-6),
                  description="W>0 gives a real catenary empty-cable state; "
                              "W=0 (CB3D13) degenerates to the 100 N bisection floor")

    return {
        "id": CASE_ID, "title": TITLE, "description": DESCRIPTION,
        "params": dict(E=cfg["E"], A=cfg["A"], LH=cfg["LH"],
                       NSEG=int(cfg["NSEG"])),
        "rtol": RTOL,
        "sag_file": SAG_FILE, "fem_csv": FEM_CSV,
        "cases": case_meta,
        "checks": checks,
        "derived": derived,
        "derived_rules": {
            "H": "mean over segments of T_k*|dx_k|/L_k (x = x_init + UX)",
            "f": "z_first - min(z), z = z_init + UZ, z_init=-amp*sin(pi*k/N)",
            "T_k": "LINK180 SMISC,1 from *_link.csv, k=1..NSEG",
            "ansys_fd": "central difference of ANSYS H over the L0 sweep",
        },
        "tangent": tangent,
        "degeneration": degeneration,
        "bridge": bridge,
        "ansys_log_expect": ANSYS_LOG_EXPECT,
        "fem_exe": "cable_sag_probe.exe",
        "fem_exe_args": ["<sag_file>", "<out_csv>"],
    }


def check(doc: dict) -> bool:
    cfg, cases = parse_sag(CASES / SAG_FILE)
    ok = True
    n = 0
    for cs in cases:
        cid = cs["id"]
        d = derive(cs, cfg)
        want = {"H": d["H"], "f": d["f"]}
        for k, t in enumerate(d["T"], 1):
            want["T_%02d" % k] = t
        for c in doc["checks"]:
            if c["case"] != cid:
                continue
            n += 1
            if c["source"] != "ansys_csv":
                continue
            got = want[c["key"]]
            if abs(got - c["ref"]) > 1e-9 + 1e-12 * abs(c["ref"]):
                print("[FAIL] {} {} got={!r} ref={!r}".format(
                    cid, c["key"], got, c["ref"]))
                ok = False
    print("[{}] {}  {} checks match frozen ANSYS CSVs".format(
        "PASS" if ok else "FAIL", doc["id"], n))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    doc = build()
    if a.check:
        return 0 if check(doc) else 1
    REF.mkdir(parents=True, exist_ok=True)
    out = REF / (doc["id"] + "_ansys.json")
    out.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print("[WRITE] {}  {} checks / {} cases / tangent rel={:.3e}".format(
        out, len(doc["checks"]), len(doc["cases"]), doc["tangent"]["rel"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
