#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS CSV 生成 CB3D16 大转动梁柱效应黄金参考值。

CB3D16 = FEM-Eigen ``fem_3d_nl.hpp`` 共转列式 CR 梁 vs ANSYS 18.0
BEAM188 + NLGEOM。8 个工况分两组：

  CIRCLE (4) : 悬臂端弯矩卷圆，lambda = 0.25/0.50/0.75/1.00
               -> 有限转动 2*pi，尖端回到根部
  PDELTA (4) : 简支梁轴压 P/Pcr = 0/0.25/0.50/0.75 + 跨中横向力
               -> 解析几何刚度，放大系数 1 -> ~3.96

每一项 ``check`` 都是真实的 ANSYS CSV 节点位移/转角值；
``compare_beam_nl_ansys_3d.py`` 用 FEM 重跑值按
``|got-ref| <= ATOL + RTOL*scale`` 逐项比对，同时校验 FEM CLI 对
冻结 ``*_fem.csv`` 的逐位复现、ANSYS 日志 0 error / 0 warning，
以及解析标尺判据（卷圆闭合、尖端转角、P-Delta 放大）。

用法：
    python build_reference_beamnl.py            # 生成 reference/*.json
    python build_reference_beamnl.py --check    # 只校验一致性
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
ATOL = 1e-8
CASE_ID = "cb3d16_beam_nl"
BEAM_FILE = CASE_ID + ".beam"
FEM_CSV = CASE_ID + "_fem.csv"
MODEL_JSON = CASE_ID + ".model.json"

TITLE = "CB3D16 large-rotation beam-column (CR kernel) vs ANSYS 18.0 BEAM188+NLGEOM"
DESCRIPTION = (
    "FEM-Eigen fem_3d_nl.hpp co-rotational beam (bisecting-vector triad, "
    "analytic geometric stiffness crBeamTangent12, Newton with load steps) "
    "is exported by beam_nl_probe.exe. ANSYS 18.0 uses 2-node BEAM188 with "
    "ASEC section, NLGEOM,ON, NSUBST/NEQIT and writes nodal UX/UY/UZ/ROTX/"
    "ROTY/ROTZ CSV. CIRCLE: cantilever end moment M=lambda*2*pi*E*Iz/L, "
    "n=8, lambda=1 rolls the beam into a full circle. PDELTA: simply "
    "supported beam-column with P=lambda*Pcr and midspan Q=2e4 N, n=64, "
    "amplification 1 to ~3.96. PDELTA group uses KEYOPT(3)=2 cubic "
    "interpolation and G raised 1000x to remove Timoshenko shear locking "
    "and shear flexibility so that BEAM188 matches the EB kernel."
)

ANSYS_LOG_EXPECT = {"errors": 0, "warnings": 0}


def read_disp_csv(path: Path) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        p = line.split()
        if not p:
            continue
        try:
            out[int(float(p[0]))] = [float(v) for v in p[1:]]
        except ValueError:
            continue
    return out


def read_fem_csv(path: Path) -> tuple[dict, dict]:
    nodal: dict[str, dict[int, list[float]]] = {}
    summ: dict[str, dict[str, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        p = line.split(",")
        if p[0] == "N":
            nodal.setdefault(p[1], {})[int(p[2])] = [float(v) for v in p[3:]]
        elif p[0] == "S":
            summ.setdefault(p[1], {})[p[2]] = float(p[3])
    return nodal, summ


def parse_beam(path: Path) -> tuple[dict[str, dict], list[dict]]:
    groups: dict[str, dict] = {}
    cases: list[dict] = []
    for line in path.read_text(encoding="ascii").splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        p = t.split()
        if p[0] == "G":
            kv = {}
            for tok in p[2:]:
                k, v = tok.split("=", 1)
                kv[k] = float(v)
            kv["NELEM"] = int(kv["NELEM"])
            groups[p[1]] = kv
        elif p[0] == "C":
            cases.append(dict(id=p[1], group=p[2], kind=p[3], lam=float(p[4])))
    return groups, cases


def exact_values(g: dict, c: dict) -> dict[str, float]:
    """与 ansys_bridge_beamnl.exact_values / beam_nl_probe.cpp 同源闭式。"""
    L, E, IZ, Q = g["L"], g["E"], g["IZ"], g["Q"]
    lam = c["lam"]
    out: dict[str, float] = {"lambda": lam}
    if c["kind"] == "CIRCLE":
        M = lam * 2.0 * math.pi * E * IZ / L
        theta = M * L / (E * IZ)
        R = L / theta if abs(theta) > 1e-14 else float("inf")
        out.update(M=M, theta=theta, xtip_ex=R * math.sin(theta),
                   ytip_ex=R * (1.0 - math.cos(theta)), rot_ex=theta)
    else:
        pcr = math.pi ** 2 * E * IZ / (L * L)
        P = lam * pcr
        ylin = Q * L ** 3 / (48.0 * E * IZ)
        if abs(P) > 1e-12:
            k = math.sqrt(P / (E * IZ))
            yex = (Q / (2.0 * P * k)) * (math.tan(k * L / 2.0) - k * L / 2.0)
        else:
            yex = ylin
        out.update(P=P, Pcr=pcr, ylin=ylin, y_pd_ex=yex,
                   amp_ex=yex / ylin if ylin else 0.0)
    return out


def derive(c: dict, g: dict) -> dict:
    """从冻结 ANSYS CSV 派生逐节点 checks 与本工况解析判据量。"""
    cid = c["id"]
    disp = read_disp_csv(CASES / (cid + "_disp.csv"))
    n = g["NELEM"]
    if set(disp) != set(range(1, n + 2)):
        raise ValueError("{}: incomplete ANSYS CSV ({} nodes, expected {})".format(
            cid, len(disp), n + 1))
    scales = [max(abs(row[k]) for row in disp.values()) for k in range(6)]
    checks = []
    for i in sorted(disp):
        for k in range(6):
            checks.append(dict(
                case=cid, key="U_{:03d}_{}".format(i, k + 1),
                ref=disp[i][k], group="dof{}".format(k + 1),
                scale=scales[k], tol=RTOL, atol=ATOL, source="ansys_csv"))
    ex = exact_values(g, c)
    tip, mid = n + 1, n // 2 + 1
    derived: dict[str, float] = {"scales": scales}
    if c["kind"] == "CIRCLE":
        rot = disp[tip][5]
        derived.update(xtip=g["L"] + disp[tip][0], ytip=disp[tip][1],
                       rot=rot, theta=ex["theta"], rot_ex=ex["rot_ex"],
                       rot_rel=abs(rot - ex["rot_ex"]) / max(abs(ex["rot_ex"]), 1.0),
                       closure=math.hypot(g["L"] + disp[tip][0], disp[tip][1]))
    else:
        ymid = -disp[mid][1]
        derived.update(ymid=ymid, y_exact=ex["y_pd_ex"],
                       amp=ymid / ex["ylin"], amp_ex=ex["amp_ex"],
                       y_rel=abs(ymid - ex["y_pd_ex"]) / abs(ex["y_pd_ex"]))
    return {"checks": checks, "derived": derived}


def build() -> dict:
    groups, cases = parse_beam(CASES / BEAM_FILE)
    nodal, summ = read_fem_csv(CASES / FEM_CSV)
    all_checks: list[dict] = []
    case_meta: list[dict] = []
    derived: dict[str, dict] = {}
    for c in cases:
        cid = c["id"]
        g = groups[c["group"]]
        d = derive(c, g)
        all_checks.extend(d["checks"])
        derived[cid] = d["derived"]
        case_meta.append(dict(
            id=cid, group=c["group"], kind=c["kind"], lam=c["lam"],
            nelem=g["NELEM"], apdl=cid + ".inp",
            csv={"disp": cid + "_disp.csv"}, out=cid + "_ansys.out"))
        if cid not in nodal or cid not in summ:
            raise ValueError("{}: missing FEM rows in {}".format(cid, FEM_CSV))

    cir = [c for c in cases if c["kind"] == "CIRCLE"]
    pd = [c for c in cases if c["kind"] == "PDELTA"]
    full = next(c for c in cir if abs(c["lam"] - 1.0) < 1e-12)
    rot_worst = max(derived[c["id"]]["rot_rel"] for c in cir)
    amps = [derived[c["id"]]["amp"] for c in pd]
    y_worst = max(derived[c["id"]]["y_rel"] for c in pd)
    criteria = {
        "full_circle_closure": dict(value=derived[full["id"]]["closure"],
                                    limit=1e-6,
                                    note="|tip-root| at M=2*pi*EI/L"),
        "max_tip_rotation_rel_error": dict(value=rot_worst, limit=1e-6,
                                           note="ANSYS ROTZ vs theta=M*L/(E*Iz)"),
        "pdelta_amplification": dict(values=amps, monotonic=True,
                                     min_final=3.0,
                                     note="ymid/ylin grows 1 -> ~3.96"),
        "pdelta_y_vs_exact_max_rel": dict(value=y_worst, limit=5e-3,
                                          note="ANSYS ymid vs Timoshenko beam-column"),
    }
    return {
        "id": CASE_ID, "title": TITLE, "description": DESCRIPTION,
        "rtol": RTOL, "atol": ATOL,
        "beam_file": BEAM_FILE, "fem_csv": FEM_CSV,
        "model_json": MODEL_JSON,
        "groups": groups, "cases": case_meta,
        "checks": all_checks,
        "derived": derived,
        "derived_rules": {
            "U_iii_k": "ANSYS *_disp.csv row (node iii, DOF k=1..6)",
            "scale": "max |ANSYS value| of that DOF over all nodes in the case",
            "criterion": "|got-ref| <= atol + rtol*scale",
            "rot_rel": "|ROTZ_tip - theta| / |theta|",
            "closure": "hypot(L + UX_tip, UY_tip)",
            "y_rel": "|ymid - y_exact| / |y_exact|",
        },
        "criteria": criteria,
        "ansys_log_expect": ANSYS_LOG_EXPECT,
        "fem_exe": "beam_nl_probe.exe",
        "fem_exe_args": ["<beam_file>", "<out_csv>"],
    }


def check(doc: dict) -> bool:
    groups, cases = parse_beam(CASES / BEAM_FILE)
    ok = True
    n = 0
    rebuilt: list[dict] = []
    for c in cases:
        rebuilt.extend(derive(c, groups[c["group"]])["checks"])
    if len(rebuilt) != len(doc["checks"]):
        print("[FAIL] check count {} vs {}".format(len(rebuilt), len(doc["checks"])))
        return False
    for a, b in zip(rebuilt, doc["checks"]):
        n += 1
        if a["case"] != b["case"] or a["key"] != b["key"] or \
           abs(a["ref"] - b["ref"]) > 1e-12 + 1e-12 * abs(b["ref"]):
            print("[FAIL] {}/{} ref {:.12g} vs {:.12g}".format(
                a["case"], a["key"], a["ref"], b["ref"]))
            ok = False
    crit = doc["criteria"]
    cur = build()
    for k, v in crit.items():
        if "value" in v and abs(cur["criteria"][k]["value"] - v["value"]) > \
                1e-12 + 1e-12 * abs(v["value"]):
            print("[FAIL] criterion {} drift".format(k))
            ok = False
    print("[{}] {} checks + {} criteria consistent with frozen ANSYS CSV".format(
        "PASS" if ok else "FAIL", n, len(crit)))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        doc = json.loads((REF / (CASE_ID + "_ansys.json")).read_text(encoding="utf-8"))
        return 0 if check(doc) else 1
    REF.mkdir(parents=True, exist_ok=True)
    doc = build()
    path = REF / (CASE_ID + "_ansys.json")
    path.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    print("[GEN ] {} -> {} checks, {} cases, {} criteria".format(
        path.name, len(doc["checks"]), len(doc["cases"]), len(doc["criteria"])))
    return 0 if check(doc) else 1


if __name__ == "__main__":
    raise SystemExit(main())
