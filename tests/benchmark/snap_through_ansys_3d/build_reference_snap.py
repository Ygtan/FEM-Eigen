#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS 逐子步 CSV 生成 CB3D17 snap-through 黄金参考值。

CB3D17 = FEM-Eigen ``stability_path_3d.hpp`` 球面弧长法 vs ANSYS 18.0
LINK180 + NLGEOM + ARCLEN,ON。独立解析锚点是 von Mises 两杆浅桁架的
``P(v)`` 闭式曲线（含极限点、下降段、零点、反向荷载段）。

每一项 ``check`` 都是真实的 ANSYS CSV 值（控制点位移 v、反力反算的
lambda、ANSYS 内部 TIME 荷载因子）。``compare_snap_through_ansys_3d.py``
离线重跑 snap_probe.exe 复现冻结 FEM CSV，再用解析锚点与二次插值完成
FEM/ANSYS 曲线互比；``--live-ansys`` 重跑 3 个 job 逐 CSV 单元 diff。

用法：
    python build_reference_snap.py            # 生成 reference/*.json
    python build_reference_snap.py --check    # 只校验一致性
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference"

CASE_ID = "cb3d17_snap_through_arc_length"
SNAP_FILE = CASE_ID + ".snap"
FEM_CSV = CASE_ID + "_fem.csv"
MODEL_JSON = CASE_ID + ".model.json"
REF_FILE = "cb3d17_snap_through_ansys.json"

TITLE = ("CB3D17 post-buckling / snap-through: FEM spherical arc-length "
         "vs ANSYS 18.0 LINK180 + ARCLEN")
DESCRIPTION = (
    "FEM-Eigen stability_path_3d.hpp implements a Crisfield spherical "
    "arc-length path follower (load factor lambda as an unknown, consistent "
    "tangent reused from fem_3d_nl.hpp co-rotational elements, quadratic "
    "constraint with correct-root selection). ANSYS 18.0 uses LINK180 with "
    "KEYOPT(2)=1 (logarithmic strain), NLGEOM,ON, ARCLEN,ON and ARCTRM,U; "
    "per-substep apex UY and support FY reactions are written to CSV. "
    "The von Mises two-bar shallow truss has an independent analytic "
    "equilibrium path P(v)=-2*N*(h-v)/L with N=EA*ln(L/L0) for the log "
    "strain measure and N=EA*(L-L0)/L0 for the FEM kernel engineering "
    "strain. The path passes the limit point, the descending branch, "
    "lambda=0 at v=h, negative loads, and lambda=0 again at v=2h. "
    "KEYOPT(2)=0 is archived as a strain-measure contrast only."
)

# 与 run_cases_snap.py 的判据同源
TOL_DP_ANSYS = 1.0e-3      # N
TOL_DP_FEM = 1.0e-3        # N
TOL_LAMBDA_TIME = 1.0e-5
TOL_EXTREMUM = 1.0e-3
TOL_CURVE = 5.0e-5

ANSYS_CASES = [
    {"id": "sn_log_fine", "keyopt2": 1, "maxarc": 0.005, "fem": "log_fine"},
    {"id": "sn_log_coarse", "keyopt2": 1, "maxarc": 0.02, "fem": "log_coarse"},
    {"id": "sn_eng", "keyopt2": 0, "maxarc": 0.005, "fem": None},
]


def parse_snap(path: Path) -> tuple[dict, list[dict]]:
    cfg: dict = {}
    fem: list[dict] = []
    for line in path.read_text(encoding="ascii").splitlines():
        p = line.split()
        if not p or p[0].startswith("#"):
            continue
        if p[0] == "G":
            for t in p[2:]:
                k, _, v = t.partition("=")
                try:
                    cfg[k] = float(v)
                except ValueError:
                    pass
        elif p[0] == "C":
            d = {"id": p[1], "group": p[2]}
            for t in p[3:]:
                k, _, v = t.partition("=")
                if k in ("DS", "PSI", "STRAIN"):
                    d[k] = float(v)
                elif k == "TAG":
                    d[k] = v
            fem.append(d)
    return cfg, fem


def read_curve_csv(path: Path, pref: float, tend: float):
    """冻结 CSV 5 列：iset, UY, R1, R3, TIME -> (v, lam_r, lam_t, iset)。"""
    out = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        p = line.split()
        if len(p) < 5:
            continue
        try:
            iset = float(p[0])
            uy, r1, r3, tnow = float(p[1]), float(p[2]), float(p[3]), float(p[4])
        except ValueError:
            continue
        out.append((-uy, (r1 + r3) / pref, tnow / tend, iset))
    return out


def analytic_p(v: float, cfg: dict, strain: int) -> float:
    B, H = cfg["B"], cfg["H"]
    L0 = math.hypot(B, H)
    L = math.hypot(B, H - v)
    N = (cfg["E"] * cfg["A"] * math.log(L / L0) if strain == 1
         else cfg["E"] * cfg["A"] * (L - L0) / L0)
    return -2.0 * N * (H - v) / L


def analytic_limit(cfg: dict, strain: int) -> tuple[float, float]:
    lo, hi = 0.0, cfg["H"]
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.3819660112501051
        m2 = lo + (hi - lo) * 0.6180339887498949
        if analytic_p(m1, cfg, strain) < analytic_p(m2, cfg, strain):
            lo = m1
        else:
            hi = m2
    v = 0.5 * (lo + hi)
    return v, analytic_p(v, cfg, strain)


def quadratic_extremum(curve, find_max: bool):
    if not curve:
        return 0.0, 0.0
    pts = sorted(curve, key=lambda x: x[0])
    vals = [p[1] for p in pts]
    i = (max if find_max else min)(range(len(vals)), key=lambda k: vals[k])
    if i <= 0 or i + 1 >= len(pts):
        return pts[i]
    x0, y0 = pts[i - 1]
    x1, y1 = pts[i]
    x2, y2 = pts[i + 1]
    s1 = (y1 - y0) / (x1 - x0)
    s2 = (y2 - y1) / (x2 - x1)
    d2 = (s2 - s1) / (x2 - x0)
    if abs(d2) < 1e-14 or abs(x1 - x0) < 1e-14:
        return x1, y1
    xs = 0.5 * (x0 + x1) - s1 / (2.0 * d2)
    return xs, y0 + s1 * (xs - x0) + d2 * (xs - x0) * (xs - x1)


def derive(job: dict, cfg: dict, tend: float) -> dict:
    cid = job["id"]
    rows = read_curve_csv(CASES / (cid + "_curve.csv"), cfg["PREF"], tend)
    if len(rows) < 10:
        raise ValueError("{}: too few ANSYS substeps ({})".format(cid, len(rows)))
    checks = []
    for i, (v, lam_r, lam_t, iset) in enumerate(rows, start=1):
        checks.append(dict(job=cid, row=i, key="V", ref=v,
                           analytic=None, tol=None, source="ansys_csv"))
        # KEYOPT(2)=0 是另一种应变度量，只作对照：不给解析期望值，
        # 避免把 ~2% 的本构差异误判成算法失败。
        lam_exact = (analytic_p(v, cfg, 1) / cfg["PREF"]
                     if job["keyopt2"] == 1 else None)
        checks.append(dict(job=cid, row=i, key="LAM", ref=lam_r,
                           analytic=lam_exact, tol=None, source="ansys_csv"))
        checks.append(dict(job=cid, row=i, key="TIME", ref=lam_t,
                           analytic=lam_r, tol=None, source="ansys_csv"))
    curve = [(v, lam_r) for v, lam_r, _, _ in rows]
    _, lmaxfit = quadratic_extremum(curve, True)
    vminfit, lminfit = quadratic_extremum(curve, False)
    worst_dp = max(abs(lam_r * cfg["PREF"] - analytic_p(v, cfg, 1))
                   for v, lam_r, _, _ in rows if 0.0 <= v <= 2 * cfg["H"])
    return {
        "checks": checks,
        "derived": {
            "nset": len(rows),
            "v_max": max(v for v, _, _, _ in rows),
            "lambda_max": max(lam for _, lam, _, _ in rows),
            "lambda_min": min(lam for _, lam, _, _ in rows),
            "lambda_max_fit": lmaxfit,
            "lambda_min_fit": lminfit,
            "v_at_lambda_min_fit": vminfit,
            "worst_dP": worst_dp,
            "max_lambda_reaction_minus_time": max(
                abs(lam_r - lam_t) for _, lam_r, lam_t, _ in rows),
        },
        "curve": [[v, lam_r] for v, lam_r, _, _ in rows],
    }


def build() -> dict:
    cfg, fem_cases = parse_snap(CASES / SNAP_FILE)
    model = json.loads((CASES / MODEL_JSON).read_text(encoding="utf-8"))
    tend = float(model["config"]["TIME_END"])
    vcr_e, pcr_e = analytic_limit(cfg, 0)
    vcr_l, pcr_l = analytic_limit(cfg, 1)
    all_checks, case_meta, derived, curves = [], [], {}, {}
    for job in ANSYS_CASES:
        d = derive(job, cfg, tend)
        all_checks.extend(d["checks"])
        derived[job["id"]] = d["derived"]
        curves[job["id"]] = d["curve"]
        case_meta.append(dict(
            id=job["id"], keyopt2=job["keyopt2"], maxarc=job["maxarc"],
            fem=job["fem"], apdl=job["id"] + ".inp",
            csv=job["id"] + "_curve.csv", out=job["id"] + "_ansys.out"))
    main = derived["sn_log_fine"]
    criteria = {
        "ansys_vs_analytic_worst_dP": dict(
            value=main["worst_dP"], limit=TOL_DP_ANSYS, unit="N",
            note="ANSYS KEYOPT(2)=1 per-substep lambda*P_ref vs analytic P_log(v)"),
        "ansys_lambda_reaction_vs_time": dict(
            value=main["max_lambda_reaction_minus_time"], limit=TOL_LAMBDA_TIME,
            note="support-reaction lambda vs ANSYS internal TIME load factor"),
        "ansys_limit_point_fit": dict(
            value=abs(main["lambda_max_fit"] - pcr_l / cfg["PREF"]),
            limit=TOL_EXTREMUM,
            note="quadratic-fit limit load vs analytic dP/dv=0"),
        "descending_branch": dict(
            lambda_min_fit=main["lambda_min_fit"], v_max=main["v_max"],
            lambda_min_limit=-0.5, v_min_limit=0.99,
            note="lambda decreases through zero into negative loads; v reaches 2h"),
        "fem_vs_analytic_worst_dP": dict(
            limit=TOL_DP_FEM, unit="N",
            note="fresh snap_probe.exe path vs analytic P_strain(v)"),
        "fem_limit_point_fit": dict(
            limit=TOL_EXTREMUM,
            note="fresh FEM quadratic-fit limit load vs analytic"),
        "fem_ansys_curve_worst": dict(
            limit=TOL_CURVE,
            note="quadratic interpolation of FEM path at ANSYS substep v"),
        "strain_measure_contrast": dict(
            value=100 * (derived["sn_eng"]["lambda_max"]
                         - derived["sn_log_fine"]["lambda_max"])
            / derived["sn_log_fine"]["lambda_max"],
            lo=1.5, hi=3.0, unit="%",
            note="KEYOPT(2)=0 vs =1 limit-load difference (constitutive, not algorithmic)"),
        "ansys_log_errors": dict(limit=0, note="all archived ANSYS logs"),
    }
    return {
        "id": CASE_ID, "title": TITLE, "description": DESCRIPTION,
        "snap_file": SNAP_FILE, "fem_csv": FEM_CSV, "model_json": MODEL_JSON,
        "config": cfg, "time_end": tend,
        "fem_cases": fem_cases,
        "cases": case_meta,
        "analytic": {
            "v_cr_eng": vcr_e, "P_cr_eng": pcr_e,
            "v_cr_log": vcr_l, "P_cr_log": pcr_l,
            "lambda_cr_eng": pcr_e / cfg["PREF"],
            "lambda_cr_log": pcr_l / cfg["PREF"],
        },
        "checks": all_checks,
        "derived": derived,
        "curves": curves,
        "tolerances": dict(
            dp_ansys=TOL_DP_ANSYS, dp_fem=TOL_DP_FEM,
            lambda_time=TOL_LAMBDA_TIME, extremum=TOL_EXTREMUM,
            curve=TOL_CURVE),
        "criteria": criteria,
        "derived_rules": {
            "V_iii": "apex displacement v = -UY from ANSYS substep iii CSV",
            "LAM_iii": "load factor from support reactions (R1+R3)/P_ref",
            "TIME_iii": "load factor from ANSYS internal TIME/TIME_end",
            "analytic": "independent P_log(v)/P_ref for LAM, LAM itself for TIME",
            "criterion": "|ref-analytic| <= limit, plus FEM reproduction and curve match",
        },
        "fem_exe": "snap_probe.exe",
        "fem_exe_args": ["<snap_file>", "<out_csv>"],
    }


def check(doc: dict) -> bool:
    cfg, _ = parse_snap(CASES / SNAP_FILE)
    ok = True
    rebuilt = []
    for job in ANSYS_CASES:
        rebuilt.extend(derive(job, cfg, doc["time_end"])["checks"])
    if len(rebuilt) != len(doc["checks"]):
        print("[FAIL] check count {} vs {}".format(len(rebuilt), len(doc["checks"])))
        return False
    for a, b in zip(rebuilt, doc["checks"]):
        if a["job"] != b["job"] or a["row"] != b["row"] or a["key"] != b["key"] \
           or abs(a["ref"] - b["ref"]) > 1e-12 + 1e-12 * abs(b["ref"]):
            print("[FAIL] {}/{}/{} ref {:.12g} vs {:.12g}".format(
                a["job"], a["row"], a["key"], a["ref"], b["ref"]))
            ok = False
    cur = build()
    for k, v in doc["criteria"].items():
        if "value" in v and "value" in cur["criteria"][k]:
            if abs(cur["criteria"][k]["value"] - v["value"]) > \
                    1e-12 + 1e-12 * abs(v["value"]):
                print("[FAIL] criterion {} drift".format(k))
                ok = False
    print("[{}] {} checks + {} criteria consistent with frozen ANSYS CSV".format(
        "PASS" if ok else "FAIL", len(rebuilt), len(doc["criteria"])))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        doc = json.loads((REF / REF_FILE).read_text(encoding="utf-8"))
        return 0 if check(doc) else 1
    REF.mkdir(parents=True, exist_ok=True)
    doc = build()
    path = REF / REF_FILE
    path.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    print("[GEN ] {} -> {} checks, {} cases, {} criteria".format(
        path.name, len(doc["checks"]), len(doc["cases"]), len(doc["criteria"])))
    return 0 if check(doc) else 1


if __name__ == "__main__":
    raise SystemExit(main())
