#!/usr/bin/env python3
"""从冻结的 cable_solver STABILITY=1 输出生成第 20 项 CTest 黄金值。

数据源只有两类：
  1. cases/*_stability.txt —— 本机 ANSYS 18.0 无关的 FEM 内核输出
     （内核口径已由 CB3D17 与 ANSYS ARCLEN 对过标）；
  2. von Mises 两杆浅桁架的**闭式解析** P(v)（独立 oracle，不读求解器）。

用法：
    python build_reference_stable.py            # 生成 reference/cb3d18_stable_path.json
    python build_reference_stable.py --check   # 只校验 JSON 与冻结输出一致
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference" / "cb3d18_stable_path.json"

E, A, B, H, PREF = 2.1e11, 0.001, 2.0, 0.5, 2.0e6
L0 = math.hypot(B, H)

# 每个 (文件名, 是否应有平衡点, 预期 ERROR 码或 None)
CASE_SPEC = [
    ("von_mises_log", True, None),
    ("von_mises_eng", True, None),
    ("von_mises_bad_ctrl", False, "invalid_control_dof"),
    ("von_mises_bad_method", False, "unsupported_path_method"),
    ("von_mises_no_load", False, "zero_reference_load"),
]

_RE_SUM = re.compile(
    r"^NSTEP=(\d+) NPOINT=(\d+) CONVERGED=([01]) "
    r"CTRLNODE=(\d+) CTRLDOF=(\d+) STRAIN=([01]) PSI=(\S+)$")
_RE_FINAL = re.compile(r"^LAMBDA_FINAL=(\S+) ERROR=(\S+)$")
_RE_POINT = re.compile(
    r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+) "
    r"ITER=(\d+) DS=(\S+) PIVOT=(\S+) CONVERGED=([01])$")


def p_of_v(v: float, strain: int) -> float:
    """独立解析锚点：两杆浅桁架竖向荷载–位移闭式解。

    L(v) = sqrt(B^2 + (H-v)^2)，N = EA*ln(L/L0)（strain=1）
    或 EA*(L-L0)/L0（strain=0），P(v) = -2*N*(H-v)/L。
    P(H)=0，P(2H)=0；极值点由黄金分割求 dP/dv=0。
    """
    length = math.hypot(B, H - v)
    n = E * A * (math.log(length / L0) if strain
                 else (length - L0) / L0)
    return -2.0 * n * (H - v) / length


def analytic_critical(strain: int) -> tuple[float, float]:
    lo, hi = 0.0, H
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.3819660112501051
        m2 = lo + (hi - lo) * 0.6180339887498949
        if p_of_v(m1, strain) < p_of_v(m2, strain):
            lo = m1
        else:
            hi = m2
    v = 0.5 * (lo + hi)
    return v, p_of_v(v, strain)


def parse_stability(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    out: dict = {"points": [], "error": "", "has_buckling_error": 0,
                 "kind": "", "nstep": 0, "npoint": 0, "converged": 0,
                 "ctrlnode": 0, "ctrldof": 0, "strain": -1, "psi": -1.0,
                 "lambda_final": None}
    for line in lines:
        s = line.strip()
        if s.startswith("[BUCKLING]") or s.startswith("ERROR=nonlinear_state_not_supported"):
            if "nonlinear_state_not_supported" in s:
                out["has_buckling_error"] = 1
        if s.startswith("KIND="):
            out["kind"] = s.split()[0].split("=", 1)[1]
        m = _RE_SUM.match(s)
        if m:
            out.update(nstep=int(m.group(1)), npoint=int(m.group(2)),
                       converged=int(m.group(3)), ctrlnode=int(m.group(4)),
                       ctrldof=int(m.group(5)), strain=int(m.group(6)),
                       psi=float(m.group(7)))
            continue
        m = _RE_FINAL.match(s)
        if m:
            out["lambda_final"] = float(m.group(1))
            out["error"] = "" if m.group(2) == "none" else m.group(2)
            continue
        if s.startswith("ERROR="):
            out["error"] = s.split("=", 1)[1]
            continue
        m = _RE_POINT.match(s)
        if m:
            out["points"].append({
                "step": int(m.group(1)), "lambda": float(m.group(2)),
                "v": -float(m.group(3)), "res": float(m.group(4)),
                "iter": int(m.group(5)), "ds": float(m.group(6)),
                "pivot": float(m.group(7)),
            })
    return out


def build() -> dict:
    checks: list[list] = []
    strings: dict[str, str] = {}
    for stem, want_points, want_error in CASE_SPEC:
        r = parse_stability(CASES / f"{stem}_stability.txt")
        if bool(r["points"]) != want_points:
            raise SystemExit(f"[FATAL] {stem}: 平衡点存在性与预期不符")
        got_error = r["error"] or None
        if got_error != want_error:
            raise SystemExit(f"[FATAL] {stem}: ERROR={got_error} != {want_error}")
        checks += [
            [f"{stem}/nstep", r["nstep"]],
            [f"{stem}/npoint", r["npoint"]],
            [f"{stem}/converged", r["converged"]],
            [f"{stem}/has_buckling_error", r["has_buckling_error"]],
            [f"{stem}/n_points_parsed", len(r["points"])],
        ]
        strings[f"{stem}/error"] = r["error"] or "none"
        if r["points"]:
            checks += [
                [f"{stem}/ctrlnode", r["ctrlnode"]],
                [f"{stem}/ctrldof", r["ctrldof"]],
                [f"{stem}/strain", r["strain"]],
                [f"{stem}/psi", r["psi"]],
                [f"{stem}/lambda_final", r["lambda_final"]],
            ]
            for p in r["points"]:
                pre = f"{stem}/p{p['step']}"
                checks += [[f"{pre}/lambda", p["lambda"]],
                           [f"{pre}/v", p["v"]],
                           [f"{pre}/res", p["res"]],
                           [f"{pre}/iter", p["iter"]],
                           [f"{pre}/ds", p["ds"]],
                           [f"{pre}/pivot", p["pivot"]]]
    criteria = {
        "E": E, "A": A, "B": B, "H": H, "PREF": PREF, "L0": L0,
        "analytic_P_at_H": p_of_v(H, 1),
        "analytic_P_at_2H": p_of_v(2.0 * H, 1),
    }
    for strain, tag in ((1, "log"), (0, "eng")):
        vcr, pcr = analytic_critical(strain)
        criteria[f"analytic_v_cr_{tag}"] = vcr
        criteria[f"analytic_P_cr_{tag}"] = pcr
        criteria[f"analytic_lambda_cr_{tag}"] = pcr / PREF
    return {
        "meta": {
            "source": "frozen cable_solver.exe STABILITY=1 output + closed-form "
                      "von Mises two-bar shallow truss (independent oracle)",
            "cases": [c[0] for c in CASE_SPEC],
            "n_checks": len(checks),
            "n_strings": len(strings),
            "n_criteria": len(criteria),
            "note": "CB3D17 snap_probe / ANSYS ARCLEN 同源；主链路 opt-in 默认全关",
        },
        "checks": checks,
        "strings": strings,
        "criteria": criteria,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="只校验现有 JSON 与冻结输出一致，不写文件")
    a = ap.parse_args()
    doc = build()
    if a.check:
        cur = json.loads(REF.read_text(encoding="utf-8"))
        ok = (cur["checks"] == doc["checks"]
              and cur["strings"] == doc["strings"]
              and all(abs(cur["criteria"][k] - v) <= 1e-12 * max(1.0, abs(v))
                      for k, v in doc["criteria"].items()))
        print(f"[{'PASS' if ok else 'FAIL'}] {doc['meta']['n_checks']} checks + "
              f"{doc['meta']['n_strings']} strings + "
              f"{doc['meta']['n_criteria']} criteria consistent with frozen output")
        return 0 if ok else 1
    REF.parent.mkdir(parents=True, exist_ok=True)
    REF.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"[GEN ] {REF.name} -> {doc['meta']['n_checks']} checks, "
          f"{doc['meta']['n_strings']} strings, "
          f"{doc['meta']['n_criteria']} criteria, "
          f"{len(doc['meta']['cases'])} cases")
    return 0


if __name__ == "__main__":
    sys.exit(main())
