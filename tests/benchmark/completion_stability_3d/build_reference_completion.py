#!/usr/bin/env python3
"""从冻结的 cable_solver.exe 输出生成第 22 项 CTest 黄金值。

数据源只有两类：
  1. cases/*_cable.txt / *_stability.txt —— 本机 ANSYS 无关的 FEM 内核输出
     （`STABASIS=COMPLETION` 成桥态基准 + 第 ① 阶段三类负例 + 默认零污染）；
  2. V 形双索单自由度体系的**独立闭式解**（不读求解器）：
         L(v)  = hypot(B, H+v)
         N(v)  = EA/L0 * (L(v) - L0)
         P(v)  = 2 N(v) (H+v) / L(v)
         F0    = P(V0)
         λ(v)  = [P(v) - F0] / PINC
     以及成桥态切线刚度
         K_yy  = 2[(EA/L0)c_y² + (N/L)(1-c_y²)],  c_y = (H+v)/L
     给出 dλ/dv = K_yy/PINC。

用法：
    python build_reference_completion.py                     # 生成 reference JSON
    python build_reference_completion.py --check            # 只校验 JSON 自洽
    python build_reference_completion.py --cases <dir>      # 指定冻结 cases 目录
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
REF = HERE / "reference" / "cb3d20_completion_stability.json"

# ---- 与 cases/*.dat 同源的几何 / 材料常量（唯一数据源在 .dat 头注释）----
B, H, E, A = 2.0, 0.6, 2.1e11, 1.0e-3
EA = E * A
L_REF = math.hypot(B, H)
L0 = 0.995 * L_REF
V0 = 0.10
PINC = 5.0e5

# (stem, 输出后缀, 是否应有平衡点, 预期 ERROR 码或 None, 预期 BASIS)
CASE_SPEC = [
    ("vshape_static", "_cable.txt", True, None, "completion_state"),
    ("vshape_tension", "_cable.txt", True, None, "completion_state"),
    ("vshape_no_load", "_cable.txt", False, "zero_incremental_load",
     "completion_state"),
    ("vshape_bad_basis", "_stability.txt", False,
     "unsupported_stability_basis", "unknown"),
    ("vshape_unconverged", "_cable.txt", False, "completion_not_converged",
     "completion_state"),
    ("vshape_default", "_cable.txt", False, "no_stability_section", ""),
]

_RE_SUM = re.compile(
    r"^NSTEP=(\d+) NPOINT=(\d+) CONVERGED=([01]) "
    r"CTRLNODE=(\d+) CTRLDOF=(\d+) STRAIN=([01]) PSI=(\S+)$")
_RE_STATE = re.compile(
    r"^STATE CTRL0=(\S+) F0_NORM=(\S+) FINC_NORM=(\S+) RES0=(\S+)$")
_RE_FINAL = re.compile(r"^LAMBDA_FINAL=(\S+) ERROR=(\S+)$")
_RE_POINT = re.compile(
    r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+) "
    r"ITER=(\d+) DS=(\S+) PIVOT=(\S+) CONVERGED=([01])$")


# ---- 独立闭式 oracle ---------------------------------------------------

def N_of_v(v: float) -> float:
    """索轴力（工程应变 N = EA(L−L0)/L0，与 fem_3d_nl.hpp 同口径）。"""
    return EA / L0 * (math.hypot(B, H + v) - L0)


def P_of_v(v: float) -> float:
    """低点竖向荷载–位移闭式解：两根对称索的竖向合力。"""
    lv = math.hypot(B, H + v)
    return 2.0 * N_of_v(v) * (H + v) / lv


F0 = P_of_v(V0)
TSTAR = N_of_v(V0)


def lambda_of_v(v: float) -> float:
    """成桥态基准：λ 只放大增量荷载，F0 是恒定不缩放的成桥态恒载。"""
    return (P_of_v(v) - F0) / PINC


def k_yy_of_v(v: float) -> float:
    """成桥态 y 向一致切线（两根轴向索，共转解析切线，独立推导）。"""
    lv = math.hypot(B, H + v)
    cy = (H + v) / lv
    return 2.0 * ((EA / L0) * cy * cy
                  + (N_of_v(v) / lv) * (1.0 - cy * cy))


def numeric_dP_dv(v: float, h: float = 1.0e-6) -> float:
    return (P_of_v(v + h) - P_of_v(v - h)) / (2.0 * h)


# ---- 冻结输出解析 -------------------------------------------------------

def parse_result(path: Path) -> dict:
    """从完整结果文件里切出 [STABILITY] 段并逐字段解析。

    段外只统计 [BUCKLING] ERROR / [NEWTON] / [LENGTH] / [U]，用于证明
    「先走完成桥态闭环、再挂稳定段」，不在这里重算任何力学量。
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    out: dict = {
        "has_section": False, "kind": "", "status": "", "basis": "",
        "pathmethod": "", "nstep": 0, "npoint": 0, "converged": 0,
        "ctrlnode": 0, "ctrldof": 0, "strain": -1, "psi": -1.0,
        "lambda_final": None, "error": "", "points": [],
        "ctrl0": None, "f0_norm": None, "finc_norm": None, "res0": None,
        "has_buckling_error": 0, "has_newton": 0, "has_length": 0,
        "u_ctrl": None,
    }
    sec = ""
    in_stab = False
    for raw in text.splitlines():
        s = raw.strip()
        if s.startswith("["):
            sec = s
            in_stab = (s == "[STABILITY]")
            if in_stab:
                out["has_section"] = True
            continue
        if s.startswith("ERROR=nonlinear_state_not_supported"):
            out["has_buckling_error"] = 1
        if sec == "[NEWTON]" and s.startswith("KIND=newton_cr"):
            out["has_newton"] = 1
        if sec == "[LENGTH]" and s.startswith("KIND=cable_length"):
            out["has_length"] = 1
        if sec == "[U]":
            p = s.split()
            if len(p) >= 3 and p[0] == "2" and p[1] == "2":
                out["u_ctrl"] = float(p[2])
        if not in_stab:
            continue
        if s.startswith("KIND="):
            for kv in s.split():
                k, _, v = kv.partition("=")
                if k == "KIND":
                    out["kind"] = v
                elif k == "STATUS":
                    out["status"] = v
                elif k == "BASIS":
                    out["basis"] = v
                elif k == "PATHMETHOD":
                    out["pathmethod"] = v
            continue
        m = _RE_SUM.match(s)
        if m:
            out.update(nstep=int(m.group(1)), npoint=int(m.group(2)),
                       converged=int(m.group(3)), ctrlnode=int(m.group(4)),
                       ctrldof=int(m.group(5)), strain=int(m.group(6)),
                       psi=float(m.group(7)))
            continue
        m = _RE_STATE.match(s)
        if m:
            out.update(ctrl0=float(m.group(1)), f0_norm=float(m.group(2)),
                       finc_norm=float(m.group(3)), res0=float(m.group(4)))
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
                "ctrl": float(m.group(3)), "res": float(m.group(4)),
                "iter": int(m.group(5)), "ds": float(m.group(6)),
                "pivot": float(m.group(7)),
            })
    if not out["has_section"]:
        out["error"] = "no_stability_section"
        out["basis"] = ""
    return out


def build(cases_dir: Path | None = None) -> dict:
    cd = cases_dir or CASES
    checks: list[list] = []
    strings: dict[str, str] = {}
    for stem, suffix, want_points, want_error, want_basis in CASE_SPEC:
        p = cd / f"{stem}{suffix}"
        if not p.is_file():
            raise SystemExit(f"[FATAL] 缺少冻结输出: {p}")
        r = parse_result(p)
        if bool(r["points"]) != want_points:
            raise SystemExit(f"[FATAL] {stem}: 平衡点存在性与预期不符")
        got_error = r["error"] or None
        if got_error != want_error:
            raise SystemExit(f"[FATAL] {stem}: ERROR={got_error} != {want_error}")
        if want_basis and r["basis"] != want_basis:
            raise SystemExit(f"[FATAL] {stem}: BASIS={r['basis']} != {want_basis}")
        pre = stem
        checks += [
            [f"{pre}/has_section", int(r["has_section"])],
            [f"{pre}/nstep", r["nstep"]],
            [f"{pre}/npoint", r["npoint"]],
            [f"{pre}/converged", r["converged"]],
            [f"{pre}/has_buckling_error", r["has_buckling_error"]],
            [f"{pre}/has_newton", r["has_newton"]],
            [f"{pre}/has_length", r["has_length"]],
            [f"{pre}/n_points_parsed", len(r["points"])],
        ]
        if r["u_ctrl"] is not None:
            checks.append([f"{pre}/u_ctrl", r["u_ctrl"]])
        strings[f"{pre}/error"] = r["error"] or "none"
        strings[f"{pre}/basis"] = r["basis"]
        strings[f"{pre}/kind"] = r["kind"]
        if r["points"]:
            checks += [
                [f"{pre}/ctrlnode", r["ctrlnode"]],
                [f"{pre}/ctrldof", r["ctrldof"]],
                [f"{pre}/strain", r["strain"]],
                [f"{pre}/psi", r["psi"]],
                [f"{pre}/lambda_final", r["lambda_final"]],
                [f"{pre}/ctrl0", r["ctrl0"]],
                [f"{pre}/f0_norm", r["f0_norm"]],
                [f"{pre}/finc_norm", r["finc_norm"]],
                [f"{pre}/res0", r["res0"]],
            ]
            for q in r["points"]:
                qp = f"{pre}/p{q['step']}"
                checks += [[f"{qp}/lambda", q["lambda"]],
                           [f"{qp}/ctrl", q["ctrl"]],
                           [f"{qp}/res", q["res"]],
                           [f"{qp}/iter", q["iter"]],
                           [f"{qp}/ds", q["ds"]],
                           [f"{qp}/pivot", q["pivot"]]]
    # 独立闭式 criteria（一个字节都不读求解器）
    criteria = {
        "B": B, "H": H, "E": E, "A": A, "EA": EA,
        "L_ref": L_REF, "L0": L0, "V0": V0, "PINC": PINC,
        "analytic_F0": F0,
        "analytic_Tstar": TSTAR,
        "analytic_N_ref": N_of_v(0.0),
        "analytic_P_at_V0": P_of_v(V0),
        "analytic_lambda_at_V0": lambda_of_v(V0),
        "analytic_lambda_at_V0_plus_0.05": lambda_of_v(V0 + 0.05),
        "analytic_k_yy_at_V0": k_yy_of_v(V0),
        "analytic_dP_dv_at_V0": numeric_dP_dv(V0),
        "analytic_dlambda_dv_at_V0": k_yy_of_v(V0) / PINC,
    }
    return {
        "meta": {
            "source": "frozen cable_solver.exe STABASIS=COMPLETION output + "
                      "closed-form V-shape twin-cable single-DOF oracle "
                      "(independent, does not read the solver)",
            "cases": [c[0] for c in CASE_SPEC],
            "n_checks": len(checks),
            "n_strings": len(strings),
            "n_criteria": len(criteria),
            "note": "第十九次 ①：成桥态基准 opt-in，默认全关、REFERENCE 逐位不变；"
                    "B3 能力矩阵仍 P 不 Y（ANSYS 两段对标属阶段 ②）",
        },
        "checks": checks,
        "strings": strings,
        "criteria": criteria,
    }


def check(doc: dict, cases_dir: Path | None = None) -> tuple[bool, str]:
    fresh = build(cases_dir)
    ok = (doc["checks"] == fresh["checks"]
          and doc["strings"] == fresh["strings"]
          and all(abs(doc["criteria"][k] - v) <= 1e-12 * max(1.0, abs(v))
                  for k, v in fresh["criteria"].items()))
    return ok, (f"{fresh['meta']['n_checks']} checks + "
                f"{fresh['meta']['n_strings']} strings + "
                f"{fresh['meta']['n_criteria']} criteria")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="只校验现有 JSON 与冻结输出一致，不写文件")
    ap.add_argument("--cases", default=None, help="冻结 cases 目录")
    ap.add_argument("--reference", default=None, help="reference JSON 目录")
    a = ap.parse_args()
    cd = Path(a.cases) if a.cases else CASES
    ref = (Path(a.reference) / REF.name) if a.reference else REF
    doc = build(cd)
    if a.check:
        cur = json.loads(ref.read_text(encoding="utf-8"))
        ok, msg = check(cur, cd)
        print(f"[{'PASS' if ok else 'FAIL'}] {msg} consistent with frozen output")
        return 0 if ok else 1
    ref.parent.mkdir(parents=True, exist_ok=True)
    ref.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"[GEN ] {ref.name} -> {doc['meta']['n_checks']} checks, "
          f"{doc['meta']['n_strings']} strings, "
          f"{doc['meta']['n_criteria']} criteria, "
          f"{len(doc['meta']['cases'])} cases")
    return 0


if __name__ == "__main__":
    sys.exit(main())
