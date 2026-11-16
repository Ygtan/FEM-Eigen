#!/usr/bin/env python3
"""从冻结的 cable_solver.exe PATHMETHOD=RIKS 输出生成第 27 项 CTest 黄金值。

数据源只有两类：
  1. cases/*_stage_stability.txt —— 用 mingw64 工具链冻结的 FEM 内核输出；
  2. von Mises 两杆浅桁架的**闭式解析** P(v)（独立 oracle，不读求解器）。

用法：
    python build_reference_riks_stage.py
    python build_reference_riks_stage.py --check
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference" / "cb3d24_riks_stage.json"

# von Mises 两杆：STAGE-1 预载 P0，*STABILITY_LOAD 给 PREF；STRAIN=1
E, A, B, H = 2.1e11, 0.001, 2.0, 0.5
P0, PREF = 1.0e6, 2.0e6
L0 = math.hypot(B, H)

CASE_SPEC = [
    ("vm_stage_riks", True, None),
    ("vm_stage_sphere", True, None),
    ("cb3d23_riks", True, None),
    ("st_bad_method", False, "unsupported_path_method"),
]
FROZEN_FILES = [f"{stem}{suffix}" for stem, _, _ in CASE_SPEC
                for suffix in (".dat", "_stage_stability.txt")]

_RE_KIND = re.compile(
    r"^KIND=(\S+) STATUS=(\S+) BASIS=(\S+) PATHMETHOD=(\S+)$")
_RE_SUM = re.compile(
    r"^NSTEP=(\d+) NPOINT=(\d+) CONVERGED=([01]) "
    r"CTRLNODE=(\d+) CTRLDOF=(\d+) STRAIN=(\d+) PSI=(\S+)$")
_RE_STATE = re.compile(
    r"^STATE CTRL0=(\S+) F0_NORM=(\S+) FINC_NORM=(\S+) RES0=(\S+)$")
_RE_STAGE = re.compile(r"^STAGE=(\d+) NELEM=(\d+) NDOF=(\d+)$")
_RE_FINAL = re.compile(r"^LAMBDA_FINAL=(\S+) ERROR=(\S+)$")
_RE_POINT = re.compile(
    r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+) "
    r"ITER=(\d+) DS=(\S+) PIVOT=(\S+) CONVERGED=([01])$")
_RE_STAGE_SUM = re.compile(
    r"^S ID=(\d+) CONVERGED=([01]) NELEM=(\d+) NDLOAD=(\d+) "
    r"NDOF=(\d+) U_NORM=(\S+)$")


def parse(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    out: dict = {
        "kind": "", "status": "", "basis": "", "pathmethod": "",
        "nstep": 0, "npoint": 0, "converged": 0,
        "ctrlnode": 0, "ctrldof": 0, "strain": -1, "psi": -1.0,
        "state": None, "stage": None, "lambda_final": None, "error": "",
        "points": [], "stages": [], "has_buckling_error": 0,
    }
    for line in lines:
        s = line.strip()
        if s.startswith("ERROR=nonlinear_state_not_supported"):
            out["has_buckling_error"] = 1
            continue   # BUCKLING 段的 ERROR 不是 [STABILITY] 的失败码
        m = _RE_KIND.match(s)
        if m:
            out.update(kind=m.group(1), status=m.group(2),
                       basis=m.group(3), pathmethod=m.group(4))
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
            out["state"] = {"ctrl0": float(m.group(1)),
                            "f0_norm": float(m.group(2)),
                            "finc_norm": float(m.group(3)),
                            "res0": float(m.group(4))}
            continue
        m = _RE_STAGE.match(s)
        if m:
            out["stage"] = {"id": int(m.group(1)), "nelem": int(m.group(2)),
                            "ndof": int(m.group(3))}
            continue
        m = _RE_FINAL.match(s)
        if m:
            out["lambda_final"] = float(m.group(1))
            out["error"] = "" if m.group(2) == "none" else m.group(2)
            continue
        if s.startswith("ERROR=") and not out["error"]:
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
            continue
        m = _RE_STAGE_SUM.match(s)
        if m:
            out["stages"].append({
                "id": int(m.group(1)), "converged": int(m.group(2)),
                "nelem": int(m.group(3)), "ndload": int(m.group(4)),
                "ndof": int(m.group(5)), "u_norm": float(m.group(6)),
            })
    return out


def p_of_v(v: float) -> float:
    """独立解析锚点：外荷载幅值 P(v) = -2*N*(H-v)/L，
    STRAIN=1 -> N = EA*ln(L/L0)。起点残差 RES0 = P(V0) - P0 也由它给出。
    """
    length = math.hypot(B, H - v)
    n = E * A * math.log(length / L0)
    return -2.0 * n * (H - v) / length


def analytic_critical() -> tuple[float, float]:
    lo, hi = 0.0, H
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.3819660112501051
        m2 = lo + (hi - lo) * 0.6180339887498949
        if p_of_v(m1) < p_of_v(m2):
            lo = m1
        else:
            hi = m2
    v = 0.5 * (lo + hi)
    return v, p_of_v(v)


def lambda_at_ctrl(points: list[dict], target: float) -> float | None:
    """在一条平衡路径上按 CTRL 线性插值 lambda（两条路径离散点不同）。"""
    pairs = sorted((p["ctrl"], p["lambda"]) for p in points if p["step"] > 0)
    for (x0, y0), (x1, y1) in zip(pairs, pairs[1:]):
        if min(x0, x1) <= target <= max(x0, x1) and x1 != x0:
            t = (target - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return None


def cross_method_worst(riks: dict, sphere: dict) -> tuple[float, int]:
    worst, n = 0.0, 0
    for p in riks["points"]:
        if p["step"] == 0:
            continue
        y = lambda_at_ctrl(sphere["points"], p["ctrl"])
        if y is None:
            continue
        worst = max(worst, abs(p["lambda"] - y))
        n += 1
    return worst, n


def closed_form_worst(riks: dict, v0: float) -> float:
    worst = 0.0
    for p in riks["points"]:
        if p["step"] == 0:
            continue
        v = v0 - p["ctrl"]
        closed = (p_of_v(v) - P0) / PREF
        worst = max(worst, abs(p["lambda"] - closed))
    return worst


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build() -> dict:
    parsed = {stem: parse(CASES / f"{stem}_stage_stability.txt")
              for stem, _, _ in CASE_SPEC}
    checks: list[list] = []
    strings: dict[str, str] = {}
    for stem, want_points, want_error in CASE_SPEC:
        r = parsed[stem]
        got_error = r["error"] or None
        if bool(r["points"]) != want_points or got_error != want_error:
            raise SystemExit(f"[FATAL] {stem}: 平衡点/ERROR 与预期不符 "
                             f"({bool(r['points'])}/{got_error})")
        checks += [
            [f"{stem}/nstep", r["nstep"]],
            [f"{stem}/npoint", r["npoint"]],
            [f"{stem}/converged", r["converged"]],
            [f"{stem}/has_buckling_error", r["has_buckling_error"]],
            [f"{stem}/n_points_parsed", len(r["points"])],
            [f"{stem}/n_stages", len(r["stages"])],
        ]
        strings[f"{stem}/error"] = r["error"] or "none"
        strings[f"{stem}/kind"] = r["kind"]
        strings[f"{stem}/status"] = r["status"]
        strings[f"{stem}/basis"] = r["basis"]
        strings[f"{stem}/pathmethod"] = r["pathmethod"]
        if r["points"]:
            checks += [
                [f"{stem}/ctrlnode", r["ctrlnode"]],
                [f"{stem}/ctrldof", r["ctrldof"]],
                [f"{stem}/strain", r["strain"]],
                [f"{stem}/psi", r["psi"]],
                [f"{stem}/lambda_final", r["lambda_final"]],
            ]
            st = r["state"]
            sg = r["stage"]
            checks += [
                [f"{stem}/state/ctrl0", st["ctrl0"]],
                [f"{stem}/state/f0_norm", st["f0_norm"]],
                [f"{stem}/state/finc_norm", st["finc_norm"]],
                [f"{stem}/state/res0", st["res0"]],
                [f"{stem}/stage/id", sg["id"]],
                [f"{stem}/stage/nelem", sg["nelem"]],
                [f"{stem}/stage/ndof", sg["ndof"]],
            ]
            for s in r["stages"]:
                pre = f"{stem}/st{s['id']}"
                checks += [[f"{pre}/converged", s["converged"]],
                           [f"{pre}/nelem", s["nelem"]],
                           [f"{pre}/ndload", s["ndload"]],
                           [f"{pre}/ndof", s["ndof"]],
                           [f"{pre}/u_norm", s["u_norm"]]]
            for p in r["points"]:
                pre = f"{stem}/p{p['step']}"
                checks += [[f"{pre}/lambda", p["lambda"]],
                           [f"{pre}/ctrl", p["ctrl"]],
                           [f"{pre}/res", p["res"]],
                           [f"{pre}/iter", p["iter"]],
                           [f"{pre}/ds", p["ds"]],
                           [f"{pre}/pivot", p["pivot"]]]

    riks = parsed["vm_stage_riks"]
    sphere = parsed["vm_stage_sphere"]
    cb = parsed["cb3d23_riks"]
    v0 = -riks["state"]["ctrl0"]
    vcr, pcr = analytic_critical()
    cross_worst, cross_n = cross_method_worst(riks, sphere)
    criteria = {
        "E": E, "A": A, "B": B, "H": H, "P0": P0, "PREF": PREF, "L0": L0,
        "V0": v0,
        "analytic_P_at_V0": p_of_v(v0),
        "analytic_RES0": p_of_v(v0) - P0,
        "analytic_v_cr": vcr,
        "analytic_P_cr": pcr,
        "analytic_lambda_cr": (pcr - P0) / PREF,
        "closed_form_worst": closed_form_worst(riks, v0),
        "closed_form_n": len(riks["points"]) - 1,
        "cross_method_worst": cross_worst,
        "cross_method_n": cross_n,
        "lambda_end_riks": riks["lambda_final"],
        "lambda_end_sphere": sphere["lambda_final"],
        "lambda_end_cb3d23_riks": cb["lambda_final"],
    }
    manifest = {name: {"sha256": sha256(CASES / name),
                       "bytes": (CASES / name).stat().st_size}
                for name in FROZEN_FILES}
    return {
        "meta": {
            "source": "frozen cable_solver.exe STAGE+RIKS output + independent "
                      "closed-form von Mises two-bar oracle",
            "cases": [c[0] for c in CASE_SPEC],
            "n_checks": len(checks),
            "n_strings": len(strings),
            "n_criteria": len(criteria),
            "note": "PATHMETHOD=RIKS is Wempner-Riks orthogonal, NOT ANSYS "
                    "ARCLEN / ABAQUS *STATIC,RIKS (spherical). STAGE only.",
        },
        "checks": checks,
        "strings": strings,
        "criteria": criteria,
        "manifest": manifest,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="只校验现有 JSON，不写文件")
    a = ap.parse_args()
    doc = build()
    if a.check:
        cur = json.loads(REF.read_text(encoding="utf-8"))
        ok = (cur["checks"] == doc["checks"]
              and cur["strings"] == doc["strings"]
              and cur["manifest"] == doc["manifest"]
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
          f"{len(doc['meta']['cases'])} cases, "
          f"{len(doc['manifest'])} frozen files")
    print(f"       closed_form_worst={doc['criteria']['closed_form_worst']:.3e} "
          f"cross_method_worst={doc['criteria']['cross_method_worst']:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
