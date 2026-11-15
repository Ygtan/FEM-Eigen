#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 CB3D21 证据生成第 23 项 CTest 黄金值。

冻结证据（cases/，全部由 BridgeStudio/cases/run_cases_completion.py 同源生成）：
  * cb3d21_completion_stability_cable.txt —— FEM 主链路 `[STABILITY]`
    `BASIS=completion_state` 逐平衡点路径 + `[U]` / `[NEWTON]`；
  * cs_default_cable.txt —— 默认零污染负例（无 `[STABILITY]`）；
  * cs_load / cs_arc 的 .inp / _curve.csv / _ansys.out —— ANSYS 两段求解；
  * .cfg.json / .model.json / _compare_report.txt —— 唯一数据源与结论。

独立锚点（criteria，纯 Python、不读求解器）：
  1. 16 根索 `L0*(1+EPS0) == 设计线形弦长`（由 model.json 节点坐标独立算）；
  2. 16 条 `INISTATE,DEFINE,...,EPEL,X,eps` 的 eps 恰等于 EPS0，
     证明 FEM 与 ANSYS 用同一份 L0 / 同一参考构形；
  3. APDL step2 的 `F,1016,FZ,F_end` 恰等于 `F_orig - LAMMAX*PINC`；
  4. λ 换算口径 `lambda=(TIME-1)*LAMMAX` 与 `DeltaUZ=UZ-UZ_step1`
     对 FEM 路径插值的逐子步偏差、初始斜率、终点交叉验证；
  5. 起点闭合 `RES0/|F0|` 与 `RES0/[NEWTON] RESLAST`。

用法：
    python build_reference_bridge.py                  # 生成 reference JSON
    python build_reference_bridge.py --check         # 只校验 JSON 自洽
    python build_reference_bridge.py --cases <dir>   # 指定冻结 cases 目录
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
REF = HERE / "reference" / "cb3d21_completion_bridge.json"

STEM = "cb3d21_completion_stability"
DEFAULT_STEM = "cs_default"
ANSYS_JOBS = ("cs_load", "cs_arc")

_RE_SUM = re.compile(
    r"^NSTEP=(\d+) NPOINT=(\d+) CONVERGED=([01]) "
    r"CTRLNODE=(\d+) CTRLDOF=(\d+) STRAIN=([01]) PSI=(\S+)$")
_RE_STATE = re.compile(
    r"^STATE CTRL0=(\S+) F0_NORM=(\S+) FINC_NORM=(\S+) RES0=(\S+)$")
_RE_FINAL = re.compile(r"^LAMBDA_FINAL=(\S+) ERROR=(\S+)$")
_RE_POINT = re.compile(
    r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+) "
    r"ITER=(\d+) DS=(\S+) PIVOT=(\S+) CONVERGED=([01])$")


def parse_stability(path: Path, ctrl_node: int, ctrl_dof: int) -> dict:
    """切出 [STABILITY] 段逐字段解析；段外只取 [U]/[NEWTON]/[BUCKLING] 标志。"""
    text = path.read_text(encoding="utf-8", errors="replace")
    out: dict = {
        "has_section": False, "kind": "", "status": "", "basis": "",
        "pathmethod": "", "nstep": 0, "npoint": 0, "converged": 0,
        "ctrlnode": 0, "ctrldof": 0, "strain": -1, "psi": -1.0,
        "lambda_final": None, "error": "", "points": [],
        "ctrl0": None, "f0_norm": None, "finc_norm": None, "res0": None,
        "has_buckling_error": 0, "has_newton": 0,
        "reslast": None, "u_ctrl": None,
    }
    sec = ""
    in_stab = False
    m_rl = re.search(r"\[NEWTON\][\s\S]*?RESLAST=(\S+)", text)
    if m_rl:
        out["reslast"] = float(m_rl.group(1))
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
        if sec == "[U]":
            p = s.split()
            if (len(p) >= 3 and p[0] == str(ctrl_node)
                    and p[1] == str(ctrl_dof)):
                out["u_ctrl"] = float(p[2])
        if not in_stab:
            continue
        if s.startswith("KIND="):
            for kv in s.split():
                k, _, v = kv.partition("=")
                if k in ("KIND", "STATUS", "BASIS", "PATHMETHOD"):
                    out[k.lower()] = v
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
                "pivot": float(m.group(7)), "converged": int(m.group(8)),
            })
    if not out["has_section"]:
        out["error"] = "no_stability_section"
        out["basis"] = ""
    return out


def parse_curve(path: Path) -> list[tuple[float, float, float, float]]:
    """ANSYS POST1 CSV：[(iset, TIME, UZ, 0)]。"""
    rows: list[tuple[float, float, float, float]] = []
    for ln in path.read_text(encoding="ascii", errors="replace").splitlines():
        q = ln.split()
        if len(q) < 3:
            continue
        try:
            rows.append((float(q[0]), float(q[1]), float(q[2]),
                         float(q[3]) if len(q) > 3 else 0.0))
        except ValueError:
            continue
    return rows


def split_curve(rows):
    s1 = [p for p in rows if p[1] <= 1.0 + 1e-9]
    s2 = [p for p in rows if p[1] > 1.0 + 1e-9]
    return (s1[-1][2] if s1 else None), s1, s2


def interp_ctrl(fem: list[tuple[float, float]], lam: float) -> float:
    if not fem:
        raise ValueError("FEM 路径为空")
    if lam <= fem[0][0]:
        return fem[0][1]
    if lam >= fem[-1][0]:
        return fem[-1][1]
    lo, hi = 0, len(fem) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if fem[mid][0] <= lam:
            lo = mid
        else:
            hi = mid
    a, b = fem[lo], fem[hi]
    return a[1] + (b[1] - a[1]) * (lam - a[0]) / (b[0] - a[0])


def parse_apdl(path: Path) -> dict:
    txt = path.read_text(encoding="ascii", errors="replace")
    loads = [float(x) for x in re.findall(r"^F,1016,FZ,(\S+)", txt, re.M)]
    eps = [float(x) for x in re.findall(
        r"^INISTATE,DEFINE,\d+,,EPEL,X,(\S+)", txt, re.M)]
    times = [float(x) for x in re.findall(r"^TIME,(\S+)", txt, re.M)]
    return {
        "f_orig": loads[0] if loads else None,
        "f_end": loads[-1] if len(loads) > 1 else None,
        "n_loads": len(loads), "eps": eps, "times": times,
        "time2": times[-1] if times else None,
        "arclen": "ARCLEN,ON" in txt,
        "n_solve": len(re.findall(r"^SOLVE\s*$", txt, re.M)),
        "n_inistate": len(eps),
    }


def parse_log(path: Path) -> tuple[int, int]:
    txt = path.read_text(encoding="ascii", errors="replace")
    def grab(label: str) -> int:
        m = re.findall(r"NUMBER OF %s\s+MESSAGES ENCOUNTERED=\s*(\d+)" % label, txt)
        return int(m[-1]) if m else -1
    return grab("ERROR"), grab("WARNING")


def build(cases_dir: Path | None = None) -> dict:
    cd = cases_dir or CASES
    cfg = json.loads((cd / (STEM + ".cfg.json")).read_text(encoding="utf-8"))
    model = json.loads((cd / (STEM + ".model.json")).read_text(encoding="utf-8"))
    EPS0, PINC, LAMMAX = float(cfg["EPS0"]), float(cfg["PINC"]), float(cfg["LAMMAX"])
    CN, CD = int(cfg["CTRL_NODE"]), int(cfg["CTRL_DOF"])

    fem = parse_stability(cd / (STEM + "_cable.txt"), CN, CD)
    dflt = parse_stability(cd / (DEFAULT_STEM + "_cable.txt"), CN, CD)
    if not fem["points"] or fem["basis"] != "completion_state":
        raise SystemExit("[FATAL] 主算例缺少 completion_state 平衡点")
    if dflt["has_section"]:
        raise SystemExit("[FATAL] 默认负例泄漏了 [STABILITY] 段")

    checks: list[list] = []
    strings: dict[str, str] = {}
    fem_path = [(q["lambda"], q["ctrl"]) for q in fem["points"]]

    # ---- FEM 段头 + 逐平衡点路径 ---------------------------------------
    pre = "fem"
    checks += [
        [f"{pre}/has_section", int(fem["has_section"])],
        [f"{pre}/nstep", fem["nstep"]],
        [f"{pre}/npoint", fem["npoint"]],
        [f"{pre}/converged", fem["converged"]],
        [f"{pre}/ctrlnode", fem["ctrlnode"]],
        [f"{pre}/ctrldof", fem["ctrldof"]],
        [f"{pre}/strain", fem["strain"]],
        [f"{pre}/psi", fem["psi"]],
        [f"{pre}/lambda_final", fem["lambda_final"]],
        [f"{pre}/ctrl0", fem["ctrl0"]],
        [f"{pre}/f0_norm", fem["f0_norm"]],
        [f"{pre}/finc_norm", fem["finc_norm"]],
        [f"{pre}/res0", fem["res0"]],
        [f"{pre}/u_ctrl", fem["u_ctrl"]],
        [f"{pre}/reslast", fem["reslast"]],
        [f"{pre}/has_newton", fem["has_newton"]],
        [f"{pre}/has_buckling_error", fem["has_buckling_error"]],
        [f"{pre}/n_points_parsed", len(fem["points"])],
    ]
    strings[f"{pre}/basis"] = fem["basis"]
    strings[f"{pre}/kind"] = fem["kind"]
    strings[f"{pre}/status"] = fem["status"]
    strings[f"{pre}/error"] = fem["error"] or "none"
    for q in fem["points"]:
        qp = f"{pre}/p{q['step']}"
        checks += [[f"{qp}/step", q["step"]],
                   [f"{qp}/lambda", q["lambda"]],
                   [f"{qp}/ctrl", q["ctrl"]],
                   [f"{qp}/res", q["res"]],
                   [f"{qp}/iter", q["iter"]],
                   [f"{qp}/ds", q["ds"]],
                   [f"{qp}/pivot", q["pivot"]],
                   [f"{qp}/converged", q["converged"]]]

    # ---- 默认零污染负例 -------------------------------------------------
    checks += [
        ["default/has_section", int(dflt["has_section"])],
        ["default/has_newton", dflt["has_newton"]],
        ["default/has_buckling_error", dflt["has_buckling_error"]],
        ["default/u_ctrl", dflt["u_ctrl"]],
    ]
    strings["default/error"] = dflt["error"] or "none"

    # ---- ANSYS 冻结曲线 + APDL + 日志 -----------------------------------
    curve_data: dict = {}
    apdl_data: dict = {}
    for job in ANSYS_JOBS:
        rows = parse_curve(cd / (job + "_curve.csv"))
        uz0, s1, s2 = split_curve(rows)
        ad = parse_apdl(cd / (job + ".inp"))
        er, wa = parse_log(cd / (job + "_ansys.out"))
        curve_data[job] = (rows, uz0, s1, s2)
        apdl_data[job] = (ad, er, wa)
        for i, r in enumerate(rows):
            checks += [[f"ansys/{job}/r{i}/iset", r[0]],
                       [f"ansys/{job}/r{i}/time", r[1]],
                       [f"ansys/{job}/r{i}/uz", r[2]],
                       [f"ansys/{job}/r{i}/pad", r[3]]]
        checks += [
            [f"ansys/{job}/n_rows", len(rows)],
            [f"ansys/{job}/n_step1", len(s1)],
            [f"ansys/{job}/n_step2", len(s2)],
            [f"ansys/{job}/uz_step1", uz0],
            [f"ansys/{job}/time_final", rows[-1][1]],
            [f"ansys/{job}/uz_final", rows[-1][2]],
            [f"apdl/{job}/f_orig", ad["f_orig"]],
            [f"apdl/{job}/f_end", ad["f_end"]],
            [f"apdl/{job}/time2", ad["time2"]],
            [f"apdl/{job}/n_inistate", ad["n_inistate"]],
            [f"apdl/{job}/n_solve", ad["n_solve"]],
            [f"apdl/{job}/n_loads", ad["n_loads"]],
            [f"apdl/{job}/arclen", int(ad["arclen"])],
            [f"log/{job}/errors", er],
            [f"log/{job}/warnings", wa],
        ]
        for i, e in enumerate(ad["eps"]):
            checks.append([f"apdl/{job}/eps{i}", e])
        strings[f"log/{job}/summary"] = "errors=%d warnings=%d" % (er, wa)
        strings[f"apdl/{job}/step2"] = (
            "ARCLEN" if ad["arclen"] else "load_control")

    # ---- 同源导出结论（strings）-----------------------------------------
    report = (cd / (STEM + "_compare_report.txt")).read_text(
        encoding="utf-8", errors="replace")
    m = re.search(r"^RESULT:\s*(.+)$", report, re.M)
    strings["report/result"] = m.group(1).strip() if m else "MISSING"
    strings["cfg/mode"] = str(cfg["COMPLETION_CFG"]["mode"])
    strings["cfg/method"] = str(cfg["COMPLETION_CFG"]["method"])

    # ---- 独立锚点：纯 Python 从 model.json + cfg.json + .inp 计算 --------
    coords = {n["id"]: (float(n["x"]), float(n["y"]), float(n["z"]))
              for n in model["nodes"]}
    kinds = {t["id"]: t["kind"] for t in model["etypes"]}
    cables = [e for e in model["elements"]
              if kinds.get(e["etype"]) == "CABLE"]
    l0map = {int(k): float(v) for k, v in cfg["cable_L0"].items()}
    chord_dev = []
    for e in cables:
        p1, p2 = coords[e["n1"]], coords[e["n2"]]
        L = math.sqrt(sum((b - a) ** 2 for a, b in zip(p1, p2)))
        chord_dev.append(abs(L - l0map[e["id"]] * (1.0 + EPS0)))
    eps_all = apdl_data["cs_load"][0]["eps"]
    eps_dev = [abs(x - EPS0) for x in eps_all]
    f_orig = apdl_data["cs_load"][0]["f_orig"]
    f_end_closed = f_orig - LAMMAX * PINC
    f_end_dev = abs(apdl_data["cs_load"][0]["f_end"] - f_end_closed)

    # ---- FEM ↔ ANSYS 交叉量（也从冻结证据重算，确定性的）---------------
    rows, uz0, s1, s2 = curve_data["cs_load"]
    worst, wlam, cnt = 0.0, 0.0, 0
    for iset, tnow, uz, _ in s2:
        lam = (tnow - 1.0) * LAMMAX
        if lam < -1e-8 or lam > fem_path[-1][0] + 1e-8:
            continue
        d = abs((uz - uz0) - interp_ctrl(fem_path, lam))
        cnt += 1
        if d > worst:
            worst, wlam = d, lam
    small = [p for p in s2 if 0.0 < (p[1] - 1.0) * LAMMAX <= 0.05]
    k_ansys = (small[0][1] - 1.0) * LAMMAX * PINC / abs(small[0][2] - uz0)
    k_fem = (fem_path[1][0] - fem_path[0][0]) * PINC / abs(
        fem_path[1][1] - fem_path[0][1])
    arc_rows = curve_data["cs_arc"][0]
    endpoint_gap = abs(arc_rows[-1][2] - rows[-1][2])

    criteria = {
        "eps0": EPS0, "pinc": PINC, "lammax": LAMMAX,
        "vmax": float(cfg["VMAX"]),
        "ctrl_node": CN, "ctrl_dof": CD,
        "f_orig": f_orig,
        "f_end_closed": f_end_closed,
        "f_end_apdl_dev": f_end_dev,
        "n_cables": len(cables),
        "n_cable_l0": len(l0map),
        "n_inistate": len(eps_all),
        "l0_min": min(l0map.values()), "l0_max": max(l0map.values()),
        "cable_chord_dev_max": max(chord_dev),
        "inistate_eps_dev_max": max(eps_dev),
        "fem_npoint": len(fem["points"]),
        "fem_lambda_final": fem["lambda_final"],
        "fem_ctrl0": fem["ctrl0"],
        "fem_res0": fem["res0"],
        "fem_f0_norm": fem["f0_norm"],
        "fem_reslast": fem["reslast"],
        "res0_over_f0": fem["res0"] / fem["f0_norm"],
        "res0_over_reslast": fem["res0"] / fem["reslast"],
        "ansys_uz_step1": uz0,
        "start_uz_gap": abs(fem["ctrl0"] - uz0),
        "ansys_n_step2": len(s2),
        "worst_duz": worst,
        "worst_duz_at_lambda": wlam,
        "n_compared_substeps": cnt,
        "initial_slope_ansys": k_ansys,
        "initial_slope_fem": k_fem,
        "initial_slope_rel": abs(k_ansys - k_fem) / k_fem,
        "endpoint_gap": endpoint_gap,
        "arc_n_rows": len(arc_rows),
    }
    return {
        "meta": {
            "source": "frozen CB3D21 two-step evidence (cable_solver.exe "
                      "BASIS=completion_state + ANSYS180 load-control / "
                      "ARCLEN) + independent same-source geometry oracle",
            "cases": [STEM, DEFAULT_STEM] + list(ANSYS_JOBS),
            "n_checks": len(checks),
            "n_strings": len(strings),
            "n_criteria": len(criteria),
            "note": "第十九次 ②：真实斜拉桥成桥态两段对标；路径单调硬化无"
                    "极限点，MODE=COMPLETION 起点不同源缺陷用 STATIC 绕开并"
                    "登记为独立后续项；② 完成后 B3 能力矩阵仍 P 不 Y",
        },
        "checks": checks,
        "strings": strings,
        "criteria": criteria,
    }


def check(doc: dict, cases_dir: Path | None = None) -> tuple[bool, str]:
    fresh = build(cases_dir)
    ok = (doc["checks"] == fresh["checks"]
          and doc["strings"] == fresh["strings"]
          and all(abs(doc["criteria"][k] - v)
                  <= 1e-10 * max(1.0, abs(v))
                  for k, v in fresh["criteria"].items()))
    return ok, (f"{fresh['meta']['n_checks']} checks + "
                f"{fresh['meta']['n_strings']} strings + "
                f"{fresh['meta']['n_criteria']} criteria")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="只校验现有 JSON 与冻结证据一致，不写文件")
    ap.add_argument("--cases", default=None, help="冻结 cases 目录")
    ap.add_argument("--reference", default=None, help="reference JSON 目录")
    a = ap.parse_args()
    cd = Path(a.cases) if a.cases else CASES
    ref = (Path(a.reference) / REF.name) if a.reference else REF
    doc = build(cd)
    if a.check:
        cur = json.loads(ref.read_text(encoding="utf-8"))
        ok, msg = check(cur, cd)
        print(f"[{'PASS' if ok else 'FAIL'}] {msg} consistent with frozen evidence")
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
