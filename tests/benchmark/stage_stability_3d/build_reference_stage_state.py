#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 26 项 CTest 黄金值生成器：从 CB3D23 冻结证据**派生**，不重跑 FEM / ANSYS。

    python build_reference_stage_state.py [--cases <dir>] [--out <reference dir>]

冻结证据由 BridgeStudio `cases/run_cases_stage_state.py --sync-to` 同步到
`cases/`。本脚本只读这些文件，输出
`reference/cb3d23_stage_state.json`；离线比对脚本
`compare_stage_state_3d.py` 会 import 本模块的解析函数做同样的派生，
再逐项与黄金值核对。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REFERENCE = HERE / "reference"
REF_NAME = "cb3d23_stage_state.json"
STEM = "cb3d23_stage_state"
BAD_STEM = "st_bad_stage"
ANSYS_CASES = ("st_load", "st_arc")
LAMMAX = 2.7
PINC = -20000.0
STAGE_ID = 2
CTRLTARGET = 0.005

# 冻结证据清单（与 BS 侧 manifest 一致；ANSYS 二进制临时文件不入库）。
FILES = (
    f"{STEM}.dat", f"{STEM}.cfg.json", f"{STEM}.model.json",
    f"{STEM}_stage_stability.txt", f"{STEM}_compare_report.txt",
    f"{STEM}_manifest.txt", f"{STEM}_run_log.txt",
    f"{BAD_STEM}.dat", f"{BAD_STEM}_stage_stability.txt",
    "st_load.inp", "st_load_curve.csv", "st_load_ansys.out",
    "st_arc.inp", "st_arc_curve.csv", "st_arc_ansys.out",
)


def normalized(txt: str) -> str:
    return txt.replace("\r\n", "\n").replace("\r", "\n")


def read_text(path: Path) -> str:
    return normalized(path.read_text(encoding="utf-8", errors="replace"))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_stages(txt: str) -> list[dict]:
    """`[STAGES]` 段的逐阶段摘要行。"""
    m = re.search(r"\[STAGES\]\n(.*?)(?=\n\[)", txt, re.S)
    if not m:
        return []
    out = []
    for ln in m.group(1).splitlines():
        kv = dict(re.findall(r"(\w+)=([^\s]+)", ln))
        if "ID" in kv:
            out.append({k: float(v) for k, v in kv.items()})
    return out


def parse_stability(txt: str) -> dict:
    """`[STABILITY]` 段：段头 / STATE / STAGE / 逐点平衡点。"""
    m = re.search(r"\[STABILITY\]\n(.*?)(?=\n\[END\])", txt, re.S)
    if not m:
        return {"present": False}
    sec = m.group(1)
    head = {}
    for k, v in re.findall(r"(\w+)=([^\s]+)", sec):
        head.setdefault(k, v)
    state, stage = {}, {}
    sm = re.search(r"^STATE (.+)$", sec, re.M)
    if sm:
        state = {k: float(v) for k, v in re.findall(r"(\w+)=([^\s]+)", sm.group(1))}
    sm = re.search(r"^STAGE=(.+)$", sec, re.M)
    if sm:
        toks = sm.group(1).split()
        stage = {"STAGE": float(toks[0])}
        stage.update({k: float(v)
                      for k, v in re.findall(r"(\w+)=([^\s]+)", sm.group(1))})
    pts = [[int(s), float(lam), float(ctrl), float(res)]
           for s, lam, ctrl, res in
           re.findall(r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+)",
                      sec, re.M)]
    return {"present": True, "head": head, "state": state,
            "stage": stage, "points": pts}


def fem_metrics(txt: str) -> dict:
    """FEM 侧派生指标。"""
    st = parse_stability(txt)
    stages = parse_stages(txt)
    if not st["present"]:
        return {"present": False, "stages": stages}
    pts = st["points"]
    f0n = st["state"].get("F0_NORM", 0.0)
    res0 = st["state"].get("RES0", float("inf"))
    return {
        "present": True,
        "head": st["head"], "state": st["state"], "stage": st["stage"],
        "npoint": len(pts),
        "nstep": int(float(st["head"].get("NSTEP", 0))),
        "res0_rel": res0 / f0n if f0n else float("inf"),
        "lambda_end": pts[-1][1] if pts else 0.0,
        "ctrl_end": pts[-1][2] if pts else 0.0,
        "lambda_max": max(p[1] for p in pts) if pts else 0.0,
        "n_descending": sum(1 for i in range(len(pts) - 1)
                            if pts[i + 1][1] < pts[i][1]),
        "stages": stages,
        "points": pts,
    }


def read_curve_csv(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    for ln in normalized(path.read_text(encoding="ascii", errors="replace")).splitlines():
        p = ln.split()
        if len(p) < 3:
            continue
        try:
            rows.append((float(p[0]), float(p[1]), float(p[2])))
        except ValueError:
            continue
    return rows


def ansys_step2_curve(rows: list[tuple[float, float, float]]) -> tuple[float, list[list[float]]]:
    """返回 (UY_step1, [[|ΔUY|, λ=TIME−1], ...])，ΔUY 向下为正。"""
    s1 = [r for r in rows if r[1] <= 1.0 + 1e-9]
    if not s1:
        raise ValueError("CSV 缺 step1 子步")
    u0 = s1[-1][2]
    out: list[list[float]] = []
    last = -1e300
    for r in rows:
        if r[1] <= 1.0 + 1e-9:
            continue
        du = -(r[2] - u0)
        if du > last + 1e-12:
            out.append([du, r[1] - 1.0])
            last = du
    return u0, out


def ansys_errors(txt: str) -> int:
    m = re.findall(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
    return int(m[-1]) if m else -1


def interp_lambda_by_ctrl(curve: list[list[float]], ctrl: float) -> float:
    if not curve:
        raise ValueError("曲线为空")
    if ctrl <= curve[0][0]:
        return curve[0][1]
    if ctrl >= curve[-1][0]:
        return curve[-1][1]
    lo, hi = 0, len(curve) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if curve[mid][0] <= ctrl:
            lo = mid
        else:
            hi = mid
    a, b = curve[lo], curve[hi]
    return a[1] + (b[1] - a[1]) * (ctrl - a[0]) / (b[0] - a[0])


def compare_curve(fem_pts: list[list[float]],
                  ansys_curve: list[list[float]]) -> dict:
    """FEM↔ANSYS 按 |ΔUY| 插值 λ 逐点对拍（只在共同区间比较）。"""
    fem = [[-p[2], p[1]] for p in fem_pts]     # 向下位移幅值
    lo = max(fem[0][0], ansys_curve[0][0])
    hi = min(fem[-1][0], ansys_curve[-1][0])
    if hi <= lo:
        return {"n": 0, "worst": float("inf"), "mean": float("inf"),
                "lo": lo, "hi": hi,
                "coverage": 0.0,
                "fem_span": fem[-1][0] - fem[0][0]}
    diffs = [abs(lam - interp_lambda_by_ctrl(ansys_curve, ctrl))
             for ctrl, lam in fem if lo <= ctrl <= hi]
    return {
        "n": len(diffs),
        "worst": max(diffs) if diffs else float("inf"),
        "mean": (sum(diffs) / len(diffs)) if diffs else float("inf"),
        "lo": lo, "hi": hi,
        "coverage": (hi - lo) / max(fem[-1][0] - fem[0][0], 1e-300),
        "fem_span": fem[-1][0] - fem[0][0],
    }


def apdl_anchor(dat_txt: str, inp_txt: str) -> dict:
    """独立锚点：FEM `.dat` 与 ANSYS `.inp` 的节点/单元/约束/阶段荷载交叉校核。

    逐行 `split(',')` + 严格 `float()`：早期用跨字段正则时 `[^,]` 会吞掉
    换行与 `!` 注释行，导致 `ValueError`（CB3D22 已踩过）。ANSYS 单元号
    是 `E,` 行的**创建顺序**，不是第一个节点号。
    """
    nodes_dat: dict[int, tuple[float, float, float]] = {}
    in_node_section = False
    for ln in dat_txt.splitlines():
        s = ln.strip()
        if s.startswith("*"):
            in_node_section = (s == "*NODE")
            continue
        if not in_node_section:
            continue
        m = re.match(r"^(\d+)\s+(-?[\d.eE+]+)\s+(-?[\d.eE+]+)\s+(-?[\d.eE+]+)\s*$", s)
        if m:
            nodes_dat[int(m.group(1))] = tuple(float(m.group(i)) for i in (2, 3, 4))

    nodes_inp: dict[int, tuple[float, float, float]] = {}
    elems_inp: list[tuple[int, int]] = []
    loads_inp: dict[int, list[float]] = {}
    fixes: list[str] = []
    n_solve = 0
    for ln in inp_txt.splitlines():
        q = [s.strip() for s in ln.split(",")]
        if not q:
            continue
        tag = q[0].upper()
        if tag == "N" and len(q) >= 5 and q[1].isdigit():
            try:
                nodes_inp[int(q[1])] = (float(q[2]), float(q[3]), float(q[4]))
            except ValueError:
                pass
        elif tag == "E" and len(q) >= 3 and q[1].isdigit() and q[2].isdigit():
            elems_inp.append((int(q[1]), int(q[2])))
        elif tag == "F" and len(q) >= 4 and q[1].isdigit() and q[2].upper() == "FY":
            try:
                loads_inp.setdefault(int(q[1]), []).append(float(q[3]))
            except ValueError:
                pass
        elif tag == "D" and len(q) >= 3:
            fixes.append("%s,%s" % (q[1], q[2].upper()))
        elif tag == "SOLVE":
            n_solve += 1

    coord_dev = 0.0
    for nid, xyz in nodes_dat.items():
        if nid not in nodes_inp:
            coord_dev = float("inf")
            break
        coord_dev = max(coord_dev, max(abs(a - b) for a, b in zip(xyz, nodes_inp[nid])))
    return {
        "n_nodes_dat": len(nodes_dat), "n_nodes_inp": len(nodes_inp),
        "coord_max_dev": coord_dev,
        "n_elems_inp": len(elems_inp), "elems_inp": [list(e) for e in elems_inp],
        "loads_inp": {str(k): v for k, v in sorted(loads_inp.items())},
        "loads_step1": {str(k): v[0] for k, v in sorted(loads_inp.items()) if v},
        "loads_step2": {str(k): v[-1] for k, v in sorted(loads_inp.items()) if v},
        "n_fixes": len(fixes), "fixes": sorted(fixes),
        "n_solve": n_solve,
        "nodes_dat": {str(k): list(v) for k, v in sorted(nodes_dat.items())},
    }


def build(cases: Path) -> dict:
    files = {}
    for name in FILES:
        p = cases / name
        if not p.is_file():
            raise FileNotFoundError(f"缺冻结证据: {p}")
        files[name] = {"sha256": sha256_file(p), "bytes": p.stat().st_size}

    fem_txt = read_text(cases / f"{STEM}_stage_stability.txt")
    fem = fem_metrics(fem_txt)
    bad_txt = read_text(cases / f"{BAD_STEM}_stage_stability.txt")
    dat_txt = read_text(cases / f"{STEM}.dat")
    cfg = json.loads(read_text(cases / f"{STEM}.cfg.json"))
    report_txt = read_text(cases / f"{STEM}_compare_report.txt")

    ansys, comparisons, anchors = {}, {}, {}
    for job in ANSYS_CASES:
        inp_txt = read_text(cases / f"{job}.inp")
        out_txt = read_text(cases / f"{job}_ansys.out")
        rows = read_curve_csv(cases / f"{job}_curve.csv")
        u0, curve = ansys_step2_curve(rows)
        ansys[job] = {
            "u0": u0, "n_rows": len(rows), "n_step2": len(curve),
            "lam_end": curve[-1][1] if curve else 0.0,
            "errors": ansys_errors(out_txt),
            "n_solve_inp": inp_txt.count("\nSOLVE\n") + (1 if inp_txt.startswith("SOLVE") else 0),
            "has_arclen": "ARCLEN,ON" in inp_txt,
            "curve": curve,
        }
        comparisons[job] = compare_curve(fem["points"], curve)
        anchors[job] = apdl_anchor(dat_txt, inp_txt)

    return {
        "case": "CB3D23",
        "title": "③b 施工阶段中间态起点的全过程稳定 (BASIS=stage_state)",
        "basis": "stage_state",
        "stage_id": STAGE_ID,
        "lammax": LAMMAX,
        "pinc": PINC,
        "ctrltarget": CTRLTARGET,
        "n_files": len(FILES),
        "files": files,
        "fem": fem,
        "ansys": ansys,
        "comparisons": comparisons,
        "anchors": anchors,
        "negative": {
            "error": "stage_not_found" if "ERROR=stage_not_found" in bad_txt else "",
            "basis": "stage_state" if "BASIS=stage_state" in bad_txt else "",
            "has_points": bool(re.search(r"^I STEP=", bad_txt, re.M)),
        },
        "config": {
            "stage_id": cfg.get("stage_id"),
            "pinc": cfg.get("pinc"),
            "stage_active": cfg.get("stage_active"),
            "stage_loads": cfg.get("stage_loads"),
            "stage_constant_loads": cfg.get("stage_constant_loads"),
            "stage_active_elements": cfg.get("stage_active_elements"),
            "fixed_x": cfg.get("fixed_x"),
            "ansys_cases": cfg.get("ansys_cases"),
        },
        "report": {
            "result": "ALL CHECKS PASSED" if "RESULT: ALL CHECKS PASSED" in report_txt else "",
            "checks": int(re.search(r"(\d+) checks, (\d+) failed", report_txt).group(1)),
            "failed": int(re.search(r"(\d+) checks, (\d+) failed", report_txt).group(2)),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default=str(CASES))
    ap.add_argument("--out", default=str(REFERENCE))
    args = ap.parse_args()
    ref = build(Path(args.cases))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    dest = out / REF_NAME
    dest.write_text(json.dumps(ref, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    print(f"reference -> {dest}")
    print("FEM  NPOINT=%d lambda_end=%.9g ctrl0=%.9g res0_rel=%.3g"
          % (ref["fem"]["npoint"], ref["fem"]["lambda_end"],
             ref["fem"]["state"]["CTRL0"], ref["fem"]["res0_rel"]))
    for job in ANSYS_CASES:
        a, c = ref["ansys"][job], ref["comparisons"][job]
        print("%s: step2=%d lam_end=%.6g errors=%d worst=%.4g mean=%.4g"
              % (job, a["n_step2"], a["lam_end"], a["errors"],
                 c["worst"], c["mean"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
