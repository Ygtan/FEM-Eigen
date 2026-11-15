#!/usr/bin/env python3
"""第 25 项 CTest 黄金值生成器（第十九次 ③a / CB3D22）。

从 BS 侧 `cases/CB3D22_limit_point_vs_ansys/` 同步来的冻结证据里，重算
所有指标并写 `reference/cb3d22_limit_point.json`。**不重跑 FEM / ANSYS**：
cases 是唯一冻结源，本脚本只做派生量计算与哈希。

用法：
    python build_reference_limit_point.py <cases_dir> [reference_dir]
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF_NAME = "cb3d22_limit_point.json"
STEM = "cb3d22_limit_point"
EPS0 = 8.0e-4
PINC = 1.0e6
CTRL = (2009, 1)
ANSYS_JOBS = ("lp_load", "lp_arc")


def normalized(txt: str) -> str:
    return re.sub(r"elapsed_ms=[0-9.eE+-]+", "elapsed_ms=*", txt)


def parse_stability(txt: str) -> dict:
    m = re.search(r"\[STABILITY\]\n(.*?)(?=\n\[END\])", txt, re.S)
    if not m:
        return {"present": False}
    sec = m.group(1)
    head = {}
    for k, v in re.findall(r"(\w+)=([^\s]+)", sec):
        head.setdefault(k, v)
    state = {}
    sm = re.search(r"^STATE (.+)$", sec, re.M)
    if sm:
        state = {k: float(v) for k, v in re.findall(r"(\w+)=([^\s]+)", sm.group(1))}
    pts = [[int(s), float(lam), float(ctrl), float(res)]
           for s, lam, ctrl, res in
           re.findall(r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+)",
                      sec, re.M)]
    return {"present": True, "head": head, "state": state, "points": pts}


def fem_metrics(txt: str) -> dict:
    st = parse_stability(txt)
    usec = re.search(r"^\[U\]\n(.*?)(?=^\[[A-Z_]+\]|\Z)", txt, re.S | re.M)
    uline = re.search(r"^2009 1 (\S+)", usec.group(1), re.M) if usec else None
    reslast = re.search(r"RESLAST=(\S+)", txt)
    return {
        "stability": st,
        "u_ctrl": float(uline.group(1)) if uline else None,
        "reslast": float(reslast.group(1)) if reslast else None,
        "sections": re.findall(r"^\[([A-Z_]+)\]", txt, re.M),
        "buckling_error_kept": "ERROR=nonlinear_state_not_supported" in txt,
    }


def read_curve(txt: str) -> list[list[float]]:
    out = []
    for ln in txt.splitlines():
        q = ln.split()
        if len(q) >= 3:
            try:
                out.append([float(q[0]), float(q[1]), float(q[2])])
            except ValueError:
                pass
    return out


def step2_curve(rows: list[list[float]]) -> tuple[float, list[list[float]]]:
    s1 = [r for r in rows if r[1] <= 1.0 + 1e-9]
    u0 = s1[-1][2] if s1 else 0.0
    out: list[list[float]] = []
    last = -1e300
    for r in rows:
        if r[1] <= 1.0 + 1e-9:
            continue
        du = r[2] - u0
        if du > last + 1e-12:
            out.append([du, r[1] - 1.0])
            last = du
    return u0, out


def limit_of(curve: list[list[float]]) -> tuple[int, float, float]:
    i = max(range(len(curve)), key=lambda k: curve[k][1])
    return i, curve[i][1], curve[i][0]


def apdl_anchor(dat: str, inp: str) -> dict:
    """独立锚点：从 ANSYS `N,`/`E,`/`INISTATE` 行与 FEM `.dat` 的 L0 交叉校核。

    解析必须逐行且严格数值化：早期用 `^N,(\\d+),([^,]+),...\\s*$` 的正则，
    `[^,]` 会跨行吞掉 `\n` 和后续 `!` 注释行，把 orientation node 的说明文字
    当成 Z 坐标（`ValueError: could not convert string to float: '25\\n! ---- ...'`）。
    ANSYS 单元号是 `E,` 行的**创建顺序**（文件头已注明 "elem no. = creation
    order"），不是该行的第一个节点号；BEAM4 的第三字段是 K 方位节点，LINK180
    索单元只有两个节点，取前两个即可。
    """
    nodes: dict[int, tuple[float, float, float]] = {}
    elems: dict[int, tuple[int, int]] = {}
    eid = 0
    for ln in inp.splitlines():
        q = [s.strip() for s in ln.split(",")]
        if not q or not q[0]:
            continue
        tag = q[0].upper()
        if tag == "N" and len(q) >= 5 and q[1].isdigit():
            try:
                nodes[int(q[1])] = (float(q[2]), float(q[3]), float(q[4]))
            except ValueError:
                pass
        elif tag == "E" and len(q) >= 3 and q[1].isdigit() and q[2].isdigit():
            eid += 1
            elems[eid] = (int(q[1]), int(q[2]))
    eps_lines = re.findall(r"INISTATE,DEFINE,(\d+),[^,\n]*,EPEL,X,([\d.eE+-]+)", inp)
    l0_by_pair = {}
    for m in re.finditer(r"^\*EN,(90\d\d),(\d+),(\d+) L0=([\d.eE+-]+)", dat, re.M):
        l0_by_pair[(int(m.group(2)), int(m.group(3)))] = float(m.group(4))
    eps_vals, chord_errs, count = [], [], 0
    for eid, eps in eps_lines:
        e = elems.get(int(eid))
        if not e or e[0] not in nodes or e[1] not in nodes:
            continue
        p1, p2 = nodes[e[0]], nodes[e[1]]
        L = sum((b - a) ** 2 for a, b in zip(p1, p2)) ** 0.5
        eps_vals.append(float(eps))
        key = e if e in l0_by_pair else (e[1], e[0])
        if key in l0_by_pair:
            chord_errs.append(abs(l0_by_pair[key] * (1.0 + float(eps)) - L))
            count += 1
    return {
        "n_nodes": len(nodes), "n_elems": len(elems),
        "n_inistate": len(eps_lines), "n_l0_matched": count,
        "eps_max_dev": max((abs(v - EPS0) for v in eps_vals), default=1e300),
        "chord_max_dev": max(chord_errs, default=1e300),
    }


def ansys_errors(out: str) -> int:
    m = re.findall(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", out)
    return int(m[-1]) if m else -1


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    cases = Path(sys.argv[1]).resolve()
    refdir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else HERE / "reference"
    refdir.mkdir(parents=True, exist_ok=True)

    ref: dict = {
        "case": "CB3D22",
        "title": "③a 真实斜拉桥越过极值点 / 下降段 —— FEM vs ANSYS 两段 ARCLEN",
        "bridge_source": "BridgeStudio/cases/CB3D22_limit_point_vs_ansys",
        "tower_scale": 0.012,
        "eps0": EPS0, "pinc": PINC,
        "ctrl": list(CTRL),
        "files": {}, "fem": {}, "ansys": {}, "anchors": {},
    }
    all_files = sorted(p.name for p in cases.iterdir() if p.is_file())
    for name in all_files:
        raw = (cases / name).read_bytes()
        ref["files"][name] = {
            "sha256": hashlib.sha256(raw).hexdigest(),
            "n_bytes": len(raw),
        }
    fem_txt = (cases / f"{STEM}_cable.txt").read_text(encoding="utf-8", errors="replace")
    dflt_txt = (cases / "lp_default_cable.txt").read_text(encoding="utf-8", errors="replace")
    dat = (cases / f"{STEM}.dat").read_text(encoding="utf-8")
    ref["fem"] = {"main": fem_metrics(fem_txt),
                  "default": {"has_stability": "[STABILITY]" in dflt_txt,
                              "n_lines": len(dflt_txt.splitlines())}}
    for job in ANSYS_JOBS:
        inp = (cases / f"{job}.inp").read_text(encoding="ascii", errors="replace")
        out = (cases / f"{job}_ansys.out").read_text(encoding="ascii", errors="replace")
        curve = read_curve((cases / f"{job}_curve.csv").read_text(errors="replace"))
        u0, s2 = step2_curve(curve)
        i, lam, ctrl = limit_of(s2) if s2 else (0, 0.0, 0.0)
        ref["ansys"][job] = {
            "n_solve": inp.count("SOLVE"),
            "n_fx": len(re.findall(r"^F,2009,FX,", inp, re.M)),
            "n_inistate": len(re.findall(r"^INISTATE,DEFINE,", inp, re.M)),
            "arclen": "ARCLEN,ON" in inp,
            "error_messages": ansys_errors(out),
            "u0": u0, "n_step2_points": len(s2),
            "limit_index": i, "limit_lambda": lam, "limit_ctrl": ctrl,
            "curve": s2,
        }
        ref["anchors"][job] = apdl_anchor(dat, inp)
    fpts = ref["fem"]["main"]["stability"]["points"]
    curve = [(p[2], p[1]) for p in fpts]
    fi, flam, fctrl = limit_of([[c, l] for c, l in curve])
    ref["fem"]["limit_index"] = fi
    ref["fem"]["limit_lambda"] = flam
    ref["fem"]["limit_ctrl"] = fctrl
    ref["fem"]["n_descending"] = sum(
        1 for k in range(fi + 1, len(curve)) if curve[k][1] < curve[k - 1][1] - 1e-9)
    ref["fem"]["lambda_end"] = curve[-1][1] if curve else 0.0
    target = refdir / REF_NAME
    target.write_text(json.dumps(ref, ensure_ascii=False, indent=2, sort_keys=True),
                      encoding="utf-8")
    print("reference ->", target)
    print(f"FEM  NPOINT={len(fpts)} lambda_max={flam:.7g} @ dUX={fctrl:.7g} "
          f"descending={ref['fem']['n_descending']}")
    for job in ANSYS_JOBS:
        a = ref["ansys"][job]
        print(f"{job}: step2={a['n_step2_points']} "
              f"lambda_max={a['limit_lambda']:.7g} @ dUX={a['limit_ctrl']:.7g} "
              f"errors={a['error_messages']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
