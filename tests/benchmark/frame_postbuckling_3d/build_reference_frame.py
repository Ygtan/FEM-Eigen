#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 21 项 CTest 黄金值生成器：CB3D19 CR 梁框架后屈曲。

从冻结的 cases/（ANSYS 逐子步 CSV + frame_probe.exe 输出 + 报告 / 日志）
生成 reference/cb3d19_frame_postbuckling.json，含：

  * checks  : ANSYS 逐子步 V/LAM/TIME、FEM 逐平衡点 6 字段、FEM S 行汇总
  * strings : 负例错误码、报告结论、ANSYS 日志摘要
  * criteria: 独立锚点 k0、极限点、下降段、曲线一致性等派生判据

本文件**自带第三套独立线性刚架刚度实现**（纯标准库、不跨仓库 import、
不调用 FEM 内核），与 BridgeStudio/ansys_bridge_frame.py 的实现、
frame_probe.exe 的 k0_fem 三方互证。

用法：
    python build_reference_frame.py            # 写 reference JSON
    python build_reference_frame.py --check    # 只校验已有 JSON 自洽
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF_FILE = HERE / "reference" / "cb3d19_frame_postbuckling.json"

CASE_ID = "cb3d19_frame_postbuckling"
PARAMS_FILE = CASE_ID + ".frame"
FEM_CSV = CASE_ID + "_fem.csv"
MODEL_JSON = CASE_ID + ".model.json"
REPORT_TXT = CASE_ID + "_compare_report.txt"
ANSYS_JOBS = ("fr_fine", "fr_coarse")
FEM_CASES = ("fine", "coarse", "psi_half")
NEG_CASE = "bad_noload"
P_KEYS = ("LAMBDA", "V", "RES", "ITER", "DS", "PIVOT")


def cfg_of(cases_dir: Path | None = None) -> dict:
    cs = Path(cases_dir) if cases_dir is not None else CASES
    return json.loads((cs / MODEL_JSON).read_text(encoding="utf-8"))["cfg"]


# ------------------------------------------------ 独立线性刚架刚度锚点
def frame_nodes(cfg: dict) -> list[tuple[float, float]]:
    n = int(cfg["NELEM"])
    B, H = cfg["B"], cfg["H"]
    pts = [(-B + (i / n) * B, (i / n) * H) for i in range(n + 1)]
    pts += [((i / n) * B, H - (i / n) * H) for i in range(1, n + 1)]
    return pts


def _local_stiffness(L: float, E: float, A: float, I: float) -> list[list[float]]:
    ea, ei = E * A / L, E * I
    return [
        [ea, 0.0, 0.0, -ea, 0.0, 0.0],
        [0.0, 12.0 * ei / L**3, 6.0 * ei / L**2,
         0.0, -12.0 * ei / L**3, 6.0 * ei / L**2],
        [0.0, 6.0 * ei / L**2, 4.0 * ei / L,
         0.0, -6.0 * ei / L**2, 2.0 * ei / L],
        [-ea, 0.0, 0.0, ea, 0.0, 0.0],
        [0.0, -12.0 * ei / L**3, -6.0 * ei / L**2,
         0.0, 12.0 * ei / L**3, -6.0 * ei / L**2],
        [0.0, 6.0 * ei / L**2, 2.0 * ei / L,
         0.0, -6.0 * ei / L**2, 4.0 * ei / L],
    ]


def _transform(c: float, s: float) -> list[list[float]]:
    return [
        [c, s, 0.0, 0.0, 0.0, 0.0],
        [-s, c, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, c, s, 0.0],
        [0.0, 0.0, 0.0, -s, c, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ]


def _spd_solve(K: list[list[float]], F: list[float]) -> list[float]:
    n = len(F)
    L = [[0.0] * n for _ in range(n)]
    D = [0.0] * n
    for i in range(n):
        for j in range(i + 1):
            s = K[i][j]
            for k in range(j):
                s -= L[i][k] * L[j][k] * D[k]
            if i == j:
                if abs(s) < 1.0e-300:
                    raise ValueError("singular linear anchor stiffness")
                D[i] = s
                L[i][i] = 1.0
            else:
                L[i][j] = s / D[j]
    z = [0.0] * n
    for i in range(n):
        z[i] = F[i] - sum(L[i][k] * z[k] for k in range(i))
    y = [z[i] / D[i] for i in range(n)]
    u = [0.0] * n
    for i in range(n - 1, -1, -1):
        u[i] = y[i] - sum(L[k][i] * u[k] for k in range(i + 1, n))
    return u


def linear_frame_stiffness(cfg: dict) -> float:
    """独立锚点：顶点初始切线刚度 k0（N/m），不依赖任何求解器。"""
    pts = frame_nodes(cfg)
    nn = len(pts)
    nd = 3 * nn
    K = [[0.0] * nd for _ in range(nd)]
    E, A, IZ = cfg["E"], cfg["A"], cfg["IZ"]
    for i in range(nn - 1):
        dx = pts[i + 1][0] - pts[i][0]
        dy = pts[i + 1][1] - pts[i][1]
        L = math.hypot(dx, dy)
        c, s = dx / L, dy / L
        kl = _local_stiffness(L, E, A, IZ)
        T = _transform(c, s)
        tmp = [[sum(kl[r][k] * T[k][q] for k in range(6)) for q in range(6)]
               for r in range(6)]
        kg = [[sum(T[k][r] * tmp[k][q] for k in range(6)) for q in range(6)]
              for r in range(6)]
        dofs = [3 * i, 3 * i + 1, 3 * i + 2,
                3 * (i + 1), 3 * (i + 1) + 1, 3 * (i + 1) + 2]
        for r in range(6):
            for q in range(6):
                K[dofs[r]][dofs[q]] += kg[r][q]
    apex = int(cfg["NELEM"])
    F = [0.0] * nd
    F[3 * apex + 1] = -1.0
    fixed = {0, 1, 3 * (nn - 1), 3 * (nn - 1) + 1}
    free = [i for i in range(nd) if i not in fixed]
    Kf = [[K[r][q] for q in free] for r in free]
    uf = _spd_solve(Kf, [F[r] for r in free])
    v = -uf[free.index(3 * apex + 1)]
    if v <= 0.0:
        raise ValueError("non-positive apex flexibility")
    return 1.0 / v


# ------------------------------------------------ 冻结文件解析
def parse_fem(path: Path):
    pts: dict[str, list[tuple]] = {}
    summ: dict[str, dict[str, float]] = {}
    errors: dict[str, str] = {}
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        p = line.split(",")
        if not p:
            continue
        try:
            if p[0] == "P" and len(p) >= 9:
                pts.setdefault(p[1], []).append((
                    int(p[2]), float(p[3]), float(p[4]), float(p[5]),
                    float(p[6]), float(p[7]), float(p[8])))
            elif p[0] == "S" and len(p) >= 4:
                summ.setdefault(p[1], {})[p[2]] = float(p[3])
            elif p[0] == "E" and len(p) >= 3:
                errors[p[1]] = p[2]
        except (ValueError, IndexError):
            continue
    return pts, summ, errors


def parse_ansys(path: Path, cfg: dict) -> list[tuple[float, float, float, float]]:
    out = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        p = line.split()
        if len(p) < 5:
            continue
        try:
            iset, uy, r1, r3, tnow = (float(x) for x in p[:5])
        except ValueError:
            continue
        out.append((-uy, (r1 + r3) / cfg["PREF"],
                    tnow / cfg["TIME_END"], iset))
    return out


def interp_lambda(curve: list[tuple[float, float]], vq: float) -> float:
    """三点二次 Lagrange 插值（与 L 阶段同口径，避免线性插值误差）。"""
    if not curve:
        return 0.0
    if len(curve) < 3:
        for (a, ya), (b, yb) in zip(curve, curve[1:]):
            if a <= vq <= b and b > a:
                return ya + (yb - ya) * (vq - a) / (b - a)
        return 0.0
    i = min(range(len(curve)), key=lambda k: abs(curve[k][0] - vq))
    j = max(0, min(i - 1, len(curve) - 3))
    xs = [curve[k][0] for k in (j, j + 1, j + 2)]
    ys = [curve[k][1] for k in (j, j + 1, j + 2)]
    if any(abs(xs[a] - xs[b]) < 1e-14 for a in range(3) for b in range(a + 1, 3)):
        return 0.0
    s = 0.0
    for a in range(3):
        t = ys[a]
        for b in range(3):
            if a != b:
                t *= (vq - xs[b]) / (xs[a] - xs[b])
        s += t
    return s


def first_limit_point(curve: list[tuple[float, float]]):
    """第一个 lambda 局部极大（不是全局最大：后段会强烈再硬化）。"""
    if len(curve) < 3:
        return 0.0, 0.0, 0.0, 0.0
    i = -1
    for k in range(1, len(curve) - 1):
        if curve[k][1] > curve[k - 1][1] and curve[k][1] >= curve[k + 1][1]:
            i = k
            break
    if i < 0:
        return 0.0, 0.0, 0.0, 0.0
    x0, y0 = curve[i - 1]
    x1, y1 = curve[i]
    x2, y2 = curve[i + 1]
    if abs(x1 - x0) < 1e-14 or abs(x2 - x1) < 1e-14:
        return x1, y1, x1, y1
    s1 = (y1 - y0) / (x1 - x0)
    s2 = (y2 - y1) / (x2 - x1)
    d2 = (s2 - s1) / (x2 - x0)
    if abs(d2) < 1e-14:
        return x1, y1, x1, y1
    vs = 0.5 * (x0 + x1) - s1 / (2.0 * d2)
    ls = y0 + s1 * (vs - x0) + d2 * (vs - x0) * (vs - x1)
    return vs, ls, x1, y1


def curve_of(pts: dict[str, list[tuple]], cid: str) -> list[tuple[float, float]]:
    return [(x[2], x[1]) for x in pts.get(cid, [])]


# ------------------------------------------------ 黄金值构建
def build(cases_dir: Path | None = None) -> dict:
    cs = Path(cases_dir) if cases_dir is not None else CASES
    cfg = cfg_of(cs)
    k0 = linear_frame_stiffness(cfg)
    pts, summ, errors = parse_fem(cs / FEM_CSV)
    ansys = {j: parse_ansys(cs / (j + "_curve.csv"), cfg)
             for j in ANSYS_JOBS}

    checks: list[dict] = []
    # ANSYS 逐子步
    for job, rows in ansys.items():
        for i, (v, lr, lt, iset) in enumerate(rows):
            for key, val in (("V", v), ("LAM", lr), ("TIME", lt)):
                checks.append({"job": job, "row": i, "key": key,
                               "ref": val, "source": "ansys_csv"})
    # FEM 逐平衡点 + S 行汇总
    for cid in FEM_CASES:
        for row in pts.get(cid, []):
            for key, val in zip(P_KEYS, row[1:]):
                checks.append({"case": cid, "row": int(row[0]), "key": key,
                               "ref": float(val), "source": "fem_csv"})
        for key, val in sorted(summ.get(cid, {}).items()):
            checks.append({"case": cid, "key": "S_" + key,
                           "ref": float(val), "source": "fem_summary"})
    for key, val in sorted(summ.get(NEG_CASE, {}).items()):
        checks.append({"case": NEG_CASE, "key": "S_" + key,
                       "ref": float(val), "source": "fem_summary"})

    report = (cs / REPORT_TXT).read_text(encoding="utf-8")
    strings = {
        "bad_noload/error": errors.get(NEG_CASE, ""),
        "report/result": "ALL CHECKS PASSED" if "ALL CHECKS PASSED" in report else "FAILED",
    }
    for job in ANSYS_JOBS:
        txt = (cs / (job + "_ansys.out")).read_text(
            encoding="ascii", errors="replace")
        import re as _re
        m = _re.findall(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
        w = _re.findall(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
        strings[job + "/log"] = "errors=%s warnings=%s" % (
            m[-1] if m else "-1", w[-1] if w else "-1")

    criteria: dict[str, float] = {
        "k0_analytic": k0,
        "k0_over_PREF": k0 / cfg["PREF"],
        "k0_rel_fem": max(
            abs(summ.get(cid, {}).get("k0_fem", 0.0) - k0) / k0
            for cid in FEM_CASES),
    }
    for cid in FEM_CASES:
        s = summ.get(cid, {})
        criteria["npoints_" + cid] = s.get("npoints", 0.0)
        criteria["lambda_max_fit_" + cid] = s.get("lambda_max_fit", 0.0)
        criteria["v_at_lambda_max_fit_" + cid] = s.get("v_at_lambda_max_fit", 0.0)
        criteria["lambda_min_" + cid] = s.get("lambda_min", 0.0)
        criteria["v_final_" + cid] = s.get("v_final", 0.0)
        criteria["max_res_rel_" + cid] = s.get("max_residual_rel", 9.9)
    for job, rows in ansys.items():
        cur = [(v, lr) for v, lr, lt, iset in rows]
        vf, lf, _, _ = first_limit_point(cur)
        criteria["lambda_cr_ansys_" + job] = lf
        criteria["v_cr_ansys_" + job] = vf
        criteria["lambda_min_ansys_" + job] = min(x[1] for x in cur)
        criteria["v_final_ansys_" + job] = cur[-1][0]
        criteria["worst_lambda_time_" + job] = max(
            abs(lr - lt) for v, lr, lt, iset in rows)
        small = [x for x in rows if 0.0 < x[0] <= 0.002 * cfg["H"]]
        if small:
            sl = small[-1][1] / small[-1][0]
            criteria["slope_ansys_" + job] = sl
            criteria["slope_rel_" + job] = abs(sl - k0 / cfg["PREF"]) / (k0 / cfg["PREF"])
        fem_curve = curve_of(pts, "fine" if job == "fr_fine" else "coarse")
        worst = 0.0
        for v, lr, lt, iset in rows:
            if 0.0 <= v <= cfg["VMAX"]:
                worst = max(worst, abs(lr - interp_lambda(fem_curve, v)))
        criteria["worst_curve_" + job] = worst

    doc = {
        "id": CASE_ID,
        "title": "CB3D19 CR beam frame post-buckling vs ANSYS BEAM188 ARCLEN",
        "description": (
            "Williams toggle 型刚接框架（两根 24 单元 CR 梁，端部铰支、"
            "顶点刚接），球面弧长法完整路径 vs ANSYS BEAM188 + NLGEOM + "
            "ARCLEN。独立锚点为纯 Python 2D 线性刚架刚度矩阵 k0。"),
        "config": cfg,
        "params_file": PARAMS_FILE,
        "fem_csv": FEM_CSV,
        "model_json": MODEL_JSON,
        "fem_exe": "frame_probe.exe",
        "fem_exe_args": [PARAMS_FILE, FEM_CSV],
        "ansys_jobs": list(ANSYS_JOBS),
        "fem_cases": list(FEM_CASES),
        "negative_case": NEG_CASE,
        "k0_analytic": k0,
        "checks": checks,
        "strings": strings,
        "criteria": criteria,
        "tolerances": {
            "k0_rel": 1.0e-8,
            "slope_rel": 1.0e-2,
            "lambda_time": 1.0e-4,
            "curve": 1.0e-3,
            "extremum": 2.0e-3,
            "res_rel": 1.0e-8,
            "reproduce": 1.0e-12,
        },
        "meta": {
            "source": "BridgeStudio/cases/CB3D19_frame_postbuckling_arc_length",
            "n_checks": len(checks),
            "n_strings": len(strings),
            "n_criteria": len(criteria),
            "cases": list(ANSYS_JOBS) + list(FEM_CASES) + [NEG_CASE],
            "note": ("黄金值完全来自冻结 cases/；--live-ansys 只写临时目录，"
                     "不覆写任何冻结 CSV。"),
        },
    }
    return doc


def check(doc: dict, cases_dir: Path | None = None) -> tuple[bool, str]:
    fresh = build(cases_dir)
    if len(fresh["checks"]) != len(doc["checks"]):
        return False, "checks 数量漂移 %d vs %d" % (
            len(fresh["checks"]), len(doc["checks"]))
    for i, (a, b) in enumerate(zip(fresh["checks"], doc["checks"])):
        if a.keys() != b.keys():
            return False, "checks[%d] 键不一致" % i
        for k in a:
            av, bv = a[k], b[k]
            if isinstance(av, float):
                if abs(av - bv) > 1e-12 + 1e-12 * abs(av):
                    return False, "checks[%d].%s 漂移 %.3e" % (i, k, av - bv)
            elif av != bv:
                return False, "checks[%d].%s 不一致" % (i, k)
    if fresh["strings"] != doc["strings"]:
        return False, "strings 漂移"
    if set(fresh["criteria"]) != set(doc["criteria"]):
        return False, "criteria 键漂移"
    for k, v in fresh["criteria"].items():
        if abs(v - doc["criteria"][k]) > 1e-12 + 1e-12 * abs(v):
            return False, "criteria.%s 漂移 %.3e" % (k, v - doc["criteria"][k])
    if abs(fresh["k0_analytic"] - doc["k0_analytic"]) > 1e-6:
        return False, "k0_analytic 漂移"
    return True, ("%d checks + %d strings + %d criteria consistent"
                  % (len(doc["checks"]), len(doc["strings"]),
                     len(doc["criteria"])))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="only verify the existing reference JSON")
    a = ap.parse_args()
    if a.check:
        doc = json.loads(REF_FILE.read_text(encoding="utf-8"))
        good, msg = check(doc)
        print(("[PASS] " if good else "[FAIL] ") + msg)
        return 0 if good else 1
    doc = build()
    good, msg = check(doc)
    if not good:
        print("[FAIL] " + msg)
        return 1
    REF_FILE.parent.mkdir(parents=True, exist_ok=True)
    REF_FILE.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    print("[GEN ] %s -> %d checks, %d strings, %d criteria, %d cases" % (
        REF_FILE.name, len(doc["checks"]), len(doc["strings"]),
        len(doc["criteria"]), len(doc["meta"]["cases"])))
    print("[PASS] " + msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
