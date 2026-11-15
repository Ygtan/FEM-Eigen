#!/usr/bin/env python3
"""第 25 项 CTest：真实斜拉桥越过极值点 / 下降段 —— FEM vs ANSYS（离线）。

③a 的实质技术收口之一：② 阶段只做到「稳定路径从成桥态起步」，尚未证明
内核能随到**越过极限点后的下降段**。本项把 BS 侧
`cases/CB3D22_limit_point_vs_ansys/` 的冻结证据固化为回归，判据来自
`cases/run_cases_limit_point.py --check-only`，但完全离线、不重跑 ANSYS。

构型是 CB3D21 同源双塔三跨斜拉桥（50 节点 / 63 单元 / 16 索），只把桥塔
截面按 `TOWER_SCALE=0.012` 缩放成稳定敏感研究构型，增量荷载改为塔顶
2009 节点 +X 水平力，靠 P-Δ 软化产生 λ_max 与下降段。

独立锚点（不读求解器结论，直接从冻结源文件重算）：

  1. `.dat` 的 16 个 `L0` ↔ `.cfg.json` 的 `l0`：同源下料长度；
  2. ANSYS `N,`/`E,` 弦长 ↔ `L0·(1+ε)`、`INISTATE eps ≡ EPS0`：
     单元号按 `E,` 行创建顺序，索单元为 48–63；
  3. `FINC_NORM` ↔ `.dat` 的 `*STABILITY_LOAD` 范数；
  4. `CTRL0` ↔ `[U](2009,1)`（`[U]` 段 10 位有效数字）；
  5. FEM `λ(ΔUX)` ↔ ANSYS 两段 `λ = TIME − 1`：因 λ 在极值点后下降，
     按单调自变量 ΔUX 插值 λ 逐点对拍，上升 / 下降 / 荷载控制分支
     分别统计最大偏差与均值。

用法：
    python compare_bridge_limit_point_3d.py --bin <build_dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0]))
from replay_tolerance import compare_replay  # noqa: E402
from build_reference_limit_point import (  # noqa: E402
    ANSYS_JOBS, CTRL, EPS0, PINC, REF_NAME, STEM,
    ansys_errors, apdl_anchor, fem_metrics, limit_of, normalized,
    read_curve, step2_curve,
)

EXE = "cable_solver.exe"
DEFAULT_ID = "lp_default"
TOL_RES_REL = 1.0e-8         # RES0 / |F0|
TOL_RES_RATIO = (0.1, 10.0)  # RES0 / [NEWTON] RESLAST
TOL_CTRL0 = 1.0e-10          # CTRL0 ≡ [U](2009,1)；[U] 段 10 位有效数字
TOL_FINC = 1.0e-9            # FINC_NORM ≡ |*STABILITY_LOAD|
TOL_L0 = 1.0e-9              # .dat L0 ↔ .cfg.json l0
TOL_EPS = 1.0e-12            # INISTATE eps ≡ EPS0
TOL_CHORD = 1.0e-6           # L0·(1+eps) ↔ ANSYS 弦长 (m)
TOL_U0 = 1.0e-9             # 两个 ANSYS job 的 step1 成桥态应一致
TOL_LAMBDA_MAX_REL = 1.0e-2  # 极值点荷载因子 FEM vs ANSYS
TOL_UX_LIMIT_REL = 1.0e-2   # 极值点位移 FEM vs ANSYS
TOL_ASCEND = 8.0e-3
TOL_DESCEND = 2.5e-2
TOL_LOAD = 1.2e-2
TOL_MEAN = 6.0e-3
TOL_RESIDUAL = 100.0        # 单个平衡点残差上限 (N)

FAILURES: list[str] = []
CHECKS = [0]


def ok(cond: bool, msg: str) -> None:
    CHECKS[0] += 1
    if not cond:
        FAILURES.append(msg)


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-30)


def interp_lambda_by_ctrl(curve: list[list[float]], ctrl: float) -> float:
    """按 ΔCTRL（单调自变量）在 (ΔCTRL, λ) 曲线上线性插值 λ。

    与 BS 侧 `run_cases_limit_point.py:interp_lambda_by_ctrl` 同口径；
    λ 本身在极值点后是非单调的，因此不能直接拿 λ 当插值自变量。
    """
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


def stability_load(dat: Path) -> tuple[float | None, list[tuple[int, int, float]]]:
    """从 .dat 独立读 *STABILITY_LOAD，返回 (范数, [(node, dof, value)])。"""
    sec = False
    rows: list[tuple[int, int, float]] = []
    for ln in dat.read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if s.startswith("*"):
            sec = s == "*STABILITY_LOAD"
            continue
        if sec and s:
            q = s.split()
            if len(q) >= 3:
                rows.append((int(q[0]), int(q[1]), float(q[2])))
    if not rows:
        return None, []
    return math.sqrt(sum(v * v for _, _, v in rows)), rows


def dat_l0(dat: Path) -> dict[int, float]:
    return {
        int(eid): float(l0)
        for eid, _, _, l0 in re.findall(
            r"^\*EN,(90\d\d),(\d+),(\d+) L0=([\d.eE+-]+)",
            dat.read_text(encoding="utf-8"), re.M)
    }


def run_exe(exe: Path, dat: Path, cwd: Path) -> int:
    p = subprocess.run([str(exe), dat.name], cwd=str(cwd),
                       capture_output=True, text=True, errors="replace", timeout=900)
    return p.returncode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    a = ap.parse_args()

    exe = Path(a.bin) / EXE
    cases = Path(a.cases)
    ref_file = Path(a.reference) / REF_NAME
    ok(exe.is_file(), f"缺少 {EXE}: {exe}")
    ok(ref_file.is_file(), f"缺少黄金值: {ref_file}")
    if FAILURES:
        print("[FAIL]", *FAILURES, sep="\n  - ")
        return 1
    ref = json.loads(ref_file.read_text(encoding="utf-8"))

    # ---- 1. 冻结证据文件完整性 ----------------------------------------
    ok(ref.get("case") == "CB3D22", "黄金值 case 标记不是 CB3D22")
    ok(ref.get("tower_scale") == 0.012, "TOWER_SCALE 漂移")
    ok(ref.get("eps0") == EPS0 and ref.get("pinc") == PINC, "eps0 / pinc 漂移")
    ok(ref.get("ctrl") == list(CTRL), "控制节点 / 自由度漂移")
    actual = sorted(p.name for p in cases.iterdir() if p.is_file())
    ok(actual == sorted(ref["files"]),
       f"冻结文件集合漂移: {set(actual) ^ set(ref['files'])}")
    texts: dict[str, str] = {}
    for name, meta in ref["files"].items():
        p = cases / name
        ok(p.is_file(), f"缺冻结文件 {name}")
        if not p.is_file():
            continue
        raw = p.read_bytes()
        ok(hashlib.sha256(raw).hexdigest() == meta["sha256"], f"{name} sha256 不符")
        ok(len(raw) == meta["n_bytes"], f"{name} 字节数不符")
        # 哈希 / 字节数按原始字节核对；文本解析统一换行，与生成器里的
        # Path.read_text()（通用换行归一化）保持同一口径，否则 CRLF 会让
        # `\[STABILITY\]\n` 与 `^...$` 一类锚点失配。
        texts[name] = raw.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    print(f"[PASS] frozen evidence integrity: {len(ref['files'])} files")

    fem_txt = texts.get(f"{STEM}_cable.txt", "")
    dflt_txt = texts.get(f"{DEFAULT_ID}_cable.txt", "")
    dat_path = cases / f"{STEM}.dat"
    cfg = json.loads(texts[f"{STEM}.cfg.json"])
    ok(cfg.get("case") == "CB3D22" and cfg.get("tower_scale") == 0.012,
       ".cfg.json 与黄金值构型不符")
    ok(cfg.get("ctrl") == list(CTRL) and abs(cfg.get("pinc", 0) - PINC) < TOL_FINC,
       ".cfg.json 的控制点 / 增量荷载不符")
    ok(abs(cfg.get("eps0", 0) - EPS0) < TOL_EPS, ".cfg.json 的 eps0 不符")

    # ---- 2. 黄金值自洽：从冻结证据重算全部派生量 ---------------------
    got_fem = fem_metrics(fem_txt)
    ok(got_fem == ref["fem"]["main"], "FEM 重算指标与黄金值不符")
    ok(ref["fem"]["default"] == {
        "has_stability": "[STABILITY]" in dflt_txt,
        "n_lines": len(dflt_txt.splitlines()),
    }, "FEM 默认算例重算指标与黄金值不符")
    fpts = got_fem["stability"]["points"]
    fem_curve = [[p[2], p[1]] for p in fpts]          # (ΔUX, λ)
    fi, flam, fctrl = limit_of(fem_curve)
    ndesc = sum(1 for k in range(fi + 1, len(fem_curve))
                if fem_curve[k][1] < fem_curve[k - 1][1] - 1e-9)
    ok(ref["fem"]["limit_index"] == fi and ref["fem"]["limit_lambda"] == flam
       and ref["fem"]["limit_ctrl"] == fctrl, "FEM 极值点重算与黄金值不符")
    ok(ref["fem"]["n_descending"] == ndesc, "FEM 下降段步数重算与黄金值不符")
    ok(ref["fem"]["lambda_end"] == (fem_curve[-1][1] if fem_curve else 0.0),
       "FEM 终点 λ 重算与黄金值不符")

    arcs: dict[str, list[list[float]]] = {}
    for job in ANSYS_JOBS:
        inp = texts[f"{job}.inp"]
        out = texts[f"{job}_ansys.out"]
        rows = read_curve(texts[f"{job}_curve.csv"])
        u0, s2 = step2_curve(rows)
        i, lam, ctrl = limit_of(s2) if s2 else (0, 0.0, 0.0)
        got = {
            "n_solve": inp.count("SOLVE"),
            "n_fx": len(re.findall(r"^F,2009,FX,", inp, re.M)),
            "n_inistate": len(re.findall(r"^INISTATE,DEFINE,", inp, re.M)),
            "arclen": "ARCLEN,ON" in inp,
            "error_messages": ansys_errors(out),
            "u0": u0, "n_step2_points": len(s2),
            "limit_index": i, "limit_lambda": lam, "limit_ctrl": ctrl,
            "curve": s2,
        }
        ok(got == ref["ansys"][job], f"{job} 重算指标与黄金值不符")
        ok(apdl_anchor(texts[f"{STEM}.dat"], inp) == ref["anchors"][job],
           f"{job} 同源锚点重算与黄金值不符")
        arcs[job] = s2
    print("[PASS] golden reference self-consistency (recomputed from frozen files)")

    # ---- 3. 同源 / 独立锚点 -------------------------------------------
    l0 = dat_l0(dat_path)
    cfg_l0 = {int(k): float(v) for k, v in cfg["l0"].items()}
    ok(len(l0) == 16 and set(l0) == set(cfg_l0),
       ".dat 与 .cfg.json 的 16 根索 L0 集合不一致")
    ok(all(abs(l0[k] - cfg_l0[k]) < TOL_L0 for k in l0 if k in cfg_l0),
       ".dat L0 与 .cfg.json 数值不符")
    finc, loads = stability_load(dat_path)
    state = got_fem["stability"]["state"]
    ok(finc is not None and abs(state["FINC_NORM"] - finc) < TOL_FINC,
       "FINC_NORM 与 *STABILITY_LOAD 独立范数不符")
    ok(len(loads) == 1 and loads[0][0] == CTRL[0] and loads[0][1] == CTRL[1],
       "*STABILITY_LOAD 控制点不是 2009/1")
    ok(loads and abs(loads[0][2] - PINC) < TOL_FINC,
       "*STABILITY_LOAD 幅值与 PINC 不同源")
    head = got_fem["stability"]["head"]
    ok(head.get("CTRLNODE") == str(CTRL[0]) and head.get("CTRLDOF") == str(CTRL[1]),
       "[STABILITY] 控制点与 .cfg.json 不符")
    ok(abs(state["CTRL0"] - got_fem["u_ctrl"]) < TOL_CTRL0,
       "CTRL0 与 [U](2009,1) 不同源（[U] 段 10 位有效数字）")
    for job in ANSYS_JOBS:
        an = ref["anchors"][job]
        ok(an["n_nodes"] == 66 and an["n_elems"] == 63,
           f"{job} ANSYS 网格节点/单元数与同源模型不符")
        ok(an["n_inistate"] == 16 and an["n_l0_matched"] == 16,
           f"{job} INISTATE / L0 匹配数不是 16")
        ok(an["eps_max_dev"] < TOL_EPS,
           f"{job} INISTATE eps 偏离 EPS0 ({an['eps_max_dev']:.3e})")
        ok(an["chord_max_dev"] < TOL_CHORD,
           f"{job} L0·(1+eps) 与 ANSYS 弦长不符 ({an['chord_max_dev']:.3e} m)")
    print(f"[PASS] homology anchors: 16 cables, eps dev="
          f"{ref['anchors']['lp_arc']['eps_max_dev']:.3e}, "
          f"chord dev={ref['anchors']['lp_arc']['chord_max_dev']:.3e} m")

    # ---- 4. FEM 物理判据：起点闭合 + 真实极限点 + 下降段 ---------------
    st = got_fem["stability"]
    ok(st.get("present") is True, "FEM 缺 [STABILITY] 段")
    ok(head.get("BASIS") == "completion_state", "FEM BASIS 不是 completion_state")
    ok(head.get("KIND") == "spherical_arc_length", "FEM 不是球面弧长法")
    ok(head.get("STATUS") == "ok" and head.get("ERROR") == "none"
       and head.get("CONVERGED") == "1", "FEM 路径未收敛")
    rr = state["RES0"] / state["F0_NORM"]
    ok(rr < TOL_RES_REL, f"FEM RES0/|F0|={rr:.3e} 未闭合")
    ratio = state["RES0"] / got_fem["reslast"] if got_fem["reslast"] else float("inf")
    ok(TOL_RES_RATIO[0] < ratio < TOL_RES_RATIO[1],
       f"FEM RES0/RESLAST={ratio:.3e} 与内层 Newton 不同量级")
    ok(int(head["NPOINT"]) == len(fpts) and int(head["NSTEP"]) == len(fpts) - 1,
       "NPOINT / NSTEP 与逐平衡点行数不一致")
    ok(abs(float(head["LAMBDA_FINAL"]) - fem_curve[-1][1]) < 1e-9,
       "LAMBDA_FINAL 与末点不一致")
    ok(len(fpts) > 400 and fpts[1][1] > 0.0, "FEM 路径点数不足或首步 λ 非正")
    ok(max(p[3] for p in fpts) < TOL_RESIDUAL,
       f"FEM 存在残差超 {TOL_RESIDUAL} N 的平衡点")
    ok(0 < fi < len(fem_curve) - 100, "FEM 极值点不在路径中段")
    ok(ndesc >= 100, f"FEM 下降段步数 {ndesc} < 100")
    ok((flam - fem_curve[-1][1]) / flam > 0.10,
       "FEM λ 降幅不足 10%，不构成下降段")
    ok(got_fem["buckling_error_kept"], "[BUCKLING] ERROR 段被删")
    for tag in ("CABLE", "U", "LENGTH", "NEWTON", "BUCKLING", "STABILITY"):
        ok(f"[{tag}]" in fem_txt, f"FEM 缺 [{tag}] 段")
    ok("[STABILITY]" not in dflt_txt, "STABILITY=0 泄漏了 [STABILITY] 段")
    du = re.search(r"^2009 1 (\S+)", dflt_txt, re.M)
    ok(du is not None and abs(float(du.group(1)) - got_fem["u_ctrl"]) < TOL_CTRL0,
       "默认算例与稳定算例的 [U](2009,1) 不一致（起点被污染）")
    print(f"[PASS] FEM: RES0/|F0|={rr:.3e}, NPOINT={len(fpts)}, "
          f"lambda_max={flam:.7g} @ dUX={fctrl:.7g}, "
          f"descending={ndesc}, lambda_end={fem_curve[-1][1]:.7g}")

    # ---- 5. ANSYS 两段求解证据 ----------------------------------------
    ok(abs(ref["ansys"]["lp_load"]["u0"] - ref["ansys"]["lp_arc"]["u0"]) < TOL_U0,
       "两个 ANSYS job 的 step1 成桥态 UX 不一致")
    for job in ANSYS_JOBS:
        an, inp, out = ref["ansys"][job], texts[f"{job}.inp"], texts[f"{job}_ansys.out"]
        ok(an["n_solve"] == 2, f"{job} 不是两段 SOLVE")
        ok(an["n_fx"] == 1 and an["n_inistate"] == 16,
           f"{job} 增量荷载 / INISTATE 条数不符")
        ok(an["error_messages"] == 0, f"{job} ANSYS 有 error 消息")
        ok("NUMBER OF ERROR   MESSAGES ENCOUNTERED=          0" in out,
           f"{job} ANSYS .out 未记录 0 error")
        ok("INISTATE,SET,DTYP,EPEL" in inp, f"{job} 缺 INISTATE,SET,DTYP,EPEL")
        ok("ARCTRM,U,1.6,2009,UX" in inp, f"{job} 缺 ARCTRM 位移终止判据")
        ok(an["n_step2_points"] > 20 and an["curve"][0][1] > 0.0,
           f"{job} step2 曲线点数不足或起点 λ 非正")
    la, arc = ref["ansys"]["lp_load"], arcs["lp_arc"]
    ok(la["arclen"] is False, "lp_load 不应启用 ARCLEN（荷载控制上升段）")
    ok(ref["ansys"]["lp_arc"]["arclen"] is True, "lp_arc 未启用 ARCLEN")
    ok("ARCLEN,ON,0.005,1.0E-8" in texts["lp_arc.inp"], "lp_arc ARCLEN 参数漂移")
    ok("TIME,3" in texts["lp_arc.inp"] and "TIME,2.5" in texts["lp_load.inp"],
       "两段 TIME 终点与 cfg 不符")
    ok(abs(la["limit_lambda"] - 1.5) < 1e-6, "lp_load 荷载控制终点 λ 不是 1.5")
    ok(la["curve"][-1][1] < flam, "lp_load 终点 λ 应小于 FEM 极值点 λ_max")
    ok(all(la["curve"][k][1] >= la["curve"][k - 1][1] - 1e-9
           for k in range(1, len(la["curve"]))),
       "lp_load 荷载控制段不应出现下降分支")
    ai, alam, actrl = limit_of(arc)
    ok(0 < ai < len(arc) - 100, "ANSYS lp_arc 极值点不在路径中段")
    ok(arc[-1][1] < alam - 0.05, "ANSYS lp_arc 未越过极值点进入下降段")
    ok("indefinite tangent matrix is detected" in texts["lp_arc_ansys.out"],
       "lp_arc 未出现不定切线刚度（极限点物理标志）")
    ok("indefinite tangent matrix is detected" not in texts["lp_load_ansys.out"],
       "lp_load 不应出现不定切线刚度（未到极值点）")
    print(f"[PASS] ANSYS: lp_load step2={la['n_step2_points']} pts (0 error), "
          f"lp_arc step2={len(arc)} pts, lambda_max={alam:.7g} @ dUX={actrl:.7g}")

    # ---- 6. FEM ↔ ANSYS 逐点对拍（按 ΔUX 插值 λ）----------------------
    ok(abs(flam - alam) / flam < TOL_LAMBDA_MAX_REL,
       f"极值点 λ 相对差 {abs(flam - alam) / flam:.3e} 超容差")
    ok(abs(fctrl - actrl) / fctrl < TOL_UX_LIMIT_REL,
       f"极值点 ΔUX 相对差 {abs(fctrl - actrl) / fctrl:.3e} 超容差")
    dv_asc: list[float] = []
    dv_desc: list[float] = []
    for ctrl, lam in fem_curve:
        if ctrl > arc[-1][0] or ctrl < arc[0][0]:
            continue
        d = abs(lam - interp_lambda_by_ctrl(arc, ctrl))
        (dv_asc if ctrl <= fctrl else dv_desc).append(d)
    wa = max(dv_asc) if dv_asc else float("inf")
    wd = max(dv_desc) if dv_desc else float("inf")
    ma = sum(dv_asc) / len(dv_asc) if dv_asc else float("inf")
    md = sum(dv_desc) / len(dv_desc) if dv_desc else float("inf")
    ok(len(dv_asc) > 50 and wa < TOL_ASCEND, f"上升段逐点 |Δλ| max={wa:.3e}")
    ok(len(dv_desc) > 100 and wd < TOL_DESCEND, f"下降段逐点 |Δλ| max={wd:.3e}")
    ok(ma < TOL_MEAN and md < TOL_MEAN,
       f"逐段平均 |Δλ| 超容差: asc={ma:.3e} desc={md:.3e}")
    load = la["curve"]
    dv_load = [abs(lam - interp_lambda_by_ctrl(load, ctrl))
               for ctrl, lam in fem_curve
               if load[0][0] <= ctrl <= load[-1][0]]
    wl = max(dv_load) if dv_load else float("inf")
    ml = sum(dv_load) / len(dv_load) if dv_load else float("inf")
    ok(len(dv_load) > 50 and wl < TOL_LOAD, f"lp_load 上升段逐点 |Δλ| max={wl:.3e}")
    ok(ml < TOL_MEAN, f"lp_load 平均 |Δλ| {ml:.3e} 超容差")
    print(f"[PASS] FEM↔ANSYS by dUX: d-lambda@limit={abs(flam-alam)/flam:.3e}, "
          f"d-UX@limit={abs(fctrl-actrl)/fctrl:.3e}; "
          f"asc {len(dv_asc)} pts worst={wa:.3e} mean={ma:.3e}; "
          f"desc {len(dv_desc)} pts worst={wd:.3e} mean={md:.3e}; "
          f"lp_load worst={wl:.3e} mean={ml:.3e}")

    # ---- 7. 当前二进制离线重跑（不重跑 ANSYS）------------------------
    with tempfile.TemporaryDirectory(prefix="cb3d22lp_") as td:
        tmp = Path(td)
        for name in (f"{STEM}.dat", f"{DEFAULT_ID}.dat"):
            shutil.copy2(cases / name, tmp / name)
        for name in (STEM, DEFAULT_ID):
            rc = run_exe(exe, tmp / f"{name}.dat", tmp)
            got = tmp / f"{name}_cable.txt"
            ok(rc == 0, f"{name} 离线重跑退出码 {rc}")
            ok(got.is_file(), f"{name} 离线重跑未产出结果文件")
            if got.is_file():
                fresh = got.read_text(encoding="utf-8", errors="replace")
                good, detail = compare_replay(texts[f"{name}_cable.txt"], fresh)
                ok(good, f"{name} 离线重跑结果与冻结输出不一致: {detail}")
                if name == STEM:
                    print(f"[..] replay {name}: {detail}")
    print("[PASS] offline FEM replay: 2 cases "
          "(cross-toolchain tolerance, see tests/benchmark/replay_tolerance.py)")

    if FAILURES:
        print(f"\n[FAIL] {len(FAILURES)}/{CHECKS[0]} 项未通过：")
        for f in FAILURES:
            print("  - " + f)
        return 1
    print(f"\n[PASS] bridge_limit_point_3d benchmark: "
          f"all {CHECKS[0]} checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
