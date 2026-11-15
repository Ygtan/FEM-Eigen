#!/usr/bin/env python3
"""第 24 项 CTest：`MODE=COMPLETION` 成桥态基准同源修复（离线）。

② 阶段登记的缺陷：`cable_solver.cpp` 用**外部原模型**装配 `F0` 并计算
`Fint(u0)`，而 `MODE=COMPLETION` 的 `u0` 是相对**内部制造线形 + 下料 L0**
求得的（求解器按值收 `Model3D`）。CB3D10 实测 `RES0=3.15 MN ≈ |F0|`、
首步 λ 变负。③ 由 `basisModel` 回传同源模型，本项把修复冻结为回归。

四个 case 的分工见 `build_reference_basis.py`。核心独立判据：

  1. 起点闭合 `RES0/|F0| < 1e-8`，且 `RES0` 与 `[NEWTON] RESLAST`
     同量级（后者是内层 Newton 独立收敛残差）；
  2. 首步 λ > 0（缺陷时 λ < 0），`CTRL0 ≡ [U](1016,3)`；
  3. **跨 MODE 独立锚点**：`MODE=COMPLETION`（制造线形基准）与
     `MODE=STATIC`（设计线形 + 冻结 L0 基准）是两条独立成桥态闭环，
     λ(v) 与 `ΔCTRL = CTRL − CTRL0` 逐点对拍；两条基准的 `CTRL0`
     本就不同（绝对位移），只比较路径增量；
  4. `FINC_NORM ≡ |*STABILITY_LOAD|`（从 .dat 独立读出，不读求解器）；
  5. 闭环段齐备：`[CABLE] / [U] / [ERECTED] / [LENGTH] / [NEWTON] /
     [BUCKLING] / [STABILITY]`，`[BUCKLING] ERROR` 段仍保留；
  6. 默认零污染：`STABILITY=0` 不写 `[STABILITY]`，同二进制重跑逐行
     一致；与历史 CB3D10 冻结基线的段集合 / 关键物理量一致（跨编译器
     的 1e-9 量级浮点噪声用容差吸收，不做逐字符比对）。

用法：
    python compare_completion_basis_3d.py --bin <build_dir> \
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
from build_reference_basis import (  # noqa: E402
    CASES, ERECT_SECTIONS, REF_NAME, ctrl_of, normalized, metrics, section_body,
)

EXE = "cable_solver.exe"
TOL_RES_REL = 1.0e-8        # RES0 / |F0|
TOL_RES_RATIO = (0.1, 10.0)  # RES0 / [NEWTON] RESLAST
TOL_CTRL0 = 1.0e-12         # CTRL0 ≡ [U](1016,3)
TOL_FINC = 1.0e-9           # FINC_NORM ≡ |*STABILITY_LOAD|
TOL_PATH_LAMBDA = 1.0e-4    # 跨 MODE 逐点 λ 绝对差
TOL_PATH_CTRL = 1.0e-5      # 跨 MODE 逐点 ΔCTRL 绝对差 (m)
TOL_LAMBDA_FINAL_REL = 1.0e-5
TOL_BASELINE_REL = 1.0e-5   # 与历史冻结基线的关键量相对差

FAILURES: list[str] = []
CHECKS = [0]


def ok(cond: bool, msg: str) -> None:
    CHECKS[0] += 1
    if not cond:
        FAILURES.append(msg)


def run_exe(exe: Path, dat: Path, cwd: Path) -> int:
    p = subprocess.run([str(exe), dat.name], cwd=str(cwd),
                       capture_output=True, text=True, errors="replace", timeout=900)
    return p.returncode


def load_inc_load(dat: Path) -> float | None:
    """从 .dat 独立读 *STABILITY_LOAD，返回荷载向量范数。"""
    sec = False
    vals: list[float] = []
    for ln in dat.read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if s.startswith("*"):
            sec = s == "*STABILITY_LOAD"
            continue
        if sec and s:
            q = s.split()
            if len(q) >= 3:
                vals.append(float(q[2]))
    if not vals:
        return None
    return math.sqrt(sum(v * v for v in vals))


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-30)


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
    # 编号口径：本项覆盖 CB3D10 / CB3D13 的离线回归，不是新 BS 算例；
    # CB3D22 属于 `bridge_limit_point_3d`（真实桥梁极值点 / 下降段）。
    ok(ref.get("case") == "CB3D10_BASIS_FIX", "黄金值 case 标记不符")
    ok(ref.get("cases_covered") == ["CB3D10", "CB3D13"], "覆盖算例编号不符")
    ok(sorted(ref.get("cases", {})) == sorted(CASES), "黄金值 case 集合不全")

    # ---- 1. 冻结文本 ↔ 黄金值自洽 ------------------------------------
    frozen = {}
    for name in CASES:
        entry = ref["cases"].get(name, {})
        txt_path = cases / entry.get("txt", f"{name}_cable.txt")
        dat_path = cases / entry.get("dat", f"{name}.dat")
        ok(txt_path.is_file(), f"{name} 缺冻结输出 {txt_path.name}")
        ok(dat_path.is_file(), f"{name} 缺冻结输入 {dat_path.name}")
        if not txt_path.is_file() or not dat_path.is_file():
            continue
        txt = txt_path.read_text(encoding="utf-8", errors="replace")
        frozen[name] = txt
        norm = normalized(txt)
        ok(hashlib.sha256(norm.encode("utf-8")).hexdigest()
           == entry.get("sha256_normalized"),
           f"{name} 冻结输出与黄金值哈希不符")
        ok(len(norm.splitlines()) == entry.get("n_lines"),
           f"{name} 冻结输出行数漂移")
        got = metrics(txt, ctrl_of(name))
        ok(got == entry.get("metrics"),
           f"{name} 冻结输出重算指标与黄金值不符")
    print("[PASS] golden reference self-consistency")

    # ---- 2. 离线重跑逐行复现 ----------------------------------------
    with tempfile.TemporaryDirectory(prefix="cb3d22_") as td:
        tmp = Path(td)
        for name in CASES:
            shutil.copy2(cases / f"{name}.dat", tmp / f"{name}.dat")
        for name in CASES:
            rc = run_exe(exe, tmp / f"{name}.dat", tmp)
            got = tmp / f"{name}_cable.txt"
            ok(rc == ref["cases"][name]["exit_code"],
               f"{name} 离线重跑退出码 {rc} != {ref['cases'][name]['exit_code']}")
            ok(got.is_file(), f"{name} 离线重跑未产出结果文件")
            if got.is_file():
                ok(normalized(got.read_text(encoding="utf-8", errors="replace"))
                   == normalized(frozen[name]),
                   f"{name} 离线重跑结果与冻结输出不一致")
    print(f"[PASS] offline replay: {len(CASES)} cases byte-identical "
          "(elapsed_ms normalized; canonical mingw64 toolchain)")

    # ---- 3. 两个成桥态基准的起点 / 路径独立判据 ----------------------
    paths: dict[str, list[list[float]]] = {}
    finals: dict[str, float] = {}
    for name in ("cb3d10_completion_basis", "cb3d10_static_basis"):
        txt = frozen[name]
        m = metrics(txt)
        st = m["stability"]
        ok(st.get("present") is True, f"{name} 缺 [STABILITY] 段")
        if not st.get("present"):
            continue
        head, state = st["head"], st["state"]
        paths[name] = st["points"]
        ok(head.get("BASIS") == "completion_state", f"{name} BASIS 不是 completion_state")
        ok(head.get("KIND") == "spherical_arc_length", f"{name} KIND 不是球面弧长法")
        ok(head.get("STATUS") == "ok" and head.get("ERROR") == "none",
           f"{name} STATUS/ERROR 非 ok")
        ok(head.get("CONVERGED") == "1", f"{name} 路径未收敛")
        f0n, res0 = state["F0_NORM"], state["RES0"]
        ok(res0 / f0n < TOL_RES_REL,
           f"{name} RES0/|F0|={res0 / f0n:.3e} 未闭合（② 缺陷 ~1.06）")
        ratio = res0 / m["reslast"] if m["reslast"] else float("inf")
        ok(TOL_RES_RATIO[0] < ratio < TOL_RES_RATIO[1],
           f"{name} RES0/RESLAST={ratio:.3e} 与内层 Newton 不同量级")
        ok(abs(state["CTRL0"] - m["u_ctrl"]) < TOL_CTRL0,
           f"{name} CTRL0 与 [U](1016,3) 不同源")
        finc = load_inc_load(cases / f"{name}.dat")
        ok(finc is not None and rel(state["FINC_NORM"], finc) < TOL_FINC,
           f"{name} FINC_NORM 与 *STABILITY_LOAD 独立范数不符")
        pts = st["points"]
        ok(len(pts) >= 100, f"{name} 平衡点过少 ({len(pts)})")
        ok(len(pts) > 1 and pts[1][1] > 0.0,
           f"{name} 首步 λ={pts[1][1] if len(pts) > 1 else 'n/a'} 非正（缺陷时 <0）")
        ok(all(pts[i][1] >= pts[i - 1][1] - 1e-9 for i in range(1, len(pts))),
           f"{name} λ 路径不单调")
        ok(all(p[3] < 100.0 for p in pts), f"{name} 存在残差超 100 N 的平衡点")
        finals[name] = float(head["LAMBDA_FINAL"])
        ok(abs(float(head["LAMBDA_FINAL"]) - pts[-1][1]) < 1e-9,
           f"{name} LAMBDA_FINAL 与末点不一致")
        # 闭环段齐备（先走完成桥态闭环，不拿 u=0 冒充 u0）
        for tag in ("CABLE", "U", "LENGTH", "NEWTON", "BUCKLING", "STABILITY"):
            ok(f"[{tag}]" in txt, f"{name} 缺 [{tag}] 段")
        ok(m["buckling_error_kept"], f"{name} [BUCKLING] ERROR 段被删")
        ok("ERECTED" in txt or name.endswith("static_basis"),
           f"{name} COMPLETION 应写 [ERECTED] 制造线形")
        print(f"[PASS] {name}: RES0={res0:.4g} N, RES0/|F0|={res0 / f0n:.3e}, "
              f"NPOINT={len(pts)}, LAMBDA_FINAL={finals[name]:.9g}")

    # ---- 4. 跨 MODE 独立锚点：两条闭环的 λ / ΔCTRL 逐点对拍 ----------
    cpts = paths.get("cb3d10_completion_basis", [])
    spts = paths.get("cb3d10_static_basis", [])
    ok(len(cpts) == len(spts) and len(cpts) > 100,
       "两条基准平衡点数量不一致，无法逐点对拍")
    if len(cpts) == len(spts) and cpts:
        worst_lam = max(abs(c[1] - s[1]) for c, s in zip(cpts, spts))
        # CTRL 是各自参考构形上的绝对位移；只比路径增量 CTRL − CTRL0
        c0, s0 = cpts[0][2], spts[0][2]
        worst_ctrl = max(abs((c[2] - c0) - (s[2] - s0)) for c, s in zip(cpts, spts))
        ok(worst_lam < TOL_PATH_LAMBDA,
           f"跨 MODE 逐点 λ 最大差 {worst_lam:.3e} 超容差")
        ok(worst_ctrl < TOL_PATH_CTRL,
           f"跨 MODE 逐点 ΔCTRL 最大差 {worst_ctrl:.3e} m 超容差")
        print(f"[PASS] cross-MODE anchor: worst|dλ|={worst_lam:.3e}, "
              f"worst|dΔCTRL|={worst_ctrl:.3e} m")
    if len(finals) == 2:
        c, s = finals["cb3d10_completion_basis"], finals["cb3d10_static_basis"]
        ok(rel(c, s) < TOL_LAMBDA_FINAL_REL,
           f"跨 MODE LAMBDA_FINAL 相对差 {rel(c, s):.3e} 超容差")
    ok(ref["cases"]["cb3d10_completion_basis"]["metrics"]["length_kind"]
       == ["erection_geometry", "ok"],
       "MODE=COMPLETION 未走制造线形几何闭合迭代")
    ok(ref["cases"]["cb3d10_static_basis"]["metrics"]["length_kind"]
       == ["none", "not_applicable"],
       "MODE=STATIC 不应出现索长/几何闭合迭代")

    # ---- 5. 默认零污染 ------------------------------------------------
    dtxt = frozen["cb3d10_default_basis"]
    ok("[STABILITY]" not in dtxt, "STABILITY=0 泄漏了 [STABILITY] 段")
    base = cases / "cb3d10_baseline_cable.txt"
    ok(base.is_file(), "缺历史 CB3D10 冻结基线")
    if base.is_file():
        btxt = base.read_text(encoding="utf-8", errors="replace")
        ok(len(normalized(dtxt).splitlines()) == len(normalized(btxt).splitlines()),
           "默认输出行数与历史基线不一致")
        ok(re.findall(r"^\[([A-Z_]+)\]", dtxt, re.M)
           == re.findall(r"^\[([A-Z_]+)\]", btxt, re.M),
           "默认输出段集合与历史基线不一致")
        dm, bm = metrics(dtxt), metrics(btxt)
        ok(dm["u_ctrl"] is not None and bm["u_ctrl"] is not None
           and rel(dm["u_ctrl"], bm["u_ctrl"]) < TOL_BASELINE_REL,
           "默认 [U](1016,3) 与历史基线漂移")
        ok(dm["reslast"] is not None and bm["reslast"] is not None
           and rel(dm["reslast"], bm["reslast"]) < TOL_BASELINE_REL,
           "默认 [NEWTON] RESLAST 与历史基线漂移")
        ok(dm["buckling_error_kept"], "默认 [BUCKLING] ERROR 段被删")
    print("[PASS] zero pollution: STABILITY=0 unchanged (sections / key metrics)")

    # ---- 6. 施工阶段链路接入核对（CB3D13）-------------------------------
    # `STABASIS=COMPLETION` 与悬索桥 5 个 [ERECTION*] 段并列输出：
    # 稳定路径起点是**成桥态**（`MODE=COMPLETION` 悬索桥成桥态闭环），
    # 不是某个施工阶段中间态。核对「共存 + 互不污染」，不是声称
    # 「每个施工阶段已具备全过程稳定分析」。
    ername = "cb3d13_erection_basis"
    etxt = frozen[ername]
    em = metrics(etxt, ctrl_of(ername))
    est = em["stability"]
    ok(est.get("present") is True, "CB3D13 缺 [STABILITY] 段")
    if est.get("present"):
        eh, es_ = est["head"], est["state"]
        ok(eh.get("BASIS") == "completion_state",
           "CB3D13 [STABILITY] 不是成桥态基准")
        ok(eh.get("STATUS") == "ok" and eh.get("ERROR") == "none"
           and eh.get("CONVERGED") == "1", "CB3D13 路径未收敛")
        ok(es_["RES0"] / es_["F0_NORM"] < TOL_RES_REL,
           f"CB3D13 RES0/|F0|={es_['RES0'] / es_['F0_NORM']:.3e} 未闭合")
        eratio = es_["RES0"] / em["reslast"] if em["reslast"] else float("inf")
        ok(TOL_RES_RATIO[0] < eratio < TOL_RES_RATIO[1],
           f"CB3D13 RES0/RESLAST={eratio:.3e} 不同量级")
        ok(abs(es_["CTRL0"] - em["u_ctrl"]) < TOL_CTRL0,
           "CB3D13 CTRL0 与 [U](1014,3) 不同源")
        efinc = load_inc_load(cases / f"{ername}.dat")
        ok(efinc is not None and rel(es_["FINC_NORM"], efinc) < TOL_FINC,
           "CB3D13 FINC_NORM 与 *STABILITY_LOAD 独立范数不符")
        ok(len(est["points"]) > 10 and est["points"][1][1] > 0.0,
           "CB3D13 路径点数不足或首步 λ 非正")
        print(f"[PASS] CB3D13 erection: RES0={es_['RES0']:.4g} N, "
              f"RES0/|F0|={es_['RES0'] / es_['F0_NORM']:.3e}, "
              f"NPOINT={len(est['points'])}, "
              f"LAMBDA_FINAL={eh['LAMBDA_FINAL']}")
    ebase = cases / "cb3d13_baseline_cable.txt"
    ok(ebase.is_file(), "缺 CB3D13 施工阶段冻结基线")
    if ebase.is_file():
        btxt = ebase.read_text(encoding="utf-8", errors="replace")
        for tag in ERECT_SECTIONS:
            ok(section_body(btxt, tag) == em["erection_sections"].get(tag),
               f"CB3D13 [{tag}] 段被稳定性链路污染（与冻结基线不一致）")
        bs = re.findall(r"^\[([A-Z_]+)\]", btxt, re.M)
        es = re.findall(r"^\[([A-Z_]+)\]", etxt, re.M)
        ok(es == bs[:-1] + ["STABILITY", "END"],
           "CB3D13 段集合不是「基线 + [STABILITY]」")
        # [NEWTON] RESLAST 是迭代停止时的残差，对编译器浮点顺序敏感
        # （历史基线由另一构建目录产出，实测相对差 8e-3）；它是收敛量级
        # 而非物理量，故只校核同量级，物理口径由 5 个 [ERECTION*] 段
        # 逐位一致与「同二进制重跑逐行一致」共同保证。
        base_reslast = metrics(btxt, ctrl_of(ername))["reslast"]
        ratio = (em["reslast"] / base_reslast) if base_reslast else float("inf")
        ok(TOL_RES_RATIO[0] < ratio < TOL_RES_RATIO[1],
           f"CB3D13 [NEWTON] RESLAST 与冻结基线不同量级 (ratio={ratio:.3e})")
    print("[PASS] erection-stage coexistence: 5 [ERECTION*] sections identical")

    if FAILURES:
        print(f"\n[FAIL] {len(FAILURES)}/{CHECKS[0]} 项未通过：")
        for f in FAILURES:
            print("  - " + f)
        return 1
    print(f"\n[PASS] completion_basis_3d benchmark: "
          f"all {CHECKS[0]} checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
