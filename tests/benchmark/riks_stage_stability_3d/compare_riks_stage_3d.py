#!/usr/bin/env python3
"""第 27 项 CTest：PATHMETHOD=RIKS(STAGE 起点) 纯离线守卫。

检查五类：
  1. 冻结证据 sha256 / 字节数（集合不得增删）；
  2. 黄金值自洽：用同一套解析从冻结 .txt 重新派生 checks / strings；
  3. 独立闭式 oracle：von Mises 两杆 P(v) 逐点核验 RIKS 路径；
  4. 同问题双方法：SPHERE 与 RIKS 按 CTRL 插值交叉比较；
  5. 物理判据 + 当前二进制逐字节重放。

不重跑 ANSYS、不重新建模、不调用 LLM。闭式 P(v) 在本文件独立实现，
不复用生成器的实现，避免 oracle 与被测代码同源。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_reference_riks_stage import (  # noqa: E402
    CASES, CASE_SPEC, FROZEN_FILES, parse, lambda_at_ctrl,
    E, A, B, H, P0, PREF, L0,
)

EXE = "cable_solver.exe"
REF_FILE = "cb3d24_riks_stage.json"
TOL_CLOSED = 1.0e-8       # 闭式逐点：机器精度量级
TOL_CROSS = 3.0e-4        # 双方法按 CTRL 线性插值
TOL_CB_ANCHOR = 1.0e-5    # 与第 26 项冻结的 SPHERE 终点对拍


def p_closed(v: float) -> float:
    """独立闭式：外荷载幅值 P(v) = 2*EA*ln(L0/L)*(H-v)/L（对数应变）。"""
    length = math.hypot(B, H - v)
    return 2.0 * E * A * math.log(L0 / length) * (H - v) / length


class Guard:
    def __init__(self) -> None:
        self.n = 0
        self.fails: list[str] = []

    def eq(self, name: str, got, want) -> None:
        self.n += 1
        if got != want:
            self.fails.append(f"{name}: got={got!r} want={want!r}")

    def ok(self, name: str, cond: bool, detail: str = "") -> None:
        self.n += 1
        if not cond:
            self.fails.append(f"{name}: FAILED {detail}")

    def close(self, name: str, got: float, want: float, tol: float) -> None:
        self.n += 1
        if not (math.isfinite(got) and abs(got - want) <= tol):
            self.fails.append(f"{name}: |{got} - {want}| > {tol}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="含 cable_solver.exe 的构建目录")
    ap.add_argument("--cases", default=str(CASES))
    ap.add_argument("--reference", default=str(HERE / "reference"))
    a = ap.parse_args()

    cases = Path(a.cases)
    refpath = Path(a.reference) / REF_FILE
    exe = Path(a.bin) / EXE
    g = Guard()

    if not exe.exists():
        print(f"[FATAL] 找不到 {exe}")
        return 2
    if not refpath.exists():
        print(f"[FATAL] 找不到 {refpath}")
        return 2

    doc = json.loads(refpath.read_text(encoding="utf-8"))
    parsed = {stem: parse(cases / f"{stem}_stage_stability.txt")
              for stem, _, _ in CASE_SPEC}

    # ---- 1. 冻结证据完整性 ----
    present = sorted(p.name for p in cases.iterdir() if p.is_file())
    g.eq("manifest/file_set", present, sorted(FROZEN_FILES))
    for name, want in doc["manifest"].items():
        p = cases / name
        g.ok(f"manifest/{name}/exists", p.exists())
        if p.exists():
            g.eq(f"manifest/{name}/sha256",
                 hashlib.sha256(p.read_bytes()).hexdigest(), want["sha256"])
            g.eq(f"manifest/{name}/bytes", p.stat().st_size, want["bytes"])

    # ---- 2. 黄金值自洽（重新派生 checks / strings）----
    checks: list[list] = []
    strings: dict[str, str] = {}
    for stem, _, _ in CASE_SPEC:
        r = parsed[stem]
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
            st, sg = r["state"], r["stage"]
            checks += [
                [f"{stem}/ctrlnode", r["ctrlnode"]],
                [f"{stem}/ctrldof", r["ctrldof"]],
                [f"{stem}/strain", r["strain"]],
                [f"{stem}/psi", r["psi"]],
                [f"{stem}/lambda_final", r["lambda_final"]],
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
    g.eq("reference/checks", checks, doc["checks"])
    g.eq("reference/strings", strings, doc["strings"])

    riks = parsed["vm_stage_riks"]
    sphere = parsed["vm_stage_sphere"]
    cb = parsed["cb3d23_riks"]
    bad = parsed["st_bad_method"]

    # ---- 3. 独立闭式 oracle ----
    v0 = -riks["state"]["ctrl0"]
    g.close("criteria/V0", v0, doc["criteria"]["V0"], 1e-14)
    g.close("closed_form/RES0", riks["state"]["res0"],
            p_closed(v0) - P0, 1e-6)
    worst_c = 0.0
    n_closed = 0
    for p in riks["points"]:
        if p["step"] == 0:
            continue
        v = v0 - p["ctrl"]
        closed = (p_closed(v) - P0) / PREF
        worst_c = max(worst_c, abs(p["lambda"] - closed))
        n_closed += 1
    g.eq("closed_form/n", n_closed, len(riks["points"]) - 1)
    g.ok("closed_form/worst", worst_c <= TOL_CLOSED,
         f"worst={worst_c:.3e} tol={TOL_CLOSED}")
    g.close("criteria/closed_form_worst", worst_c,
            doc["criteria"]["closed_form_worst"], 1e-14)

    # ---- 4. 同问题双方法交叉 ----
    worst_x, n_cross = 0.0, 0
    for p in riks["points"]:
        if p["step"] == 0:
            continue
        y = lambda_at_ctrl(sphere["points"], p["ctrl"])
        if y is None:
            continue
        worst_x = max(worst_x, abs(p["lambda"] - y))
        n_cross += 1
    g.ok("cross_method/n", n_cross >= 150, f"n={n_cross}")
    g.ok("cross_method/worst", worst_x <= TOL_CROSS,
         f"worst={worst_x:.3e} tol={TOL_CROSS}")
    g.close("criteria/cross_method_worst", worst_x,
            doc["criteria"]["cross_method_worst"], 1e-14)

    # ---- 5. 物理判据 ----
    g.eq("riks/kind", riks["kind"], "wempner_riks_orthogonal")
    g.eq("riks/pathmethod", riks["pathmethod"], "RIKS")
    g.eq("riks/basis", riks["basis"], "stage_state")
    g.eq("sphere/kind", sphere["kind"], "spherical_arc_length")
    g.eq("sphere/pathmethod", sphere["pathmethod"], "SPHERE")
    g.eq("cb/kind", cb["kind"], "wempner_riks_orthogonal")
    g.eq("cb/pathmethod", cb["pathmethod"], "RIKS")
    g.eq("bad/error", bad["error"], "unsupported_path_method")
    # ERROR 段走 stabilityErrorSection3D，没有完整 KIND/PATHMETHOD 头
    g.eq("bad/kind", bad["kind"], "")
    g.eq("riks/stage", (riks["stage"]["id"], riks["stage"]["nelem"],
                        riks["stage"]["ndof"]), (1, 2, 9))
    g.eq("cb/stage", (cb["stage"]["id"], cb["stage"]["nelem"],
                      cb["stage"]["ndof"]), (2, 3, 12))
    g.ok("riks/nonzero_start", abs(v0) > 1e-3, f"v0={v0}")
    g.ok("riks/snap_through",
         riks["lambda_final"] < 0.0 and riks["converged"] == 1,
         f"lam_end={riks['lambda_final']}")
    g.ok("cb/res0_tiny",
         abs(cb["state"]["res0"]) / cb["state"]["f0_norm"] < 1e-9)
    # 第 26 项冻结的 SPHERE 终点（同源模型，独立外部锚点）
    g.close("cb/anchor_to_test26", cb["lambda_final"], 2.57577583,
            TOL_CB_ANCHOR)
    g.eq("criteria/cases", [c[0] for c in CASE_SPEC], doc["meta"]["cases"])

    # ---- 6. 当前二进制逐字节重放 ----
    with tempfile.TemporaryDirectory(prefix="riks_ctest_") as td:
        work = Path(td)
        for stem, want_points, _ in CASE_SPEC:
            src = cases / f"{stem}.dat"
            dat = work / src.name
            shutil.copy2(src, dat)
            cp = subprocess.run([str(exe), dat.name], cwd=str(work),
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, timeout=180)
            out = work / f"{stem}_stage_stability.txt"
            g.ok(f"replay/{stem}/output_exists", out.exists())
            if out.exists():
                frozen = (cases / f"{stem}_stage_stability.txt").read_bytes()
                g.eq(f"replay/{stem}/bytes_identical", out.read_bytes(), frozen)
            want_rc = 0 if want_points else 3
            g.eq(f"replay/{stem}/exit_code", cp.returncode, want_rc)

    total = g.n
    if g.fails:
        print(f"[FAIL] {len(g.fails)}/{total} checks failed")
        for line in g.fails[:60]:
            print("  -", line)
        if len(g.fails) > 60:
            print(f"  ... and {len(g.fails) - 60} more")
        return 1
    print(f"[PASS] {total}/{total} checks: "
          f"closed_form_worst={worst_c:.3e}, cross_method_worst={worst_x:.3e}, "
          f"cb_lambda_end={cb['lambda_final']:.10f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
