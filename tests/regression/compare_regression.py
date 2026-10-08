#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FEM-Eigen .dat 逐字节回归比对（CTest 用）。

对每个 verify_*.dat 算例，用当前构建出的 example_*.exe 重算，把输出与
tests/regression/baselines/ 下的“黄金基线”比对：
  1) 逐字节相同；
  2) 否则把数值为 0 的 token 规整成 "0"（消除 -0 / 0.0 / 0e+00 这类写法差异）后逐字节相同；
  3) 否则做数值容差比对（相对 1e-6、绝对 1e-8，且双方 |v|<1e-6 视为等价噪声）。
三步都通过才判 PASS，否则 FAIL 并给出最大偏差。

生成/刷新基线：
    python compare_regression.py --freeze --bin <build> --src <source> --baseline <dir>
比对模式（CTest 调用）：
    python compare_regression.py --bin <build> --src <source> --baseline <dir>
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RTOL = 1e-6
ATOL = 1e-8
ZERO = 1e-6
NUMBER_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

# (exe, dat, extra_args, [输出文件名])
CASES = [
    ("example_mixed.exe", "verify_ce.dat", [], ["verify_ce.result"]),
    ("example_mixed.exe", "verify_prestress_equiv.dat", [], ["verify_prestress_equiv.result"]),
    ("example_mixed.exe", "verify_prestress_tendon.dat", [], ["verify_prestress_tendon.result"]),
    ("example_mixed.exe", "verify_prestress_ecc_tendon.dat", [], ["verify_prestress_ecc_tendon.result"]),
    ("example_mixed.exe", "verify_selfweight.dat", [], ["verify_selfweight.result"]),
    ("example_mixed.exe", "verify_creep_shrinkage.dat", [], ["verify_creep_shrinkage.result"]),
    ("example_mixed.exe", "verify_cable_nonlinear.dat", [], ["verify_cable_nonlinear.result"]),
    ("example_mixed.exe", "verify_cable_slack.dat", [], ["verify_cable_slack.result"]),
    ("example_mixed.exe", "verify_cable_tension.dat", [], ["verify_cable_tension.result"]),
    ("example_mixed.exe", "verify_stress.dat", [], ["verify_stress.result"]),
    ("example_mixed.exe", "verify_strut.dat", [], ["verify_strut.result"]),
    ("example_mixed.exe", "verify_temp_gradient.dat", [], ["verify_temp_gradient.result"]),
    ("example_mixed.exe", "verify_temp_profile.dat", [], ["verify_temp_profile.result"]),
    ("example_mixed.exe", "verify_temp_uniform.dat", [], ["verify_temp_uniform.result"]),
    ("example_mixed.exe", "verify_udl.dat", [], ["verify_udl.result"]),
    ("example_influence.exe", "verify_influence.dat", [],
     ["verify_influence_influence_lane1.txt", "verify_influence_envelope_lane1.txt"]),
    ("example_influence.exe", "verify_lane_nonlinear.dat", [],
     ["verify_lane_nonlinear_influence_lane1.txt", "verify_lane_nonlinear_envelope_lane1.txt"]),
    ("example_stage.exe", "verify_stage.dat", [],
     ["verify_stage_stages.txt", "verify_stage_camber.txt"]),
    ("example_stage.exe", "verify_stage.dat", ["--backward"],
     ["verify_stage_stages_backward.txt"]),
    ("example_stage.exe", "verify_stage.dat", ["--camber-iterate"],
     ["verify_stage_stages_camber_iterated.txt"]),
    ("example_stage.exe", "verify_stage_advanced.dat", [], ["verify_stage_advanced_stages.txt"]),
    ("example_stage.exe", "verify_stage_ceb90.dat", [], ["verify_stage_ceb90_stages.txt"]),
    ("example_stage.exe", "verify_stage_curves.dat", [], ["verify_stage_curves_stages.txt"]),
    ("example_stage.exe", "verify_stage_iterate.dat", ["--iterate"],
     ["verify_stage_iterate_stages_iterated.txt"]),
    ("example_stage.exe", "verify_stage_multi_iterate.dat", ["--iterate"],
     ["verify_stage_multi_iterate_stages_iterated.txt"]),
    ("example_stage_envelope.exe", "verify_stage.dat", [], ["verify_stage_stage_envelope.txt"]),
    ("example_3d_beam_live.exe", "verify_3d_beam_live.dat", [],
     ["verify_3d_beam_live_live_envelope.csv"]),
    ("example_3d_influence_surface.exe", "verify_3d_influence_surface.dat", [],
     ["verify_3d_influence_surface_surface.csv"]),
]


def normalize_zeros(text):
    """把数值恰好为 0 的浮点 token 统一替换成 "0"。"""
    return NUMBER_RE.sub(lambda m: "0" if float(m.group()) == 0.0 else m.group(), text)


def extract_numbers(text):
    return [float(m.group()) for m in NUMBER_RE.finditer(text)]


def compare(a, b):
    """返回 (level, detail)。level: identical / normalized / numeric / fail。"""
    if a == b:
        return "identical", ""
    if normalize_zeros(a) == normalize_zeros(b):
        return "normalized", ""
    an = extract_numbers(a)
    bn = extract_numbers(b)
    if len(an) != len(bn):
        return "fail", f"数字个数不一致 base={len(an)} mod={len(bn)}"
    max_rel = 0.0
    first = -1
    for i, (x, y) in enumerate(zip(an, bn)):
        if abs(x) < ZERO and abs(y) < ZERO:
            continue
        d = abs(x - y)
        rel = d / max(abs(x), abs(y), 1e-300)
        if rel > max_rel:
            max_rel = rel
        if d > ATOL + RTOL * max(abs(x), abs(y)):
            if first < 0:
                first = i
    if first < 0:
        return "numeric", f"max_rel={max_rel:.3g}"
    return "fail", f"max_rel={max_rel:.3g} first_idx={first}"


def run_case(exe, dat, args, scratch, bin_dir, src_dir):
    exe_path = bin_dir / exe
    dat_path = src_dir / dat
    if not exe_path.exists():
        return {"case": f"{exe} {dat}", "exit": "EXE缺失", "outputs": []}
    if not dat_path.exists():
        return {"case": f"{exe} {dat}", "exit": "DAT缺失", "outputs": []}
    shutil.copy2(dat_path, scratch / dat)
    env = dict(os.environ)
    for msys in (r"C:\msys64\mingw64\bin", r"C:\msys64\usr\bin"):
        if os.path.isdir(msys):
            env["PATH"] = msys + ";" + env.get("PATH", "")
    proc = subprocess.run([str(exe_path), str(scratch / dat)] + args,
                          cwd=str(scratch), env=env, capture_output=True,
                          text=True, encoding="utf-8", errors="replace", timeout=600)
    return {"case": f"{exe} {dat}", "exit": proc.returncode,
            "stderr": (proc.stderr or "").strip()[-400:], "outputs": []}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="构建目录（含 example_*.exe）")
    ap.add_argument("--src", required=True, help="源码目录（含 verify_*.dat）")
    ap.add_argument("--baseline", required=True, help="基线目录")
    ap.add_argument("--freeze", action="store_true", help="用当前输出刷新基线")
    args = ap.parse_args()

    bin_dir = Path(args.bin)
    src_dir = Path(args.src)
    baseline_dir = Path(args.baseline)
    baseline_dir.mkdir(parents=True, exist_ok=True)

    failed = 0
    total = 0
    with tempfile.TemporaryDirectory(prefix="fem_reg_") as td:
        scratch = Path(td)
        for exe, dat, extra, outs in CASES:
            r = run_case(exe, dat, extra, scratch, bin_dir, src_dir)
            if r["exit"] != 0:
                print(f"[FAIL] {r['case']} 运行失败 exit={r['exit']}")
                if r.get("stderr"):
                    print(f"       {r['stderr']}")
                failed += 1
                continue
            for out in outs:
                total += 1
                produced = scratch / out
                baseline = baseline_dir / out
                if args.freeze:
                    if produced.exists():
                        shutil.copy2(produced, baseline)
                        print(f"[FREEZE] {out}")
                    else:
                        print(f"[FREEZE-MISS] {out}（未生成）")
                    continue
                if not baseline.exists():
                    print(f"[FAIL] {out} 基线缺失（可先 --freeze 生成）")
                    failed += 1
                    continue
                if not produced.exists():
                    print(f"[FAIL] {out} 输出缺失")
                    failed += 1
                    continue
                btxt = baseline.read_text(encoding="utf-8", errors="replace")
                mtxt = produced.read_text(encoding="utf-8", errors="replace")
                level, detail = compare(btxt, mtxt)
                if level == "fail":
                    print(f"[FAIL] {out}  {detail}")
                    failed += 1
                else:
                    print(f"[PASS] {out}  ({level})")

    if args.freeze:
        print(f"\n基线已刷新到: {baseline_dir}")
        return 0
    print(f"\n{'=' * 60}\n结果: {total - failed}/{total} 通过, {failed} 失败")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
