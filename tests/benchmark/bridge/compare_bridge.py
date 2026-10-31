#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""真实桥梁 ANSYS 对标基线（CTest 用）。

用当前构建的 example_mixed.exe 重算 bridge_tied_arch_80m.dat，解析 .result，
把节点位移（ux/uy）与单元内力（N、M1、M2）和 ANSYS 结果（冻结在
tests/benchmark/bridge/reference/）比对。

容差：位移相对 0.5%，轴力 5%，弯矩 5%（小量值处相对误差偏大）。

用法：
    python compare_bridge.py --bin <build> --src <source> --reference <reference_dir>
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# 控制台编码兜底：Windows 默认 cp936 里既没有 ✅ 也没有 mm² 的 ²，
# 一旦输出被重定向到管道/文件就抛 UnicodeEncodeError、半路崩掉（退出码 1，
# 看起来像用例失败）。统一成 UTF-8，不可编码的字符替换为 ?，绝不崩。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass
del _s



def parse_result(text):
    disp = {}    # node -> [u, v, theta]
    force = {}   # elem -> [N1, N2, V1, M1, V2, M2]
    section = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("节点位移"):
            section = "disp"
            continue
        if line.startswith("单元内力"):
            section = "force"
            continue
        if line.startswith("单元应力") or line.startswith("支座反力"):
            section = None
            continue
        if section is None:
            continue
        toks = line.split()
        try:
            first = int(toks[0])
        except (ValueError, IndexError):
            continue
        if section == "disp" and len(toks) >= 4:
            disp[first] = [float(toks[1]), float(toks[2]), float(toks[3])]
        elif section == "force" and len(toks) >= 8:
            force[first] = [float(t) for t in toks[2:8]]
    return disp, force


def rel_err(got, ref):
    return abs(got - ref) / max(abs(ref), 1e-9)


def within(got, ref, rtol, atol):
    return abs(got - ref) <= atol + rtol * abs(ref)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--src", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()

    bin_dir = Path(args.bin)
    src_dir = Path(args.src)
    refs = sorted(Path(args.reference).glob("*.json"))
    if not refs:
        print("[FAIL] 参考目录为空: %s" % args.reference)
        return 1

    failed = total = 0
    for ref_file in refs:
        case = json.loads(ref_file.read_text(encoding="utf-8"))
        tol = case.get("tolerances", {})
        d_tol = tol.get("disp_rtol", 0.005)
        d_atol = tol.get("disp_atol", 1e-5)
        f_tol = tol.get("force_rtol", 0.05)
        m_tol = tol.get("moment_rtol", 0.05)

        dat = src_dir / case["dat"]
        exe = bin_dir / case["exe"]
        if not exe.exists():
            print("[FAIL] %s  EXE缺失" % case["id"])
            failed += 1
            continue
        scratch = Path(tempfile.mkdtemp(prefix="fem_bridge_"))
        try:
            shutil.copy2(dat, scratch / case["dat"])
            env = dict(os.environ)
            for msys in (r"C:\msys64\mingw64\bin", r"C:\msys64\usr\bin"):
                if os.path.isdir(msys):
                    env["PATH"] = msys + ";" + env.get("PATH", "")
            proc = subprocess.run([str(exe), str(scratch / case["dat"])],
                                  cwd=str(scratch), env=env, capture_output=True,
                                  text=True, encoding="utf-8", errors="replace", timeout=600)
            if proc.returncode != 0:
                print("[FAIL] %s  运行失败 exit=%d" % (case["id"], proc.returncode))
                failed += 1
                continue
            res = scratch / (Path(case["dat"]).stem + ".result")
            disp, force = parse_result(res.read_text(encoding="utf-8", errors="replace"))

            for nid, ref in case["nodes"].items():
                nid_i = int(nid)
                row = disp.get(nid_i)
                for key, ridx in (("ux", 0), ("uy", 1)):
                    got = row[ridx] if row else None
                    total += 1
                    if got is None:
                        print("[FAIL] %s:N%s.%s 缺失" % (case["id"], nid, key))
                        failed += 1
                    elif not within(got, ref[key], d_tol, d_atol):
                        print("[FAIL] %s:N%s.%s got=%.6g ref=%.6g" % (case["id"], nid, key, got, ref[key]))
                        failed += 1
                    else:
                        print("[PASS] %s:N%s.%s %.6g" % (case["id"], nid, key, got))

            for eid, ref in case["elements"].items():
                eid_i = int(eid)
                row = force.get(eid_i)
                for key, ridx in (("N", 0), ("M1", 3), ("M2", 5)):
                    if key not in ref:
                        continue
                    got = row[ridx] if row else None
                    total += 1
                    tt = m_tol if key.startswith("M") else f_tol
                    if got is None:
                        print("[FAIL] %s:E%s.%s 缺失" % (case["id"], eid, key))
                        failed += 1
                    elif not within(got, ref[key], tt, 0.0):
                        print("[FAIL] %s:E%s.%s got=%.6g ref=%.6g" % (case["id"], eid, key, got, ref[key]))
                        failed += 1
                    else:
                        print("[PASS] %s:E%s.%s %.6g" % (case["id"], eid, key, got))
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

    print("\n" + "=" * 60)
    print("结果: %d/%d 通过, %d 失败" % (total - failed, total, failed))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
