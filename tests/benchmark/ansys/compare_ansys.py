#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ANSYS 对标基线比对（CTest 用）。

对 tests/benchmark/ansys/reference/ 下每个 JSON 参考算例：
  1) 用当前构建的 exe 重算 .dat；
  2) 解析 .result；
  3) 抽取参考 JSON 中列出的特征量（跨中挠度/支座反力/跨中弯矩/跨中应力）；
  4) 与黄金参考值做相对容差比对：|got-ref| <= atol + rtol*|ref|。

说明：
  - 当前仓库没有 ANSYS 数值结果文件，参考值来自 ANSYS BEAM4 对应的
    欧拉-伯努利闭式解（BEAM4 默认不计剪切变形，与该闭式解一致）。
  - 若要接入真实 ANSYS 结果，只需把 JSON 中的 ref/rtol 改成 ANSYS 输出值，
    解析器无需改动。

用法：
    python compare_ansys.py --bin <build> --src <source> --reference <reference_dir>
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

FORCE_COLS = ["N1", "N2", "V1", "M1", "V2", "M2"]
STRESS_COLS = ["s1", "s2", "stop1", "sbot1", "stop2", "sbot2"]


def parse_number(token):
    try:
        return float(token)
    except ValueError:
        return None


def parse_result(text):
    """解析 .result 文本，返回 {disp, force, stress, reaction}。"""
    disp = {}      # node -> [u, v, theta]
    force = {}     # elem -> {col: value}
    stress = {}    # elem -> {col: value}
    reaction = []  # [(node, dof, value)]
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
        if line.startswith("单元应力"):
            section = "stress"
            continue
        if line.startswith("支座反力"):
            section = "reaction"
            continue
        if section is None:
            continue
        toks = line.split()
        if not toks:
            continue
        try:
            first = int(toks[0])
        except ValueError:
            continue  # 表头行
        if section == "disp":
            if len(toks) >= 4:
                disp[first] = [parse_number(toks[1]), parse_number(toks[2]), parse_number(toks[3])]
        elif section == "force":
            if len(toks) >= 8:
                force[first] = dict(zip(FORCE_COLS, [parse_number(t) for t in toks[2:8]]))
        elif section == "stress":
            if len(toks) >= 8:
                stress[first] = dict(zip(STRESS_COLS, [parse_number(t) for t in toks[2:8]]))
        elif section == "reaction":
            if len(toks) >= 3:
                reaction.append((first, int(float(toks[1])), parse_number(toks[2])))
    return {"disp": disp, "force": force, "stress": stress, "reaction": reaction}


def extract(parsed, check):
    kind = check["kind"]
    if kind == "nodal_disp":
        node = int(check["node"])
        dof = int(check["dof"])
        row = parsed["disp"].get(node)
        if row is None or row[dof - 1] is None:
            return None
        return row[dof - 1]
    if kind == "reaction":
        node = int(check["node"])
        dof = int(check["dof"])
        for (n, d, v) in parsed["reaction"]:
            if n == node and d == dof:
                return v
        return None
    if kind == "force":
        row = parsed["force"].get(int(check["elem"]))
        return None if row is None else row.get(check["col"])
    if kind == "stress":
        row = parsed["stress"].get(int(check["elem"]))
        return None if row is None else row.get(check["col"])
    raise ValueError("未知 kind: " + kind)


def run_case(exe, dat, bin_dir, src_dir):
    exe_path = bin_dir / exe
    dat_path = src_dir / dat
    if not exe_path.exists():
        return None, "EXE缺失: " + exe
    if not dat_path.exists():
        return None, "DAT缺失: " + dat
    scratch = Path(tempfile.mkdtemp(prefix="fem_ansys_"))
    try:
        shutil.copy2(dat_path, scratch / dat)
        env = dict(os.environ)
        for msys in (r"C:\msys64\mingw64\bin", r"C:\msys64\usr\bin"):
            if os.path.isdir(msys):
                env["PATH"] = msys + ";" + env.get("PATH", "")
        proc = subprocess.run([str(exe_path), str(scratch / dat)],
                              cwd=str(scratch), env=env, capture_output=True,
                              text=True, encoding="utf-8", errors="replace", timeout=600)
        if proc.returncode != 0:
            return None, "运行失败 exit=%d: %s" % (proc.returncode, (proc.stderr or "").strip()[-400:])
        result_path = scratch / (Path(dat).stem + ".result")
        if not result_path.exists():
            return None, "未生成 .result"
        return parse_result(result_path.read_text(encoding="utf-8", errors="replace")), None
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="构建目录（含 example_*.exe）")
    ap.add_argument("--src", required=True, help="源码目录（含 .dat）")
    ap.add_argument("--reference", required=True, help="参考 JSON 目录")
    args = ap.parse_args()

    bin_dir = Path(args.bin)
    src_dir = Path(args.src)
    ref_dir = Path(args.reference)
    refs = sorted(ref_dir.glob("*.json"))
    if not refs:
        print("[FAIL] 参考目录为空: %s" % ref_dir)
        return 1

    failed = 0
    total = 0
    for ref_file in refs:
        case = json.loads(ref_file.read_text(encoding="utf-8"))
        parsed, err = run_case(case["exe"], case["dat"], bin_dir, src_dir)
        if err:
            print("[FAIL] %s  %s" % (case["id"], err))
            failed += 1
            continue
        for check in case["checks"]:
            total += 1
            got = extract(parsed, check)
            name = "%s:%s" % (case["id"], check["name"])
            if got is None:
                print("[FAIL] %s  未能从结果中提取（kind=%s）" % (name, check["kind"]))
                failed += 1
                continue
            ref = check["ref"]
            rtol = check.get("rtol", 0.001)
            atol = check.get("atol", 0.0)
            d = abs(got - ref)
            tol = atol + rtol * abs(ref)
            rel = d / max(abs(ref), 1e-300)
            if d <= tol:
                print("[PASS] %s  got=%.8g ref=%.8g rel=%.3g" % (name, got, ref, rel))
            else:
                print("[FAIL] %s  got=%.8g ref=%.8g rel=%.3g (rtol=%g)" % (name, got, ref, rel, rtol))
                failed += 1

    print("\n" + "=" * 60 + "\n结果: %d/%d 通过, %d 失败" % (total - failed, total, failed))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
