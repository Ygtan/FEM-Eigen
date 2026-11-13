#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""连续梁 ANSYS 对标 benchmark（CTest / run_tests.ps1 入口）。

默认（离线，CTest 用，不需要 ANSYS 许可）：
  对 reference/*.json 每个算例
    1) 用当前构建的 exe 重算 cases/<dat>，解析 <stem>.result；
    2) 逐个特征量与冻结的 ANSYS 18.0 / BEAM3 实测值比对，
       判据 |got-ref| <= atol + rtol*|ref|，rtol = 1e-3；
    3) 顺带校验归档完整性：.inp 存在、ANSYS CSV 行数与 .dat 节点/单元数一致。
  参考值由 build_reference.py 从 cases/*_ansys_*.csv 抽取，不是 FEM-Eigen 自测值。

可选（--live-ansys，需本机 ANSYS 许可，不进 CTest）：
  在同一临时目录里执行归档的 APDL 命令流，实时把 ANSYS 导出的
  disp/reac/elem CSV 与本次 FEM-Eigen 结果全量比对（位移/反力/弯矩/剪力幅值）。
  这样既有"日常无需许可的回归守护"，也能随时重做端到端商业软件校核。

运行输出一律 ASCII：CTest 日志、PowerShell 控制台、重定向到文件在任何
宿主编码下都可读；中文说明见 README.md 与 reference/*.json。

用法：
    python compare_continuous_beam.py --bin <build> --cases <cases_dir> --reference <ref_dir>
    python compare_continuous_beam.py --bin build --cases ... --reference ... --live-ansys
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
# 输出被重定向到管道/文件时可能抛 UnicodeEncodeError 半路崩掉。统一 UTF-8 +
# errors=replace，且运行输出只用 ASCII，任何宿主都不会乱码。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass
del _s

# 与 tests/benchmark/ansys 共用同一套 .result 解析器，避免两个 benchmark
# 对同一份结果文件的列映射各写一份、日后悄悄分叉。
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ansys"))
from compare_ansys import parse_result, extract  # noqa: E402

DEFAULT_ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
MSYS_DIRS = (r"C:\msys64\mingw64\bin", r"C:\msys64\usr\bin")


def msys_env():
    env = dict(os.environ)
    for msys in MSYS_DIRS:
        if os.path.isdir(msys):
            env["PATH"] = msys + ";" + env.get("PATH", "")
    return env


def read_csv_table(path):
    """读 ANSYS *VWRITE CSV -> {首列 int: 整行 float 列表}。"""
    rows = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        toks = raw.replace(",", " ").split()
        if len(toks) < 2:
            continue
        try:
            key = int(round(float(toks[0])))
            rows[key] = [float(t) for t in toks]
        except ValueError:
            continue
    return rows


def count_dat_entities(text):
    """从 .dat 统计节点数 / 单元数（*NODE 段数字行、*EN, 行）。"""
    nodes = elems = 0
    in_node = False
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        up = line.upper()
        if up.startswith("*EN,"):
            elems += 1
            continue
        if up.startswith("*NODE"):
            in_node = True
            continue
        if line.startswith("*"):
            in_node = False
            continue
        if in_node and line[0].isdigit():
            nodes += 1
    return nodes, elems


def check_archive(case, cases_dir):
    """校验归档文件齐全、ANSYS CSV 行数与 .dat 规模一致。返回错误列表。"""
    errs = []
    stem = case["id"]
    dat_path = cases_dir / case["dat"]
    if not dat_path.exists():
        return ["missing preprocessor output .dat: %s" % case["dat"]]
    if not (cases_dir / case["inp"]).exists():
        errs.append("missing same-source ANSYS APDL: %s" % case["inp"])

    nnode, nelem = count_dat_entities(
        dat_path.read_text(encoding="utf-8", errors="replace"))
    if nnode == 0 or nelem == 0:
        errs.append(".dat has no nodes/elements (got %d nodes, %d elems)" %
                    (nnode, nelem))
        return errs

    expected = {
        "%s_ansys_disp.csv" % stem: nnode,
        "%s_ansys_reac.csv" % stem: nnode,
        "%s_ansys_elem.csv" % stem: nelem,
    }
    for name, want in expected.items():
        path = cases_dir / name
        if not path.exists():
            errs.append("missing frozen ANSYS result: %s" % name)
            continue
        got = len(read_csv_table(path))
        if got != want:
            errs.append("%s has %d rows, .dat expects %d" % (name, got, want))

    if any(c.get("kind") == "spring_force" for c in case.get("checks", [])):
        sp = cases_dir / ("%s_ansys_spring.csv" % stem)
        if not sp.exists():
            errs.append("missing frozen ANSYS spring result: %s" % sp.name)
        elif not read_csv_table(sp):
            errs.append("%s has no rows" % sp.name)
    return errs


def run_fem(exe, dat, bin_dir, cases_dir, scratch):
    exe_path = bin_dir / exe
    dat_src = cases_dir / dat
    if not exe_path.exists():
        return None, "missing exe: %s" % exe_path
    shutil.copy2(dat_src, scratch / dat)
    proc = subprocess.run([str(exe_path), str(scratch / dat)],
                          cwd=str(scratch), env=msys_env(), capture_output=True,
                          text=True, encoding="utf-8", errors="replace", timeout=600)
    if proc.returncode != 0:
        return None, "FEM-Eigen failed exit=%d: %s" % (
            proc.returncode, (proc.stderr or "").strip()[-400:])
    result_path = scratch / (Path(dat).stem + ".result")
    if not result_path.exists():
        return None, "no .result produced"
    return parse_result(result_path.read_text(encoding="utf-8", errors="replace")), None


def compare_checks(case, parsed):
    """与冻结的 ANSYS 值比对，返回 (通过数, 失败数, 输出行列表)。"""
    lines, failed = [], 0
    for chk in case["checks"]:
        name = "%s:%s" % (case["id"], chk["name"])
        if chk["kind"] == "spring_force":
            # FEM-Eigen 把弹性支承当自由自由度，.result 反力段不列弹簧力，
            # 只能由 F_spring = -K*(UY - spring_base) 反算。
            row = parsed["disp"].get(int(chk["node"]))
            u = row[1] if row and row[1] is not None else None
            got = None if u is None else -float(chk["k"]) * (
                u - float(chk.get("value", 0.0)))
        else:
            got = extract(parsed, chk)
        if got is not None and chk.get("abs"):
            got = abs(got)
        if got is None:
            lines.append("[FAIL] %-36s cannot extract from .result (kind=%s)" %
                         (name, chk["kind"]))
            failed += 1
            continue
        ref = chk["ref"]
        rtol = chk.get("rtol", 0.001)
        atol = chk.get("atol", 0.0)
        d = abs(got - ref)
        tol = atol + rtol * abs(ref)
        # ref=0 时用 atol 作分母，输出"占用了多少容差"，避免除零得到 1e290 之类
        denom = max(abs(ref), atol, 1e-300)
        ratio = d / denom
        tag = "PASS" if d <= tol else "FAIL"
        lines.append("[%s] %-36s got=%-15.8g ref=%-15.8g ratio=%-9.3g%s" %
                     (tag, name, got, ref, ratio,
                      "" if d <= tol else " (tol=%g)" % tol))
        if d > tol:
            failed += 1
    return len(case["checks"]) - failed, failed, lines


def live_ansys(ansys_exe, case, cases_dir, parsed, scratch, timeout=600):
    """可选：实时跑归档 APDL，全量比对本轮 FEM-Eigen 结果。返回 (错误列表, 输出行)。"""
    stem = case["id"]
    inp_name = Path(case["inp"]).name
    lines, errs = [], []
    shutil.copy2(cases_dir / case["inp"], scratch / inp_name)
    out_name = stem + "_live_ansys.out"
    cmd = [str(ansys_exe), "-b", "-j", stem, "-i", inp_name, "-o", out_name]
    try:
        proc = subprocess.run(cmd, cwd=str(scratch), env=msys_env(),
                              capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as e:
        return ["ANSYS launch failed/timeout: %s" % e], lines
    if proc.returncode != 0:
        errs.append("ANSYS exit code %d" % proc.returncode)

    csv_names = {
        "disp": "%s_disp.csv" % stem,
        "reac": "%s_reac.csv" % stem,
        "elem": "%s_elem.csv" % stem,
    }
    for kind, name in csv_names.items():
        if not (scratch / name).exists():
            errs.append("ANSYS did not export %s (see %s)" % (name, out_name))
    if errs:
        return errs, lines
    disp = read_csv_table(scratch / csv_names["disp"])
    reac = read_csv_table(scratch / csv_names["reac"])
    elem = read_csv_table(scratch / csv_names["elem"])
    rtol = 1e-3

    # 位移：全部节点 UY，按最大位移量级归一（对 0 值用绝对容差）
    pairs = [(n, row[1], disp[n][2]) for n, row in parsed["disp"].items()
             if n in disp and len(disp[n]) > 2 and row[1] is not None]
    if pairs:
        scale = max(abs(a) for _, a, _ in pairs) or 1.0
        dmax = max(abs(a - b) for _, a, b in pairs)
        ok = dmax <= rtol * scale + 1e-9
        lines.append("  %s LIVE UY (%d nodes): maxdiff %.4g m / scale %.4g m = %.3g  %s" %
                     (stem, len(pairs), dmax, scale, dmax / scale,
                      "PASS" if ok else "FAIL"))
        if not ok:
            errs.append("%s LIVE UY out of tolerance" % stem)

    # 支座反力 FY
    rerr = 0.0
    cnt = 0
    for (node, dof, val) in parsed["reaction"]:
        if dof != 2 or node not in reac or len(reac[node]) <= 2:
            continue
        ref = reac[node][2]
        rerr = max(rerr, abs(val - ref) / max(abs(ref), 1.0))
        cnt += 1
    lines.append("  %s LIVE reaction FY (%d): max rel diff %.3g  %s" %
                 (stem, cnt, rerr, "PASS" if rerr <= rtol else "FAIL"))
    if cnt == 0 or rerr > rtol:
        errs.append("%s LIVE reaction out of tolerance or missing" % stem)

    # 弯矩 M1/M2（符号约定相同）、剪力 V1/V2（幅值）
    mscale = vscale = mdiff = vdiff = 0.0
    mcnt = vcnt = 0
    fem_force = parsed["force"]
    for e, row in elem.items():
        f = fem_force.get(e)
        if f is None:
            continue
        for key, col, is_moment in (("M1", 3, True), ("M2", 6, True),
                                    ("V1", 2, False), ("V2", 5, False)):
            if col >= len(row) or f.get(key) is None:
                continue
            ref = row[col] if is_moment else abs(row[col])
            got = f[key] if is_moment else abs(f[key])
            if is_moment:
                mscale = max(mscale, abs(ref))
                mdiff = max(mdiff, abs(got - ref))
                mcnt += 1
            else:
                vscale = max(vscale, abs(ref))
                vdiff = max(vdiff, abs(got - ref))
                vcnt += 1
    for label, diff, scale, cnt in (("moment", mdiff, mscale, mcnt),
                                    ("shear |V|", vdiff, vscale, vcnt)):
        ratio = diff / (scale or 1.0)
        ok = cnt > 0 and ratio <= rtol
        lines.append("  %s LIVE %s (%d items): maxdiff %.6g / scale %.6g = %.3g  %s" %
                     (stem, label, cnt, diff, scale, ratio, "PASS" if ok else "FAIL"))
        if not ok:
            errs.append("%s LIVE %s out of tolerance" % (stem, label))

    spring_checks = [c for c in case.get("checks", [])
                     if c.get("kind") == "spring_force"]
    if spring_checks:
        sp = scratch / ("%s_spring.csv" % stem)
        if not sp.exists():
            errs.append("%s LIVE: ANSYS did not export spring CSV" % stem)
        else:
            stab = read_csv_table(sp)
            worst = 0.0
            for c in spring_checks:
                nd = int(c["node"])
                row = stab.get(nd)
                fem_row = parsed["disp"].get(nd)
                u = fem_row[1] if fem_row else None
                if row is None or len(row) < 2 or u is None:
                    errs.append("%s LIVE spring node %d missing" % (stem, nd))
                    continue
                f_fem = -float(c["k"]) * (u - float(c.get("value", 0.0)))
                worst = max(worst, abs(f_fem - abs(row[1])) / max(abs(f_fem), 1.0))
            ok = worst <= rtol
            lines.append("  %s LIVE spring force (%d): max rel diff %.3g  %s" %
                         (stem, len(spring_checks), worst, "PASS" if ok else "FAIL"))
            if not ok:
                errs.append("%s LIVE spring force out of tolerance" % stem)
    return errs, lines


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="build dir (contains example_mixed.exe)")
    ap.add_argument("--cases", default=str(here / "cases"), help="case inputs / frozen ANSYS CSV")
    ap.add_argument("--reference", default=str(here / "reference"), help="reference JSON dir")
    ap.add_argument("--live-ansys", nargs="?", const=DEFAULT_ANSYS, default=None,
                    metavar="ANSYS_EXE", help="also run archived APDL live (needs license)")
    ap.add_argument("--ansys-timeout", type=int, default=600)
    args = ap.parse_args()

    bin_dir, cases_dir, ref_dir = Path(args.bin), Path(args.cases), Path(args.reference)
    refs = sorted(ref_dir.glob("*.json"))
    if not refs:
        print("[FAIL] empty reference dir: %s" % ref_dir)
        return 1

    total = failed = 0
    for ref_file in refs:
        case = json.loads(ref_file.read_text(encoding="utf-8"))
        stem = case["id"]
        print("\n== %s ==" % stem)

        errs = check_archive(case, cases_dir)
        if errs:
            for e in errs:
                print("[FAIL] %s archive: %s" % (stem, e))
            failed += len(errs)
            continue
        print("[PASS] %s archive: .dat/.inp/ANSYS CSV present, row counts match" % stem)

        scratch = Path(tempfile.mkdtemp(prefix="fem_cb_"))
        try:
            parsed, err = run_fem(case["exe"], case["dat"], bin_dir, cases_dir, scratch)
            if err:
                print("[FAIL] %s %s" % (stem, err))
                failed += 1
                continue
            passed, bad, lines = compare_checks(case, parsed)
            for line in lines:
                print(line)
            total += len(case["checks"])
            failed += bad

            if args.live_ansys:
                live_errs, live_lines = live_ansys(
                    args.live_ansys, case, cases_dir, parsed, scratch,
                    timeout=args.ansys_timeout)
                for line in live_lines:
                    print(line)
                for e in live_errs:
                    print("[FAIL] %s %s" % (stem, e))
                failed += len(live_errs)
                print("  %s LIVE ANSYS cross-check: %s" %
                      (stem, "PASS" if not live_errs else "FAIL"))
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

    print("\n" + "=" * 68)
    print("Result: %d/%d passed, %d failed" % (total - failed, total, failed))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
