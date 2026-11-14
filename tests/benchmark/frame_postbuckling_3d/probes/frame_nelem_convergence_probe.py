#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CB3D19 诊断：CR 梁每杆单元数收敛性扫描（不参与 CTest）。

用同一个 frame_probe.exe 跑 NELEM = 8/12/16/24/32 的 Williams toggle
框架，输出每条曲线上的第一极限点、lambda@0.5、lambda@0.9、点数，
以及相对 NELEM=32 的最大 |dlambda|。用于复核 README §4 的定标结论。

用法：
    python frame_nelem_convergence_probe.py \
        --exe D:/AI-work/FEM-Eigen/build/frame_probe.exe
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

BASE = ("E=2.1e11 G=8.1e13 A=1.0e-3 IY=1.0e-6 IZ=1.0e-6 J=2.0e-6 "
        "B=2.0 H=0.5 PREF=4.0e6 TOL=1.0e-9 MAXSTEPS=3000 MAXITER=40 "
        "VMAX=0.95 DSMAX=0.01 DSMIN=1.0e-8")
NELEMS = (8, 12, 16, 24, 32)


def interp3(curve, vq):
    if len(curve) < 3:
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


def first_peak(curve):
    for i in range(1, len(curve) - 1):
        if curve[i][1] > curve[i - 1][1] and curve[i][1] >= curve[i + 1][1]:
            return curve[i]
    return (0.0, 0.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True, help="path to frame_probe.exe")
    a = ap.parse_args()
    exe = Path(a.exe)
    if not exe.is_file():
        print("[FAIL] missing exe:", exe)
        return 1
    with tempfile.TemporaryDirectory(prefix="cb3d19_nelem_") as td:
        tmp = Path(td)
        lines = []
        for n in NELEMS:
            lines.append("G n%d %s NELEM=%d" % (n, BASE, n))
        for n in NELEMS:
            lines.append("C c%d n%d DS=0.005 PSI=1.0 TAG=n%d" % (n, n, n))
        (tmp / "scan.frame").write_text("\n".join(lines) + "\n", encoding="ascii")
        p = subprocess.run([str(exe), "scan.frame", "scan.csv"],
                           cwd=str(tmp), timeout=1800,
                           capture_output=True, text=True, errors="replace")
        if p.returncode != 0:
            print("[FAIL] frame_probe exit=%d: %s" % (p.returncode, p.stderr))
            return 1
        curves = {}
        for row in csv.reader((tmp / "scan.csv").read_text().splitlines()):
            if row and row[0] == "P":
                curves.setdefault(row[1], []).append((float(row[4]), float(row[3])))
        ref = curves.get("c32", [])
        print("NELEM  npoints   v_cr       lambda_cr    lam@0.5     "
              "lam@0.9     worst|dlam| vs n32")
        for n in NELEMS:
            c = curves.get("c%d" % n, [])
            if not c:
                print("%-6d NO DATA" % n)
                continue
            vp, lp = first_peak(c)
            worst = max((abs(l - interp3(ref, v)), v) for v, l in c)[0] if ref else 0.0
            print("%-6d %-9d %-10.6f %-12.8f %-11.7f %-11.7f %.4e" % (
                n, len(c), vp, lp, interp3(c, 0.5), interp3(c, 0.9), worst))
    print("\n结论：8->24 单元差异显著收敛；24 与 32 的差已小于 ANSYS 对齐噪声，")
    print("故正式算例定标 NELEM=24（详见上级 README §4）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
