#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BEAM188 第二轮：网格细化收敛性 + SECCONTROL 剪切开关 + PREsOL 内力标签。

悬臂梁 L=10m 沿 +X，仅 FZ=-10000。
Euler-Bernoulli: UZ = -PL^3/(3EIyy) = -3.90097e-5, ROTY = -PL^2/(2EIyy)
Timoshenko(As=A): UZ = EB + PL/(G*A) = -4.1660e-5
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 2.06e11, 7.9e10
A, IYY, IZZ, J = 0.478, 0.4148, 0.2, 0.3
L, PZ = 10.0, -10000.0
EB = PZ * L ** 3 / (3 * E * IYY)
TM = EB + PZ * L / (G * A)

TEMPLATE = r"""/BATCH
/TITLE,BEAM188 mesh/shear probe {tag}
/PREP7
ET,1,BEAM188
KEYOPT,1,1,0
MP,EX,1,{E}
MP,PRXY,1,0.3
MP,GXY,1,{G}
SECTYPE,1,BEAM,ASEC
SECDATA,{A},{IYY},0,{IZZ},0,{J},0,0,0,0
{seccontrol}
SLIST,1
N,1,0,0,0
{nodes}
TYPE,1
MAT,1
SECNUM,1
{elems}
D,1,ALL,0
F,{tip},FZ,{PZ}
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,3
*GET,PD(1,1),NODE,{tip},U,Z
*GET,PD(1,2),NODE,{tip},ROT,Y
*GET,PD(1,3),NODE,1,RF,MY
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3)
 (1X,3(1X,E20.12))
*CFCLOS
PRESOL,F
PRESOL,M
FINISH
/EXIT,NOSAVE
"""

# 试各种 SECCONTROL 写法，观察 SLIST 里 Shear Correction 是否变化
SECCONTROLS = {
    "none": "",
    "c1": "SECCONTROL,0,0",
    "c2": "SECCONTROL,,,0,0",
    "c3": "SECCONTROL,1e20,1e20",
    "c4": "SECCONTROL,,,1e20,1e20",
    "c5": "SECCONTROL,0.001,0.001",
}


def make(tag, n, sec):
    nodes = "\n".join("N,%d,%.6g,0,0" % (i + 1, L * i / n)
                      for i in range(1, n + 1))
    elems = "\n".join("E,%d,%d" % (i, i + 1) for i in range(1, n + 1))
    return TEMPLATE.format(tag=tag, E="%.10g" % E, G="%.10g" % G,
                           A=A, IYY=IYY, IZZ=IZZ, J=J, PZ=PZ,
                           nodes=nodes, elems=elems, tip=n + 1,
                           seccontrol=sec)


def run(tag, text):
    (HERE / (tag + ".inp")).write_text(text, encoding="ascii")
    proc = subprocess.run([ANSYS, "-b", "-j", tag, "-i", tag + ".inp",
                           "-o", tag + ".out"], cwd=str(HERE),
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=300)
    p = HERE / (tag + ".csv")
    vals = None
    if p.exists():
        try:
            vals = [float(v) for v in p.read_text(encoding="ascii").split()]
        except ValueError:
            pass
    out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace") \
        if (HERE / (tag + ".out")).exists() else ""
    shear = [l.strip() for l in out.splitlines() if "Shear Correction-yy" in l]
    errs = [l.strip() for l in out.splitlines() if "*** ERROR ***" in l][:2]
    return vals, (shear[0] if shear else ""), errs, proc.returncode


def main():
    print("EB UZ = %.9g   Timoshenko(As=A) UZ = %.9g" % (EB, TM))
    print("\n--- mesh refinement (no SECCONTROL) ---")
    for n in (1, 2, 5, 10, 20, 40):
        tag = "b2m%02d" % n
        vals, shear, errs, rc = run(tag, make(tag, n, ""))
        if vals:
            print("  n=%2d UZ=%-16.9g UZ/EB=%.6f ROTY=%-14.9g RF_MY=%-14.9g rc=%d %s" %
                  (n, vals[0], vals[0] / EB, vals[1], vals[2], rc,
                   " ".join(errs)))
        else:
            print("  n=%2d NO RESULT rc=%d %s" % (n, rc, " ".join(errs)))
    print("\n--- SECCONTROL variants (n=10) ---")
    for key, sec in SECCONTROLS.items():
        tag = "b2s_" + key
        vals, shear, errs, rc = run(tag, make(tag, 10, sec))
        print("  %-5s %-28s UZ/EB=%s  %s  %s" %
              (key, sec or "(none)",
               ("%.6f" % (vals[0] / EB)) if vals else "NO RESULT",
               shear, " ".join(errs)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
