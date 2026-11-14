#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BEAM188 第四轮：用 MP,GXY 消除 Timoshenko 剪切，检验是否精确退化为 EB。

梁仅受 FZ 与 FY 弯曲（无扭矩），此时 G 在 FEM-Eigen 3D 内核中只用于扭转、
不参与求解，因此把 ANSYS 的 GXY 放大不改变可比物理量，却能让 G*As -> inf。
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 2.06e11, 7.9e10
A, IYY, IZZ, J = 0.478, 0.4148, 0.2, 0.3
L, PZ, PY = 10.0, -10000.0, -5000.0
EB_Z = PZ * L ** 3 / (3 * E * IYY)
EB_Y = PY * L ** 3 / (3 * E * IZZ)

TEMPLATE = r"""/BATCH
/TITLE,BEAM188 GXY probe {tag}
/PREP7
ET,1,BEAM188
KEYOPT,1,1,0
MP,EX,1,{E}
MP,PRXY,1,0.3
MP,GXY,1,{gxy}
SECTYPE,1,BEAM,ASEC
SECDATA,{A},{IYY},0,{IZZ},0,{J},0,0,0,0
N,1,0,0,0
{nodes}
TYPE,1
MAT,1
SECNUM,1
{elems}
D,1,ALL,0
F,{tip},FZ,{PZ}
F,{tip},FY,{PY}
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,4
*GET,PD(1,1),NODE,{tip},U,Z
*GET,PD(1,2),NODE,{tip},U,Y
*GET,PD(1,3),NODE,{tip},ROT,Y
*GET,PD(1,4),NODE,{tip},ROT,Z
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3),PD(1,4)
 (1X,4(1X,E20.12))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""


def make(tag, n, gxy):
    nodes = "\n".join("N,%d,%.10g,0,0" % (i + 1, L * i / n)
                      for i in range(1, n + 1))
    elems = "\n".join("E,%d,%d" % (i, i + 1) for i in range(1, n + 1))
    return TEMPLATE.format(tag=tag, E="%.10g" % E, gxy="%.10g" % gxy,
                           A=A, IYY=IYY, IZZ=IZZ, J=J, PZ=PZ, PY=PY,
                           nodes=nodes, elems=elems, tip=n + 1)


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
    errs = [l.strip() for l in out.splitlines() if "*** ERROR ***" in l][:2]
    return vals, errs


def main():
    print("EB: UZ=%.10g  UY=%.10g" % (EB_Z, EB_Y))
    for n in (1, 10, 30):
        for mult in (1, 1e3, 1e6, 1e9, 1e12):
            tag = "b4n%02dm%02d" % (n, len(str(int(mult))) - 1)
            vals, errs = run(tag, make(tag, n, G * mult))
            if vals:
                print("  n=%2d G*%-6.0e UZ/EB=%.9f UY/EB=%.9f  rotY=%-14.8g %s" %
                      (n, mult, vals[0] / EB_Z, vals[1] / EB_Y, vals[2],
                       " ".join(errs)))
            else:
                print("  n=%2d G*%-6.0e NO RESULT %s" % (n, mult, " ".join(errs)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
