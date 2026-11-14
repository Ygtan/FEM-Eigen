#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BEAM188 最后一轮：确定正式算例配置。

A. SECCONTROL 量级 x 网格：找与 Euler-Bernoulli 差异 < 1e-3 的稳定配置
B. 同时加扭矩 MX，验证 SECCONTROL 不影响 G*J 扭转刚度
C. 用 ETABLE + *GET(...,ETAB,...) 标定 i/j 端内力（对比直接 SMISC）
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 2.06e11, 7.9e10
A, IYY, IZZ, J = 0.478, 0.4148, 0.2, 0.3
L, PZ, PY, TX = 10.0, -10000.0, -5000.0, 1000.0
EB_Z = PZ * L ** 3 / (3 * E * IYY)
EB_Y = PY * L ** 3 / (3 * E * IZZ)
RX = TX * L / (G * J)

TEMPLATE = r"""/BATCH
/TITLE,BEAM188 final probe {tag}
/PREP7
ET,1,BEAM188
KEYOPT,1,1,0
MP,EX,1,{E}
MP,PRXY,1,0.3
MP,GXY,1,{G}
SECTYPE,1,BEAM,ASEC
SECDATA,{A},{IYY},0,{IZZ},0,{J},0,0,0,0
{seccontrol}
N,1,0,0,0
{nodes}
TYPE,1
MAT,1
SECNUM,1
{elems}
D,1,ALL,0
F,{tip},FZ,{PZ}
F,{tip},FY,{PY}
F,{tip},MX,{TX}
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,6
*GET,PD(1,1),NODE,{tip},U,Z
*GET,PD(1,2),NODE,{tip},U,Y
*GET,PD(1,3),NODE,{tip},ROT,X
*GET,PD(1,4),NODE,{tip},ROT,Y
*GET,PD(1,5),NODE,{tip},ROT,Z
*GET,PD(1,6),NODE,1,RF,MY
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3),PD(1,4),PD(1,5),PD(1,6)
 (1X,6(1X,E20.12))
*CFCLOS
ETABLE,EN_I,SMISC,1
ETABLE,EMY_I,SMISC,2
ETABLE,EMZ_I,SMISC,3
ETABLE,ETQ_I,SMISC,4
ETABLE,ESFZ_I,SMISC,5
ETABLE,ESFY_I,SMISC,6
ETABLE,EN_J,SMISC,14
ETABLE,EMY_J,SMISC,15
ETABLE,EMZ_J,SMISC,16
ETABLE,ETQ_J,SMISC,17
ETABLE,ESFZ_J,SMISC,18
ETABLE,ESFY_J,SMISC,19
*DIM,PE,ARRAY,1,12
*GET,PE(1,1),ELEM,1,ETAB,EN_I
*GET,PE(1,2),ELEM,1,ETAB,EMY_I
*GET,PE(1,3),ELEM,1,ETAB,EMZ_I
*GET,PE(1,4),ELEM,1,ETAB,ETQ_I
*GET,PE(1,5),ELEM,1,ETAB,ESFZ_I
*GET,PE(1,6),ELEM,1,ETAB,ESFY_I
*GET,PE(1,7),ELEM,1,ETAB,EN_J
*GET,PE(1,8),ELEM,1,ETAB,EMY_J
*GET,PE(1,9),ELEM,1,ETAB,EMZ_J
*GET,PE(1,10),ELEM,1,ETAB,ETQ_J
*GET,PE(1,11),ELEM,1,ETAB,ESFZ_J
*GET,PE(1,12),ELEM,1,ETAB,ESFY_J
*CFOPEN,{tag}_et,csv
*VWRITE,PE(1,1),PE(1,2),PE(1,3),PE(1,4),PE(1,5),PE(1,6),PE(1,7),PE(1,8),PE(1,9),PE(1,10),PE(1,11),PE(1,12)
 (1X,12(1X,E20.12))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""


def make(tag, n, sec):
    nodes = "\n".join("N,%d,%.12g,0,0" % (i + 1, L * i / n)
                      for i in range(1, n + 1))
    elems = "\n".join("E,%d,%d" % (i, i + 1) for i in range(1, n + 1))
    return TEMPLATE.format(tag=tag, E="%.10g" % E, G="%.10g" % G,
                           A=A, IYY=IYY, IZZ=IZZ, J=J, PZ=PZ, PY=PY, TX=TX,
                           seccontrol=sec, nodes=nodes, elems=elems, tip=n + 1)


def run(tag, text):
    (HERE / (tag + ".inp")).write_text(text, encoding="ascii")
    subprocess.run([ANSYS, "-b", "-j", tag, "-i", tag + ".inp",
                    "-o", tag + ".out"], cwd=str(HERE),
                   capture_output=True, text=True, encoding="utf-8",
                   errors="replace", timeout=300)
    def vec(name):
        p = HERE / name
        if not p.exists():
            return None
        try:
            return [float(v) for v in p.read_text(encoding="ascii").split()]
        except ValueError:
            return None
    return vec(tag + ".csv"), vec(tag + "_et.csv")


def main():
    print("EB: UZ=%.10g UY=%.10g ROTX=%.10g" % (EB_Z, EB_Y, RX))
    print("\n--- SECCONTROL x mesh ---")
    for val in ("1e12", "1e16", "1e20"):
        for n in (10, 30, 60):
            tag = "b5v%02dn%02d" % (int(val[-2:]), n)
            d, et = run(tag, make(tag, n, "SECCONTROL,%s,%s" % (val, val)))
            if d:
                print("  sec=%-5s n=%2d UZ/EB=%.9f UY/EB=%.9f ROTX/EB=%.9f" %
                      (val, n, d[0] / EB_Z, d[1] / EB_Y, d[2] / RX))
            else:
                print("  sec=%-5s n=%2d NO RESULT" % (val, n))
    tag = "b5def"
    d, et = run(tag, make(tag, 30, ""))
    if d:
        print("  default(Timoshenko) n=30 UZ/EB=%.9f UY/EB=%.9f ROTX/EB=%.9f" %
              (d[0] / EB_Z, d[1] / EB_Y, d[2] / RX))
    print("\n--- ETABLE end forces (elem1, n=30, SECCONTROL 1e20) ---")
    tag = "b5v20n30"
    d, et = run(tag, make(tag, 30, "SECCONTROL,1e20,1e20"))
    if et:
        labels = ["N_I", "My_I", "Mz_I", "T_I", "SFz_I", "SFy_I",
                  "N_J", "My_J", "Mz_J", "T_J", "SFz_J", "SFy_J"]
        for lab, v in zip(labels, et):
            print("   %-6s = %-16.8g" % (lab, v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
