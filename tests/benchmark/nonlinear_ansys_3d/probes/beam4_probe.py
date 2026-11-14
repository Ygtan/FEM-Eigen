#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""验证 ANSYS BEAM4（3D Euler-Bernoulli，SHEARZ/SHEARY=0）与 FEM-Eigen 3D 严格一致。

FEM-Eigen 3D cant3d.report.txt 基准：
  UZ_tip = -3.9009715e-05, ROTY_tip = 5.8514572e-06
  elem1 i 端: N=0 Vy=0 Vz=10000 T=0 My=100000 Mz=0
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 2.06e11, 7.9e10
A, IYY, IZZ, J = 0.478, 0.4148, 0.2, 0.3   # IYY 用于 z 向弯曲, IZZ 用于 y 向弯曲
L, PZ, PY, TX = 10.0, -10000.0, -5000.0, 1000.0
EB_Z = PZ * L ** 3 / (3 * E * IYY)
EB_Y = PY * L ** 3 / (3 * E * IZZ)
RX = TX * L / (G * J)

TEMPLATE = r"""/BATCH
/TITLE,BEAM4 strict EB cross-check {tag}
/PREP7
ET,1,BEAM4
MP,EX,1,{E}
MP,PRXY,1,0.3
MP,GXY,1,{G}
R,1,{A},{IZZ},{IYY},2.0,8.0,0.0
RMORE,0.0,{J},0.0,0.0,0.0,0.0
N,1,0,0,0
{nodes}
TYPE,1
MAT,1
REAL,1
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
*DIM,PD,ARRAY,1,7
*GET,PD(1,1),NODE,{tip},U,Z
*GET,PD(1,2),NODE,{tip},U,Y
*GET,PD(1,3),NODE,{tip},ROT,X
*GET,PD(1,4),NODE,{tip},ROT,Y
*GET,PD(1,5),NODE,{tip},ROT,Z
*GET,PD(1,6),NODE,1,RF,MY
*GET,PD(1,7),NODE,1,RF,MZ
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3),PD(1,4),PD(1,5),PD(1,6),PD(1,7)
 (1X,7(1X,E20.12))
*CFCLOS
ETABLE,EN_I,SMISC,1
ETABLE,EMY_I,SMISC,6
ETABLE,EMZ_I,SMISC,5
ETABLE,ESFZ_I,SMISC,2
ETABLE,ESFY_I,SMISC,3
ETABLE,ETQ_I,SMISC,4
ETABLE,EN_J,SMISC,7
ETABLE,EMY_J,SMISC,12
ETABLE,EMZ_J,SMISC,11
ETABLE,ESFZ_J,SMISC,8
ETABLE,ESFY_J,SMISC,9
ETABLE,ETQ_J,SMISC,10
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
PRETAB,EN_I,EMY_I,EMZ_I,ETQ_I,ESFZ_I,ESFY_I
*CFOPEN,{tag}_et,csv
*VWRITE,PE(1,1),PE(1,2),PE(1,3),PE(1,4),PE(1,5),PE(1,6),PE(1,7),PE(1,8),PE(1,9),PE(1,10),PE(1,11),PE(1,12)
 (1X,12(1X,E20.12))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""


def make(tag, n):
    nodes = "\n".join("N,%d,%.12g,0,0" % (i + 1, L * i / n)
                      for i in range(1, n + 1))
    elems = "\n".join("E,%d,%d" % (i, i + 1) for i in range(1, n + 1))
    return TEMPLATE.format(tag=tag, E="%.10g" % E, G="%.10g" % G,
                           A=A, IYY=IYY, IZZ=IZZ, J=J, PZ=PZ, PY=PY, TX=TX,
                           nodes=nodes, elems=elems, tip=n + 1)


def vec(name):
    p = HERE / name
    if not p.exists():
        return None
    try:
        return [float(v) for v in p.read_text(encoding="ascii").split()]
    except ValueError:
        return None


def main():
    print("EB: UZ=%.10g UY=%.10g ROTX=%.10g" % (EB_Z, EB_Y, RX))
    tag = "b4eb10"
    (HERE / (tag + ".inp")).write_text(make(tag, 10), encoding="ascii")
    subprocess.run([ANSYS, "-b", "-j", tag, "-i", tag + ".inp",
                    "-o", tag + ".out"], cwd=str(HERE),
                   capture_output=True, text=True, encoding="utf-8",
                   errors="replace", timeout=300)
    d, et = vec(tag + ".csv"), vec(tag + "_et.csv")
    print("disp:", d)
    if d:
        print("  UZ/EB=%.9f UY/EB=%.9f ROTX/EB=%.9f ROTY=%g ROTZ=%g RF_MY=%g RF_MZ=%g" %
              (d[0] / EB_Z, d[1] / EB_Y, d[2] / RX, d[3], d[4], d[5], d[6]))
    print("etab:", et)
    out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace")
    print("\nPREtab section:")
    for i, l in enumerate(out.splitlines()):
        if "ELEM" in l and "EN_I" in l:
            print("\n".join(out.splitlines()[i:i + 5]))
            break
    print("\nerrors:", [l.strip() for l in out.splitlines()
                        if "*** ERROR ***" in l][:3])
    return 0


if __name__ == "__main__":
    sys.exit(main())
