#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibrate BEAM4 SFBEAM LKEY=2 (local y / global Y) sign.

Same simply-supported beam, but lateral (Y) DOFs are released so a lateral
UDL produces a real deflection. Closed form uses IZZ (bending about z).
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 200e9, 7.9e10
A, IYY, IZZ, J = 0.01, 1e-4, 1e-4, 1e-4
L, Q = 4.0, 10000.0

TEMPLATE = r"""/BATCH
/TITLE,lateral SFBEAM calibration {tag}
/PREP7
ET,1,BEAM4
MP,EX,1,{E}
MP,PRXY,1,0.3
MP,GXY,1,{G}
R,1,{A},{IZZ},{IYY},0.02,0.02,0.0
RMORE,0.0,{J},0.0,0.0,0.0,0.0
N,1,0,0,0
N,2,{x2},0,0
N,3,{x3},0,0
TYPE,1
MAT,1
REAL,1
E,1,2
E,2,3
! vertical plane (UZ/ROTY) fully restrained -> no Z displacement
D,1,ALL,0
D,2,UY,0
D,2,UZ,0
D,2,ROTX,0
D,2,ROTY,0
D,3,UY,0
D,3,UZ,0
D,3,ROTX,0
D,3,ROTY,0
! release lateral rotation at 2/3 is automatic (ROTZ free)
{sfbeam}
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,6
*GET,PD(1,1),NODE,2,U,Y
*GET,PD(1,2),NODE,1,ROT,Z
*GET,PD(1,3),NODE,1,RF,FY
*GET,PD(1,4),NODE,3,RF,FY
*GET,PD(1,5),NODE,2,U,Z
*GET,PD(1,6),NODE,1,ROT,Y
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3),PD(1,4),PD(1,5),PD(1,6)
 (1X,6(1X,E20.12))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""


def main():
    uy = -5 * Q * L ** 4 / (384 * E * IZZ)
    rotz = Q * L ** 3 / (24 * E * IZZ)
    print("closed form (q=%.0f, IZZ=%.3g): |UY_mid|=%.9g |ROTZ_n1|=%.9g |RF_FY|=%.9g" %
          (Q, IZZ, abs(uy), abs(rotz), Q * L / 2))
    print()
    print("%-16s %14s %14s %14s %14s %6s" %
          ("job", "UY_mid", "ROTZ_n1", "RF_FY_n1", "RF_FY_n3", "err"))
    for key in (1, 2, 3):
        for sgn, val in (("pos", Q), ("neg", -Q)):
            tag = "lat_k%d_%s" % (key, sgn)
            sf = "\n".join("SFBEAM,%d,%d,PRES,%.12g,%.12g" % (e, key, val, val)
                           for e in (1, 2))
            txt = TEMPLATE.format(tag=tag, E="%.10g" % E, G="%.10g" % G,
                                  A=A, IYY=IYY, IZZ=IZZ, J=J,
                                  x2=L / 2, x3=L, sfbeam=sf)
            (HERE / (tag + ".inp")).write_text(txt, encoding="ascii")
            subprocess.run([ANSYS, "-b", "-j", tag, "-i", tag + ".inp",
                            "-o", tag + ".out"], cwd=str(HERE),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=300)
            p = HERE / (tag + ".csv")
            if not p.exists():
                print("%-16s  (no CSV)" % tag)
                continue
            d = [float(v) for v in p.read_text(encoding="ascii").split()]
            out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace")
            ne = sum(1 for ln in out.splitlines() if "*** ERROR ***" in ln)
            print("%-16s %14.7g %14.7g %14.7g %14.7g %6d" %
                  (tag, d[0], d[1], d[2], d[3], ne))
    return 0


if __name__ == "__main__":
    sys.exit(main())
