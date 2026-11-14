#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibrate ANSYS BEAM4 SFBEAM LKEY/sign and ACEL selfweight for 3D Z loads.

Simply supported 3D beam, X longitudinal, Z vertical.
FEM-Eigen 3D uses local z = global Z for horizontal X beams (IY -> IYY).
Closed form (q downward, L=4 m, E*IYY=2e7):
    mid UZ   = -5 q L^4 / (384 E I)
    end ROTY =  q L^3 / (24 E I)   (sign to be read)
    support RF_FZ = q L / 2
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 200e9, 7.9e10
A, IYY, IZZ, J = 0.01, 1e-4, 1e-4, 1e-4
L, Q = 4.0, 10000.0
RHO, GACC = 7850.0, 9.81
QSW = RHO * A * GACC          # 770.085 N/m

TEMPLATE = r"""/BATCH
/TITLE,SFBEAM/selfweight calibration {tag}
/PREP7
ET,1,BEAM4
MP,EX,1,{E}
MP,PRXY,1,0.3
MP,GXY,1,{G}
{dens}
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
D,1,UX,0
D,1,UY,0
D,1,UZ,0
D,1,ROTX,0
D,1,ROTZ,0
D,2,UY,0
D,2,ROTX,0
D,2,ROTZ,0
D,3,UY,0
D,3,UZ,0
D,3,ROTX,0
D,3,ROTZ,0
{sfbeam}
FINISH
/SOLU
ANTYPE,STATIC
{acel}
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,6
*GET,PD(1,1),NODE,2,U,Z
*GET,PD(1,2),NODE,1,ROT,Y
*GET,PD(1,3),NODE,1,RF,FZ
*GET,PD(1,4),NODE,3,RF,FZ
*GET,PD(1,5),NODE,2,U,Y
*GET,PD(1,6),NODE,2,ROT,Z
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3),PD(1,4),PD(1,5),PD(1,6)
 (1X,6(1X,E20.12))
*CFCLOS
ETABLE,ESFZ_I,SMISC,2
ETABLE,EMY_I,SMISC,6
ETABLE,ESFY_I,SMISC,3
ETABLE,EMZ_I,SMISC,5
*DIM,PE,ARRAY,1,4
*GET,PE(1,1),ELEM,1,ETAB,ESFZ_I
*GET,PE(1,2),ELEM,1,ETAB,EMY_I
*GET,PE(1,3),ELEM,1,ETAB,ESFY_I
*GET,PE(1,4),ELEM,1,ETAB,EMZ_I
*CFOPEN,{tag}_et,csv
*VWRITE,PE(1,1),PE(1,2),PE(1,3),PE(1,4)
 (1X,4(1X,E20.12))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""


def make(tag, sfbeam, dens, acel):
    return TEMPLATE.format(tag=tag, E="%.10g" % E, G="%.10g" % G,
                           A=A, IYY=IYY, IZZ=IZZ, J=J,
                           x2=L / 2, x3=L, dens=dens,
                           sfbeam=sfbeam, acel=acel)


def vec(name):
    p = HERE / name
    if not p.exists():
        return None
    try:
        return [float(v) for v in p.read_text(encoding="ascii").split()]
    except ValueError:
        return None


def run(tag, sfbeam, dens="", acel=""):
    (HERE / (tag + ".inp")).write_text(make(tag, sfbeam, dens, acel), encoding="ascii")
    subprocess.run([ANSYS, "-b", "-j", tag, "-i", tag + ".inp", "-o", tag + ".out"],
                   cwd=str(HERE), capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=300)
    out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace")
    nerr = sum(1 for ln in out.splitlines() if "*** ERROR ***" in ln)
    nwarn = sum(1 for ln in out.splitlines() if "*** WARNING ***" in ln)
    return vec(tag + ".csv"), vec(tag + "_et.csv"), nerr, nwarn


def reference(q):
    return (-5 * q * L ** 4 / (384 * E * IYY),
            q * L ** 3 / (24 * E * IYY),
            q * L / 2)


def main():
    jobs = []
    for key in (1, 2, 3):
        for sgn, val in (("pos", Q), ("neg", -Q)):
            tag = "sw_k%d_%s" % (key, sgn)
            sf = "\n".join("SFBEAM,%d,%d,PRES,%.12g,%.12g" % (e, key, val, val)
                           for e in (1, 2))
            jobs.append((tag, sf, "", ""))
    # selfweight via ACEL + MP,DENS only
    jobs.append(("sw_acel", "", "MP,DENS,1,%.12g" % RHO,
                 "ACEL,0.0,0.0,%.12g" % GACC))
    # same magnitude via SFBEAM, using each candidate key
    for key in (1, 2, 3):
        tag = "sw_ref_k%d" % key
        sf = "\n".join("SFBEAM,%d,%d,PRES,%.12g,%.12g" % (e, key, QSW, QSW)
                       for e in (1, 2))
        jobs.append((tag, sf, "", ""))

    print("closed form q=%.3f: UZ=%.9g ROTY=%.9g RF=%.9g" % ((Q,) + reference(Q)))
    print("selfweight line load q=%.6f N/m" % QSW)
    print("closed form q=QSW : UZ=%.9g ROTY=%.9g RF=%.9g" % reference(QSW))
    print()
    hdr = "%-14s %14s %14s %14s %14s %14s %6s %6s" % (
        "job", "UZ_mid", "ROTY_n1", "RF_FZ_n1", "RF_FZ_n3", "UY_mid", "err", "warn")
    print(hdr)
    print("-" * len(hdr))
    for tag, sf, dens, acel in jobs:
        d, et, ne, nw = run(tag, sf, dens, acel)
        if not d:
            print("%-14s  (no CSV)" % tag)
            continue
        print("%-14s %14.7g %14.7g %14.7g %14.7g %14.7g %6d %6d" %
              (tag, d[0], d[1], d[2], d[3], d[4], ne, nw))
        if et:
            print("%-14s   elem1 i: SFZ=%g MY=%g SFY=%g MZ=%g" %
                  ("", et[0], et[1], et[2], et[3]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
