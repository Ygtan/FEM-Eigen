#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""探测 ANSYS 18.0 COMBIN14 与 BEAM3 混用时的 2D 竖向接地弹簧写法。

基准：悬臂梁 L=5m, E=2e11, A=0.02, I=1e-4 -> k_beam = 3EI/L^3 = 4.8e5 N/m
端部接地弹簧 K=1e6 N/m，端部竖向力 P=-10000 N（向下）。
并联总刚度 1.48e6 N/m -> v_tip = -10000/1.48e6 = -0.0067567568 m
（与 FEM-Eigen tests/regression/baselines/verify_elastic_support.result 一致）。

逐个测 KEYOPT(2)=0/1/2/3，找出哪个给出正确的 UY 弹簧行为。
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
V_REF = -10000.0 / (3.0 * 2e11 * 1e-4 / 5.0 ** 3 + 1e6)

TEMPLATE = """/BATCH
/TITLE,COMBIN14 2D vertical spring probe KOP2={kop2}
/PREP7
ET,1,BEAM3
ET,2,COMBIN14
KEYOPT,2,2,{kop2}
R,1,0.02,1e-4,1.0
R,2,1e6
MP,EX,1,2e11
MP,PRXY,1,0.3
! beam node 1 (fixed) and node 2 (free tip)
N,1,0,0,0
N,2,5,0,0
! ground node coincident with tip
N,3,5,0,0
TYPE,1
MAT,1
REAL,1
E,1,2
TYPE,2
REAL,2
E,3,2
D,1,UX,0
D,1,UY,0
D,1,ROTZ,0
D,3,ALL,0
F,2,FY,-10000
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PV,ARRAY,1,2
*GET,PV(1,1),NODE,2,U,Y
*GET,PV(1,2),ELEM,2,SMISC,1
*CFOPEN,sp{kop2}_res,csv
*VWRITE,PV(1,1),PV(1,2)
 (1X,2(1X,E18.10))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""


def main():
    print("theoretical v_tip = %.10g m" % V_REF)
    for kop2 in (0, 1, 2, 3):
        job = "sp%d" % kop2
        inp = HERE / (job + ".inp")
        inp.write_text(TEMPLATE.format(kop2=kop2), encoding="ascii")
        proc = subprocess.run(
            [ANSYS, "-b", "-j", job, "-i", inp.name, "-o", job + ".out"],
            cwd=str(HERE), capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=300)
        res = HERE / (job + "_res.csv")
        out = HERE / (job + ".out")
        vals = None
        if res.exists():
            try:
                vals = [float(v) for v in res.read_text().split()]
            except ValueError:
                pass
        errs = ""
        if out.exists():
            lines = out.read_text(encoding="utf-8", errors="replace").splitlines()
            errs = " | ".join(l.strip() for l in lines
                              if "ERROR" in l.upper() or "***" in l)[:220]
        print("KOP2=%d exit=%d vals=%s  %s" %
              (kop2, proc.returncode, vals, errs))
        if vals and len(vals) >= 1 and abs(vals[0] - V_REF) < 1e-6:
            print("  ==> MATCH: KEYOPT(2)=%d gives the correct vertical spring" % kop2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
