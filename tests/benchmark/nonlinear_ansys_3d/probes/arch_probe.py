#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Probe: inclined 3D BEAM (X-Z plane) + vertical TRUSS hanger.

FEM-Eigen beam3DTransformation: ez_ref = global Z (or global Y if ex || Z),
ey = ez_ref x ex, ez = ex x ey.  Check whether ANSYS BEAM4's default
orientation agrees for inclined members (arch rib / inclined tower leg).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, r"D:\AI-work\BridgeStudio")
from model_spec_3d import (  # noqa: E402
    Boundary3D, Element3D, ElemTypeDef3D, Material3D, ModelSpec3D,
    Node3D, NodalLoad3D, RealConst3D,
)

FEM_EXE = Path(r"D:\AI-work\FEM-Eigen\build\example_3d_report.exe")
ANSYS_EXE = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
MSYS_BIN = r"C:\msys64\mingw64\bin"

E, A, IY, IZ, JTOR, G = 2.06e11, 0.02, 1e-3, 4e-4, 2e-3, 7.9e10
A_TRUSS = 0.001


def n(x):
    return "{:.12g}".format(float(x))


def main():
    # node 1 (0,0,0) fixed; node 2 (5,0,4) arch crown-ish; node 3 (10,0,0)
    # node 4 (5,0,0) deck point under node 2 -> vertical hanger TRUSS
    nodes = [Node3D(id=1, x=0.0, y=0.0, z=0.0),
             Node3D(id=2, x=5.0, y=0.0, z=4.0),
             Node3D(id=3, x=10.0, y=0.0, z=0.0),
             Node3D(id=4, x=5.0, y=0.0, z=0.0)]
    els = [Element3D(id=1, etype=1, n1=1, n2=2),      # inclined BEAM
           Element3D(id=2, etype=1, n1=2, n2=3),      # inclined BEAM
           Element3D(id=3, etype=1, n1=1, n2=4),       # horizontal BEAM (deck)
           Element3D(id=4, etype=1, n1=4, n2=3),       # horizontal BEAM (deck)
           Element3D(id=5, etype=2, n1=2, n2=4)]        # vertical TRUSS hanger
    model = ModelSpec3D(
        title="inclined beam probe", nodes=nodes,
        materials=[Material3D(id=1, E=E)],
        reals=[RealConst3D(id=1, A=A, IY=IY, IZ=IZ, J=JTOR, G=G),
               RealConst3D(id=2, A=A_TRUSS)],
        etypes=[ElemTypeDef3D(id=1, kind="BEAM", mat=1, real=1),
                ElemTypeDef3D(id=2, kind="TRUSS", mat=1, real=2)],
        elements=els,
        boundaries=[Boundary3D(node=1, dofs=(1, 6), value=0.0),
                    Boundary3D(node=3, dofs=(2, 3), value=0.0)],
        loads=[NodalLoad3D(node=4, dof=3, value=-100000.0)],
    )
    (HERE / "archp.dat").write_text(model.to_dat(compact=False), encoding="utf-8")
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    p = subprocess.run([str(FEM_EXE), "archp.dat"], cwd=str(HERE), env=env,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    print("FEM rc=", p.returncode)
    print((HERE / "archp.report.txt").read_text(encoding="utf-8", errors="replace"))

    L = []
    a = L.append
    a("/BATCH")
    a("/TITLE,inclined beam probe")
    a("/PREP7")
    a("ET,1,BEAM4")
    a("ET,2,LINK180")
    a("R,1,{},{},{},1,1,0.0".format(n(A), n(IZ), n(IY)))
    a("RMORE,0.0,{},0.0,0.0,0.0,0.0".format(n(JTOR)))
    a("R,2,{},0.0".format(n(A_TRUSS)))
    a("MP,EX,1,{}".format(n(E)))
    a("MP,PRXY,1,0.3")
    a("MP,GXY,1,{}".format(n(G)))
    for nd in nodes:
        a("N,{},{},{},{}".format(nd.id, n(nd.x), n(nd.y), n(nd.z)))
    a("TYPE,1")
    a("MAT,1")
    a("REAL,1")
    for e in els[:4]:
        a("E,{},{}".format(e.n1, e.n2))
    a("TYPE,2")
    a("MAT,1")
    a("REAL,2")
    a("E,2,4")
    for d in ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"):
        a("D,1,{},0.0".format(d))
    for d in ("UY", "UZ"):
        a("D,3,{},0.0".format(d))
    a("F,4,FZ,-100000.0")
    a("FINISH")
    a("/SOLU")
    a("ANTYPE,STATIC")
    a("SOLVE")
    a("FINISH")
    a("/POST1")
    a("SET,LAST")
    a("*DIM,ADSP,ARRAY,4,7")
    a("*DO,I,1,4")
    a(" ADSP(I,1)=I")
    for k, c in enumerate(("X", "Y", "Z"), start=2):
        a(" *GET,TMP,NODE,I,U,{}".format(c))
        a(" ADSP(I,{})=TMP".format(k))
    for k, c in enumerate(("X", "Y", "Z"), start=5):
        a(" *GET,TMP,NODE,I,ROT,{}".format(c))
        a(" ADSP(I,{})=TMP".format(k))
    a("*ENDDO")
    a("*CFOPEN,archp_disp,csv")
    a("*VWRITE," + ",".join("ADSP(1,{})".format(k) for k in range(1, 8)))
    a(" (1X,F8.0,6(1X,E18.10))")
    a("*CFCLOS")
    a("ESEL,S,TYPE,,1")
    for k in range(1, 7):
        a("ETABLE,BI{},SMISC,{}".format(k, k))
    a("*DIM,ABM,ARRAY,4,7")
    a("*DO,I,1,4")
    a(" ABM(I,1)=I")
    for k in range(1, 7):
        a(" *GET,TMP,ELEM,I,ETAB,BI{}".format(k))
        a(" ABM(I,{})=TMP".format(k + 1))
    a("*ENDDO")
    a("*CFOPEN,archp_beam,csv")
    a("*VWRITE," + ",".join("ABM(1,{})".format(k) for k in range(1, 8)))
    a(" (1X,F8.0,6(1X,E18.10))")
    a("*CFCLOS")
    a("ESEL,S,TYPE,,2")
    a("ETABLE,LKN,SMISC,1")
    a("*DIM,ALK,ARRAY,1,2")
    a(" ALK(1,1)=5")
    a(" *GET,TMP,ELEM,5,ETAB,LKN")
    a(" ALK(1,2)=TMP")
    a("*CFOPEN,archp_link,csv")
    a("*VWRITE,ALK(1,1),ALK(1,2)")
    a(" (1X,F8.0,1(1X,E18.10))")
    a("*CFCLOS")
    a("ESEL,ALL")
    a("FINISH")
    a("/EXIT,NOSAVE")
    (HERE / "archp.inp").write_text("\n".join(L) + "\n", encoding="ascii")
    q = subprocess.run([str(ANSYS_EXE), "-b", "-j", "archp", "-i", "archp.inp",
                        "-o", "archp.out"], cwd=str(HERE), env=env,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    print("ANSYS rc=", q.returncode)
    for s in ("_disp", "_beam", "_link"):
        f = HERE / ("archp" + s + ".csv")
        print("---", f.name, "---")
        print(f.read_text(encoding="ascii", errors="replace").strip() if f.exists() else "MISSING")
    txt = (HERE / "archp.out").read_text(encoding="utf-8", errors="replace")
    import re
    print("warnings:", re.findall(r"NUMBER OF WARNING MESSAGES ENCOUNTERED=\s*\d+", txt))


if __name__ == "__main__":
    main()
