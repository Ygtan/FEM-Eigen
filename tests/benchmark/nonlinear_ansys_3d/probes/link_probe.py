#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Probe: ANSYS 3D axial element (LINK8 / LINK180) vs FEM-Eigen 3D TRUSS.

Minimal hybrid model: beam + tower column (BEAM) + one stay cable (TRUSS).
All nodes are attached to a BEAM -> 6 DOF everywhere, so both solvers
write the same displacement vector shape.
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
    Node3D, NodalLoad3D, RealConst3D, validate_model_3d,
)

FEM_EXE = Path(r"D:\AI-work\FEM-Eigen\build\example_3d_report.exe")
ANSYS_EXE = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
MSYS_BIN = r"C:\msys64\mingw64\bin"

E_BEAM, E_CABLE = 2.06e11, 1.95e11
A_BEAM, IY, IZ, JTOR, G = 0.02, 1e-3, 1e-3, 2e-3, 7.9e10
A_CABLE = 0.002


def build_model() -> ModelSpec3D:
    nodes = [Node3D(id=1, x=0.0, y=0.0, z=0.0),
             Node3D(id=2, x=10.0, y=0.0, z=0.0),
             Node3D(id=3, x=10.0, y=0.0, z=8.0)]
    elements = [Element3D(id=1, etype=1, n1=1, n2=2),
                Element3D(id=2, etype=1, n1=2, n2=3),
                Element3D(id=3, etype=2, n1=3, n2=1)]
    model = ModelSpec3D(
        title="hybrid beam+truss probe",
        nodes=nodes,
        materials=[Material3D(id=1, E=E_BEAM), Material3D(id=2, E=E_CABLE)],
        reals=[RealConst3D(id=1, A=A_BEAM, IY=IY, IZ=IZ, J=JTOR, G=G),
               RealConst3D(id=2, A=A_CABLE)],
        etypes=[ElemTypeDef3D(id=1, kind="BEAM", mat=1, real=1),
                ElemTypeDef3D(id=2, kind="TRUSS", mat=2, real=2)],
        elements=elements,
        boundaries=[Boundary3D(node=1, dofs=(1, 6), value=0.0)],
        loads=[NodalLoad3D(node=2, dof=3, value=-100000.0)],
    )
    errs, warns = validate_model_3d(model)
    print("validate:", errs, warns)
    return model


def n(x):
    return "{:.12g}".format(float(x))


def apdl(kind: str, job: str) -> str:
    L = []
    a = L.append
    a("/BATCH")
    a("/TITLE,hybrid probe - " + kind)
    a("/PREP7")
    a("ET,1,BEAM4")
    a("ET,2," + kind)
    if kind == "LINK180":
        a("KEYOPT,2,2,0")
    a("R,1,{},{},{},{},{},0.0".format(n(A_BEAM), n(IZ), n(IY), n(1.0), n(1.0)))
    a("RMORE,0.0,{},0.0,0.0,0.0,0.0".format(n(JTOR)))
    a("R,2,{},0.0".format(n(A_CABLE)))
    a("MP,EX,1,{}".format(n(E_BEAM)))
    a("MP,PRXY,1,0.3")
    a("MP,GXY,1,{}".format(n(G)))
    a("MP,EX,2,{}".format(n(E_CABLE)))
    a("MP,PRXY,2,0.3")
    a("N,1,0.0,0.0,0.0")
    a("N,2,10.0,0.0,0.0")
    a("N,3,10.0,0.0,8.0")
    a("TYPE,1")
    a("MAT,1")
    a("REAL,1")
    a("E,1,2")
    a("E,2,3")
    a("TYPE,2")
    a("MAT,2")
    a("REAL,2")
    a("E,3,1")
    a("D,1,UX,0.0")
    a("D,1,UY,0.0")
    a("D,1,UZ,0.0")
    a("D,1,ROTX,0.0")
    a("D,1,ROTY,0.0")
    a("D,1,ROTZ,0.0")
    a("F,2,FZ,-100000.0")
    a("FINISH")
    a("/SOLU")
    a("ANTYPE,STATIC")
    a("SOLVE")
    a("FINISH")
    # ---- post: displacements ----
    a("/POST1")
    a("SET,LAST")
    a("*DIM,ADSP,ARRAY,3,7")
    a("*DO,I,1,3")
    a(" ADSP(I,1)=I")
    a(" *GET,TMP,NODE,I,U,X")
    a(" ADSP(I,2)=TMP")
    a(" *GET,TMP,NODE,I,U,Y")
    a(" ADSP(I,3)=TMP")
    a(" *GET,TMP,NODE,I,U,Z")
    a(" ADSP(I,4)=TMP")
    a(" *GET,TMP,NODE,I,ROT,X")
    a(" ADSP(I,5)=TMP")
    a(" *GET,TMP,NODE,I,ROT,Y")
    a(" ADSP(I,6)=TMP")
    a(" *GET,TMP,NODE,I,ROT,Z")
    a(" ADSP(I,7)=TMP")
    a("*ENDDO")
    a("*CFOPEN,{}_disp,csv".format(job))
    a("*VWRITE," + ",".join("ADSP(1,{})".format(k) for k in range(1, 8)))
    a(" (1X,F8.0,6(1X,E18.10))")
    a("*CFCLOS")
    # ---- post: element forces, selected by TYPE ----
    a("! BEAM4 SMISC 1..6 = i-end N/Vy/Vz/T/My/Mz")
    a("ESEL,S,TYPE,,1")
    for k in range(1, 7):
        a("ETABLE,BI{},SMISC,{}".format(k, k))
    a("*DIM,ABM,ARRAY,2,7")
    a("*DO,I,1,2")
    a(" ABM(I,1)=I")
    for k in range(1, 7):
        a(" *GET,TMP,ELEM,I,ETAB,BI{}".format(k))
        a(" ABM(I,{})=TMP".format(k + 1))
    a("*ENDDO")
    a("*CFOPEN,{}_beam,csv".format(job))
    a("*VWRITE," + ",".join("ABM(1,{})".format(k) for k in range(1, 8)))
    a(" (1X,F8.0,6(1X,E18.10))")
    a("*CFCLOS")
    a("! LINK axial force: probe SMISC 1..6")
    a("ESEL,S,TYPE,,2")
    for k in range(1, 7):
        a("ETABLE,LK{},SMISC,{}".format(k, k))
    a("*DIM,ALK,ARRAY,1,7")
    a(" ALK(1,1)=3")
    for k in range(1, 7):
        a(" *GET,TMP,ELEM,3,ETAB,LK{}".format(k))
        a(" ALK(1,{})=TMP".format(k + 1))
    a("*CFOPEN,{}_link,csv".format(job))
    a("*VWRITE," + ",".join("ALK(1,{})".format(k) for k in range(1, 8)))
    a(" (1X,F8.0,6(1X,E18.10))")
    a("*CFCLOS")
    a("ESEL,ALL")
    a("FINISH")
    a("/EXIT,NOSAVE")
    return "\n".join(L) + "\n"


def run(cmd, cwd):
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    p = subprocess.run([str(c) for c in cmd], cwd=str(cwd), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=300)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main():
    model = build_model()
    (HERE / "hyb.dat").write_text(model.to_dat(compact=False), encoding="utf-8")
    rc, out = run([FEM_EXE, "hyb.dat"], HERE)
    print("FEM rc=", rc)
    print("--- hyb.report.txt ---")
    print((HERE / "hyb.report.txt").read_text(encoding="utf-8", errors="replace"))

    for kind in ("LINK8", "LINK180"):
        job = "hyb_" + kind.lower()
        (HERE / (job + ".inp")).write_text(apdl(kind, job), encoding="ascii")
        rc, out = run([ANSYS_EXE, "-b", "-j", job, "-i", job + ".inp",
                       "-o", job + ".out"], HERE)
        print("=== {} rc={} ===".format(kind, rc))
        outf = HERE / (job + ".out")
        if outf.exists():
            txt = outf.read_text(encoding="utf-8", errors="replace")
            import re
            print("  errors/warnings:", re.findall(r"NUMBER OF (?:ERROR|WARNING) MESSAGES[^\n]*", txt)[:4])
            for line in txt.splitlines():
                if "*** ERROR ***" in line or "*** WARNING ***" in line:
                    print("  ", line.strip()[:160])
        for suffix in ("_disp", "_beam", "_link"):
            f = HERE / (job + suffix + ".csv")
            print("  {}: {}".format(f.name, f.read_text(encoding="ascii", errors="replace").strip()
                                    if f.exists() else "MISSING"))


if __name__ == "__main__":
    sys.exit(main())
