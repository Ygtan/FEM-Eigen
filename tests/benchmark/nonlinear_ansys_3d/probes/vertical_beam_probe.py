#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F 阶段探测：竖直 BEAM4 局部坐标轴与方向节点 K 的语义标定。

FEM-Eigen 3D 的 beam3DTransformation 对「轴平行全局 Z」的单元取
    ez_ref = global Y  ->  ey = Y x Z = +X,  ez = Z x X = +Y
即 local y = +X, local z = +Y。局部刚度里
    v-theta_z (local y 向弯曲) 用 IZ，w-theta_y (local z 向弯曲) 用 IY，
所以全局 FX 使 UX = P L^3/(3 E IZ)，FY 使 UY = P L^3/(3 E IY)。

若 ANSYS BEAM4 对竖直单元默认 local y = global Y，则 IY/IZ 互换，
UX/UY 与 FEM 对调。本脚本用 IY != IZ 的截面同时施加 FX、FY，比较
「默认 / 方向节点 K 在 +Y / K 在 +X」，找出与 FEM 严格一致的写法。
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
os.environ["PATH"] = MSYS_BIN + os.pathsep + os.environ.get("PATH", "")

E, A = 2.06e11, 0.10
IY, IZ = 0.010, 0.040          # 故意取 4 倍差，便于一眼看出互换
JTOR, GS = 0.020, 7.9e10
L, P = 10.0, 1.0e5


def n(x: float) -> str:
    return "{:.12g}".format(float(x))


def build_fem() -> ModelSpec3D:
    return ModelSpec3D(
        title="vertical BEAM orientation probe",
        nodes=[Node3D(id=1, x=0.0, y=0.0, z=0.0),
               Node3D(id=2, x=0.0, y=0.0, z=L)],
        materials=[Material3D(id=1, E=E)],
        reals=[RealConst3D(id=1, A=A, IY=IY, IZ=IZ, J=JTOR, G=GS)],
        etypes=[ElemTypeDef3D(id=1, kind="BEAM", mat=1, real=1)],
        elements=[Element3D(id=1, etype=1, n1=1, n2=2)],
        boundaries=[Boundary3D(node=1, dofs=(1, 6), value=0.0)],
        loads=[NodalLoad3D(node=2, dof=1, value=P),
               NodalLoad3D(node=2, dof=2, value=P)],
    )


def apdl(job: str, case: str) -> str:
    out: list[str] = []
    a = out.append
    a("/BATCH")
    a("/TITLE,vertical BEAM4 orientation probe - " + case)
    a("/PREP7")
    a("ET,1,BEAM4")
    # ANSYS: R,id,A,IZZ,IYY,TKZ,TKY,...（与 spec_mixed_to_apdl 一致）
    a("R,1,{},{},{},{},{},0.0".format(n(A), n(IZ), n(IY), n(1.0), n(1.0)))
    a("RMORE,0.0,{},0.0,0.0,0.0,0.0".format(n(JTOR)))
    a("MP,EX,1,{}".format(n(E)))
    a("MP,PRXY,1,0.3")
    a("MP,GXY,1,{}".format(n(GS)))
    a("N,1,0.0,0.0,0.0")
    a("N,2,0.0,0.0,{}".format(n(L)))
    if case == "kY":
        a("N,3,0.0,1.0,0.0")            # K 在全局 +Y
    elif case == "kX":
        a("N,3,1.0,0.0,0.0")            # K 在全局 +X
    elif case == "kY_fixed":
        a("N,3,0.0,1.0,0.0")
    a("TYPE,1")
    a("MAT,1")
    a("REAL,1")
    a("E,1,2{}".format(",3" if case.startswith("k") else ""))
    a("D,1,UX,0.0")
    a("D,1,UY,0.0")
    a("D,1,UZ,0.0")
    a("D,1,ROTX,0.0")
    a("D,1,ROTY,0.0")
    a("D,1,ROTZ,0.0")
    if case == "kY_fixed":
        a("D,3,UX,0.0")
        a("D,3,UY,0.0")
        a("D,3,UZ,0.0")
    a("F,2,FX,{}".format(n(P)))
    a("F,2,FY,{}".format(n(P)))
    a("FINISH")
    a("/SOLU")
    a("ANTYPE,STATIC")
    a("SOLVE")
    a("FINISH")
    a("/POST1")
    a("SET,LAST")
    a("*DIM,ARR,ARRAY,1,3")
    a(" *GET,TMP,NODE,2,U,X")
    a(" ARR(1,1)=TMP")
    a(" *GET,TMP,NODE,2,U,Y")
    a(" ARR(1,2)=TMP")
    a(" *GET,TMP,NODE,2,U,Z")
    a(" ARR(1,3)=TMP")
    a("*CFOPEN,{}_disp,csv".format(job))
    a("*VWRITE,ARR(1,1),ARR(1,2),ARR(1,3)")
    a(" (3(1X,E18.10))")
    a("*CFCLOS")
    a("FINISH")
    a("/EXIT,NOSAVE")
    return "\n".join(out) + "\n"


def run(cmd: list[str]) -> int:
    p = subprocess.run([str(c) for c in cmd], cwd=str(HERE),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=300)
    return p.returncode


def main() -> int:
    # FEM-Eigen 参考解
    model = build_fem()
    (HERE / "vbeam.dat").write_text(model.to_dat(compact=False), encoding="utf-8")
    rc = run([FEM_EXE, "vbeam.dat"])
    print("FEM-Eigen exit =", rc)
    fem = {}
    for line in (HERE / "vbeam.report.txt").read_text(
            encoding="utf-8", errors="replace").splitlines():
        t = line.replace(",", " ").split()
        if len(t) >= 7 and t[0] == "2":
            try:
                fem = {"UX": float(t[1]), "UY": float(t[2])}
            except ValueError:
                pass
            break
    print("FEM  UX={:.8e}  UY={:.8e}".format(fem.get("UX", float("nan")),
                                              fem.get("UY", float("nan"))))
    print("理论 UX(FX,IZ)={:.8e}  UY(FY,IY)={:.8e}".format(
        P * L ** 3 / (3 * E * IZ), P * L ** 3 / (3 * E * IY)))
    print()

    cases = [("default", "vbd"), ("kY", "vby"), ("kX", "vbx"),
             ("kY_fixed", "vbyf")]
    for case, job in cases:
        (HERE / (job + ".inp")).write_text(apdl(job, case), encoding="ascii")
        rc = run([ANSYS_EXE, "-b", "-j", job, "-i", job + ".inp",
                  "-o", job + ".out"])
        f = HERE / (job + "_disp.csv")
        vals = []
        if f.exists():
            vals = [float(x) for x in f.read_text().split()]
        ux = vals[0] if len(vals) > 0 else float("nan")
        uy = vals[1] if len(vals) > 1 else float("nan")
        log = (HERE / (job + ".out")).read_text(
            encoding="utf-8", errors="replace") if (HERE / (job + ".out")).exists() else ""
        warn = [s.strip() for s in log.splitlines() if "WARNING" in s.upper()][:2]
        ferr = fem.get("UX", float("nan"))
        ok = (abs(ux - ferr) <= 1e-6 * max(abs(ferr), 1e-30)
              and abs(uy - fem.get("UY", float("nan")))
              <= 1e-6 * max(abs(fem.get("UY", float("nan"))), 1e-30))
        print("{:<10s} rc={}  UX={:.8e}  UY={:.8e}  {}  {}".format(
            case, rc, ux, uy, "MATCH_FEM" if ok else "DIFF", "; ".join(warn)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
