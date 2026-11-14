#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G 阶段探测：LINK180 的 ISTRN / 只受拉 / NLGEOM 应变度量标定。

目的：把 FEM-Eigen 共转列式索单元（fem_3d_nl.hpp）
    N = EA/L0 * (Ln - L0),  CABLE 且 w=0 时只受拉
精确映射到 ANSYS 18.0 LINK180 的
    R,id,A,ISTRN  +  KEYOPT(3)  +  NLGEOM,ON
必须实测三件事：
  1) ISTRN 的符号与基准长度（零位移时的初始索力）；
  2) 大位移下 ANSYS 用工程应变 (Ln-Lref)/Lref 还是对数应变 ln(Ln/Lref)；
  3) KEYOPT(3) 哪一个值才是 tension-only，以及压缩时是否真的退出工作。
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANSYS_EXE = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
MSYS_BIN = r"C:\msys64\mingw64\bin"
os.environ["PATH"] = MSYS_BIN + os.pathsep + os.environ.get("PATH", "")

E, A, L = 2.0e11, 0.001, 10.0
EA = E * A


def n(x: float) -> str:
    return "{:.12g}".format(float(x))


def apdl(job: str, *, istrn: float, kop3: int, kop2: int,
         nlgeom: bool, delta: float) -> str:
    o: list[str] = []
    a = o.append
    a("/BATCH")
    a("/TITLE,LINK180 ISTRN/NLGEOM probe - " + job)
    a("/PREP7")
    a("ET,1,LINK180")
    a("KEYOPT,1,2,{}".format(kop2))
    a("KEYOPT,1,3,{}".format(kop3))
    a("R,1,{},{}".format(n(A), n(istrn)))
    a("MP,EX,1,{}".format(n(E)))
    a("MP,PRXY,1,0.3")
    a("MP,DENS,1,0.0")
    a("N,1,0.0,0.0,0.0")
    a("N,2,{},0.0,0.0".format(n(L)))
    a("TYPE,1")
    a("MAT,1")
    a("REAL,1")
    a("E,1,2")
    a("D,1,UX,0.0")
    a("D,1,UY,0.0")
    a("D,1,UZ,0.0")
    a("D,2,UY,0.0")
    a("D,2,UZ,0.0")
    a("D,2,UX,{}".format(n(delta)))
    a("FINISH")
    a("/SOLU")
    a("ANTYPE,STATIC")
    if nlgeom:
        a("NLGEOM,ON")
        a("NSUBST,20,50,5")
    else:
        a("NLGEOM,OFF")
    a("OUTRES,ALL,LAST")
    a("SOLVE")
    a("FINISH")
    a("/POST1")
    a("SET,LAST")
    a("ETABLE,LKN,SMISC,1")
    a("*DIM,ARR,ARRAY,1,3")
    a(" *GET,TMP,ELEM,1,ETAB,LKN")
    a(" ARR(1,1)=TMP")
    a(" *GET,TMP,NODE,2,RF,FX")
    a(" ARR(1,2)=TMP")
    a(" *GET,TMP,NODE,2,U,X")
    a(" ARR(1,3)=TMP")
    a("*CFOPEN,{}_res,csv".format(job))
    a("*VWRITE,ARR(1,1),ARR(1,2),ARR(1,3)")
    a(" (3(1X,E18.10))")
    a("*CFCLOS")
    a("FINISH")
    a("/EXIT,NOSAVE")
    return "\n".join(o) + "\n"


def run(job: str, **kw) -> tuple[int, list[float], str]:
    (HERE / (job + ".inp")).write_text(apdl(job, **kw), encoding="ascii")
    p = subprocess.run([str(ANSYS_EXE), "-b", "-j", job,
                        "-i", job + ".inp", "-o", job + ".out"],
                       cwd=str(HERE), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=300)
    f = HERE / (job + "_res.csv")
    vals = [float(x) for x in f.read_text().split()] if f.exists() else []
    log = (HERE / (job + ".out")).read_text(
        encoding="utf-8", errors="replace") if (HERE / (job + ".out")).exists() else ""
    notes = [s.strip() for s in log.splitlines()
             if "WARNING" in s.upper() or "ERROR" in s.upper()]
    return p.returncode, vals, " | ".join(notes[:2])


def main() -> int:
    cases = [
        # job,           istrn, kop3, kop2, nlgeom, delta
        ("ln_base",      0.0,   0,    0,    False, 0.01),
        ("ln_eng_small", 0.0,   0,    0,    True,  0.01),
        ("ln_eng_big",   0.0,   0,    0,    True,  0.50),
        ("ln_log_big",   0.0,   0,    1,    True,  0.50),
        ("ln_is_pos",    0.01,  0,    0,    True,  0.00),
        ("ln_is_neg",   -0.01,  0,    0,    True,  0.00),
        ("ln_is_pos_d",  0.01,  0,    0,    True,  0.01),
        ("ln_to1_t",     0.0,   1,    0,    True,  0.01),
        ("ln_to1_c",     0.0,   1,    0,    True, -0.01),
        ("ln_to2_c",     0.0,   2,    0,    True, -0.01),
    ]
    print("{:<14s} {:>5s} {:>7s} {:>7s} {:>7s}  {:>16s} {:>16s} {:>16s}  {}".format(
        "job", "rc", "istrn", "kop3", "kop2", "N(SMISC1)", "RF_FX", "UX", "note"))
    for job, istrn, k3, k2, nl, d in cases:
        rc, v, note = run(job, istrn=istrn, kop3=k3, kop2=k2,
                          nlgeom=nl, delta=d)
        vals = (v + [float("nan")] * 3)[:3]
        print("{:<14s} {:>5d} {:>7.4f} {:>7d} {:>7d}  {:>16.8e} {:>16.8e} {:>16.8e}  {}".format(
            job, rc, istrn, k3, k2, vals[0], vals[1], vals[2], note[:70]))

    print()
    print("=== 解析参考（EA = {:.6g} N）===".format(EA))
    for d in (0.01, 0.50):
        Ln = L + d
        print("  delta={:.3f}  Ln={:.4f}  工程 EA*(Ln-L)/L   = {:.8e}".format(
            d, Ln, EA * (Ln - L) / L))
        print("  delta={:.3f}  Ln={:.4f}  对数 EA*ln(Ln/L)    = {:.8e}".format(
            d, Ln, EA * math.log(Ln / L)))
        print("  delta={:.3f}  Ln={:.4f}  FEM EA/L0*(Ln-L0)   = {:.8e} (L0=L)".format(
            d, Ln, EA / L * (Ln - L)))
    print("  ISTRN=+0.01 若为零位移初拉力:  EA*ISTRN = {:.8e}".format(EA * 0.01))
    print("  ISTRN=+0.01 对应 L0 = L/(1+ISTRN) = {:.8f} m".format(L / 1.01))
    return 0


if __name__ == "__main__":
    sys.exit(main())
