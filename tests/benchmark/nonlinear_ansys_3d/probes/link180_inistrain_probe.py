#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G 阶段探测 2：LINK180 施加初始索力/初应变的可行机制。

实测已知：R,id,A,ISTRN 的第 3 个参数对 LINK180 无效（ISTRN 被忽略）。
本脚本比较三种替代机制在「零节点位移」下能否产生预期初始拉力
    N0 = EA * 0.01 = 2.0e+6 N  (对应 L0 = L/(1+0.01))
  A) INISTATE,SET,DTYP,EPEL + INISTATE,DEFINE
  B) 温度荷载 BFE,ELEM,TEMP  + MP,ALPX / TREF
  C) 温度荷载 BF,NODE,TEMP   + MP,ALPX / TREF
并分别看 NLGEOM OFF / ON 下的表现。
"""
from __future__ import annotations

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
EPS0 = 0.01
N0 = EA * EPS0


def n(x: float) -> str:
    return "{:.12g}".format(float(x))


def apdl(job: str, mechanism: str, nlgeom: bool, kop2: int = 0) -> str:
    o: list[str] = []
    a = o.append
    a("/BATCH")
    a("/TITLE,LINK180 initial strain mechanism - " + job)
    a("/PREP7")
    a("ET,1,LINK180")
    a("KEYOPT,1,2,{}".format(kop2))
    a("KEYOPT,1,3,0")
    a("R,1,{},0.0".format(n(A)))
    a("MP,EX,1,{}".format(n(E)))
    a("MP,PRXY,1,0.3")
    a("MP,ALPX,1,1.0")
    a("N,1,0.0,0.0,0.0")
    a("N,2,{},0.0,0.0".format(n(L)))
    a("TYPE,1")
    a("MAT,1")
    a("REAL,1")
    a("E,1,2")
    if mechanism == "inistate":
        a("INISTATE,SET,DTYP,EPEL")
        a("INISTATE,DEFINE,1,,EPEL,X,{}".format(n(EPS0)))
    elif mechanism == "inistate_all":
        a("INISTATE,SET,DTYP,EPEL")
        a("INISTATE,DEFINE,ALL,,EPEL,X,{}".format(n(EPS0)))
    elif mechanism == "bfe_temp":
        a("TREF,0.0")
        a("BFE,1,TEMP,1,{}".format(n(-EPS0)))
    elif mechanism == "bf_temp":
        a("TREF,0.0")
        a("BF,1,TEMP,{}".format(n(-EPS0)))
        a("BF,2,TEMP,{}".format(n(-EPS0)))
    a("D,1,UX,0.0")
    a("D,1,UY,0.0")
    a("D,1,UZ,0.0")
    a("D,2,UX,0.0")
    a("D,2,UY,0.0")
    a("D,2,UZ,0.0")
    a("FINISH")
    a("/SOLU")
    a("ANTYPE,STATIC")
    a("NLGEOM,{}".format("ON" if nlgeom else "OFF"))
    if nlgeom:
        a("NSUBST,10,20,5")
    a("OUTRES,ALL,LAST")
    a("SOLVE")
    a("FINISH")
    a("/POST1")
    a("SET,LAST")
    a("ETABLE,LKN,SMISC,1")
    a("*DIM,ARR,ARRAY,1,2")
    a(" *GET,TMP,ELEM,1,ETAB,LKN")
    a(" ARR(1,1)=TMP")
    a(" *GET,TMP,NODE,2,RF,FX")
    a(" ARR(1,2)=TMP")
    a("*CFOPEN,{}_res,csv".format(job))
    a("*VWRITE,ARR(1,1),ARR(1,2)")
    a(" (2(1X,E18.10))")
    a("*CFCLOS")
    a("FINISH")
    a("/EXIT,NOSAVE")
    return "\n".join(o) + "\n"


def run(job: str, mech: str, nl: bool, kop2: int = 0):
    (HERE / (job + ".inp")).write_text(apdl(job, mech, nl, kop2), encoding="ascii")
    p = subprocess.run([str(ANSYS_EXE), "-b", "-j", job,
                        "-i", job + ".inp", "-o", job + ".out"],
                       cwd=str(HERE), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=300)
    f = HERE / (job + "_res.csv")
    v = [float(x) for x in f.read_text().split()] if f.exists() else []
    log = (HERE / (job + ".out")).read_text(
        encoding="utf-8", errors="replace") if (HERE / (job + ".out")).exists() else ""
    errs = [s.strip() for s in log.splitlines()
            if "*** ERROR ***" in s or "INISTATE" in s.upper()
            or "UNRECOGNIZED" in s.upper() or "NOT VALID" in s.upper()]
    return p.returncode, (v + [float("nan")] * 2)[:2], " | ".join(errs[:2])[:110]


def main() -> int:
    print("目标初始拉力 N0 = EA*eps0 = {:.8e} N\n".format(N0))
    jobs = [
        ("ii_none",     "none",         False, 0),
        ("ii_inist",    "inistate",     False, 0),
        ("ii_inistall", "inistate_all", False, 0),
        ("ii_bfe",      "bfe_temp",     False, 0),
        ("ii_bf",       "bf_temp",      False, 0),
        ("ii_bfe_nl",   "bfe_temp",     True,  0),
        ("ii_bf_nl",    "bf_temp",      True,  0),
        ("ii_inist_nl", "inistate",     True,  0),
        ("ii_bfe_k1",   "bfe_temp",     True,  1),
    ]
    print("{:<14s} {:>5s} {:>16s} {:>16s}  {}".format(
        "job", "rc", "N(SMISC1)", "RF_FX", "note"))
    for job, mech, nl, k2 in jobs:
        rc, v, note = run(job, mech, nl, k2)
        print("{:<14s} {:>5d} {:>16.8e} {:>16.8e}  {}".format(
            job, rc, v[0], v[1], note))
    return 0


if __name__ == "__main__":
    sys.exit(main())
