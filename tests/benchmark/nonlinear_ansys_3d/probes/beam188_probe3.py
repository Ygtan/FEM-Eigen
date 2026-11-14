#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BEAM188 第三轮：找与 Euler-Bernoulli 严格对齐的配置。

变量：KEYOPT(3) 形函数阶次、SECCONTROL 剪切面积量级、网格数 n。
纯 FZ=-10000 悬臂梁，EB 解析 UZ=-3.9009715e-05。
另标定 SMISC 1..19（两列数组一次 *VWRITE，避免参数个数限制）。
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent
E, G = 2.06e11, 7.9e10
A, IYY, IZZ, J = 0.478, 0.4148, 0.2, 0.3
L, PZ = 10.0, -10000.0
EB = PZ * L ** 3 / (3 * E * IYY)

TEMPLATE = r"""/BATCH
/TITLE,BEAM188 align probe {tag}
/PREP7
ET,1,BEAM188
KEYOPT,1,1,0
KEYOPT,1,3,{kop3}
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
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,3
*GET,PD(1,1),NODE,{tip},U,Z
*GET,PD(1,2),NODE,{tip},ROT,Y
*GET,PD(1,3),NODE,1,RF,MY
*CFOPEN,{tag},csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3)
 (1X,3(1X,E20.12))
*CFCLOS
{smisc_block}
FINISH
/EXIT,NOSAVE
"""

SMISC_BLOCK = r"""*DIM,PT,ARRAY,19,2
*DO,I,1,19
 PT(I,1)=I
*ENDDO
*GET,PT(1,2),ELEM,1,SMISC,1
*GET,PT(2,2),ELEM,1,SMISC,2
*GET,PT(3,2),ELEM,1,SMISC,3
*GET,PT(4,2),ELEM,1,SMISC,4
*GET,PT(5,2),ELEM,1,SMISC,5
*GET,PT(6,2),ELEM,1,SMISC,6
*GET,PT(7,2),ELEM,1,SMISC,7
*GET,PT(8,2),ELEM,1,SMISC,8
*GET,PT(9,2),ELEM,1,SMISC,9
*GET,PT(10,2),ELEM,1,SMISC,10
*GET,PT(11,2),ELEM,1,SMISC,11
*GET,PT(12,2),ELEM,1,SMISC,12
*GET,PT(13,2),ELEM,1,SMISC,13
*GET,PT(14,2),ELEM,1,SMISC,14
*GET,PT(15,2),ELEM,1,SMISC,15
*GET,PT(16,2),ELEM,1,SMISC,16
*GET,PT(17,2),ELEM,1,SMISC,17
*GET,PT(18,2),ELEM,1,SMISC,18
*GET,PT(19,2),ELEM,1,SMISC,19
*CFOPEN,{tag}_sm,txt
*VWRITE,PT(1,1),PT(1,2)
 (1X,F4.0,1X,E20.12)
*CFCLOS
"""


def make(tag, n, kop3, sec, with_smisc=False):
    nodes = "\n".join("N,%d,%.10g,0,0" % (i + 1, L * i / n)
                      for i in range(1, n + 1))
    elems = "\n".join("E,%d,%d" % (i, i + 1) for i in range(1, n + 1))
    return TEMPLATE.format(tag=tag, E="%.10g" % E, G="%.10g" % G,
                           A=A, IYY=IYY, IZZ=IZZ, J=J, PZ=PZ,
                           kop3=kop3, seccontrol=sec, nodes=nodes,
                           elems=elems, tip=n + 1,
                           smisc_block=SMISC_BLOCK.format(tag=tag)
                           if with_smisc else "")


def run(tag, text):
    (HERE / (tag + ".inp")).write_text(text, encoding="ascii")
    proc = subprocess.run([ANSYS, "-b", "-j", tag, "-i", tag + ".inp",
                           "-o", tag + ".out"], cwd=str(HERE),
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=300)
    p = HERE / (tag + ".csv")
    vals = None
    if p.exists():
        try:
            vals = [float(v) for v in p.read_text(encoding="ascii").split()]
        except ValueError:
            pass
    out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace") \
        if (HERE / (tag + ".out")).exists() else ""
    errs = [l.strip() for l in out.splitlines() if "*** ERROR ***" in l][:2]
    return vals, errs, proc.returncode


def main():
    print("EB UZ = %.10g" % EB)
    print("\n--- KEYOPT(3) x SECCONTROL x n ---")
    cases = [
        ("k0s0n01", 1, 0, ""),
        ("k0s0n10", 10, 0, ""),
        ("k0s0n40", 40, 0, ""),
        ("k2s0n10", 10, 2, ""),
        ("k3s0n01", 1, 3, ""),
        ("k3s0n10", 10, 3, ""),
        ("k0sAn10", 10, 0, "SECCONTROL,1e12,1e12"),
        ("k0sAn40", 40, 0, "SECCONTROL,1e12,1e12"),
        ("k3sAn01", 1, 3, "SECCONTROL,1e12,1e12"),
        ("k3sAn10", 10, 3, "SECCONTROL,1e12,1e12"),
        ("k0sBn10", 10, 0, "SECCONTROL,1e8,1e8"),
        ("k0sCn10", 10, 0, "SECCONTROL,1e16,1e16"),
    ]
    for tag, n, kop3, sec in cases:
        vals, errs, rc = run(tag, make(tag, n, kop3, sec))
        if vals:
            print("  %-9s n=%2d KOP3=%d sec=%-24s UZ/EB=%.9f  err=%s" %
                  (tag, n, kop3, sec or "(none)", vals[0] / EB,
                   " ".join(errs)))
        else:
            print("  %-9s NO RESULT rc=%d %s" % (tag, rc, " ".join(errs)))

    tag = "b3smisc"
    vals, errs, rc = run(tag, make(tag, 10, 0, "", with_smisc=True))
    p = HERE / (tag + "_sm.txt")
    print("\n--- SMISC 1..19 (elem 1 of n=10 pure FZ) ---")
    if p.exists():
        for line in p.read_text(encoding="ascii", errors="replace").splitlines():
            t = line.split()
            if len(t) >= 2 and abs(float(t[1])) > 1e-8:
                print("   ", line.strip())
    else:
        print("   NO SMISC OUTPUT", errs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
