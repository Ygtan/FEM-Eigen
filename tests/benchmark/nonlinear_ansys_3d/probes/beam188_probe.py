#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ANSYS 18.0 BEAM188 语义探测。

悬臂梁 L=10m 沿 +X，端点同时受 FY / FZ / MX，用于一次测定：
  1. SECTYPE,BEAM,ASEC + SECDATA 参数顺序（尤其 J 的位置）
  2. 无方向节点时局部 y/z 轴朝向；显式方向节点 K 的效果
  3. 是否包含 Timoshenko 剪切变形（与 Euler-Bernoulli 解析解比较）
  4. BEAM188 SMISC 1..30 与 i/j 端 N/Vy/Vz/T/My/Mz 的映射
  5. MP,GXY 是否覆盖 G（扭转刚度）
  6. SECCONTROL 能否关闭横向剪切
"""
import subprocess
import sys
from pathlib import Path

ANSYS = r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe"
HERE = Path(__file__).resolve().parent

E, G = 2.06e11, 7.9e10
A, IYY, IZZ, J = 0.478, 0.4148, 0.2, 0.3
L = 10.0
PY, PZ, TX = -5000.0, -10000.0, 1000.0

W_EB = PZ * L ** 3 / (3.0 * E * IYY)     # UZ 端部挠度（Euler-Bernoulli）
V_EB = PY * L ** 3 / (3.0 * E * IZZ)     # UY 端部挠度
RX_EB = TX * L / (G * J)                  # ROTX 端部转角
W_SH = PZ * L / (0.85 * G * A)            # 若剪切面积 = 0.85A 时的附加挠度量级

TEMPLATE = r"""/BATCH
/TITLE,BEAM188 probe {tag}
/PREP7
ET,1,BEAM188
KEYOPT,1,1,0
MP,EX,1,{E}
MP,PRXY,1,0.3
{gxy}
SECTYPE,1,BEAM,ASEC
SECDATA,{A},{IYY},0,{IZZ},0,{J},0,0,0,0
{seccontrol}
SLIST,1
N,1,0,0,0
N,2,{L},0,0
{k_node}
TYPE,1
MAT,1
SECNUM,1
E,1,2{k_arg}
D,1,ALL,0
F,2,FY,{PY}
F,2,FZ,{PZ}
F,2,MX,{TX}
FINISH
/SOLU
ANTYPE,STATIC
SOLVE
FINISH
/POST1
SET,LAST
*DIM,PD,ARRAY,1,6
*GET,PD(1,1),NODE,2,U,X
*GET,PD(1,2),NODE,2,U,Y
*GET,PD(1,3),NODE,2,U,Z
*GET,PD(1,4),NODE,2,ROT,X
*GET,PD(1,5),NODE,2,ROT,Y
*GET,PD(1,6),NODE,2,ROT,Z
*CFOPEN,{tag}_disp,csv
*VWRITE,PD(1,1),PD(1,2),PD(1,3),PD(1,4),PD(1,5),PD(1,6)
 (1X,6(1X,E20.12))
*CFCLOS
*DIM,PS,ARRAY,1,30
{smisc_gets}
*CFOPEN,{tag}_smisc,csv
*VWRITE,{smisc_names}
 (1X,30(1X,E20.12))
*CFCLOS
*DIM,PR,ARRAY,1,6
*GET,PR(1,1),NODE,1,RF,FX
*GET,PR(1,2),NODE,1,RF,FY
*GET,PR(1,3),NODE,1,RF,FZ
*GET,PR(1,4),NODE,1,RF,MX
*GET,PR(1,5),NODE,1,RF,MY
*GET,PR(1,6),NODE,1,RF,MZ
*CFOPEN,{tag}_reac,csv
*VWRITE,PR(1,1),PR(1,2),PR(1,3),PR(1,4),PR(1,5),PR(1,6)
 (1X,6(1X,E20.12))
*CFCLOS
FINISH
/EXIT,NOSAVE
"""

VARIANTS = [
    # tag, GXY line, SECCONTROL line, K node line, E 命令的 K 参数
    ("b188a", "", "", "", ""),
    ("b188b", "MP,GXY,1,{G}".format(G=G), "", "", ""),
    ("b188c", "MP,GXY,1,{G}".format(G=G), "", "N,3,0,0,1", ",3"),
    ("b188d", "MP,GXY,1,{G}".format(G=G), "SECCONTROL,,,0,0", "", ""),
]


def build(tag, gxy, seccontrol, k_node, k_arg):
    gets, names = [], []
    for i in range(1, 31):
        gets.append("*GET,PS(1,%d),ELEM,1,SMISC,%d" % (i, i))
        names.append("PS(1,%d)" % i)
    return TEMPLATE.format(
        tag=tag, E="%.10g" % E, A=A, IYY=IYY, IZZ=IZZ, J=J, L=L,
        PY=PY, PZ=PZ, TX=TX, gxy=gxy, seccontrol=seccontrol,
        k_node=k_node, k_arg=k_arg,
        smisc_gets="\n".join(gets), smisc_names=",".join(names))


def read_vec(path):
    if not path.exists():
        return None
    try:
        return [float(v) for v in path.read_text(encoding="ascii",
                                                 errors="replace").split()]
    except ValueError:
        return None


def main():
    print("Euler-Bernoulli: UZ=%.9g  UY=%.9g  ROTX=%.9g" % (W_EB, V_EB, RX_EB))
    print("Shear-only estimate (if As=0.85A): UZ extra ~%.9g" % W_SH)
    for tag, gxy, seccontrol, k_node, k_arg in VARIANTS:
        inp = HERE / (tag + ".inp")
        inp.write_text(build(tag, gxy, seccontrol, k_node, k_arg), encoding="ascii")
        proc = subprocess.run(
            [ANSYS, "-b", "-j", tag, "-i", inp.name, "-o", tag + ".out"],
            cwd=str(HERE), capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=300)
        disp = read_vec(HERE / (tag + "_disp.csv"))
        smisc = read_vec(HERE / (tag + "_smisc.csv"))
        reac = read_vec(HERE / (tag + "_reac.csv"))
        out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace") \
            if (HERE / (tag + ".out")).exists() else ""
        errs = [l.strip() for l in out.splitlines()
                if "*** ERROR ***" in l or "ERROR:" in l][:3]
        print("\n--- %s exit=%d ---" % (tag, proc.returncode))
        print("  disp  =", disp)
        print("  reac  =", reac)
        if smisc:
            nz = [(i + 1, v) for i, v in enumerate(smisc) if abs(v) > 1e-6]
            print("  SMISC nonzero:", nz)
        if errs:
            print("  ERRORS:", " | ".join(errs))
        if disp:
            print("  UZ/EB=%.6f  UY/EB=%.6f  ROTX/EB=%.6f" %
                  (disp[2] / W_EB, disp[1] / V_EB, disp[3] / RX_EB))
    return 0


if __name__ == "__main__":
    sys.exit(main())
