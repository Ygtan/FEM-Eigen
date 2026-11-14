"""K 阶段可行性探针：ANSYS BEAM188 + NLGEOM 大转动 / 梁柱效应。

只读原型，不进入仓库。运行后与解析解对比：
  CIRCLE : 悬臂端弯矩，纯弯曲 elastica -> 圆弧闭式解
  PDELTA : 简支梁轴压 + 跨中横向力 -> Timoshenko 梁柱精确解
"""
from __future__ import annotations

import math
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANSYS = Path(r"C:/Program Files/ANSYS Inc/v180/ANSYS/bin/winx64/ANSYS180.exe")

# 与 test_3d_nl.cpp:testCRCantileverCircle / testCRBeamColumn 同源
E = 2.1e11
G = 8.1e10
NU = 0.3


def num(v: float) -> str:
    return "{:.12E}".format(v)


def read_csv(path: Path) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        p = line.split()
        if not p:
            continue
        try:
            out[int(float(p[0]))] = [float(x) for x in p[1:]]
        except ValueError:
            continue
    return out


def make_apdl(*, kind: str, tag: str, n: int, L: float,
              A: float, Iy: float, Iz: float, J: float,
              M: float = 0.0, P: float = 0.0, Q: float = 0.0,
              nsub: int = 200) -> str:
    Ls: list[str] = []
    a = Ls.append
    a("/BATCH")
    a("/TITLE,BEAM188 NLGEOM probe - " + tag)
    a("/PREP7")
    a("ET,1,BEAM188")
    a("KEYOPT,1,1,0")
    a("SECTYPE,1,BEAM,ASEC")
    # SECDATA: A, IYY, IYZ, IZZ, IW, J, CGY, CGZ, SHY, SHZ
    a("SECDATA,{},{},0,{},0,{},0,0,0,0".format(num(A), num(Iy), num(Iz), num(J)))
    a("SECNUM,1")
    a("MP,EX,1,{}".format(num(E)))
    a("MP,PRXY,1,{}".format(num(NU)))
    a("MP,GXY,1,{}".format(num(G)))
    a("! ---- nodes ----")
    for i in range(1, n + 2):
        a("N,{},{},0,0".format(i, num((i - 1) * L / n)))
    a("TYPE,1")
    a("MAT,1")
    for i in range(1, n + 1):
        a("E,{},{}".format(i, i + 1))
    a("! ---- BC ----")
    if kind == "CIRCLE":
        a("D,1,ALL,0")
        a("F,{},MZ,{}".format(n + 1, num(M)))
    elif kind == "PDELTA":
        for d in ("UX", "UY", "UZ", "ROTX"):
            a("D,1,{},0".format(d))
        for d in ("UY", "UZ"):
            a("D,{},{},0".format(n + 1, d))
        a("F,{},FX,{}".format(n + 1, num(-P)))
        a("F,{},FY,{}".format(n // 2 + 1, num(-Q)))
    else:
        raise ValueError(kind)
    a("FINISH")
    a("/SOLU")
    a("ANTYPE,STATIC")
    a("NLGEOM,ON")
    a("NSUBST,{},{},{}".format(nsub, nsub * 5, max(10, nsub // 10)))
    a("NEQIT,200")
    a("OUTRES,ALL,LAST")
    a("SOLVE")
    a("FINISH")
    a("! ---- post ----")
    a("/POST1")
    a("SET,LAST")
    a("*DIM,ADSP,ARRAY,{},7".format(n + 1))
    a("*DO,I,1,{}".format(n + 1))
    a(" ADSP(I,1)=I")
    for k, c in enumerate(("X", "Y", "Z"), start=2):
        a(" *GET,TMP,NODE,I,U,{}".format(c))
        a(" ADSP(I,{})=TMP".format(k))
    for k, c in enumerate(("X", "Y", "Z"), start=5):
        a(" *GET,TMP,NODE,I,ROT,{}".format(c))
        a(" ADSP(I,{})=TMP".format(k))
    a("*ENDDO")
    a("*CFOPEN,{}_disp,csv".format(tag))
    a("*VWRITE," + ",".join("ADSP(1,{})".format(k) for k in range(1, 8)))
    a(" (1X,F8.0,6(1X,E18.10))")
    a("*CFCLOS")
    a("FINISH")
    a("/EXIT,NOSAVE")
    return "\n".join(Ls) + "\n"


def run(tag: str, apdl: str) -> dict[int, list[float]]:
    inp = HERE / (tag + ".inp")
    inp.write_text(apdl, encoding="ascii", newline="\n")
    p = subprocess.run(
        [str(ANSYS), "-b", "-j", tag, "-i", tag + ".inp", "-o", tag + ".out"],
        cwd=str(HERE), capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=1200)
    out = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace")
    ne = out.count("*** ERROR ***")
    nw = out.count("*** WARNING ***")
    csv = HERE / (tag + "_disp.csv")
    print("  job {} rc={} errors={} warnings={} csv={}".format(
        tag, p.returncode, ne, nw, csv.exists()))
    if not csv.exists():
        return {}
    return read_csv(csv)


def circle_exact(L: float, Iz: float, M: float) -> tuple[float, float, float, float]:
    """返回 (theta, x_tip, y_tip, rot_tip)。"""
    theta = M * L / (E * Iz)
    R = L / theta if theta else float("inf")
    return theta, R * math.sin(theta), R * (1.0 - math.cos(theta)), theta


def main() -> None:
    L, A, Iy, Iz, J, n = 10.0, 0.02, 5e-2, 2e-4, 1e-2, 8
    Mfull = 2.0 * math.pi * E * Iz / L
    print("=== CIRCLE: cantilever end moment (E=%.3g, Iz=%.3g, L=%g, n=%d) ===" % (E, Iz, L, n))
    print("  M_full=2*pi*EI/L = %.9g N*m" % Mfull)
    for alpha in (0.25, 0.50, 0.75, 1.00):
        tag = "cir%03d" % int(alpha * 100)
        M = alpha * Mfull
        apdl = make_apdl(kind="CIRCLE", tag=tag, n=n, L=L, A=A, Iy=Iy, Iz=Iz,
                         J=J, M=M, nsub=100)
        d = run(tag, apdl)
        if not d:
            continue
        tip = d[n + 1]
        ux, uy, rotz = tip[0], tip[1], tip[5]
        th, xex, yex, rex = circle_exact(L, Iz, M)
        print("  alpha=%.2f theta=%.6f | ANSYS ux=%+.8g uy=%+.8g rotz=%+.8g" %
              (alpha, th, ux, uy, rotz))
        print("             exact dx=%+.8g dy=%+.8g rot=%+.8g | dpos=%.3e drot=%.3e" %
              (xex - L, yex, rex,
               math.hypot((L + ux) - xex, uy - yex), (rotz - rex) / max(rex, 1)))

    print()
    print("=== PDELTA: simply supported beam-column ===")
    L2, A2, Iy2, Iz2, J2, n2 = 10.0, 0.2, 5e-2, 2e-4, 1e-3, 32
    Pcr = math.pi ** 2 * E * Iz2 / (L2 * L2)
    Q = 2.0e4
    ylin = Q * L2 ** 3 / (48.0 * E * Iz2)
    print("  Pcr=%.9g N, Q=%g N, y_lin=%.9g m" % (Pcr, Q, ylin))
    for ratio in (0.0, 0.25, 0.50, 0.75):
        tag = "pd%03d" % int(ratio * 100)
        P = ratio * Pcr
        if ratio == 0:
            yex = ylin
        else:
            k = math.sqrt(P / (E * Iz2))
            yex = (Q / (2.0 * P * k)) * (math.tan(k * L2 / 2.0) - k * L2 / 2.0)
        apdl = make_apdl(kind="PDELTA", tag=tag, n=n2, L=L2, A=A2, Iy=Iy2,
                         Iz=Iz2, J=J2, P=P, Q=Q, nsub=40)
        d = run(tag, apdl)
        if not d:
            continue
        mid = n2 // 2 + 1
        wy = d[mid][1]
        print("  P/Pcr=%.2f  ANSYS wy=%+.9g  exact=-%.9g  rel=%+.3e  amp=%.5f" %
              (ratio, wy, yex, (abs(wy) - yex) / yex, abs(wy) / ylin))


if __name__ == "__main__":
    main()
