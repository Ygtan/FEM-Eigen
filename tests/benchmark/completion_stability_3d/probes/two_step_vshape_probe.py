"""② 阶段可行性探针：ANSYS 两段求解（成桥态静力 + 增量 ARCLEN）。

目的（不改任何产品代码）：
1. 单步静力（参考构形 + 初应变 + F0）能否收敛到 FEM 成桥态 u0 = -0.1；
2. 第二个荷载步打开 ARCLEN 后，能否继续跟踪增量荷载路径；
3. 用反力反推总荷载 F_total = -(RF1+RF3)，换算
   lambda_FEM = (F_total - F0)/PINC，与 V 形双索闭式 lambda(v) 对拍，
   从而绕开 ANSYS 内部 lambda / TIME 语义歧义。
"""
import math
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANSYS = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")

B, H, E, A = 2.0, 0.6, 2.1e11, 1.0e-3
EA = E * A
LREF = math.hypot(B, H)
L0 = 0.995 * LREF
EPS0 = (LREF - L0) / L0
V0 = 0.10
PINC = 5.0e5


def N_of(v):
    return EA / L0 * (math.hypot(B, H + v) - L0)


def P_of(v):
    L = math.hypot(B, H + v)
    return 2.0 * N_of(v) * (H + v) / L


F0 = P_of(V0)


def lam_closed(v):
    return (P_of(v) - F0) / PINC


def apdl(job, time2, keyopt2=0, step2_load=None, arclen=True):
    if step2_load is None:
        step2_load = -(F0 + PINC)
    a = [
        f"/TITLE,stage2 two-step probe - {job}",
        "/PREP7",
        "ET,1,LINK180",
        f"KEYOPT,1,2,{keyopt2}",
        "KEYOPT,1,3,0",
        f"R,1,{A:.12g}",
        f"MP,EX,1,{E:.12g}",
        "MP,PRXY,1,0.3",
        f"N,1,{-B:.12g},{H:.12g},0",
        "N,2,0,0,0",
        f"N,3,{B:.12g},{H:.12g},0",
        "E,1,2", "E,2,3",
        "INISTATE,SET,DTYP,EPEL",
        f"INISTATE,DEFINE,1,,EPEL,X,{EPS0:.15g}",
        f"INISTATE,DEFINE,2,,EPEL,X,{EPS0:.15g}",
        "D,1,ALL,0", "D,3,ALL,0", "D,2,UX,0", "D,2,UZ,0",
        "FINISH", "/SOLU", "ANTYPE,STATIC", "NLGEOM,ON", "AUTOTS,OFF",
        "! ---- load step 1: completion state (static Newton) ----",
        "TIME,1",
        f"F,2,FY,{-F0:.15g}",
        "NSUBST,20,100,10", "NEQIT,50", "OUTRES,ALL,ALL", "SOLVE",
        "! ---- load step 2: incremental load with arc-length ----",
    ]
    if arclen:
        a += ["ARCLEN,ON,0.01,1.0E-8"]
    a += [
        f"F,2,FY,{step2_load:.15g}",
        "NSUBST,400,2000,20", "NEQIT,60",
        f"TIME,{time2:.12g}",
        "ARCTRM,U,0.15,2,UY",
        "OUTRES,ALL,ALL", "SOLVE", "FINISH",
        "/POST1",
        "*GET,NSET,ACTIVE,0,SET,NSET",
        "*DIM,CURV,ARRAY,NSET,5",
        "*DO,I,1,NSET",
        " SET,,,,,,,I",
        " *GET,VY,NODE,2,U,Y",
        " *GET,R1,NODE,1,RF,FY",
        " *GET,R3,NODE,3,RF,FY",
        " *GET,TNOW,ACTIVE,0,SET,TIME",
        " CURV(I,1)=I", " CURV(I,2)=TNOW", " CURV(I,3)=VY",
        " CURV(I,4)=R1", " CURV(I,5)=R3",
        "*ENDDO",
        f"*CFOPEN,{job}_curve,csv",
        "*VWRITE,CURV(1,1),CURV(1,2),CURV(1,3),CURV(1,4),CURV(1,5)",
        " (1X,F8.0,4(1X,E20.12))",
        "*CFCLOS", "FINISH", "/EXIT,NOSAVE",
    ]
    return "\n".join(a) + "\n"


def run(job, **kw):
    (HERE / (job + ".inp")).write_text(apdl(job, **kw), encoding="ascii")
    p = subprocess.run([str(ANSYS), "-b", "-j", job, "-i", job + ".inp",
                        "-o", job + ".out"], cwd=str(HERE),
                       capture_output=True, text=True, timeout=600)
    return p.returncode


def parse(job):
    path = HERE / (job + "_curve.csv")
    if not path.exists():
        return []
    out = []
    for line in path.read_text(errors="replace").splitlines():
        q = line.split()
        if len(q) < 5:
            continue
        try:
            out.append(tuple(float(x) for x in q[:5]))
        except ValueError:
            continue
    return out  # (iset, time, uy, r1, r3)


def report(job):
    pts = parse(job)
    print(f"\n=== {job}: {len(pts)} sets ===")
    if not pts:
        return
    # step 1 末尾（time 接近 1 的最后一个点）
    s1 = [p for p in pts if p[1] <= 1.0 + 1e-9]
    print(f"  step1 last: time={s1[-1][1]:.8g} uy={s1[-1][2]:.10g} "
          f"(FEM u0 = {-V0})  err={s1[-1][2]+V0:.3g}")
    s2 = [p for p in pts if p[1] > 1.0 + 1e-9]
    if not s2:
        print("  !! no step-2 sets (ARCLEN did not advance)")
        return
    worst, at = 0.0, 0.0
    for _, t, uy, r1, r3 in s2:
        v = -uy
        ftot = r1 + r3          # 向上反力之和 = 向下外荷载量级（正值）
        lam = (ftot - F0) / PINC
        d = abs(lam - lam_closed(v))
        if d > worst:
            worst, at = d, v
    print(f"  worst |lambda_reaction - lambda_closed| = {worst:.6g} at v={at:.6g}")
    print(f"  closed lambda at v=0.15 = {lam_closed(0.15):.8g}")
    print(f"  step2 last: lambda_reaction={((s2[-1][3]+s2[-1][4])-F0)/PINC:.8g}"
          f"  lambda_time={s2[-1][1]-1:.8g}")


if __name__ == "__main__":
    print(f"F0={F0:.12g}  L0={L0:.12g}  EPS0={EPS0:.12g}  PINC={PINC:.12g}")
    cases = [
        ("pD_keyopt1", dict(time2=2.0, keyopt2=1)),
    ]
    for job, kw in cases:
        rc = run(job, **kw)
        print(f"[ANSYS] {job} exit={rc}")
        report(job)
    print("\n=== recheck pA_time2 with corrected reaction sign ===")
    report("pA_time2")
