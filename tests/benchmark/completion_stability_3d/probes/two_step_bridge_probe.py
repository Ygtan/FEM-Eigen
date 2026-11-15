"""② 探针 v2：ANSYS 两段求解的方向/口径定标。
v2 = step2 重新施加全部 50 个 F 荷载 + ARCLEN + TIME=6
v3 = step2 荷载控制 Newton（不打开 ARCLEN），TIME=2
"""
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = Path(r"D:\AI-work\BridgeStudio\cases\CB3D10_completion_vs_ansys_nl\cb3d10_completion_nl.inp")
ANSYS = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
PINC = 1.0e6
CTRL = 1016
F_ORIG = -266666.666667
F0_TOTAL = 19.8383009627e6


def make(job, arclen, time_end, reapply, f_end=None, vmax=0.0545):
    if f_end is None:
        f_end = F_ORIG - PINC
    lines = SRC.read_text(encoding="ascii", errors="replace").splitlines()
    cnd, loads = [], []
    for ln in lines:
        s = ln.strip()
        m = re.match(r"^D,(\d+),", s)
        if m and int(m.group(1)) not in cnd:
            cnd.append(int(m.group(1)))
        if s.startswith("F,"):
            loads.append(s)
    print(f"[{job}] constrained={len(cnd)} loads={len(loads)}")

    out, inserted = [], False
    for ln in lines:
        s = ln.strip()
        if s == "/POST1":
            out += ["/POST1", "*GET,NSET,ACTIVE,0,SET,NSET",
                    "*DIM,CURV,ARRAY,NSET,5", "*DO,I,1,NSET",
                    " SET,,,,,,,I", " *GET,VZ,NODE,%d,U,Z" % CTRL,
                    " *GET,TNOW,ACTIVE,0,SET,TIME", " FT=0"]
            for nd in cnd:
                out += [" *GET,RFJ,NODE,%d,RF,FZ" % nd, " FT=FT+RFJ"]
            out += [" CURV(I,1)=I", " CURV(I,2)=TNOW", " CURV(I,3)=VZ",
                    " CURV(I,4)=FT", " CURV(I,5)=0", "*ENDDO",
                    "*CFOPEN,%s_curve,csv" % job,
                    "*VWRITE,CURV(1,1),CURV(1,2),CURV(1,3),CURV(1,4),CURV(1,5)",
                    " (1X,F8.0,4(1X,E20.12))", "*CFCLOS",
                    "FINISH", "/EXIT,NOSAVE"]
            break
        if s == "OUTRES,ALL,LAST":
            out.append("OUTRES,ALL,ALL")
            continue
        out.append(ln)
        if s == "SOLVE" and not inserted:
            inserted = True
            if reapply:
                for f in loads:
                    q = f.split(",")
                    if len(q) >= 4 and int(q[1]) == CTRL and q[2].upper() == "FZ":
                        out.append("F,%d,FZ,%.12g" % (CTRL, f_end))
                    else:
                        out.append(f)
            else:
                out.append("F,%d,FZ,%.12g" % (CTRL, f_end))
            if arclen:
                out.append("ARCLEN,ON,0.01,1.0E-8")
            out += ["NSUBST,600,3000,20", "NEQIT,60",
                    "TIME,%.12g" % time_end,
                    "ARCTRM,U,%.8g,%d,UZ" % (vmax, CTRL),
                    "OUTRES,ALL,ALL", "SOLVE"]
    (HERE / (job + ".inp")).write_text("\n".join(out) + "\n", encoding="ascii")


def run(job):
    p = subprocess.run([str(ANSYS), "-b", "-j", job, "-i", job + ".inp",
                        "-o", job + ".out"], cwd=str(HERE),
                       capture_output=True, text=True, timeout=900)
    return p.returncode


def report(job):
    f = HERE / (job + "_curve.csv")
    if not f.exists():
        print(f"[{job}] no csv")
        return
    pts = []
    for ln in f.read_text(errors="replace").splitlines():
        q = ln.split()
        if len(q) < 4:
            continue
        try:
            pts.append(tuple(float(x) for x in q[:4]))
        except ValueError:
            continue
    s1 = [p for p in pts if p[1] <= 1.0 + 1e-9]
    s2 = [p for p in pts if p[1] > 1.0 + 1e-9]
    print(f"[{job}] sets={len(pts)} step1={len(s1)} step2={len(s2)}")
    if s1:
        print(f"   step1 last: UZ={s1[-1][2]:.9g} FT={s1[-1][3]:.10g}")
    if s2:
        print(f"   step2 first: t={s2[0][1]:.8g} UZ={s2[0][2]:.9g} FT={s2[0][3]:.10g}")
        print(f"   step2 last : t={s2[-1][1]:.8g} UZ={s2[-1][2]:.9g} FT={s2[-1][3]:.10g}")
        # FEM 的 lambda = (FT - F0_total)/PINC（若 FT 是总外荷载量级）
        lam0 = (s2[0][3] - F0_TOTAL) / PINC
        lam1 = (s2[-1][3] - F0_TOTAL) / PINC
        print(f"   lambda_from_reaction: {lam0:.6g} -> {lam1:.6g}")
    o = (HERE / (job + ".out")).read_text(errors="replace")
    er = re.findall(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", o)
    wa = re.findall(r"NUMBER OF WARNING\s+MESSAGES ENCOUNTERED=\s*(\d+)", o)
    print(f"   errors={er[-1] if er else '?'} warnings={wa[-1] if wa else '?'}")


if __name__ == "__main__":
    LAMMAX = 4.6
    cases = [
        ("b4_loadfull", False, 2.0, False, F_ORIG - LAMMAX * PINC, 0.0530),
        ("b5_arcfull", True, 6.0, False, F_ORIG - LAMMAX * PINC, 0.0530),
    ]
    for job, arc, tend, reap, fend, vm in cases:
        make(job, arc, tend, reap, fend, vm)
        rc = run(job)
        print(f"[{job}] ANSYS exit={rc}")
        report(job)
    # ---- FEM 对拍：lambda = (TIME-1)*LAMMAX，ctrl = UZ - UZ_step1 ----
    import bisect
    txt = Path("s_static_cable.txt").read_text(errors="replace")
    fem = [(float(a), float(b)) for a, b in re.findall(
        r"^I STEP=\d+ LAMBDA=(\S+) CTRL=(\S+)", txt, re.M)]
    fl = [p[0] for p in fem]
    fc = [p[1] for p in fem]

    def interp(lq):
        if lq <= fl[0]:
            return fc[0]
        if lq >= fl[-1]:
            return fc[-1]
        j = bisect.bisect_left(fl, lq)
        a, b = fl[j - 1], fl[j]
        return fc[j - 1] + (fc[j] - fc[j - 1]) * (lq - a) / (b - a)

    for job in ("b4_loadfull", "b5_arcfull"):
        f = HERE / (job + "_curve.csv")
        if not f.exists():
            continue
        pts = []
        for ln in f.read_text(errors="replace").splitlines():
            q = ln.split()
            if len(q) < 4:
                continue
            try:
                pts.append(tuple(float(x) for x in q[:4]))
            except ValueError:
                continue
        s1 = [p for p in pts if p[1] <= 1.0 + 1e-9]
        s2 = [p for p in pts if p[1] > 1.0 + 1e-9]
        if not s1 or not s2:
            continue
        uz0 = s1[-1][2]
        worst, at, n = 0.0, 0.0, 0
        for _, t, uz, ft in s2:
            lam = (t - 1.0) * LAMMAX
            if lam > 4.5530691289834 + 1e-6:
                continue
            got = uz - uz0
            want = interp(lam)
            d = abs(got - want)
            n += 1
            if d > worst:
                worst, at = d, lam
        print(f"[{job}] FEM vs ANSYS: {n} points, worst |dUZ|={worst:.6g}"
              f" at lambda={at:.6g}, last UZ={s2[-1][2]:.9g}")
