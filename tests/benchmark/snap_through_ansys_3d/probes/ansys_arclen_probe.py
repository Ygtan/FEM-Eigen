"""L-stage feasibility probe: von Mises two-bar shallow truss snap-through.

ANSYS: LINK180 + NLGEOM + ARCLEN,ON + ARCTRM.  For every converged
substep read apex UY and both support FY reactions; lambda is recovered
from the reactions because the only applied load is at the apex.

Independent analytic anchor:
    L(v) = sqrt(b^2 + (h-v)^2),  N(v) = EA*(L-L0)/L0
    P(v) = -2*N*(h-v)/L
"""
import math
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANSYS = Path(r"C:/Program Files/ANSYS Inc/v180/ANSYS/bin/winx64/ANSYS180.exe")

E = 2.1e11
A = 1.0e-3
B = 2.0
H = 0.5
L0 = math.hypot(B, H)
PREF = 1.0e6

TAGS = ("a0", "a1", "a2")


def n(x):
    return "%.12E" % x


def p_of_v(v):
    L = math.hypot(B, H - v)
    N = E * A * (L - L0) / L0
    return -2.0 * N * (H - v) / L


def analytic_limit():
    """Golden-section search for the first limit point dP/dv = 0."""
    lo, hi = 0.0, H
    for _ in range(200):
        m1 = lo + (hi - lo) * 0.381966
        m2 = lo + (hi - lo) * 0.618034
        if p_of_v(m1) < p_of_v(m2):
            lo = m1
        else:
            hi = m2
    v = 0.5 * (lo + hi)
    return v, p_of_v(v)


def make_apdl(tag, maxarc=0.05, nsub=400, vmax=1.0, tend=4.0):
    a = []
    a.append("/TITLE,von Mises truss snap-through probe - " + tag)
    a.append("/PREP7")
    a.append("ET,1,LINK180")
    a.append("KEYOPT,1,2,0")
    a.append("KEYOPT,1,3,0")
    a.append("R,1," + n(A))
    a.append("MP,EX,1," + n(E))
    a.append("MP,PRXY,1,0.3")
    a.append("MP,DENS,1,0")
    a.append("N,1," + n(-B) + "," + n(0.0) + ",0")
    a.append("N,2," + n(0.0) + "," + n(H) + ",0")
    a.append("N,3," + n(B) + "," + n(0.0) + ",0")
    a.append("E,1,2")
    a.append("E,2,3")
    a.append("D,1,ALL,0")
    a.append("D,3,ALL,0")
    a.append("D,2,UZ,0")
    a.append("FINISH")
    a.append("/SOLU")
    a.append("ANTYPE,STATIC")
    a.append("NLGEOM,ON")
    a.append("AUTOTS,OFF")
    a.append("ARCLEN,ON," + n(maxarc) + ",1.0E-7")
    a.append("NSUBST,%d,%d,%d" % (nsub, nsub * 5, 20))
    a.append("NEQIT,60")
    a.append("TIME," + n(tend))
    a.append("OUTRES,ALL,ALL")
    a.append("F,2,FY," + n(-PREF))
    a.append("ARCTRM,U," + n(vmax) + ",2,UY")
    a.append("SOLVE")
    a.append("FINISH")
    a.append("/POST1")
    a.append("*GET,NSET,ACTIVE,0,SET,NSET")
    a.append("*DIM,CURV,ARRAY,NSET,5")
    a.append("*DO,I,1,NSET")
    a.append(" SET,,,,,,,I")
    a.append(" *GET,VY,NODE,2,U,Y")
    a.append(" *GET,R1,NODE,1,RF,FY")
    a.append(" *GET,R3,NODE,3,RF,FY")
    a.append(" *GET,TNOW,ACTIVE,0,SET,TIME")
    a.append(" CURV(I,1)=I")
    a.append(" CURV(I,2)=VY")
    a.append(" CURV(I,3)=R1")
    a.append(" CURV(I,4)=R3")
    a.append(" CURV(I,5)=TNOW")
    a.append("*ENDDO")
    a.append("*CFOPEN," + tag + "_curve,csv")
    a.append("*VWRITE,CURV(1,1),CURV(1,2),CURV(1,3),CURV(1,4),CURV(1,5)")
    a.append(" (1X,F8.0,4(1X,E20.12))")
    a.append("*CFCLOS")
    a.append("FINISH")
    a.append("/EXIT,NOSAVE")
    return "\n".join(a) + "\n"


def read_curve(path):
    out = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        p = line.split()
        if len(p) < 5:
            continue
        try:
            out.append((int(float(p[0])), float(p[1]), float(p[2]),
                        float(p[3]), float(p[4])))
        except ValueError:
            continue
    return out


if __name__ == "__main__":
    vcr, pcr = analytic_limit()
    print("analytic: L0=%.9f  v_cr=%.9f  P_cr=%.6f N  lambda_cr=%.6f"
          % (L0, vcr, pcr, pcr / PREF))
    print("analytic: P(h)=%.6f  P(2h)=%.6f  P(0)=%.6f"
          % (p_of_v(H), p_of_v(2 * H), p_of_v(0.0)))
    cfgs = [("b0", 0.05, 400, 4.0), ("b1", 0.02, 800, 4.0),
            ("b2", 0.10, 300, 4.0)]
    for tag, mx, ns, tend in cfgs:
        inp = HERE / (tag + ".inp")
        inp.write_text(make_apdl(tag, mx, ns, 1.0, tend), encoding="ascii", newline="\n")
        p = subprocess.run(
            [str(ANSYS), "-b", "-j", tag, "-i", tag + ".inp", "-o", tag + ".out"],
            cwd=str(HERE), capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=1200)
        out = HERE / (tag + ".out")
        csv = HERE / (tag + "_curve.csv")
        txt = out.read_text(encoding="utf-8", errors="replace") if out.exists() else ""
        ne = txt.count("*** ERROR ***")
        nw = txt.count("*** WARNING ***")
        if not csv.exists():
            print("%s maxarc=%.3f nsub=%d rc=%d CSV-MISSING %dE/%dW"
                  % (tag, mx, ns, p.returncode, ne, nw))
            continue
        c = read_curve(csv)
        if not c:
            print("%s empty curve %dE/%dW" % (tag, ne, nw))
            continue
        # lambda = (R1+R3)/PREF ; v = -UY
        pts = [(-r[1], (r[2] + r[3]) / PREF, r[4]) for r in c]
        vmax = max(v for v, _, _ in pts)
        lmax = max(l for _, l, _ in pts)
        lmin = min(l for _, l, _ in pts)
        errs = []
        for v, lam, _ in pts:
            if 0.0 <= v <= 2 * H:
                errs.append(abs(lam * PREF - p_of_v(v)))
        worst = max(errs) if errs else float("nan")
        print("%s maxarc=%.3f nsub=%d TIME=%.1f  nset=%d  v_max=%.6f  "
              "lambda=[%+.5f,%+.5f]  worst|dP|=%.4g  %dE/%dW"
              % (tag, mx, ns, tend, len(c), vmax, lmin, lmax, worst, ne, nw))
        if tag == "b1":
            print("   sample (v,lambda,TIME): " + "  ".join(
                "%+.4f/%+.4f/%.3f" % x for x in pts[::max(1, len(pts) // 8)][:9]))
