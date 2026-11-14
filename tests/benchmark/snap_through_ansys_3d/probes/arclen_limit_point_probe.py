"""Diagnose ANSYS ARCLEN behaviour past the limit point.

Hypothesis: with F_ref = P_ref and default ramping, ANSYS ties the load
factor to TIME/TIME_end and stops at lambda=1, i.e. still on the rising
branch (v=0.12 < v_cr=0.214).  If the reference load exceeds the analytic
limit load, a true arc-length solver must let lambda decrease after the
limit point and follow the descending branch.
"""
import subprocess
from pathlib import Path
import probe

HERE = Path(__file__).resolve().parent
PNEW = 2.0e6  # > analytic P_cr = 1.1884e6

CFGS = [
    ("c0", 0.020, 500, 1.0, False),
    ("c1", 0.005, 1200, 1.0, False),
    ("c2", 0.020, 500, 2.0, False),
    ("c3", 0.020, 500, 1.0, True),
]


def run(tag, maxarc, nsub, tend, kbc):
    apdl = probe.make_apdl(tag, maxarc, nsub, 1.0, tend)
    apdl = apdl.replace("F,2,FY," + probe.n(-probe.PREF),
                        "F,2,FY," + probe.n(-PNEW))
    if kbc:
        apdl = apdl.replace("NLGEOM,ON", "NLGEOM,ON\nKBC,1")
    (HERE / (tag + ".inp")).write_text(apdl, encoding="ascii", newline="\n")
    p = subprocess.run(
        [str(probe.ANSYS), "-b", "-j", tag, "-i", tag + ".inp",
         "-o", tag + ".out"],
        cwd=str(HERE), capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=1800)
    txt = (HERE / (tag + ".out")).read_text(
        encoding="utf-8", errors="replace")
    csv = HERE / (tag + "_curve.csv")
    ne, nw = txt.count("*** ERROR ***"), txt.count("*** WARNING ***")
    if not csv.exists():
        print("%s F=%.1fMN TIME=%.1f KBC=%s  rc=%d CSV-MISSING %dE/%dW"
              % (tag, PNEW / 1e6, tend, kbc, p.returncode, ne, nw))
        return
    c = probe.read_curve(csv)
    pts = [(-r[1], (r[2] + r[3]) / PNEW, r[4]) for r in c]
    vmax = max(v for v, _, _ in pts)
    lmax = max(l for _, l, _ in pts)
    lmin = min(l for _, l, _ in pts)
    tmax = max(t for _, _, t in pts)
    errs = [abs(l * PNEW - probe.p_of_v(v)) for v, l, _ in pts
            if 0.0 <= v <= 2 * probe.H]
    print("%s F=%.1fMN TIME=%.1f KBC=%s  nset=%3d  v_max=%.6f  "
          "lambda=[%+.5f,%+.5f]  TIME_max=%.3f  worst|dP|=%.4g  %dE/%dW"
          % (tag, PNEW / 1e6, tend, kbc, len(c), vmax, lmin, lmax,
             tmax, max(errs) if errs else float("nan"), ne, nw))
    print("    sample: " + "  ".join("%+.3f/%+.4f" % (v, l)
          for v, l, _ in pts[::max(1, len(pts) // 10)][:11]))
    mntr = sorted(HERE.glob(tag + "*.mntr"))
    if mntr:
        lines = mntr[0].read_text(encoding="utf-8", errors="replace").splitlines()
        print("    mntr head: " + " | ".join(x.strip() for x in lines[:6]))


if __name__ == "__main__":
    vcr, pcr = probe.analytic_limit()
    print("analytic limit: v_cr=%.9f  P_cr=%.6f N  lambda_cr(P_ref=2MN)=%.6f"
          % (vcr, pcr, pcr / PNEW))
    for cfg in CFGS:
        run(*cfg)
