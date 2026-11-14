"""Locate the ~3% systematic offset between ANSYS LINK180 and the
engineering-strain truss constitutive used by fem_3d_nl.hpp.

Candidates:
  (a) strain measure (engineering vs logarithmic vs Almansi),
  (b) LINK180 KEYOPT(2),
  (c) area definition (R vs SECTYPE/SECDATA).
"""
import math
import subprocess
from pathlib import Path
import probe

HERE = Path(__file__).resolve().parent
PNEW = 2.0e6
E, A, B, H, L0 = probe.E, probe.A, probe.B, probe.H, probe.L0


def p_limit(constitutive):
    def p_of_v2(v):
        L = math.hypot(B, H - v)
        if constitutive == "eng":
            N = E * A * (L - L0) / L0
        elif constitutive == "log":
            N = E * A * math.log(L / L0)
        elif constitutive == "alm":
            N = E * A * (L - L0) / L
        elif constitutive == "green":
            N = E * A * (L * L - L0 * L0) / (2 * L0 * L0) * (L0 / L)
        return -2.0 * N * (H - v) / L
    lo, hi = 0.0, H
    for _ in range(200):
        m1, m2 = lo + (hi - lo) * 0.381966, lo + (hi - lo) * 0.618034
        if p_of_v2(m1) < p_of_v2(m2):
            lo = m1
        else:
            hi = m2
    v = 0.5 * (lo + hi)
    return v, p_of_v2(v), p_of_v2


print("=== analytic limit load under alternative strain measures ===")
for c in ("eng", "log", "alm", "green"):
    v, p, _ = p_limit(c)
    print("  %-6s  v_cr=%.9f  P_cr=%.4f N  lambda(P_ref=2MN)=%.9f"
          % (c, v, p, p / PNEW))
print("  ANSYS measured lambda_max ~ 0.61257  ->  P=%.4f N" % (0.61257 * PNEW))

CFG = [
    ("d0", "KEYOPT,1,2,0", "R"),
    ("d1", "KEYOPT,1,2,1", "R"),
    ("d2", "KEYOPT,1,2,0", "SEC"),
    ("d3", "KEYOPT,1,2,1", "SEC"),
]


def run(tag, keyopt2, area_mode):
    apdl = probe.make_apdl(tag, 0.005, 1200, 1.0, 1.0)
    apdl = apdl.replace("F,2,FY," + probe.n(-probe.PREF),
                        "F,2,FY," + probe.n(-PNEW))
    apdl = apdl.replace("KEYOPT,1,2,0", keyopt2)
    if area_mode == "SEC":
        apdl = apdl.replace("R,1," + n(A), "SECTYPE,1,LINK\nSECDATA," + n(A))
    (HERE / (tag + ".inp")).write_text(apdl, encoding="ascii", newline="\n")
    p = subprocess.run(
        [str(probe.ANSYS), "-b", "-j", tag, "-i", tag + ".inp",
         "-o", tag + ".out"],
        cwd=str(HERE), capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=1800)
    csv = HERE / (tag + "_curve.csv")
    txt = (HERE / (tag + ".out")).read_text(encoding="utf-8", errors="replace")
    ne, nw = txt.count("*** ERROR ***"), txt.count("*** WARNING ***")
    if not csv.exists():
        print("  %s %s/%s  CSV-MISSING %dE/%dW" % (tag, keyopt2, area_mode, ne, nw))
        return
    c = probe.read_curve(csv)
    pts = [(-r[1], (r[2] + r[3]) / PNEW) for r in c]
    lmax = max(l for _, l in pts)
    errs = [abs(l * PNEW - probe.p_of_v(v)) for v, l in pts if 0 <= v <= 2 * H]
    print("  %s %s/%s  nset=%3d  v_max=%.6f  lambda_max=%+.6f  "
          "worst|dP|=%.4g  %dE/%dW"
          % (tag, keyopt2, area_mode, len(c), max(v for v, _ in pts),
             lmax, max(errs), ne, nw))


def n(x):
    return "%.12E" % x


if __name__ == "__main__":
    print("=== ANSYS LINK180 variants ===")
    for cfg in CFG:
        run(*cfg)
    print("=== c1.out warnings ===")
    for line in (HERE / "c1.out").read_text(
            encoding="utf-8", errors="replace").splitlines():
        if "WARNING" in line:
            print("  " + line.strip())
