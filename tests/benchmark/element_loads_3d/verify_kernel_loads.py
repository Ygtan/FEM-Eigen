#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Kernel-level validation of 3D *ELEMENT_LOAD / *GRAVITY / *SELF_WEIGHT.

Every case is checked against an INDEPENDENT Euler-Bernoulli closed form,
not against ANSYS. ANSYS cross-checks of the production bridge cases live in
compare_element_loads_3d.py.

Coverage:
  * GZ / GY / LZ directions (both 3D bending planes, IY vs IZ)
  * trapezoidal q1 != q2 via a strict linear-superposition identity
  * automatic selfweight (DENS), ADDMASS, *SELF_WEIGHT RHO/MASS
  * selfweight of a vertical TRUSS (axial element, separate report path)

Console output is pure ASCII on purpose (CTest log / PowerShell / redirect).
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

E, G = 200e9, 7.9e10
A, IY, IZ, J = 0.01, 2e-4, 1e-4, 1e-4   # IY (z-plane) != IZ (y-plane) on purpose
L = 4.0
DENS, GACC = 7850.0, 9.81
MSYS_BIN = r"C:\msys64\mingw64\bin"

HEAD = """# 3D distributed-load kernel check: {name}
*ET
1 {etype}
*MP
1 E={E} {dens}
*R
1 A={A} IY={IY} IZ={IZ} J={J} G={G} {addmass}
*NODE
1 0.0 0.0 0.0
2 2.0 0.0 0.0
3 4.0 0.0 0.0
*TYPE,1
*MAT,1
*REAL,1
*EN,1,1,2
*EN,2,2,3
*BOUNDARY
{constraints}
{loadblocks}
*END
"""

# Horizontal beam, simply supported in BOTH bending planes:
#   vertical (UZ / ROTY) and lateral (UY / ROTZ) rotations are free,
#   all rigid-body / torsional DOFs restrained.
BC_H = "\n".join([
    "1 1 0.0", "1 2 0.0", "1 3 0.0", "1 4 0.0",
    "3 2 0.0", "3 3 0.0", "3 4 0.0",
])

TRUSS_DAT = """# vertical TRUSS under selfweight (top fixed, bottom free)
*ET
1 TRUSS
*MP
1 E=200e9 DENS=7850
*R
1 A=0.01
*NODE
1 0.0 0.0 1.0
2 0.0 0.0 0.0
*TYPE,1
*MAT,1
*REAL,1
*EN,1,1,2
*BOUNDARY
1 1-3 0.0
2 1 0.0
2 2 0.0
*GRAVITY
0 0 -9.81
*END
"""


def beam_dat(name, loadblocks, *, dens=False, addmass=0.0, etype="BEAM", bc=BC_H):
    return HEAD.format(name=name, etype=etype, E="%.10g" % E,
                       dens=("DENS=%.10g" % DENS) if dens else "",
                       A=A, IY="%.10g" % IY, IZ="%.10g" % IZ,
                       J="%.10g" % J, G="%.10g" % G,
                       addmass=("ADDMASS=%.10g" % addmass) if addmass else "",
                       constraints=bc, loadblocks=loadblocks)


def beam_expect(q, EI, *, plane="z"):
    """Exact simply-supported EB solution for a downward UDL q [N/m]."""
    vmid = -5.0 * q * L ** 4 / (384.0 * EI)
    theta = q * L ** 3 / (24.0 * EI)
    R1 = q * L / 2.0
    if plane == "z":
        return {
            "UZ_2": vmid,
            "ROTY_1": theta,
            "Vz_e1": (R1, True),
            "My_e1": (0.0, True),
            "My_e2": (q * L * L / 8.0, True),
        }
    return {
        "UY_2": vmid,
        "ROTZ_1": -theta,
        "Vy_e1": (R1, True),
        "Mz_e1": (0.0, True),
    }


CASES = [
    dict(name="udl_gz", dat=beam_dat("GZ", "*ELEMENT_LOAD\n1-2 GZ -10000 -10000"),
         expect=beam_expect(10000.0, E * IY, plane="z")),
    dict(name="udl_gy", dat=beam_dat("GY", "*ELEMENT_LOAD\n1-2 GY -10000 -10000"),
         expect=beam_expect(10000.0, E * IZ, plane="y")),
    dict(name="udl_lz", dat=beam_dat("LZ", "*ELEMENT_LOAD\n1-2 LZ -10000 -10000"),
         expect=beam_expect(10000.0, E * IY, plane="z")),
    dict(name="selfweight_dens",
         dat=beam_dat("DENS", "*GRAVITY\n0 0 -9.81", dens=True),
         expect=beam_expect(DENS * A * GACC, E * IY, plane="z")),
    dict(name="selfweight_addmass",
         dat=beam_dat("ADDMASS", "*GRAVITY\n0 0 -9.81", dens=True, addmass=50.0),
         expect=beam_expect((DENS * A + 50.0) * GACC, E * IY, plane="z")),
    dict(name="selfweight_mass",
         dat=beam_dat("SELF_WEIGHT MASS",
                      "*GRAVITY\n0 0 -9.81\n*SELF_WEIGHT\nALL MASS 100"),
         expect=beam_expect(100.0 * GACC, E * IY, plane="z")),
    dict(name="selfweight_rho",
         dat=beam_dat("SELF_WEIGHT RHO",
                      "*GRAVITY\n0 0 -9.81\n*SELF_WEIGHT\nALL RHO 7850"),
         expect=beam_expect(7850.0 * A * GACC, E * IY, plane="z")),
]


def parse_report(text):
    disp, force, section = {}, {}, None
    for ln in text.splitlines():
        if "节点位移" in ln:
            section = "disp"
            continue
        if "单元内力" in ln:
            section = "force"
            continue
        if section == "disp":
            t = ln.split()
            if t and re.match(r"^\d+$", t[0]):
                disp[int(t[0])] = [float(v) for v in t[1:]]
        elif section == "force":
            t = ln.split()
            if len(t) >= 2 and re.match(r"^\d+$", t[0]):
                kv = {}
                for tok in t[2:]:
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        try:
                            kv[k] = float(v)
                        except ValueError:
                            pass
                force[int(t[0])] = kv
    return disp, force


def read_fem(exe, dat_text, work, stem):
    p = work / (stem + ".dat")
    p.write_text(dat_text, encoding="utf-8")
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    r = subprocess.run([str(exe), p.name], cwd=str(work), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=120)
    if r.returncode != 0:
        raise RuntimeError("solver rc=%d stderr=%s" % (r.returncode, r.stderr[-300:]))
    rp = work / (stem + ".report.txt")
    if not rp.exists():
        raise RuntimeError("no report.txt produced")
    return parse_report(rp.read_text(encoding="utf-8", errors="replace"))


def got_scalars(disp, force):
    d1 = disp.get(1, [0.0] * 6)
    d2 = disp.get(2, [0.0] * 6)
    f1 = force.get(1, {})
    f2 = force.get(2, {})
    return {
        "UZ_2": d2[2] if len(d2) > 2 else 0.0,
        "UY_2": d2[1] if len(d2) > 1 else 0.0,
        "ROTY_1": d1[4] if len(d1) > 4 else 0.0,
        "ROTZ_1": d1[5] if len(d1) > 5 else 0.0,
        "Vz_e1": f1.get("Vz", 0.0),
        "My_e1": f1.get("My", 0.0),
        "Vy_e1": f1.get("Vy", 0.0),
        "Mz_e1": f1.get("Mz", 0.0),
        "My_e2": f2.get("My", 0.0),
    }


def evaluate(case, disp, force):
    got = got_scalars(disp, force)
    rows = []
    for key, ref in case["expect"].items():
        if isinstance(ref, tuple):
            want = ref[0]
            g, w = abs(got[key]), abs(want)
        else:
            want = ref
            g, w = got[key], want
        scale = max(abs(w), 1e-12)
        ok = abs(g - w) <= 1e-6 * scale + 1e-9
        rows.append((key, got[key], want, ok))
    return rows


def superposition_case(exe, work):
    """fe(q1,q2) = fe(q1,0) + fe(0,q2), strictly (FEF is linear in q1/q2).

    Note it is NOT (fe(q1,q1)+fe(q2,q2))/2: the consistent load vector uses
    L/6(2q1+q2), L/20(7q1+3q2), L^2/60(3q1+2q2), so the two end values
    carry different weights. The triangular decomposition verifies BOTH
    coefficients separately and needs no analytic trapezoid solution.
    """
    def one(stem, q1, q2):
        dat = beam_dat(stem, "*ELEMENT_LOAD\n1-2 GZ -%.10g -%.10g" % (q1, q2))
        disp, force = read_fem(exe, dat, work, stem)
        return got_scalars(disp, force), disp, force

    tri_i, di, fi = one("trap_i", 10000.0, 0.0)
    tri_j, dj, fj = one("trap_j", 0.0, 20000.0)
    mix, dm, fm = one("trap_mix", 10000.0, 20000.0)
    # report.txt carries 8 significant digits -> relative tolerance 1e-6.
    tol = lambda w: 1e-6 * max(abs(w), 1e-3) + 1e-6   # noqa: E731
    rows = []
    for node in sorted(di):
        for k in range(len(di[node])):
            key = "U_node%d_dof%d" % (node, k + 1)
            want = di[node][k] + dj[node][k]
            g = dm[node][k]
            rows.append((key, g, want, abs(g - want) <= tol(want)))
    for elem in sorted(fi):
        for key in sorted(fi[elem]):
            want = fi[elem][key] + fj[elem][key]
            g = fm[elem][key]
            rows.append(("F_e%d_%s" % (elem, key), g, want,
                         abs(g - want) <= tol(want)))
    # q1 != q2 must actually be different from a uniform q1 load
    uni, _, _ = one("trap_uniform", 10000.0, 10000.0)
    delta = abs(mix["Vz_e1"] - uni["Vz_e1"])
    rows.append(("trapezoid_is_active", delta, 8333.3333,
                 abs(delta - 8333.3333) <= 1e-2))
    return rows


def truss_case(exe, work):
    disp, force = read_fem(exe, TRUSS_DAT, work, "truss_sw")
    w = DENS * A * GACC
    u2 = disp.get(2, [0.0, 0.0, 0.0])[2]
    N1 = force.get(1, {}).get("N", 0.0)
    ref_u = -w * 1.0 ** 2 / (2.0 * E * A)
    ref_N = w * 1.0
    return [("UZ_free_end", u2, ref_u,
             abs(u2 - ref_u) <= 1e-6 * abs(ref_u) + 1e-12),
            ("N_fixed_end", N1, ref_N,
             abs(N1 - ref_N) <= 1e-6 * abs(ref_N) + 1e-9)]


def report_block(label, rows, tally):
    bad = [r for r in rows if not r[3]]
    print("-- %-24s %d checks %s" % (label, len(rows), "PASS" if not bad else "FAIL"))
    for key, got, ref, ok in rows:
        if not ok:
            print("   [FAIL] %-22s got=%-15.9g ref=%-15.9g" % (key, got, ref))
            tally[1] += 1
        else:
            tally[0] += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="directory containing example_3d_report.exe")
    a = ap.parse_args()
    exe = Path(a.bin) / "example_3d_report.exe"
    if not exe.exists():
        print("[FAIL] solver not found: %s" % exe)
        return 1
    print("== 3D distributed-load kernel validation ==")
    print("solver: %s" % exe)
    tally = [0, 0]     # pass, fail
    with tempfile.TemporaryDirectory(prefix="fem3d_loads_") as td:
        work = Path(td)
        for case in CASES:
            try:
                disp, force = read_fem(exe, case["dat"], work, case["name"])
                rows = evaluate(case, disp, force)
            except Exception as exc:                       # noqa: BLE001
                print("[FAIL] %-24s exception: %s" % (case["name"], exc))
                tally[1] += 1
                continue
            report_block(case["name"], rows, tally)
        for label, fn in (("trapezoid_superposition", superposition_case),
                          ("truss_selfweight", truss_case)):
            try:
                rows = fn(exe, work)
            except Exception as exc:                       # noqa: BLE001
                print("[FAIL] %-24s exception: %s" % (label, exc))
                tally[1] += 1
            else:
                report_block(label, rows, tally)
    print("Result: %d/%d passed" % (tally[0], tally[0] + tally[1]))
    return 0 if tally[1] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
