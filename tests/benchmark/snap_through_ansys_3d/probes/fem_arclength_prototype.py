"""Python prototype of the spherical arc-length (Crisfield) path follower.

Same augmented system that will be implemented in stability_path_3d.hpp:
    R(u,lambda) = lambda*Fref - Fint(u) = 0
    du^T du + psi^2*dlambda^2 = ds^2
    K*duR = R,  K*duF = Fref,  du = duR + dlambda*duF
    -> quadratic a*dlambda^2 + b*dlambda + c = 0
Root selection keeps the increment pointing in the same direction as the
previous converged increment (prevents backtracking at limit points).
"""
import math

E = 2.1e11
A = 1.0e-3
B = 2.0
H = 0.5
L0 = math.hypot(B, H)
PREF = 1.0e6
ND = 9
CDOF = [0, 1, 2, 5, 6, 7, 8]
FDOF = [3, 4]

P = [(-B, 0.0, 0.0), (0.0, H, 0.0), (B, 0.0, 0.0)]
ELEM = [(0, 1), (1, 2)]

FREF = [0.0] * ND
FREF[4] = -PREF


def assemble(u):
    Fint = [0.0] * ND
    K = [[0.0] * ND for _ in range(ND)]
    for (i, j) in ELEM:
        l = [P[j][k] + u[3 * j + k] - P[i][k] - u[3 * i + k] for k in range(3)]
        Ln = math.sqrt(sum(x * x for x in l))
        c = [x / Ln for x in l]
        N = E * A / L0 * (Ln - L0)
        for k in range(3):
            Fint[3 * i + k] += -N * c[k]
            Fint[3 * j + k] += N * c[k]
        G = [[(E * A / L0) * c[p] * c[q]
              + (N / Ln) * ((1.0 if p == q else 0.0) - c[p] * c[q])
              for q in range(3)] for p in range(3)]
        for p in range(3):
            for q in range(3):
                K[3 * i + p][3 * i + q] += G[p][q]
                K[3 * j + p][3 * j + q] += G[p][q]
                K[3 * i + p][3 * j + q] += -G[p][q]
                K[3 * j + p][3 * i + q] += -G[p][q]
    return Fint, K


def solve(K, rhs):
    """Gaussian elimination with partial pivoting; Dirichlet rows pre-applied."""
    n = len(rhs)
    M = [row[:] + [rhs[r]] for r, row in enumerate(K)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-30:
            raise RuntimeError("singular tangent at column %d" % col)
        M[col], M[piv] = M[piv], M[col]
        for r in range(col + 1, n):
            f = M[r][col] / M[col][col]
            if f:
                for c in range(col, n + 1):
                    M[r][c] -= f * M[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][c] * x[c] for c in range(r + 1, n))) / M[r][r]
    return x


def apply_dirichlet(K, rhs):
    K = [row[:] for row in K]
    rhs = rhs[:]
    for d in CDOF:
        for c in range(ND):
            K[d][c] = 0.0
        K[d][d] = 1.0
        rhs[d] = 0.0
    return K, rhs


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def p_of_v(v):
    L = math.hypot(B, H - v)
    N = E * A * (L - L0) / L0
    return -2.0 * N * (H - v) / L


def run(ds0=0.02, psi=1.0, vmax=1.0, maxsteps=400, verbose=False):
    u = [0.0] * ND
    lam = 0.0
    ds = ds0
    prev_du = None
    curve = [(0.0, 0.0)]
    for step in range(maxsteps):
        # ---- predictor ----
        Fint, Kt = assemble(u)
        Kd, _ = apply_dirichlet(Kt, [0.0] * ND)
        duF = solve(Kd, FREF)
        den = math.sqrt(dot(duF, duF) + psi * psi)
        sign = 1.0
        if prev_du is not None:
            # Crisfield: sign of lambda increment from previous tangent
            sign = 1.0 if dot(prev_du, duF) > 0 else -1.0
        dlam = sign * ds / den
        du0 = [dlam * x for x in duF]
        # trial state
        u_t = [u[k] + du0[k] for k in range(ND)]
        lam_t = lam + dlam
        du_acc = du0[:]
        dlam_acc = dlam
        # ---- corrector ----
        converged = False
        for it in range(30):
            Fint, Kt = assemble(u_t)
            R = [lam_t * FREF[k] - Fint[k] for k in range(ND)]
            Kd, Rd = apply_dirichlet(Kt, R)
            duR = solve(Kd, Rd)
            duF = solve(Kd, FREF)
            base = [du_acc[k] + duR[k] for k in range(ND)]
            a = dot(duF, duF) + psi * psi
            b = 2.0 * dot(base, duF) + 2.0 * psi * psi * dlam_acc
            c = dot(base, base) + psi * psi * dlam_acc * dlam_acc - ds * ds
            disc = b * b - 4.0 * a * c
            if disc < 0.0:
                disc = 0.0
            sq = math.sqrt(disc)
            r1 = (-b + sq) / (2.0 * a)
            r2 = (-b - sq) / (2.0 * a)
            cand1 = [base[k] + r1 * duF[k] for k in range(ND)]
            cand2 = [base[k] + r2 * duF[k] for k in range(ND)]
            # choose root with larger projection onto previous converged increment
            ref = prev_du if prev_du is not None else du_acc
            s1, s2 = dot(cand1, ref), dot(cand2, ref)
            dd = r1 if s1 >= s2 else r2
            du_new = cand1 if s1 >= s2 else cand2
            du = [du_new[k] - du_acc[k] for k in range(ND)]
            u_t = [u_t[k] + du[k] for k in range(ND)]
            lam_t += dd
            du_acc = du_new
            dlam_acc += dd
            if math.sqrt(dot(Rd, Rd)) < 1e-7 * PREF:
                converged = True
                break
        if not converged:
            if verbose:
                print("  step %d NOT CONVERGED, halving ds" % step)
            ds *= 0.5
            if ds < 1e-7:
                break
            continue
        u, lam = u_t, lam_t
        prev_du = du_acc
        v = -u[4]
        curve.append((v, lam))
        if verbose and step % 25 == 0:
            print("  step %3d  v=%+.6f  lambda=%+.6f  P_ex=%.4f  dP=%+.3g  ds=%.2e"
                  % (step, v, lam, p_of_v(v), lam * PREF - p_of_v(v), ds))
        if v >= vmax:
            break
        # mild step adaptation
        ds = min(ds0 * 2.0, max(ds0 * 0.25, ds * 1.05))
    return curve


if __name__ == "__main__":
    # analytic limit point by golden section
    lo, hi = 0.0, H
    for _ in range(200):
        m1, m2 = lo + (hi - lo) * 0.381966, lo + (hi - lo) * 0.618034
        if p_of_v(m1) < p_of_v(m2):
            lo = m1
        else:
            hi = m2
    vcr, pcr = 0.5 * (lo + hi), p_of_v(0.5 * (lo + hi))
    print("analytic limit: v_cr=%.9f  P_cr=%.6f N  lambda_cr=%.9f"
          % (vcr, pcr, pcr / PREF))
    print("analytic P(h)=%.6f  P(2h)=%.6f" % (p_of_v(H), p_of_v(2 * H)))
    for ds in (0.05, 0.02, 0.01):
        c = run(ds0=ds, vmax=1.0, verbose=(ds == 0.02))
        vmax = max(v for v, _ in c)
        lmax = max(l for _, l in c)
        lmin = min(l for _, l in c)
        errs = [abs(l * PREF - p_of_v(v)) for v, l in c if 0.0 <= v <= 2 * H]
        print("ds=%-6.3g npoints=%4d  v_max=%.6f  lambda=[%+.6f,%+.6f]  "
              "worst|dP|=%.4g N (%.3g%%)"
              % (ds, len(c), vmax, lmin, lmax, max(errs),
                 100 * max(errs) / pcr))
