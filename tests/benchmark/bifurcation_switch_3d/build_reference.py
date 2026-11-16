"""CB3D28 分岔 benchmark 的独立闭式锚点生成器。

只读冻结输入/输出算 sha256；所有物理量由本文件独立解析推导，
不读取 FEM CSV 里的任何物理数值，保证锚点与实现解耦。

解析模型（高拱 von Mises 两杆桁架，B=1 H=3）：
  L0 = hypot(B,H),  L(v) = hypot(B,H-v),  h = H-v
  N = EA(L-L0)/L0,  P(v) = -2 N h / L

两个临界条件（各自可由切线刚度为零严格推出）：
  对称极值点：  L^3 = B^2 * L0
  反对称分岔：  L^3 = L0 * h^2   →  (B^2+h^2)^1.5 = L0 h^2
后者在 h>0 上有两个根；较小 v（较大 h）的根先发生，本几何下
v_bif ≈ 0.446 < v_sym ≈ 1.926，故反对称分岔先于对称极限点。
"""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE / "cases" / "cb3d28_bifurcation_switch.snap"
FROZEN = HERE / "cases" / "cb3d28_bifurcation_switch_fem.csv"
OUT = HERE / "reference" / "cb3d28_bifurcation_switch.json"

E = 2.1e11
A = 1.0e-3
B = 1.0
H = 3.0
PREF = 2.0e6
EA = E * A
L0 = math.hypot(B, H)


def p_of_v(v):
    h = H - v
    L = math.hypot(B, h)
    N = EA * (L - L0) / L0
    return -2.0 * N * h / L


def f_bif(h):
    return (B * B + h * h) ** 1.5 - L0 * h * h


def positive_roots(f, hi, n=4000):
    """扫描 + 二分，返回 f 在 (0,hi] 上的全部单根（升序）。"""
    roots = []
    step = hi / n
    xs = [step * i for i in range(1, n + 1)]
    prev_x, prev_y = xs[0], f(xs[0])
    for x in xs[1:]:
        y = f(x)
        if prev_y == 0.0:
            roots.append(prev_x)
        elif prev_y * y < 0.0:
            a, b = prev_x, x
            for _ in range(200):
                m = 0.5 * (a + b)
                if f(a) * f(m) <= 0.0:
                    b = m
                else:
                    a = m
            roots.append(0.5 * (a + b))
        prev_x, prev_y = x, y
    return roots


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    roots = positive_roots(f_bif, H)
    if len(roots) < 2:
        raise RuntimeError(f"expected 2 bifurcation roots, got {roots}")
    # 较大 h → 较小 v → 路径上先遇到；取它作为被检测的分岔点。
    h_bif = max(roots)
    v_bif = H - h_bif
    p_bif = p_of_v(v_bif)

    L_sym = (B * B * L0) ** (1.0 / 3.0)
    h_sym = math.sqrt(max(0.0, L_sym * L_sym - B * B))
    v_sym = H - h_sym
    p_sym = p_of_v(v_sym)

    ref = {
        "schema": "cb3d28_bifurcation_reference_v1",
        "sha256": {
            "input": sha256(CASE),
            "frozen_csv": sha256(FROZEN),
        },
        "params": {
            "E": E, "A": A, "B": B, "H": H, "PREF": PREF,
            "EA": EA, "L0": L0,
        },
        "closed_form": {
            "bifurcation_roots_h": roots,
            "h_bif": h_bif,
            "v_bif": v_bif,
            "P_bif": p_bif,
            "lambda_bif": p_bif / PREF,
            "h_sym": h_sym,
            "v_sym": v_sym,
            "P_sym": p_sym,
            "lambda_sym": p_sym / PREF,
            "bif_before_limit_point": v_bif < v_sym,
            "condition_bif": "(B^2+h^2)^1.5 = L0*h^2",
            "condition_sym": "L^3 = B^2*L0",
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(ref, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"  bifurcation roots h = {[round(r, 10) for r in roots]}")
    print(f"  v_bif={v_bif:.10f}  lambda_bif={p_bif / PREF:.10f}")
    print(f"  v_sym={v_sym:.10f}  lambda_sym={p_sym / PREF:.10f}")


if __name__ == "__main__":
    main()
