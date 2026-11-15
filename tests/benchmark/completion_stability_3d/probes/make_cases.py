#!/usr/bin/env python3
"""a 阶段 ①：生成第 22 项 CTest 的冻结 cases（一次性，不进 CTest）。

V 形双索单自由度体系：两根 CABLE 从 (±B, H) 连到低点 (0,0)，
低点只留 dof 2 (y)。参考构形带预拉（L0 < L_ref），成桥态在恒载 F0
下低点下移 V0。独立闭式解（不读求解器）：

    L(v)  = hypot(B, H+v)
    N(v)  = EA/L0 * (L(v) - L0)
    P(v)  = 2 N(v) (H+v) / L(v)
    F0    = P(V0)                     （成桥态平衡）
    T*    = N(V0)                     （MODE=TENSION 的目标索力）
    λ(v)  = [P(v) - F0] / PINC        （成桥态基准的增量荷载因子）

用法： python make_cases.py <exe> <out_dir>
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

B, H, E, A = 2.0, 0.6, 2.1e11, 1.0e-3
EA = E * A
L_REF = math.hypot(B, H)
L0 = 0.995 * L_REF
V0 = 0.10
PINC = 5.0e5


def N_of_v(v: float) -> float:
    return EA / L0 * (math.hypot(B, H + v) - L0)


def P_of_v(v: float) -> float:
    lv = math.hypot(B, H + v)
    return 2.0 * N_of_v(v) * (H + v) / lv


F0 = P_of_v(V0)
TSTAR = N_of_v(V0)

STAB_KEYS = (
    "STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
    "CTRLNODE=2 CTRLDOF=2 CTRLTARGET=0.05\n"
    "DS0=0.002 DSMAX=0.006 DSMIN=1.0e-9 PSI=0.5\n"
    "MAXPATHSTEP=1500 STABITER=40 STABTOL=1.0e-9 STRAIN=0\n"
)

HEAD = f"""# 第十九次 ① / 第 22 项 CTest：V 形双索单自由度成桥态全过程稳定
# B={B} H={H} E={E:.6g} A={A:.6g} L0={L0:.17g} V0={V0} PINC={PINC:.6g}
*NODE
1 -{B:.10g} {H:.10g} 0
2 0 0 0
3 {B:.10g} {H:.10g} 0
*ET
1 CABLE
*MP
1 E={E:.10g}
*R
1 A={A:.10g} L0={L0:.17g}
*TYPE,1
*MAT,1
*REAL,1
*EN,1,1,2
*EN,2,2,3
*BOUNDARY
1 1-3 0.0
3 1-3 0.0
2 1 0.0
2 3 0.0
*CLOAD
2 2 -{F0:.17g}
"""


def bridge(mode: str, stab: str) -> str:
    cab = f"{TSTAR:.17g}" if mode in ("TENSION", "COMPLETION") else "0"
    return (
        "*CABLE_BRIDGE\n"
        f"MODE={mode}\n"
        "METHOD=MINENERGY\n"
        "NLSTEPS=5\nNEWTON=30\nMAXITER=20\nTOL=1e-10\nTOLT=1e-09\nVERBOSE=0\n"
        + stab
        + "*CABLE\n"
        f"1 0 {cab}\n2 0 {cab}\n"
    )


def inc_load() -> str:
    return f"*STABILITY_LOAD\n2 2 -{PINC:.17g}\n"


CASES: dict[str, str] = {
    # 成桥态基准 / 两种闭环
    "vshape_static": HEAD + bridge("STATIC", STAB_KEYS) + inc_load() + "*END\n",
    "vshape_tension": HEAD + bridge("TENSION", STAB_KEYS) + inc_load() + "*END\n",
    # 负例：无增量荷载 / 未知基准 / 成桥态未收敛
    "vshape_no_load": HEAD + bridge("STATIC", STAB_KEYS) + "*END\n",
    "vshape_bad_basis": HEAD + bridge(
        "STATIC",
        STAB_KEYS.replace("STABASIS=COMPLETION", "STABASIS=FOO")
    ) + inc_load() + "*END\n",
    "vshape_unconverged": HEAD + bridge("COMPLETION", STAB_KEYS)
                          + inc_load() + "*END\n",
    # 默认零污染对照：有 *STABILITY_LOAD、但完全没有 STABILITY 键。
    # 成桥态位移必须仍是 -V0 —— 证明增量荷载段没有混进 m.loads。
    "vshape_default": HEAD + bridge("STATIC", "") + inc_load() + "*END\n",
}


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: python make_cases.py <cable_solver.exe> <out_dir>")
        return 2
    exe, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    print(f"L_ref={L_REF:.12g} L0={L0:.12g} T*={TSTAR:.12g} F0={F0:.12g}")
    print(f"N_ref={N_of_v(0.0):.9g} N (pre-tension > 0)")
    print(f"P(V0)={P_of_v(V0):.9g} == F0")
    for stem, txt in CASES.items():
        (out / f"{stem}.dat").write_text(txt, encoding="utf-8")
        p = subprocess.run([str(exe), f"{stem}.dat"], cwd=str(out),
                           capture_output=True, text=True, errors="replace")
        produced = sorted(
            q.name for q in out.glob(f"{stem}_*.txt"))
        print(f"{stem:22s} exit={p.returncode:<3d} -> {produced}")
        if (p.stderr or "").strip():
            print("   stderr:", (p.stderr or "").strip()[:200])
        for q in out.glob(f"{stem}_*.txt"):
            body = q.read_text(encoding="utf-8", errors="replace")
            for ln in body.splitlines():
                if ln.startswith(("KIND=", "STATE ", "LAMBDA_FINAL=", "ERROR=")):
                    print("   ", ln[:130])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
