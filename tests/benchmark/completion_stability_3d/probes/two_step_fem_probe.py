"""② 阶段 FEM 可行性探针：在 CB3D10 完整双塔三跨斜拉桥上跑
STABASIS=COMPLETION（成桥态起点全过程稳定）。不修改任何产品文件。
"""
import re
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = Path(r"D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\cases\cb3d10_completion.dat")
EXE = Path(r"D:\AI-work\FEM-Eigen\build\cable_solver.exe")

STAB_KEYS = (
    "STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
    "CTRLNODE=1016 CTRLDOF=3 CTRLTARGET=0.05\n"
    "DS0=0.001 DSMAX=0.004 DSMIN=1.0e-9 PSI=0.5\n"
    "MAXPATHSTEP=2000 STABITER=40 STABTOL=1.0e-9 STRAIN=0\n"
)
INC_LOAD = 1.0e6
CTRL_NODE, CTRL_DOF = 1016, 3


def build(tag, stab_keys=STAB_KEYS, inc=True):
    txt = SRC.read_text(encoding="utf-8")
    # 在 *CABLE_BRIDGE 段末尾（*CABLE 之前）插入稳定性键
    txt = txt.replace("\n*CABLE\n", "\n" + stab_keys + "*CABLE\n", 1)
    if inc:
        txt = txt.replace("\n*END\n",
                          f"\n*STABILITY_LOAD\n{CTRL_NODE} {CTRL_DOF} {-INC_LOAD:.12g}\n*END\n", 1)
    dst = HERE / (tag + ".dat")
    dst.write_text(txt, encoding="utf-8", newline="\n")
    return dst


def run(tag):
    p = subprocess.run([str(EXE), tag + ".dat"], cwd=str(HERE),
                       capture_output=True, text=True, timeout=900)
    out = HERE / (tag + "_cable.txt")
    print(f"\n=== {tag} exit={p.returncode} ===")
    print((p.stdout or "").strip()[:400])
    if (p.stderr or "").strip():
        print("STDERR:", (p.stderr or "").strip()[:300])
    if not out.exists():
        print("!! no output file")
        return
    t = out.read_text(errors="replace")
    m = re.search(r"\[STABILITY\]\n(.*?)(?=\n\[END\])", t, re.S)
    if not m:
        print("!! no [STABILITY] section")
        print("sections:", re.findall(r"^\[[A-Z_]+\]", t, re.M)[-8:])
        return
    sec = m.group(1)
    for ln in sec.splitlines()[:5]:
        print("  ", ln)
    pts = re.findall(r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+)", sec, re.M)
    print(f"  points={len(pts)}  first={pts[:2]}  last={pts[-1] if pts else None}")
    print("  BUCKLING ERROR kept:", "ERROR=nonlinear_state_not_supported" in t)
    print("  U(1016,3) =", re.search(r"^1016 3 (\S+)", t, re.M).group(1)
          if re.search(r"^1016 3 (\S+)", t, re.M) else "n/a")


if __name__ == "__main__":
    variants = {
        "t1": ("STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
               "CTRLNODE=1016 CTRLDOF=3 CTRLTARGET=0.05\n"
               "DS0=0.01 DSMAX=0.05 DSMIN=1.0e-9 PSI=0.1\n"
               "MAXPATHSTEP=2000 STABITER=40 STABTOL=1.0e-5 STRAIN=0\n"),
        "t2": ("STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
               "CTRLNODE=1016 CTRLDOF=3 CTRLTARGET=0.05\n"
               "DS0=0.005 DSMAX=0.02 DSMIN=1.0e-9 PSI=0.5\n"
               "MAXPATHSTEP=2000 STABITER=40 STABTOL=1.0e-4 STRAIN=0\n"),
        "t3": ("STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
               "CTRLNODE=1016 CTRLDOF=3 CTRLTARGET=0.05\n"
               "DS0=0.002 DSMAX=0.01 DSMIN=1.0e-9 PSI=0.5\n"
               "MAXPATHSTEP=2000 STABITER=40 STABTOL=1.0e-6 STRAIN=0\n"),
    }
    for tag, keys in variants.items():
        build(tag, keys)
        run(tag)
