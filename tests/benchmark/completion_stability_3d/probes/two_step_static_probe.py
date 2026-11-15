"""诊断：CB3D10 改用 MODE=STATIC + 显式 L0，检查 COMPLETION 分支
改写 m 是否是起点不闭合（RES0 大）的根因。
"""
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = Path(r"D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\cases\cb3d10_completion.dat")
FROZEN = Path(r"D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\cases\cb3d10_completion_cable.txt")
EXE = Path(r"D:\AI-work\FEM-Eigen\build\cable_solver.exe")

# 从冻结 [CABLE] 段取每根索的 L0（第 4 列）
L0 = {}
sec = False
for ln in FROZEN.read_text(errors="replace").splitlines():
    if ln.startswith("["):
        sec = (ln.strip() == "[CABLE]")
        continue
    if sec and ln.strip():
        q = ln.split()
        if len(q) >= 4:
            try:
                L0[int(q[0])] = float(q[3])
            except ValueError:
                pass
print("L0 count =", len(L0))

txt = SRC.read_text(encoding="utf-8")
txt = txt.replace("MODE=COMPLETION", "MODE=STATIC", 1)
# 把 16 条 *EN,9001..9016 行改成带 L0= 覆盖
def add_l0(m):
    eid = int(m.group(1))
    if eid not in L0:
        return m.group(0)
    return m.group(0).rstrip() + " L0=%.15g" % L0[eid]
txt = re.sub(r"^\*EN,(90\d\d),[^\n]*$", add_l0, txt, flags=re.M)
# 展开的 *DO 段（9006~9008）需要单独处理：改写为三行显式 *EN
do_block = re.search(r"\*DO I 9006 9008\n\*EN,I,I-6999,I-7994\n\*ENDDO", txt)
if do_block:
    rows = []
    for i in (9006, 9007, 9008):
        rows.append("*EN,%d,%d,%d L0=%.15g" % (i, i - 6999, i - 7994, L0[i]))
    txt = txt.replace(do_block.group(0), "\n".join(rows), 1)
# 追加稳定性键与增量荷载
txt = txt.replace("\n*CABLE\n", "\n" +
    "STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
    "CTRLNODE=1016 CTRLDOF=3 CTRLTARGET=0.05\n"
    "DS0=0.002 DSMAX=0.01 DSMIN=1.0e-9 PSI=0.5\n"
    "MAXPATHSTEP=2000 STABITER=40 STABTOL=1.0e-7 STRAIN=0\n*CABLE\n", 1)
txt = txt.replace("\n*END\n", "\n*STABILITY_LOAD\n1016 3 -1000000\n*END\n", 1)
(HERE / "s_static.dat").write_text(txt, encoding="utf-8", newline="\n")

p = subprocess.run([str(EXE), "s_static.dat"], cwd=str(HERE),
                   capture_output=True, text=True, timeout=900)
print("exit", p.returncode)
print((p.stdout or "").strip()[:300])
out = (HERE / "s_static_cable.txt").read_text(errors="replace")
m = re.search(r"\[STABILITY\]\n(.*?)(?=\n\[END\])", out, re.S)
if m:
    for ln in m.group(1).splitlines()[:5]:
        print("  ", ln)
    pts = re.findall(r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+)", m.group(1), re.M)
    print("  points", len(pts), "last", pts[-1] if pts else None)
print("NEWTON:", re.search(r"\[NEWTON\]\n(.*?)\nS 1", out, re.S).group(1).strip()
      if re.search(r"\[NEWTON\]\n(.*?)\nS 1", out, re.S) else "n/a")
