#!/usr/bin/env python3
"""第 24 项 CTest 冻结 cases / 黄金值的唯一生成器（第十九次 ③）。

③ 收口了 ② 登记的缺陷：`solveCompletionState3D()` /
`iterateCableUnstrainedLengths()` 按值接收 `Model3D`，内部会把索 `L0`
改成下料长度（COMPLETION 还把参考构形改成制造线形 `X0+ΔX`），而主链路
仍拿外部原模型去装配 `F0` 并计算 `Fint(u0)`，于是成桥态稳定路径起点
不闭合（CB3D10 实测 `RES0=3.15 MN ≈ F0_NORM`、首步 λ 变负）。
修复后由 `CompletionStateResult3D::basisModel` /
`CableLengthResult3D::basisModel` 回传同源模型，本项把它冻结成回归。

四个 case：

| case | MODE | STABASIS | 用途 |
|---|---|---|---|
| `cb3d10_completion_basis` | COMPLETION | COMPLETION | 修复主体：制造线形基准、下料 L0 同源 |
| `cb3d10_static_basis` | STATIC | COMPLETION | 独立物理锚点：设计线形 + 冻结 `[CABLE]` L0 |
| `cb3d10_default_basis` | COMPLETION | (STABILITY=0) | 默认零污染：与历史 CB3D10 逐行一致 |
| `cb3d13_erection_basis` | COMPLETION（悬索桥） | COMPLETION | ③ 施工阶段链路接入核对：`[STABILITY]` 与 5 个 `[ERECTION*]` 段共存、互不污染 |

用法：
    python build_reference_basis.py <cable_solver.exe> <cases_dir> [reference_dir]

正常回归**不需要**重跑本脚本（会覆写冻结输出）。
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "tests/benchmark/nonlinear_ansys_3d/cases/cb3d10_completion.dat"
FROZEN_BASELINE = ROOT / "tests/benchmark/nonlinear_ansys_3d/cases/cb3d10_completion_cable.txt"
SRC_ERECT = ROOT / "tests/benchmark/suspension_erection_3d/cases/cb3d13_balanced_erection.dat"
FROZEN_ERECT = ROOT / "tests/benchmark/suspension_erection_3d/cases/cb3d13_balanced_erection_cable.txt"
# 编号口径：本项不是一个新的 BridgeStudio 算例，而是对既有 CB3D10 /
# CB3D13 冻结输出的离线回归，因此不占用 CB3D22（该编号属于 BS 侧
# `cases/CB3D22_limit_point_vs_ansys/` 真实桥梁极限点算例）。
REF_NAME = "completion_basis_3d.json"

# 与 ② 阶段 CB3D10 诊断探针同一组弧长参数；MAXPATHSTEP 收窄到 400，
# 实测 129 个平衡点即可到达 CTRLTARGET=0.05，CTest 不需要 2000 步。
STAB_KEYS = (
    "STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
    "CTRLNODE=1016 CTRLDOF=3 CTRLTARGET=0.05\n"
    "DS0=0.005 DSMAX=0.02 DSMIN=1.0e-9 PSI=0.5\n"
    "MAXPATHSTEP=400 STABITER=40 STABTOL=1.0e-4 STRAIN=0\n"
)
INC_SECTION = "*STABILITY_LOAD\n1016 3 -1000000\n"
CTRL_NODE, CTRL_DOF, PINC = 1016, 3, 1.0e6

# 悬索桥施工阶段算例（CB3D13）：控制点取主跨跨中梁节点 1014，
# 目标位移 0.02 m，其余弧长参数与 CB3D10 同组。
ERECT_STAB_KEYS = (
    "STABILITY=1 STABASIS=COMPLETION PATHMETHOD=SPHERE "
    "CTRLNODE=1014 CTRLDOF=3 CTRLTARGET=0.02\n"
    "DS0=0.005 DSMAX=0.02 DSMIN=1.0e-9 PSI=0.5\n"
    "MAXPATHSTEP=400 STABITER=40 STABTOL=1.0e-4 STRAIN=0\n"
)
ERECT_INC_SECTION = "*STABILITY_LOAD\n1014 3 -1000000\n"
ERECT_CTRL = (1014, 3)

CASES = (
    "cb3d10_completion_basis",
    "cb3d10_static_basis",
    "cb3d10_default_basis",
    "cb3d13_erection_basis",
)
ERECT_SECTIONS = ("ERECTION", "ERECTION_SHAPE", "ERECTION_CHAIN",
                  "ERECTION_LOAD", "ERECTION_FORCE")


def ctrl_of(name: str) -> tuple[int, int]:
    return ERECT_CTRL if name == "cb3d13_erection_basis" else (CTRL_NODE, CTRL_DOF)


def normalized(txt: str) -> str:
    """耗时不是物理量；其余逐字符保留。"""
    return re.sub(r"elapsed_ms=[0-9.eE+-]+", "elapsed_ms=*", txt)


def baseline_L0() -> dict[int, float]:
    """从冻结 [CABLE] 段取 16 根索的下料长度（第 4 列）。"""
    out: dict[int, float] = {}
    sec = False
    for ln in FROZEN_BASELINE.read_text(errors="replace").splitlines():
        if ln.startswith("["):
            sec = ln.strip() == "[CABLE]"
            continue
        if not sec or not ln.strip():
            continue
        q = ln.split()
        if len(q) >= 4:
            try:
                out[int(q[0])] = float(q[3])
            except ValueError:
                pass
    if len(out) != 16:
        raise RuntimeError(f"冻结 [CABLE] 段应有 16 根索，实得 {len(out)}")
    return out


def make_dat(name: str) -> str:
    if name == "cb3d13_erection_basis":
        # 稳定性键必须落在 *CABLE_BRIDGE 段内（*SUSP_BRIDGE 之前），
        # 不能插到 *ERECTION 段后的 *CABLE 前——那样会被当成施工阶段参数。
        er = SRC_ERECT.read_text(encoding="utf-8")
        er = er.replace("\n*SUSP_BRIDGE\n",
                        "\n" + ERECT_STAB_KEYS + "*SUSP_BRIDGE\n", 1)
        er = er.replace("\n*END\n", "\n" + ERECT_INC_SECTION + "*END\n", 1)
        return er
    txt = SRC.read_text(encoding="utf-8")
    if name == "cb3d10_static_basis":
        txt = txt.replace("MODE=COMPLETION", "MODE=STATIC", 1)
        L0 = baseline_L0()

        def add_l0(m: re.Match[str]) -> str:
            eid = int(m.group(1))
            return (m.group(0).rstrip() + f" L0={L0[eid]:.15g}") if eid in L0 else m.group(0)

        txt = re.sub(r"^\*EN,(90\d\d),[^\n]*$", add_l0, txt, flags=re.M)
        do_block = re.search(
            r"\*DO I 9006 9008\n\*EN,I,I-6999,I-7994\n\*ENDDO", txt)
        if do_block:
            rows = ["*EN,%d,%d,%d L0=%.15g" % (i, i - 6999, i - 7994, L0[i])
                    for i in (9006, 9007, 9008)]
            txt = txt.replace(do_block.group(0), "\n".join(rows), 1)
    if name != "cb3d10_default_basis":
        txt = txt.replace("\n*CABLE\n", "\n" + STAB_KEYS + "*CABLE\n", 1)
        txt = txt.replace("\n*END\n", "\n" + INC_SECTION + "*END\n", 1)
    return txt


def parse_stability(txt: str) -> dict:
    m = re.search(r"\[STABILITY\]\n(.*?)(?=\n\[END\])", txt, re.S)
    if not m:
        return {"present": False}
    sec = m.group(1)
    head = {}
    for ln in sec.splitlines():
        if "=" in ln and not ln.startswith("I "):
            for k, v in re.findall(r"(\w+)=([^\s]+)", ln):
                head.setdefault(k, v)
    state = {}
    sm = re.search(r"^STATE (.+)$", sec, re.M)
    if sm:
        for k, v in re.findall(r"(\w+)=([^\s]+)", sm.group(1)):
            state[k] = float(v)
    pts = re.findall(r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+)", sec, re.M)
    return {
        "present": True,
        "head": head,
        "state": state,
        "points": [[int(s), float(lam), float(ctrl), float(res)]
                   for s, lam, ctrl, res in pts],
    }


def section_body(txt: str, tag: str) -> str | None:
    m = re.search(r"^\[" + tag + r"\]\n(.*?)(?=^\[[A-Z_]+\]|\Z)",
                  txt, re.S | re.M)
    return m.group(1) if m else None


def metrics(txt: str, ctrl: tuple[int, int] = (CTRL_NODE, CTRL_DOF)) -> dict:
    st = parse_stability(txt)
    reslast = re.search(r"RESLAST=(\S+)", txt)
    u_ctrl = re.search(r"^%d %d (\S+)" % ctrl, txt, re.M)
    sections = re.findall(r"^\[([A-Z_]+)\]", txt, re.M)
    length_kind = re.search(r"^\[LENGTH\]\nKIND=(\S+) STATUS=(\S+)", txt, re.M)
    return {
        "stability": st,
        "reslast": float(reslast.group(1)) if reslast else None,
        "u_ctrl": float(u_ctrl.group(1)) if u_ctrl else None,
        "sections": sections,
        "length_kind": list(length_kind.groups()) if length_kind else None,
        "buckling_error_kept": "ERROR=nonlinear_state_not_supported" in txt,
        "finc_norm_expected": PINC,
        "erection_sections": {tag: section_body(txt, tag) for tag in ERECT_SECTIONS},
    }


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    exe = Path(sys.argv[1]).resolve()
    cases = Path(sys.argv[2]).resolve()
    refdir = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else HERE / "reference"
    cases.mkdir(parents=True, exist_ok=True)
    refdir.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(FROZEN_BASELINE, cases / "cb3d10_baseline_cable.txt")
    shutil.copyfile(FROZEN_ERECT, cases / "cb3d13_baseline_cable.txt")
    ref: dict = {
        "case": "CB3D10_BASIS_FIX",
        "cases_covered": ["CB3D10", "CB3D13"],
        "title": "第十九次 ③：MODE=COMPLETION 成桥态基准同源修复 + 施工阶段链路共存核对",
        "source_dat": "tests/benchmark/nonlinear_ansys_3d/cases/cb3d10_completion.dat",
        "baseline_txt": "cb3d10_baseline_cable.txt",
        "erection_source_dat": "tests/benchmark/suspension_erection_3d/cases/cb3d13_balanced_erection.dat",
        "erection_baseline_txt": "cb3d13_baseline_cable.txt",
        "erection_sections": list(ERECT_SECTIONS),
        "stab_keys": STAB_KEYS.strip(),
        "erect_stab_keys": ERECT_STAB_KEYS.strip(),
        "incremental_load_section": INC_SECTION.strip(),
        "erect_incremental_load_section": ERECT_INC_SECTION.strip(),
        "ctrl": [CTRL_NODE, CTRL_DOF],
        "erect_ctrl": list(ERECT_CTRL),
        "cases": {},
    }
    for name in CASES:
        dat = cases / f"{name}.dat"
        dat.write_text(make_dat(name), encoding="utf-8", newline="\n")
        p = subprocess.run([str(exe), dat.name], cwd=str(cases),
                           capture_output=True, text=True, errors="replace", timeout=900)
        out = cases / f"{name}_cable.txt"
        if not out.exists():
            raise RuntimeError(f"{name} 未产出结果文件：{p.stdout} {p.stderr}")
        txt = out.read_text(errors="replace")
        norm = normalized(txt)
        ref["cases"][name] = {
            "exit_code": p.returncode,
            "dat": dat.name,
            "txt": out.name,
            "sha256_normalized": hashlib.sha256(norm.encode("utf-8")).hexdigest(),
            "n_lines": len(norm.splitlines()),
            "ctrl": list(ctrl_of(name)),
            "metrics": metrics(txt, ctrl_of(name)),
        }
        st = ref["cases"][name]["metrics"]["stability"]
        print(f"{name}: exit={p.returncode} lines={len(norm.splitlines())} "
              f"stability={st.get('present')} "
              f"RES0={st.get('state', {}).get('RES0')} "
              f"NPOINT={len(st.get('points', []))}")
    target = refdir / REF_NAME
    target.write_text(json.dumps(ref, ensure_ascii=False, indent=2, sort_keys=True),
                      encoding="utf-8")
    print("reference ->", target)
    return 0


if __name__ == "__main__":
    sys.exit(main())
