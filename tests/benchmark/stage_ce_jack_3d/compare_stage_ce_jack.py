#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 42 项 CTest：分级顶推 CE 施工阶段 benchmark（.dat 解析 → 内核端到端）。

离线（无需 ANSYS / 网络 / 许可）：
  * 冻结分级顶推 .dat 的 sha256 与 reference 对齐；
  * 当前 example_3d_stage_nl.exe 重放：三阶段全部收敛，*STAGE CE= 的
    增量 δ 经非线性 stage 内核累加为累计 Δ（UX(C)-UX(T) 逐阶段命中）；
  * 分隔符兼容性：把 CE=1:0.01 改成 CE=1;0.01 的等价写法仍命中（';' 分隔）；
  * 拒绝路径（动态生成）：CE= 引用不存在的 CE id 被显式拒绝；
  * 源码契约：内核两处复制 m.ces 并应用 ceConstants、constrainDisconnected3D
    对 CE 引用节点豁免、*STAGE 解析 CE= 关键分支在位。

用法：
    python compare_stage_ce_jack.py --bin <dir> \
        --cases <cases_dir> --reference <reference_dir>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FE_ROOT = HERE.parents[2]
RESULTS: list[tuple[bool, str, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((bool(ok), name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_nlstage(text: str) -> dict[int, dict]:
    """解析 example_3d_stage_nl 的 .nlstage.txt：每阶段 [U] 的 (node,dof)->u。"""
    stages: dict[int, dict] = {}
    cur: int | None = None
    section = ""
    for line in text.splitlines():
        m = re.match(r"^STAGE\s+(\d+)\s+converged=(\d+)", line)
        if m:
            cur = int(m.group(1))
            stages[cur] = {"converged": int(m.group(2)) == 1, "U": {}}
            section = ""
            continue
        if line == "[U]":
            section = "U"
            continue
        if line.startswith("[") or line.startswith("#") or line.startswith("STAGE="):
            if line.startswith("[") and not line.startswith("[U]"):
                section = ""
            continue
        if section == "U" and cur is not None:
            parts = line.split()
            if len(parts) == 3:
                try:
                    stages[cur]["U"][(int(parts[0]), int(parts[1]))] = float(parts[2])
                except ValueError:
                    pass
    return stages


def slip(stage: dict) -> float:
    u = stage["U"]
    return u.get((3, 1), 0.0) - u.get((2, 1), 0.0)


def run_case(exe: Path, dat_text: str, workdir: Path, tag: str):
    dat = workdir / f"{tag}.dat"
    dat.write_text(dat_text, encoding="utf-8")
    r = subprocess.run([str(exe), str(dat)], cwd=str(workdir),
                       capture_output=True, text=True, timeout=120)
    out = dat.with_suffix(dat.suffix + ".nlstage.txt")
    return r, (out.read_text(encoding="utf-8", errors="replace") if out.exists() else "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()

    bin_dir = Path(args.bin)
    cases = Path(args.cases)
    ref_dir = Path(args.reference)
    ref = json.loads((ref_dir / "stage_ce_jack.json").read_text(encoding="utf-8"))
    exe = bin_dir / "example_3d_stage_nl.exe"
    check("example_3d_stage_nl.exe 存在", exe.exists(), str(exe))

    dat_path = cases / ref["case"]["file"]
    base = dat_path.read_text(encoding="utf-8")
    check("冻结 .dat sha256 对齐",
          sha256_file(dat_path) == ref["case"]["sha256"],
          sha256_file(dat_path))

    atol = float(ref["slip"]["atol"])
    targets = [(int(s["id"]), float(s["target_cumulative_slip"])) for s in ref["stages"]]

    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        # 主算例重放
        r, text = run_case(exe, base, wd, "jack")
        check("重放 returncode=0", r.returncode == 0, f"rc={r.returncode}")
        st = parse_nlstage(text)
        check("阶段数=3", len(st) == 3, f"got={len(st)}")
        for sid, target in targets:
            s = st.get(sid)
            check(f"阶段 {sid} 收敛", bool(s and s["converged"]))
            if s:
                got = slip(s)
                check(f"阶段 {sid} 累计滑移=Δ {target}", abs(got - target) <= atol,
                      f"got={got:.6g} target={target} atol={atol}")

        # 分隔符兼容：';' 是多组之间的分隔符（':' 才是 id 与值的分隔）。
        # 在每组后追加尾随 ';'：归一化后产生空段，应被跳过而结果不变。
        alt = base.replace("CE=1:0.01", "CE=1:0.01;")
        check("尾随 ';' 写法已构造", alt != base)
        r2, text2 = run_case(exe, alt, wd, "jack_sep")
        check("尾随 ';' 写法重放 returncode=0", r2.returncode == 0, f"rc={r2.returncode}")
        st2 = parse_nlstage(text2)
        if 3 in st2:
            check("尾随 ';' 写法阶段3累计滑移=0.03",
                  abs(slip(st2[3]) - 0.03) <= atol,
                  f"got={slip(st2[3]):.6g}")

        # 拒绝路径：CE= 引用不存在的 CE id 999
        bad = base.replace("CE=1:0.01", "CE=999:0.01")
        check("不存在 CE id 负例已构造", bad != base)
        r3, text3 = run_case(exe, bad, wd, "jack_bad")
        check("不存在 CE id 被拒绝（rc≠0）", r3.returncode != 0,
              f"rc={r3.returncode} err={r3.stderr.strip()[:120]}")

    # 源码契约（避免 benchmark 通过但内核契约被悄悄改掉）
    hsrc = (FE_ROOT / "fem_3d.hpp").read_text(encoding="utf-8", errors="replace")
    nsrc = (FE_ROOT / "fem_3d_nl.hpp").read_text(encoding="utf-8", errors="replace")
    check("fem_3d.hpp Stage3D 有 ceConstants", "ceConstants" in hsrc)
    check("fem_3d.hpp *STAGE 解析 CE=", 'kv.rfind("CE=",0)==0' in hsrc)
    check("fem_3d.hpp 线性 stage 复制 m.ces", "mc.ces = m.ces;" in hsrc)
    check("fem_3d.hpp constrainDisconnected3D 豁免 CE 节点", "ceRef.count(id)" in hsrc)
    check("fem_3d_nl.hpp 非线性 stage 累加 δ 为 Δ", "ceB[cid] += d" in nsrc)

    ok = all(x[0] for x in RESULTS)
    print(f"\n结果: {'全部通过' if ok else str(sum(not x[0] for x in RESULTS)) + ' 项失败'}"
          f"  ({sum(x[0] for x in RESULTS)}/{len(RESULTS)})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
