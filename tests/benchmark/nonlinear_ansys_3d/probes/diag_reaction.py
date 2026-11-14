#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断：CB3D10 中 ANSYS 支座反力总和与外荷载不闭合的来源。"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = Path(r"D:\AI-work\BridgeStudio\cases\CB3D10_completion_vs_ansys_nl")
SRC = (CASE / "cb3d10_completion_nl.inp").read_text(encoding="ascii")
ANSYS = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
MSYS = r"C:\msys64\mingw64\bin"
os.environ["PATH"] = MSYS + os.pathsep + os.environ.get("PATH", "")


def read_csv(path: Path) -> dict[int, list[float]]:
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        t = line.split()
        if len(t) >= 2:
            try:
                out[int(round(float(t[0])))] = [float(x) for x in t[1:]]
            except ValueError:
                pass
    return out


def run(name: str, transform) -> None:
    work = HERE / ("diag_" + name)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    lines = transform(SRC.splitlines())
    txt = "\n".join(lines).replace("cb3d10_completion_nl", name) + "\n"
    (work / (name + ".inp")).write_text(txt, encoding="ascii")
    p = subprocess.run([str(ANSYS), "-b", "-j", name, "-i", name + ".inp",
                        "-o", name + ".out"], cwd=str(work),
                       capture_output=True, text=True, timeout=600)
    reac = read_csv(work / (name + "_reac.csv"))
    link = read_csv(work / (name + "_link.csv"))
    disp = read_csv(work / (name + "_disp.csv"))
    sup = {1001, 1009, 1024, 1032, 2001, 3001}
    s = sum(reac[n][2] for n in sup if n in reac and len(reac[n]) >= 3)
    allsum = sum(r[2] for r in reac.values() if len(r) >= 3)
    tf = sum(v[0] for k, v in link.items() if k >= 9000) if link else 0.0
    dmax = max((abs(v[2]) for v in disp.values() if len(v) >= 3),
               default=0.0)
    log = (work / (name + ".out")).read_text(
        encoding="utf-8", errors="replace")
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", log)
    nerr = int(m.group(1)) if m else -1
    print("{:<16s} rc={} err={} support_FZ={:>14.1f} all_FZ={:>14.1f} "
          "sum_T={:>12.1f} max|UZ|={:.6g} nlink={}".format(
              name, p.returncode, nerr, s, allsum, tf, dmax, len(link)))


def main() -> int:
    print("external FZ = -18036250.0 N; FEM support sum = +18030800 N")
    run("A_base", lambda L: L)
    run("B_noini", lambda L: [x for x in L if not x.startswith("INISTATE")])
    run("C_kop0", lambda L: [x.replace("KEYOPT,3,3,1", "KEYOPT,3,3,0")
                             for x in L])
    def no_k(lines):
        out = []
        for x in lines:
            if re.match(r"^N,30(1[0-9]|2[0-5]),", x):
                continue                      # drop orientation nodes
            if x.startswith("E,"):
                x = re.sub(r",30\d\d$", "", x)
            out.append(x)
        return out

    run("D_noini_kop0", lambda L: [x.replace("KEYOPT,3,3,1", "KEYOPT,3,3,0")
                                   for x in L if not x.startswith("INISTATE")])
    run("E_noK", no_k)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
