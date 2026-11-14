#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G 阶段可行性试验：小规模 COMPLETION 成桥态，确认输出段可用。"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, r"D:\AI-work\BridgeStudio")
os.environ["PATH"] = (r"C:\msys64\mingw64\bin" + os.pathsep
                      + os.environ.get("PATH", ""))

from cable_stayed_3d import build_cable_stayed_3d  # noqa: E402
from fem_eigen_cable import solve_cable_bridge  # noqa: E402


def main() -> int:
    model, info = build_cable_stayed_3d(
        main_span=80.0, side_span=30.0, tower_height=25.0,
        n_cable=4, panel=10.0,
        B=12.0, h=2.5, tf=0.28, tw=0.45,
        tower_b=3.0, tower_h=3.0,
        q_deck=50e3, cable_area_mm2=3000.0,
        cable_dens=0.0,                 # 关闭索自重垂度 / Ernst，便于 ANSYS 对标
        mode="COMPLETION", method="DECK",
    )
    print("nodes {}  elements {}  cables {}".format(
        len(model.nodes), len(model.elements), info["n_cable_actual"]))
    print("node id range {} .. {}".format(
        min(n.id for n in model.nodes), max(n.id for n in model.nodes)))
    print("elem id range {} .. {}".format(
        min(e.id for e in model.elements), max(e.id for e in model.elements)))
    print("etypes", [(t.id, t.kind, t.mat, t.real) for t in model.etypes])
    print("reals", [(r.id, round(r.A, 6), round(r.IY, 6), round(r.W, 6))
                    for r in model.reals])
    print("boundaries", [(b.node, b.dofs, b.value) for b in model.boundaries])
    print("loads", len(model.loads), "sum FZ",
          sum(ld.value for ld in model.loads if ld.dof == 3))
    (HERE / "trial.dat").write_text(model.to_dat(), encoding="utf-8")
    dat_head = (HERE / "trial.dat").read_text(encoding="utf-8").splitlines()
    print("--- .dat CABLE_BRIDGE section ---")
    for i, line in enumerate(dat_head):
        if "*CABLE_BRIDGE" in line:
            print("\n".join(dat_head[i:i + 12]))
            break

    res = solve_cable_bridge(model, out_dir=HERE, dat_name="trial.dat",
                             timeout=600)
    md = res["meta"]
    print("--- solver result ---")
    print("converged", res["converged"], "rc", res["returncode"],
          "iters", md.get("iterations"),
          "closure_err", md.get("closure_err"),
          "tension_err", md.get("tension_err"))
    print("tension kN", md.get("tension_min_kN"), md.get("tension_max_kN"))
    print("cables {}  erected {}  u {}  reaction {}".format(
        len(res["cables"]), len(res["erected"]), len(res["u"]),
        len(res["reaction"])))
    if res["cables"]:
        print("L0 range {:.6f} .. {:.6f}".format(
            min(c["L0"] for c in res["cables"]),
            max(c["L0"] for c in res["cables"])))
        print("T range {:.1f} .. {:.1f} N".format(
            min(c["T"] for c in res["cables"]),
            max(c["T"] for c in res["cables"])))
        print("first cable", res["cables"][0])
    if res["u"]:
        print("max |U| {:.6e}".format(max(abs(v) for v in res["u"].values())))
    if res["camber"]:
        print("max |CAMBER| {:.6e}".format(
            max(max(abs(x) for x in v) for v in res["camber"].values())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
