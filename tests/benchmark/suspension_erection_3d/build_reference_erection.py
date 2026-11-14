#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS CSV + FEM [_cable.txt] 生成 CB3D13 黄金参考值。

CB3D13 = 自锚式悬索桥「自动平衡边跨 + 施工阶段找形」，三层验证：

  [A] 全桥成桥态（BEAM4 + LINK180 + NLGEOM + INISTATE）
      位移 / 索力 / 反力，与 suspension_ansys_3d 同口径。
  [B] 自动平衡边跨（纯 FEM 判据，不依赖 ANSYS CSV）
      梁端竖向反力 |Rz| / (q_deck*L_main) <= 1e-3。
  [C] 施工阶段 5 个有荷载阶段，每阶段一份 LINK180 索链 APDL，
      从 0.20 m 扰动几何独立 NLGEOM 求解；逐阶段比较
      主缆线形 z / 主缆段轴力 T / 水平力 H。

每一项 ``ref`` 都来自真实 ANSYS 18.0 批处理导出的 CSV；[C] 的线形与
水平力是从 CSV 派生量（z_init + UZ、T 的水平分量），派生口径写在
``derived`` 字段中，并由 compare_suspension_erection_3d.py 重算校验。

容差口径与 BridgeStudio/ansys_bridge_susp.py::compare_susp() 以及
ansys_bridge_erection.py::compare_erection() 一致：
    每个分组取 scale = max(|FEM|, |ANSYS|)，
    逐项 |got-ref| <= atol + rtol * scale，rtol = 2e-3。

用法：
    python build_reference_erection.py            # 生成 reference/*.json
    python build_reference_erection.py --check    # 只校验一致性
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference"

RTOL = 2e-3
CLOSURE_TOL = 1e-2
TENSION_TOL = 1e-2
BALANCE_TOL = 1e-3

CASE_ID = "cb3d13_balanced_erection"
APDL_STEM = CASE_ID + "_nl_ansys"
FEM_RESULT = CASE_ID + "_cable.txt"
STAGES = [1, 2, 3, 4, 5]
SKIPPED_STAGES = [0]

TITLE = "CB3D13 self-anchored suspension, balanced side span + erection stages vs ANSYS"
DESCRIPTION = (
    "60 nodes / 83 elements, 26 main-cable segments + 23 hangers; "
    "side_span=0 -> auto balanced L_bal=sqrt(2*H*dz_top/q)=71.414284 m. "
    "[A] Full-bridge completion: design geometry D=[ERECTED]+[U], "
    "per-cable absolute L0 -> INISTATE. [B] Anchor vertical component V~0 "
    "so girder ends need no uplift restraint. [C] cable_solver *ERECTION "
    "6 stages (5 with load); per stage a LINK180 main-span cable chain is "
    "solved from a 0.20 m perturbed geometry with stage-specific unstrained "
    "lengths L0=L_stage/(1+T_stage/EA) and stage hanger loads."
)

CSV_PATTERN = {
    "disp": "{stem}_disp.csv",
    "reac": "{stem}_reac.csv",
    "link": "{stem}_link.csv",
}
CHECK_SOURCE_CSV = dict(CSV_PATTERN, force=CSV_PATTERN["link"])
CSV_LABEL = {
    "disp": "UX/UY/UZ/ROTX/ROTY/ROTZ (relative to design geometry D)",
    "reac": "FX/FY/FZ/MX/MY/MZ",
    "link": "LINK180 axial force SMISC 1 (full bridge, tension only)",
    "stage_disp": "UX/UY/UZ of the main-span main-cable chain",
    "stage_link": "LINK180 axial force SMISC 1 (main-span cable segments)",
}

#: 全桥模型有两个 CABLE 单元类型，ANSYS 18.0 各给一条 undocumented 警告；
#: 阶段索链用 KEYOPT(3)=0（各阶段均受拉），0 warning。
ANSYS_LOG_EXPECT = {
    "full": {"errors": 0, "warnings": 2,
             "warning_signature": "KEYOPT(3) = 1 is an undocumented"},
    "stage": {"errors": 0, "warnings": 0, "warning_signature": None},
}


def read_csv_table(path: Path) -> dict[int, list[float]]:
    rows: dict[int, list[float]] = {}
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        toks = raw.replace(",", " ").split()
        if len(toks) < 2:
            continue
        try:
            vals = [float(t) for t in toks]
        except ValueError:
            continue
        rows[int(round(vals[0]))] = vals[1:]
    return rows


def _kv(tok: str) -> tuple[str, float] | None:
    k, sep, v = tok.partition("=")
    if not sep:
        return None
    try:
        return k, float(v)
    except ValueError:
        return None


def parse_fem_result(path: Path) -> dict:
    """解析 <stem>_cable.txt，含 I 阶段新增的 [ERECTION_*] 段。"""
    u: dict[tuple[int, int], float] = {}
    cables: dict[int, dict] = {}
    reaction: dict[tuple[int, int], float] = {}
    meta: dict[str, float] = {}
    erection: dict[int, dict] = {}
    shape: dict[int, list[float]] = {}
    chain: list[tuple[int, float]] = []
    load: dict[int, list[float]] = {}
    force: dict[int, list[float]] = {}
    sec = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("["):
            sec = line.strip("[]").upper()
            continue
        parts = line.split()
        if sec == "":
            for tok in parts:
                kv = _kv(tok)
                if kv:
                    meta[kv[0]] = kv[1]
            continue
        try:
            if sec == "U" and len(parts) >= 3:
                u[(int(parts[0]), int(parts[1]))] = float(parts[2])
            elif sec == "CABLE" and len(parts) >= 7:
                cables[int(parts[0])] = {
                    "anchor": int(parts[1]), "T": float(parts[2]),
                    "L0": float(parts[3]), "chord": float(parts[4]),
                }
            elif sec == "REACTION" and len(parts) >= 3:
                reaction[(int(parts[0]), int(parts[1]))] = float(parts[2])
            elif sec == "ERECTION" and len(parts) >= 3:
                row: dict = {"stage": int(parts[0]), "name": parts[1]}
                for tok in parts[2:]:
                    kv = _kv(tok)
                    if kv:
                        row[kv[0]] = kv[1]
                erection[int(parts[0])] = row
            elif sec == "ERECTION_SHAPE" and len(parts) >= 2:
                shape[int(parts[0])] = [float(v) for v in parts[1:]]
            elif sec == "ERECTION_CHAIN" and len(parts) >= 2:
                chain.append((int(parts[0]), float(parts[1])))
            elif sec == "ERECTION_LOAD" and len(parts) >= 2:
                load[int(parts[0])] = [float(v) for v in parts[1:]]
            elif sec == "ERECTION_FORCE" and len(parts) >= 2:
                force[int(parts[0])] = [float(v) for v in parts[1:]]
        except ValueError:
            continue
    return {"u": u, "cables": cables, "reaction": reaction, "meta": meta,
            "erection": erection, "shape": shape, "chain": chain,
            "load": load, "force": force}


def stage_apdl_zinit(path: Path) -> list[float]:
    """从归档 APDL 的 ``N,<i>,<x>,0,<z>`` 行读出扰动初始几何。"""
    z: list[float] = []
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        if raw.startswith("N,"):
            z.append(float(raw.split(",")[4]))
    return z


def stage_ansys(disp: dict[int, list[float]], link: dict[int, list[float]],
                *, xs: list[float], z_init: list[float]) -> dict:
    """由冻结 ANSYS CSV 还原阶段构形 / 索力 / 水平力（与算例脚本同口径）。"""
    n = len(xs)
    x = [xs[i - 1] + disp[i][0] for i in range(1, n + 1)]
    z = [z_init[i - 1] + disp[i][2] for i in range(1, n + 1)]
    T = [link[j][0] for j in range(1, n)]
    Hs = []
    for k in range(n - 1):
        dx, dz = x[k + 1] - x[k], z[k + 1] - z[k]
        L = math.hypot(dx, dz)
        Hs.append(T[k] * abs(dx) / L if L > 0 else 0.0)
    return {"x": x, "z": z, "T": T, "H": sum(Hs) / len(Hs)}


def build_one() -> dict:
    fem = parse_fem_result(CASES / FEM_RESULT)
    tables = {k: read_csv_table(CASES / v.format(stem=APDL_STEM))
              for k, v in CSV_PATTERN.items()}

    # ---------------- [A] 全桥成桥态 ----------------
    disp_pairs, cable_pairs, reac_pairs = [], [], []
    for (nd, dof), fv in sorted(fem["u"].items()):
        row = tables["disp"].get(nd)
        if row is None or len(row) < dof:
            raise ValueError("disp 缺节点 %d dof %d" % (nd, dof))
        disp_pairs.append((nd, dof, fv, row[dof - 1]))
    for eid in sorted(fem["cables"]):
        row = tables["link"].get(eid)
        if row is None or not row:
            raise ValueError("link 缺索单元 %d" % eid)
        cable_pairs.append((eid, fem["cables"][eid]["T"], abs(row[0])))
    for (nd, dof), fv in sorted(fem["reaction"].items()):
        row = tables["reac"].get(nd)
        if row is None or len(row) < dof:
            raise ValueError("reac 缺节点 %d dof %d" % (nd, dof))
        reac_pairs.append((nd, dof, fv, row[dof - 1]))

    # ---------------- [C] 施工阶段 ----------------
    xs = [x for _, x in fem["chain"]]
    n = len(xs)
    stage_data: dict[int, dict] = {}
    stage_files: dict[str, dict] = {}
    for st in STAGES:
        s = CASE_ID + "_stage%d" % st
        dname, lname, iname = s + "_disp.csv", s + "_link.csv", s + ".inp"
        for f in (dname, lname, iname):
            if not (CASES / f).exists():
                raise FileNotFoundError("缺少冻结阶段文件: " + f)
        zi = stage_apdl_zinit(CASES / iname)
        if len(zi) != n:
            raise ValueError("stage %d APDL 节点数 %d != CHAIN %d" % (st, len(zi), n))
        d, lk = read_csv_table(CASES / dname), read_csv_table(CASES / lname)
        stage_data[st] = stage_ansys(d, lk, xs=xs, z_init=zi)
        stage_files[str(st)] = {"disp": dname, "link": lname,
                                "apdl": iname, "out": s + "_ansys.out",
                                "z_init": zi}

    # ---------------- 分组量级 ----------------
    def scale_of(pairs) -> float:
        return max([max(abs(v) for v in p[2:]) for p in pairs] + [0.0])

    scales = {
        "disp": scale_of(disp_pairs),
        "force": scale_of(cable_pairs),
        "reac": scale_of(reac_pairs),
        "stage_shape": max(
            max(abs(a), abs(b))
            for st in STAGES
            for a, b in zip(fem["shape"][st], stage_data[st]["z"])),
        "stage_force": max(
            max(abs(a), abs(b))
            for st in STAGES
            for a, b in zip(fem["force"][st], stage_data[st]["T"])),
        "stage_H": max(
            max(abs(fem["erection"][st]["H"]), abs(stage_data[st]["H"]))
            for st in STAGES),
    }

    checks: list[dict] = []

    def mk(name: str, kind: str, *, key: int | None = None,
          col: int | None = None, ref: float, dof: int | None = None,
          stage: int | None = None, abs_: bool = False,
          atol: float = 1e-6, source: str = "") -> None:
        item = {"name": name, "kind": kind, "ref": ref,
                "rtol": RTOL, "atol": atol, "scale": scales[kind],
                "abs": abs_, "source": source}
        if key is not None:
            item["key"] = key
        if col is not None:
            item["col"] = col
        if dof is not None:
            item["dof"] = dof
        if stage is not None:
            item["stage"] = stage
        checks.append(item)

    for nd, dof, _, av in disp_pairs:
        mk("disp_n%d_d%d" % (nd, dof), "disp", key=nd, col=dof - 1,
           ref=av, dof=dof, atol=1e-9,
           source="%s:%d:col%d" % (
               CHECK_SOURCE_CSV["disp"].format(stem=APDL_STEM), nd, dof))
    for eid, _, av in cable_pairs:
        mk("cable_e%d_T" % eid, "force", key=eid, col=0, ref=av, abs_=True,
           source="%s:%d:col1" % (
               CHECK_SOURCE_CSV["force"].format(stem=APDL_STEM), eid))
    for nd, dof, _, av in reac_pairs:
        mk("reac_n%d_d%d" % (nd, dof), "reac", key=nd, col=dof - 1,
           ref=av, dof=dof,
           source="%s:%d:col%d" % (
               CHECK_SOURCE_CSV["reac"].format(stem=APDL_STEM), nd, dof))
    for st in STAGES:
        sd = stage_data[st]
        for k, av in enumerate(sd["z"]):
            mk("stage%d_z%d" % (st, k), "stage_shape", key=k, ref=av,
               stage=st, atol=1e-6,
               source="%s:%d:col4+APDL_N_z" % (
                   stage_files[str(st)]["disp"], k + 1))
        for k, av in enumerate(sd["T"]):
            mk("stage%d_T%d" % (st, k), "stage_force", key=k, ref=abs(av),
               stage=st, abs_=True,
               source="%s:%d:col1" % (stage_files[str(st)]["link"], k + 1))
        mk("stage%d_H" % st, "stage_H", ref=abs(sd["H"]), stage=st,
           source="%s+%s:derived(mean T*dx/L)" % (
               stage_files[str(st)]["link"], stage_files[str(st)]["disp"]))

    # ---------------- [B] 自动平衡边跨（FEM 判据）----------------
    # 梁端 = 桥面高程 z≈0 且 x 最小/最大的支座节点。不能简单地对
    # [REACTION] 的 dof=3 节点取 min/max id：塔顶节点 (3001/4001) 也有
    # 竖向反力，会混进来，把 4682 kN 当成「梁端抗拔力」。
    model = json.loads((CASES / (CASE_ID + ".model.json")).read_text(
        encoding="utf-8"))
    nodes_xyz = {n["id"]: (float(n["x"]), float(n["y"]), float(n["z"]))
                 for n in model["nodes"]}
    deck_reac = [nd for (nd, dof) in fem["reaction"]
                if dof == 3 and nd in nodes_xyz and abs(nodes_xyz[nd][2]) < 1e-9]
    if len(deck_reac) < 2:
        raise ValueError("找不到桥面两端支座节点")
    end_nodes = [min(deck_reac, key=lambda n: nodes_xyz[n][0]),
                max(deck_reac, key=lambda n: nodes_xyz[n][0])]
    end_rz = [fem["reaction"].get((nd, 3), 0.0) for nd in end_nodes]
    q_deck = 20e3
    L_main = 120.0
    balance_ratio = max(abs(v) for v in end_rz) / (q_deck * L_main)
    balance_expect = {
        "end_nodes": end_nodes,
        "end_rz_N": end_rz,
        "q_deck_N_per_m": q_deck,
        "main_span_m": L_main,
        "max_end_rz_ratio": balance_ratio,
        "tol": BALANCE_TOL,
        "side_span_balanced_m": math.sqrt(5100.0),
        "note": ("L_bal = sqrt(2*H0*dz_top/q); side_span=0 asks the "
                 "generator to use it. V=0 at the girder-end anchor means "
                 "no uplift restraint is needed; CB3D12 with a manual 40 m "
                 "side span needs -872.621 kN (uplift) at each end."),
    }

    csv_names = {k: v.format(stem=APDL_STEM) for k, v in CSV_PATTERN.items()}
    archive = [
        CASE_ID + ".dat", CASE_ID + ".model.json", FEM_RESULT,
        CASE_ID + "_nl.inp", APDL_STEM + ".out",
        *[csv_names[k] for k in ("disp", "reac", "link")],
        CASE_ID + "_nl_compare_report.txt",
        CASE_ID + "_balance_report.txt",
        CASE_ID + "_erection_compare_report.txt",
        CASE_ID + "_manifest.txt", CASE_ID + "_run_log.txt",
    ]
    for st in STAGES:
        sf = stage_files[str(st)]
        archive += [sf["apdl"], sf["out"], sf["disp"], sf["link"]]

    md = fem["meta"]
    return {
        "id": CASE_ID,
        "title": TITLE,
        "description": DESCRIPTION,
        "dat": CASE_ID + ".dat",
        "apdl": CASE_ID + "_nl.inp",
        "apdl_stem": APDL_STEM,
        "exe": "cable_solver.exe",
        "fem_result": FEM_RESULT,
        "reference_source": (
            "ANSYS 18.0 batch results frozen in cases/*_ansys_*.csv and "
            "cases/*_stage*_*.csv; generated by BridgeStudio/"
            "ansys_bridge_susp.py, ansys_bridge_nl.py and "
            "ansys_bridge_erection.py from the same ModelSpec3D + FEM "
            "form-finding result"
        ),
        "geometry_note": (
            "Full bridge: ANSYS initial geometry = design geometry "
            "D = [ERECTED] + [U], per-cable absolute L0 -> INISTATE. "
            "Erection stages: ANSYS initial geometry = FEM stage shape "
            "+ 0.20 m sinusoidal perturbation; stage unstrained length "
            "L0_stage = L_stage/(1+T_stage/EA) injected per LINK180 segment."
        ),
        "csv": csv_names,
        "csv_label": CSV_LABEL,
        "stage_files": stage_files,
        "stage_csv_label": {"disp": CSV_LABEL["stage_disp"],
                            "link": CSV_LABEL["stage_link"]},
        "stages": STAGES,
        "skipped_stages": SKIPPED_STAGES,
        "skipped_reason": (
            "stage 0 (empty cable) hits the H bisection lower bound 100 N "
            "with sag=0; a Load-free LINK180 chain is a mechanism, so it "
            "is recorded but not cross-checked against ANSYS."
        ),
        "archive": archive,
        "expected_rows": {k: len(v) for k, v in tables.items()},
        "fem_rows": {
            "u": len(fem["u"]), "cables": len(fem["cables"]),
            "reaction": len(fem["reaction"]),
            "erection": len(fem["erection"]),
            "shape": {str(k): len(v) for k, v in fem["shape"].items()},
            "chain": len(fem["chain"]),
            "load": {str(k): len(v) for k, v in fem["load"].items()},
            "force": {str(k): len(v) for k, v in fem["force"].items()},
        },
        "chain_x": xs,
        "check_counts": {
            k: sum(1 for c in checks if c["kind"] == k) for k in scales
        },
        "group_scale": scales,
        "tolerance": (
            "|got-ref| <= atol + rtol*scale, rtol=2e-3; scale is the "
            "per-group max(|FEM|,|ANSYS|), identical to compare_susp() "
            "and compare_erection()"
        ),
        "balance_expect": balance_expect,
        "ansys_log_expect": ANSYS_LOG_EXPECT,
        "converged": int(md.get("converged", 0)),
        "closure_err": md.get("closure_err"),
        "tension_err": md.get("tension_err"),
        "closure_tol": CLOSURE_TOL,
        "tension_tol": TENSION_TOL,
        "derived": {
            "stage_shape": "z_ans[k] = APDL N-node z_init[k] + disp_csv[k+1][2] (UZ)",
            "stage_force": "T_ans[k] = link_csv[k+1][0] (SMISC 1)",
            "stage_H": "mean over k of T[k]*|x[k+1]-x[k]|/L[k], "
                       "with x = chain_x + disp_csv[i][0] (UX)",
        },
        "checks": checks,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="只校验已生成 JSON 与冻结文件是否一致")
    a = ap.parse_args()
    REF.mkdir(parents=True, exist_ok=True)
    data = build_one()
    out = REF / (CASE_ID + "_ansys.json")
    if a.check:
        if not out.exists():
            print("[FAIL] missing %s" % out.name)
            return 1
        old = json.loads(out.read_text(encoding="utf-8"))
        if old == data:
            print("[PASS] %s  %d checks match frozen files" %
                  (out.name, len(data["checks"])))
            return 0
        print("[FAIL] %s out of sync with frozen files" % out.name)
        return 1
    out.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print("[WRITE] %s  %d checks  scales=%s" %
          (out.name, len(data["checks"]), data["group_scale"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
