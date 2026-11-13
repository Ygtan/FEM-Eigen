#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS CSV + FEM [_cable.txt] 生成 nonlinear_ansys_3d 黄金参考值。

每一项检查的 ref 都来自真实 ANSYS 18.0 (NLGEOM + LINK180 + INISTATE)
批处理导出的 CSV；检查项清单（哪些节点 / 自由度 / 索单元）来自 FEM
成桥态结果，保证 CTest 重跑 cable_solver 后逐项有对应。

容差口径与 BridgeStudio/ansys_bridge_nl.py::compare_nl() 一致：
    对每个分组（位移 / 索力 / 反力）取 scale = max(|FEM|, |ANSYS|)，
    逐项判据 |got-ref| <= atol + rtol * scale，rtol = 2e-3。
不能退化成 rtol * |ref|：转动自由度量级只有 1e-5 rad，与 ~1e-2 m 的
平动共用分组量级，逐项相对容差会把 ANSYS/FEM 的一阶差异误判为 FAIL。

用法：
    python build_reference_nl.py            # 生成 reference/*.json
    python build_reference_nl.py --check    # 只校验一致性，不写文件
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference"

CASE_ID = "cb3d10_completion"
APDL_STEM = "cb3d10_completion_nl"
CSV_PATTERN = {
    "disp": "{stem}_ansys_disp.csv",
    "reac": "{stem}_ansys_reac.csv",
    "link": "{stem}_ansys_link.csv",
}
CSV_LABEL = {
    "disp": "UX/UY/UZ/ROTX/ROTY/ROTZ (relative to [ERECTED] geometry)",
    "reac": "FX/FY/FZ/MX/MY/MZ",
    "link": "LINK180 axial force SMISC 1 (tension only)",
}
# CSV_PATTERN 描述的是冻结表（同时驱动 doc["csv"]、归档行数审计和
# --live-ansys 的输出文件名）。检查项 kind 是独立命名空间：索力检查从
# FEM 结果取值，但 source 指向 link 表，所以需要这张「仅用于溯源」的映射。
CHECK_SOURCE_CSV = dict(CSV_PATTERN, force=CSV_PATTERN["link"])
FEM_RESULT = CASE_ID + "_cable.txt"
RTOL = 2e-3

ARCHIVE = [
    CASE_ID + ".dat",
    CASE_ID + ".model.json",
    FEM_RESULT,
    APDL_STEM + ".inp",
    APDL_STEM + "_ansys.out",
    APDL_STEM + "_ansys_disp.csv",
    APDL_STEM + "_ansys_reac.csv",
    APDL_STEM + "_ansys_link.csv",
    APDL_STEM + "_compare_report.txt",
    APDL_STEM + "_manifest.txt",
    APDL_STEM + "_run_log.txt",
]
# ANSYS 18.0 对 LINK180 KEYOPT(3)=1 会如实给出一条
# "is an undocumented option" 警告；该警告不隐藏，作为已知基线记录。
ANSYS_LOG_EXPECT = {
    "errors": 0,
    "warnings": 1,
    "warning_signature": "KEYOPT(3) = 1 is an undocumented",
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


def parse_fem_result(path: Path) -> dict:
    """自包含解析 <stem>_cable.txt 的 [U] / [CABLE] / [REACTION] / 头部。"""
    u: dict[tuple[int, int], float] = {}
    cables: dict[int, dict] = {}
    reaction: dict[tuple[int, int], float] = {}
    meta: dict[str, float] = {}
    sec = ""
    for raw in path.read_text(encoding="utf-8",
                              errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("["):
            sec = line.strip("[]").upper()
            continue
        parts = line.split()
        if sec == "":
            for tok in parts:
                if "=" in tok:
                    k, _, v = tok.partition("=")
                    try:
                        meta[k] = float(v)
                    except ValueError:
                        pass
            continue
        if sec == "U" and len(parts) >= 3:
            u[(int(parts[0]), int(parts[1]))] = float(parts[2])
        elif sec == "CABLE" and len(parts) >= 7:
            cables[int(parts[0])] = {
                "anchor": int(parts[1]), "T": float(parts[2]),
                "L0": float(parts[3]), "chord": float(parts[4]),
            }
        elif sec == "REACTION" and len(parts) >= 3:
            reaction[(int(parts[0]), int(parts[1]))] = float(parts[2])
    return {"u": u, "cables": cables, "reaction": reaction, "meta": meta}


def build_one() -> dict:
    fem = parse_fem_result(CASES / FEM_RESULT)
    tables: dict[str, dict[int, list[float]]] = {}
    for kind, pat in CSV_PATTERN.items():
        p = CASES / pat.format(stem=APDL_STEM)
        if not p.exists():
            raise FileNotFoundError("缺少冻结 ANSYS CSV: %s" % p)
        tables[kind] = read_csv_table(p)

    # ---- 先取出三组配对值，计算与 compare_nl() 相同的分组量级 ----
    disp_pairs: list[tuple[int, int, float, float]] = []
    for (nd, dof), fv in sorted(fem["u"].items()):
        row = tables["disp"].get(nd)
        if row is None or len(row) < dof:
            raise ValueError("disp 缺节点 %d dof %d" % (nd, dof))
        disp_pairs.append((nd, dof, fv, row[dof - 1]))

    cable_pairs: list[tuple[int, float, float]] = []
    for eid in sorted(fem["cables"]):
        row = tables["link"].get(eid)
        if row is None or not row:
            raise ValueError("link 缺索单元 %d" % eid)
        cable_pairs.append((eid, fem["cables"][eid]["T"], abs(row[0])))

    reac_pairs: list[tuple[int, int, float, float]] = []
    for (nd, dof), fv in sorted(fem["reaction"].items()):
        row = tables["reac"].get(nd)
        if row is None or len(row) < dof:
            raise ValueError("reac 缺节点 %d dof %d" % (nd, dof))
        reac_pairs.append((nd, dof, fv, row[dof - 1]))

    def scale_of(pairs) -> float:
        return max([max(abs(v) for v in p[2:]) for p in pairs] + [0.0])

    scales = {
        "disp": scale_of(disp_pairs),
        "force": scale_of(cable_pairs),
        "reac": scale_of(reac_pairs),
    }

    checks: list[dict] = []

    def mk(name: str, kind: str, key: int, col: int, ref: float,
           *, dof: int | None = None, abs_: bool = False,
           atol: float = 1e-6) -> None:
        item = {
            "name": name, "kind": kind, "key": key, "col": col,
            "ref": ref, "rtol": RTOL, "atol": atol,
            "scale": scales[kind], "abs": abs_,
            "source": "%s:%s:col%d" % (
                CHECK_SOURCE_CSV[kind].format(stem=APDL_STEM),
                key, col + 1),
        }
        if dof is not None:
            item["dof"] = dof
        checks.append(item)

    for nd, dof, _, av in disp_pairs:
        mk("disp_n%d_d%d" % (nd, dof), "disp", nd, dof - 1,
           av, dof=dof, atol=1e-9)
    for eid, _, av in cable_pairs:
        mk("cable_e%d_T" % eid, "force", eid, 0, av, abs_=True)
    for nd, dof, _, av in reac_pairs:
        mk("reac_n%d_d%d" % (nd, dof), "reac", nd, dof - 1,
           av, dof=dof)

    md = fem["meta"]
    return {
        "id": CASE_ID,
        "title": "CB3D10 FEM-Eigen nonlinear completion state vs ANSYS NLGEOM + LINK180",
        "description": (
            "50-node / 63-element twin-pylon cable-stayed bridge, 16 CABLE stays, "
            "MODE=COMPLETION METHOD=DECK, cable self weight W=0 (no sag/Ernst). "
            "ANSYS initial geometry = FEM [ERECTED] manufactured shape; "
            "per-cable unstrained length L0 -> INISTATE,DEFINE EPEL,X eps=(L-L0)/L0. "
            "BEAM4 + LINK180(KEYOPT(3)=1 tension only) + NLGEOM,ON. "
            "Comparison covers 300 displacement DOF, 16 cable forces, 9 reactions."
        ),
        "dat": CASE_ID + ".dat",
        "apdl": APDL_STEM + ".inp",
        "apdl_stem": APDL_STEM,
        "exe": "cable_solver.exe",
        "fem_result": FEM_RESULT,
        "reference_source": (
            "ANSYS 18.0 (BEAM4 + LINK180 + NLGEOM + INISTATE) batch results "
            "frozen in cases/%s_ansys_{disp,reac,link}.csv; "
            "generated by BridgeStudio/ansys_bridge_nl.py from the same "
            "ModelSpec3D + FEM completion state" % APDL_STEM
        ),
        "csv": {k: v.format(stem=APDL_STEM) for k, v in CSV_PATTERN.items()},
        "csv_label": CSV_LABEL,
        "archive": ARCHIVE,
        "expected_rows": {k: len(v) for k, v in tables.items()},
        "fem_rows": {
            "u": len(fem["u"]),
            "cables": len(fem["cables"]),
            "reaction": len(fem["reaction"]),
        },
        "check_counts": {
            "disp": len(disp_pairs),
            "force": len(cable_pairs),
            "reac": len(reac_pairs),
        },
        "group_scale": scales,
        "tolerance": (
            "|got-ref| <= atol + rtol*scale, rtol=2e-3; scale is the "
            "per-group max(|FEM|,|ANSYS|), identical to compare_nl()"
        ),
        "ansys_log_expect": ANSYS_LOG_EXPECT,
        "converged": int(md.get("converged", 0)),
        "closure_err": md.get("closure_err"),
        "tension_err": md.get("tension_err"),
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
