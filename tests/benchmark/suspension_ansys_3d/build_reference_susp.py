#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS CSV + FEM [_cable.txt] 生成 suspension_ansys_3d 黄金参考值。

两个算例：CB3D11 地锚式（EARTH）、CB3D12 自锚式（SELF），参数完全相同，
只改体系。每一项检查的 ref 都来自真实 ANSYS 18.0
(BEAM4 + LINK180 + NLGEOM + INISTATE) 批处理导出的 CSV；检查项清单来自
FEM 找形结果，保证 CTest 重跑 cable_solver 后逐项有对应。

容差口径与 BridgeStudio/ansys_bridge_susp.py::compare_susp() 一致：
    对每个分组（位移 / 索力 / 反力）取 scale = max(|FEM|, |ANSYS|)，
    逐项判据 |got-ref| <= atol + rtol * scale，rtol = 2e-3。

用法：
    python build_reference_susp.py            # 生成 reference/*.json
    python build_reference_susp.py --check    # 只校验一致性，不写文件
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference"

RTOL = 2e-3
#: 悬索桥 FormFinding 的 closure_err = max|[U]|（设计几何到平衡构形的
#: 位移量级）、tension_err = 索力达标偏差，两者都不是「迭代残差」，
#: 不能沿用斜拉桥 COMPLETION 的 1e-6。重跑结果另由
#: compare_suspension_ansys_3d.py 的 frozen_reproduction 逐位对拍，
#: 宽松阈值不会掩盖任何漂移。
CLOSURE_TOL = 1e-2
TENSION_TOL = 1e-2

CASE_SPECS = [
    {
        "id": "cb3d11_susp_earth",
        "system": "EARTH",
        "title": "CB3D11 earth-anchored suspension bridge form-finding vs ANSYS NLGEOM",
        "description": (
            "50 nodes / 65 elements, 20 main-cable segments + 17 hangers; "
            "main cable anchored at external anchorages, H taken by ground. "
            "FEM cable_solver FormFinding branch; ANSYS initial geometry = "
            "design geometry D = [ERECTED] + [U], per-cable absolute L0 -> "
            "INISTATE,DEFINE EPEL,X eps=(L_D-L0)/L0. "
            "BEAM4 + LINK180(KEYOPT(3)=1 tension only) + NLGEOM,ON."
        ),
    },
    {
        "id": "cb3d12_susp_self",
        "system": "SELF",
        "title": "CB3D12 self-anchored suspension bridge form-finding vs ANSYS NLGEOM",
        "description": (
            "48 nodes / 65 elements, 20 main-cable segments + 17 hangers; "
            "main cable anchored at girder ends, H compresses the girder; "
            "longitudinal restraint at midspan. Same geometry/load parameters "
            "as CB3D11, only the structural system differs."
        ),
    },
]

CSV_PATTERN = {
    "disp": "{stem}_disp.csv",
    "reac": "{stem}_reac.csv",
    "link": "{stem}_link.csv",
}
CSV_LABEL = {
    "disp": "UX/UY/UZ/ROTX/ROTY/ROTZ (relative to design geometry D)",
    "reac": "FX/FY/FZ/MX/MY/MZ",
    "link": "LINK180 axial force SMISC 1 (tension only)",
}
CHECK_SOURCE_CSV = dict(CSV_PATTERN, force=CSV_PATTERN["link"])

#: ANSYS 18.0 对两个 CABLE 单元类型各给一条
#: "KEYOPT(3) = 1 is an undocumented option" 警告；不隐藏，作为已知基线。
ANSYS_LOG_EXPECT = {
    "errors": 0,
    "warnings": 2,
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
                if "=" in tok:
                    k, _, v = tok.partition("=")
                    try:
                        meta[k] = float(v)
                    except ValueError:
                        pass
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
        except ValueError:
            continue
    return {"u": u, "cables": cables, "reaction": reaction, "meta": meta}


def build_one(spec: dict) -> dict:
    case_id = spec["id"]
    apdl_stem = case_id + "_nl_ansys"       # ANSYS jobname == CSV prefix
    fem_result = case_id + "_cable.txt"
    fem = parse_fem_result(CASES / fem_result)
    tables: dict[str, dict[int, list[float]]] = {}
    for kind, pat in CSV_PATTERN.items():
        p = CASES / pat.format(stem=apdl_stem)
        if not p.exists():
            raise FileNotFoundError("缺少冻结 ANSYS CSV: %s" % p)
        tables[kind] = read_csv_table(p)

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
                CHECK_SOURCE_CSV[kind].format(stem=apdl_stem), key, col + 1),
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

    csv_names = {k: v.format(stem=apdl_stem) for k, v in CSV_PATTERN.items()}
    archive = [
        case_id + ".dat",
        case_id + ".model.json",
        fem_result,
        case_id + "_nl.inp",
        apdl_stem + ".out",
        *[csv_names[k] for k in ("disp", "reac", "link")],
        case_id + "_nl_compare_report.txt",
        case_id + "_nl_manifest.txt",
        case_id + "_nl_run_log.txt",
    ]
    md = fem["meta"]
    return {
        "id": case_id,
        "title": spec["title"],
        "description": spec["description"],
        "system": spec["system"],
        "dat": case_id + ".dat",
        "apdl": case_id + "_nl.inp",
        "apdl_stem": apdl_stem,
        "exe": "cable_solver.exe",
        "fem_result": fem_result,
        "reference_source": (
            "ANSYS 18.0 (BEAM4 + LINK180 + NLGEOM + INISTATE) batch results "
            "frozen in cases/%s_{disp,reac,link}.csv; generated by "
            "BridgeStudio/ansys_bridge_susp.py + ansys_bridge_nl.py from the "
            "same ModelSpec3D + FEM form-finding result" % apdl_stem
        ),
        "geometry_note": (
            "ANSYS initial geometry = design geometry D = [ERECTED] + [U], "
            "NOT [ERECTED]. [CAMBER]/[ERECTED] are camber suggestions only; "
            "L0 is back-computed from the equilibrated chord length. "
            "Using [ERECTED] directly gives a 0.25 relative displacement error."
        ),
        "csv": csv_names,
        "csv_label": CSV_LABEL,
        "archive": archive,
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
            "per-group max(|FEM|,|ANSYS|), identical to compare_susp()"
        ),
        "ansys_log_expect": dict(ANSYS_LOG_EXPECT),
        "converged": int(md.get("converged", 0)),
        "closure_err": md.get("closure_err"),
        "tension_err": md.get("tension_err"),
        "closure_tol": CLOSURE_TOL,
        "tension_tol": TENSION_TOL,
        "checks": checks,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="只校验已生成 JSON 与冻结文件是否一致")
    a = ap.parse_args()
    REF.mkdir(parents=True, exist_ok=True)
    bad = 0
    for spec in CASE_SPECS:
        data = build_one(spec)
        out = REF / (spec["id"] + "_ansys.json")
        if a.check:
            if not out.exists():
                print("[FAIL] missing %s" % out.name)
                bad += 1
                continue
            old = json.loads(out.read_text(encoding="utf-8"))
            if old == data:
                print("[PASS] %s  %d checks match frozen files" %
                      (out.name, len(data["checks"])))
            else:
                print("[FAIL] %s out of sync with frozen files" % out.name)
                bad += 1
            continue
        out.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                       encoding="utf-8")
        print("[WRITE] %s  %d checks  scales=%s" %
              (out.name, len(data["checks"]), data["group_scale"]))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
