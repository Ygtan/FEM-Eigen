#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS BEAM4 CSV 生成 3D 连续箱梁黄金参考值。

BEAM4 是 3D Euler-Bernoulli 梁（SHEARZ/SHEARY=0），与 FEM-Eigen 3D
内核（12x12 EB 空间梁）逐位一致，因此作为严格基准。BEAM188 是
Timoshenko 单元，差异单列在 compare_report.txt 的 [4] 节量化，
不进入黄金参考值。

用法：
    python build_reference_3d.py            # 生成 reference/*.json
    python build_reference_3d.py --check    # 校验 JSON 与冻结 CSV 是否脱节
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REFERENCE = HERE / "reference"

CSV = {
    "disp": "{stem}_b4_disp.csv",
    "reac": "{stem}_b4_reac.csv",
    "elem": "{stem}_b4_elem.csv",
}
# 列索引（0-based，首列 id 已剥离）
DISP_COL = {"UX": 0, "UY": 1, "UZ": 2, "ROTX": 3, "ROTY": 4, "ROTZ": 5}
REAC_COL = {"FX": 0, "FY": 1, "FZ": 2, "MX": 3, "MY": 4, "MZ": 5}
# 单元 CSV 统一列 [N Vy Vz T My Mz]_I + [N Vy Vz T My Mz]_J
ELEM_COL = {"N": 0, "Vy": 1, "Vz": 2, "T": 3, "My": 4, "Mz": 5}

SPECS = {
    "cb3d01_3span": dict(
        title="CB3D01 3D three-span box beam 3x30m, nodalized UDL 20kN/m",
        description="三跨连续箱梁 3D 单梁线，31 节点 / 30 单元，纯竖向弯曲",
        supports=(1, 11, 21, 31),
        mid_nodes=(6, 16, 26),
        key_elems=(1, 5, 10, 11, 15, 20, 21, 25, 30),
        disp_dofs=("UZ", "ROTY"),
        force_cols=("Vz", "My"),
    ),
    "cb3d02_space": dict(
        title="CB3D02 3D three-span box beam, lateral force + torsion",
        description="同 CB3D01 几何 + 中跨跨中横向力与扭矩，激活 UY/ROTX/Vy/Mz/T",
        supports=(1, 11, 21, 31),
        mid_nodes=(6, 16, 26),
        key_elems=(1, 5, 10, 11, 15, 20, 21, 25, 30),
        disp_dofs=("UZ", "UY", "ROTX", "ROTY"),
        force_cols=("Vz", "Vy", "My", "Mz", "T"),
    ),
}


def read_table(path: Path) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        t = raw.replace(",", " ").split()
        if len(t) < 2:
            continue
        try:
            v = [float(x) for x in t]
        except ValueError:
            continue
        out[int(round(v[0]))] = v[1:]
    return out


def load_tables(stem: str) -> dict[str, dict[int, list[float]]]:
    return {kind: read_table(CASES / pat.format(stem=stem))
            for kind, pat in CSV.items()}


def build_checks(stem: str, spec: dict, tabs: dict) -> list[dict]:
    checks: list[dict] = []

    def mk(name: str, kind: str, ref: float, source: str,
           *, node: int | None = None, dof: int | None = None,
           elem: int | None = None, col: str | None = None,
           abs_: bool = False, note: str = "", atol: float = 1e-9) -> None:
        c = {"name": name, "kind": kind, "ref": ref,
             "rtol": 1e-3, "atol": atol, "source": source}
        if node is not None:
            c["node"] = node
        if dof is not None:
            c["dof"] = dof
        if elem is not None:
            c["elem"] = elem
        if col is not None:
            c["col"] = col
        if abs_:
            c["abs"] = True
        if note:
            c["note"] = note
        order = ["name", "kind", "node", "dof", "elem", "col", "abs",
                 "ref", "rtol", "atol", "source", "note"]
        checks.append({k: c[k] for k in order if k in c})

    # 位移：支座节点（含零位移检查）+ 每跨跨中
    for nd in sorted(set(spec["supports"]) | set(spec["mid_nodes"])):
        row = tabs["disp"][nd]
        for dof_name in spec["disp_dofs"]:
            ci = DISP_COL[dof_name]
            mk("disp_n%d_%s" % (nd, dof_name), "disp", row[ci],
               "%s:node=%d:col=%d" % (CSV["disp"].format(stem=stem), nd, ci + 2),
               node=nd, dof=ci + 1)
    # 单元 i 端内力（幅值；FEM-Eigen 与 ANSYS 端力符号约定不同）
    for e in spec["key_elems"]:
        row = tabs["elem"][e]
        for col_name in spec["force_cols"]:
            ci = ELEM_COL[col_name]
            mk("force_e%d_%s" % (e, col_name), "force", abs(row[ci]),
               "%s:elem=%d:col=%d" % (CSV["elem"].format(stem=stem), e, ci + 2),
               elem=e, col=col_name, abs_=True, atol=1e-6,
               note="FEM-Eigen 3D report i-end vs ANSYS BEAM4 SMISC, compare magnitude; "
                    "atol=1e-6 N*m absorbs roundoff on theoretical zero components")
    # 支座反力：FEM-Eigen 3D 不输出反力段，此项仅做冻结 CSV 自校验（防篡改）
    for nd in spec["supports"]:
        row = tabs["reac"][nd]
        ci = REAC_COL["FZ"]
        mk("reaction_n%d_FZ" % nd, "reaction", row[ci],
           "%s:node=%d:col=%d" % (CSV["reac"].format(stem=stem), nd, ci + 2),
           node=nd, dof=ci + 1,
           note="ANSYS BEAM4 only; FEM-Eigen 3D report has no reaction section")
    return checks


def build_one(stem: str, spec: dict) -> dict:
    tabs = load_tables(stem)
    model = json.loads((CASES / (stem + ".model.json")).read_text(encoding="utf-8"))
    return {
        "id": stem,
        "dimension": "3D",
        "title": spec["title"],
        "description": spec["description"],
        "dat": stem + ".dat",
        "exe": "example_3d_report.exe",
        "report": stem + ".report.txt",
        "n_nodes": len(model["nodes"]),
        "n_elements": len(model["elements"]),
        "reference_source": (
            "ANSYS 18.0 / BEAM4 (Euler-Bernoulli, SHEARZ=SHEARY=0), "
            "frozen *_b4_disp.csv / *_b4_reac.csv / *_b4_elem.csv; "
            "each check records its source file, row key and column"),
        "theory_note": (
            "BEAM188 (Timoshenko, linear shape function) intentionally NOT used "
            "as reference: it differs from FEM-Eigen 3D EB by transverse shear. "
            "Its results are archived and quantified in *_compare_report.txt [4]."),
        "rtol": 1e-3,
        "checks": build_checks(stem, spec, tabs),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="verify existing JSON still matches frozen CSV")
    a = ap.parse_args()
    REFERENCE.mkdir(parents=True, exist_ok=True)
    bad = 0
    for stem, spec in SPECS.items():
        doc = build_one(stem, spec)
        path = REFERENCE / (stem + "_ansys.json")
        if a.check:
            if not path.exists():
                print("[FAIL] missing %s" % path.name)
                bad += 1
                continue
            old = json.loads(path.read_text(encoding="utf-8"))
            if old.get("checks") != doc["checks"]:
                print("[FAIL] %s out of sync with frozen CSV" % path.name)
                bad += 1
            else:
                print("[OK]   %s (%d checks)" % (path.name, len(doc["checks"])))
        else:
            path.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")
            print("WRITE %s  %d checks" % (path.name, len(doc["checks"])))
    if a.check:
        print("Result: %s" % ("FAIL" if bad else "PASS"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
