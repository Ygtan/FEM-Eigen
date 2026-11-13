#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Freeze ANSYS 18.0 BEAM4 results as golden references for E-stage cases.

All reference values come from REAL ANSYS batch CSVs (nothing typed by hand).
Each check carries ``source = <csv file>:<row key>:<column>`` so ``--check``
can detect any drift between the JSON and the frozen CSVs.

    python build_reference_loads.py            # write reference/*.json
    python build_reference_loads.py --check    # verify consistency only
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
REF = HERE / "reference"

CSV_PATTERN = {
    "disp": "{stem}_ansys_disp.csv",
    "reac": "{stem}_ansys_reac.csv",
    "elem": "{stem}_ansys_elem.csv",
}
CSV_LABEL = {
    "disp": "UX/UY/UZ/ROTX/ROTY/ROTZ",
    "reac": "FX/FY/FZ/MX/MY/MZ",
    "elem": "BEAM4 i-end N/Vy/Vz/T/My/Mz then j-end",
}

# CSV_PATTERN describes the frozen tables (also drives doc["csv"], the archive
# row-count audit and the --live-ansys file names).  Check "kind" is a separate
# namespace: "force" checks read the same elem table but are extracted from the
# FEM report rather than from a CSV, hence this extra source-only mapping.
CHECK_SOURCE_CSV = dict(CSV_PATTERN, force=CSV_PATTERN["elem"])

CASES_SPEC = [
    {
        "id": "cb3d06_udl",
        "title": "CB3D06 3x30m box beam, REAL *ELEMENT_LOAD q=20 kN/m",
        "description": "3D BEAM, 31 nodes / 30 elements, global -Z UDL via "
                       "*ELEMENT_LOAD (consistent nodal loads incl. fixed-end "
                       "moments); ANSYS side uses SFBEAM LKEY=1",
        "nnode": 31, "nelem": 30, "support_nodes": [1, 11, 21, 31],
    },
    {
        "id": "cb3d07_selfweight",
        "title": "CB3D07 distributed UDL + DENS selfweight + SELF_WEIGHT + point load",
        "description": "stacked load channels: *ELEMENT_LOAD 20 kN/m, "
                       "DENS=2500 with *GRAVITY 0 0 -9.81, "
                       "*SELF_WEIGHT ALL MASS 100, node 16 FZ=-100 kN",
        "nnode": 31, "nelem": 30, "support_nodes": [1, 11, 21, 31],
    },
    {
        "id": "cb3d08_trapezoid",
        "title": "CB3D08 trapezoidal *ELEMENT_LOAD q1=20 / q2=40 kN/m",
        "description": "every BEAM element carries a trapezoid; verifies the "
                       "q1/q2 fixed-end-force coefficients L/20(7q1+3q2) etc.",
        "nnode": 31, "nelem": 30, "support_nodes": [1, 11, 21, 31],
    },
]

DISP_NAMES = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")
BEAM_NAMES = ("N", "Vy", "Vz", "T", "My", "Mz")


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


def build_one(spec: dict) -> dict:
    stem = spec["id"]
    tables = {}
    for kind, pat in CSV_PATTERN.items():
        p = CASES / pat.format(stem=stem)
        if not p.exists():
            raise FileNotFoundError("missing frozen ANSYS CSV: %s" % p)
        tables[kind] = read_csv_table(p)

    checks: list[dict] = []

    def mk(name, kind, key, col, ref, *, abs_=False, atol=1e-9):
        checks.append({
            "name": name, "kind": kind, "key": key, "col": col,
            "ref": ref, "rtol": 1e-3, "atol": atol,
            "abs": abs_,
            "source": "%s:%s:col%d" % (
                CHECK_SOURCE_CSV[kind].format(stem=stem), key, col + 1),
        })

    # 1) every nodal DOF (FEM report has 6 DOFs per BEAM node)
    for nd in range(1, spec["nnode"] + 1):
        row = tables["disp"].get(nd)
        if row is None or len(row) < 6:
            raise ValueError("%s disp missing node %d" % (stem, nd))
        for c, nm in enumerate(DISP_NAMES):
            mk("disp_n%d_%s" % (nd, nm), "disp", nd, c, row[c], atol=1e-9)

    # 2) every BEAM i-end force component (magnitudes: sign conventions differ)
    for eid in range(1, spec["nelem"] + 1):
        row = tables["elem"].get(eid)
        if row is None or len(row) < 6:
            raise ValueError("%s elem missing element %d" % (stem, eid))
        for c, nm in enumerate(BEAM_NAMES):
            mk("beam_e%d_%s" % (eid, nm), "force", eid, c, abs(row[c]),
               abs_=True, atol=1e-6)

    # 3) support reactions (used by --live-ansys and by the force balance note)
    reactions = []
    for nd in spec["support_nodes"]:
        row = tables["reac"].get(nd, [0.0] * 6)
        reactions.append({"node": nd, "fz": row[2]})

    return {
        "id": stem,
        "title": spec["title"],
        "description": spec["description"],
        "dat": stem + ".dat",
        "apdl": stem + ".inp",
        "exe": "example_3d_report.exe",
        "reference_source": (
            "ANSYS 18.0 BEAM4 batch results frozen in "
            "cases/%s_ansys_{disp,reac,elem}.csv; generated by "
            "BridgeStudio/ansys_bridge_3d_loads.py from the same ModelSpec3D"
            % stem
        ),
        "csv": {k: v.format(stem=stem) for k, v in CSV_PATTERN.items()},
        "csv_label": CSV_LABEL,
        "expected_rows": {
            "disp": spec["nnode"], "reac": spec["nnode"],
            "elem": spec["nelem"],
        },
        "support_nodes": spec["support_nodes"],
        "support_fz_sum": sum(r["fz"] for r in reactions),
        "checks": checks,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="only verify JSON vs frozen CSV, do not write")
    a = ap.parse_args()
    REF.mkdir(parents=True, exist_ok=True)
    bad = 0
    for spec in CASES_SPEC:
        data = build_one(spec)
        out = REF / (spec["id"] + "_ansys.json")
        if a.check:
            if not out.exists():
                print("[FAIL] missing %s" % out.name)
                bad += 1
                continue
            old = json.loads(out.read_text(encoding="utf-8"))
            if old == data:
                print("[PASS] %s  %d checks match frozen CSV" %
                      (out.name, len(data["checks"])))
            else:
                print("[FAIL] %s out of sync with frozen CSV" % out.name)
                bad += 1
        else:
            out.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")
            print("[WRITE] %s  %d checks" % (out.name, len(data["checks"])))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
