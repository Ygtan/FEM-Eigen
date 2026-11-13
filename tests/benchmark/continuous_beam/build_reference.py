#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从冻结的 ANSYS 批处理 CSV 生成连续梁 benchmark 的黄金参考 JSON。

与 tests/benchmark/ansys 的差别：那里的 ref 是"用闭式解代表 ANSYS"（仓库当时
没有 ANSYS 数值结果文件）；这里的 ref 全部取自本机 ANSYS 18.0 / BEAM3 批处理
真实导出的 CSV（cases/*_ansys_*.csv），是外部软件实测值，不是解析解、也不是
FEM-Eigen 自己的输出。

输入（cases/ 目录，由 BridgeStudio/cases/run_cases.py 产出并归档）：
  <stem>.dat                BridgeStudio 前处理输出，FEM-Eigen new format
  <stem>.inp                与 .dat 同源的 ANSYS MAPDL 命令流
  <stem>_ansys_disp.csv     ANSYS: node UX UY ROTZ
  <stem>_ansys_reac.csv     ANSYS: node FX FY
  <stem>_ansys_elem.csv     ANSYS: elem NI VI MI NJ VJ MJ
输出：
  reference/<stem>_ansys.json

用法：
    python build_reference.py                        # 用脚本所在目录的 cases/ 与 reference/
    python build_reference.py --cases DIR --reference DIR
    python build_reference.py --check                # 只校验现有 JSON 与 CSV 是否一致，不写文件
"""
import argparse
import json
import sys
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass
del _s


# CSV 列（整行 0-based）：
#   disp: [node, UX, UY, ROTZ]
#   reac: [node, FX, FY]
#   elem: [elem, NI, VI, MI, NJ, VJ, MJ]
# FEM-Eigen .result 单元内力列（parse_result 的 FORCE_COLS）：
#   [N1, N2, V1, M1, V2, M2]，对应 ANSYS [NI, NJ, VI, MI, VJ, MJ]
#
# src = (csv 种类, 行键(node/elem), 列索引)；force_ref_zero/force_ref_abs 为后处理。
#   force_ref_zero : ANSYS 里的 ~1e-9 舍入噪声，参考值取精确 0，用 atol 判定
#   force_ref_abs  : FEM-Eigen 剪力与 ANSYS SMISC 2/8 符号约定相反，比幅值
CASES = {
    "cb01_2span": {
        "id": "cb01_2span",
        "dat": "cb01_2span.dat",
        "inp": "cb01_2span.inp",
        "exe": "example_mixed.exe",
        "description": (
            "两跨等截面连续梁 2x20 m（总长 40 m），矩形截面 0.6x1.2 m，"
            "A=0.72 m2 I=0.0864 m4，E=2.06e11 Pa，20 个 2D BEAM 单元 / 21 节点，"
            "全桥均布荷载 q=20 kN/m 向下，自重关闭。支座节点 1/11/21。"
        ),
        "ansys": {"release": "18.0", "element": "BEAM3", "mode": "batch -b"},
        "checks": [
            # 竖向挠度（m，向下为负）
            ("left_span_quarter_uy", "nodal_disp", 6, 2, ("disp", 6, 2), {}),
            ("left_span_max_uy", "nodal_disp", 5, 2, ("disp", 5, 2), {}),
            ("right_span_quarter_uy", "nodal_disp", 16, 2, ("disp", 16, 2), {}),
            ("right_span_max_uy", "nodal_disp", 17, 2, ("disp", 17, 2), {}),
            ("end_support_uy", "nodal_disp", 1, 2, ("disp", 1, 2), {"atol": 1e-9}),
            ("mid_support_uy", "nodal_disp", 11, 2, ("disp", 11, 2), {"atol": 1e-9}),
            ("far_support_uy", "nodal_disp", 21, 2, ("disp", 21, 2), {"atol": 1e-9}),
            # 支座竖向反力（N，向上为正）
            ("reaction_R1", "reaction", 1, 2, ("reac", 1, 2), {}),
            ("reaction_R2", "reaction", 11, 2, ("reac", 11, 2), {}),
            ("reaction_R3", "reaction", 21, 2, ("reac", 21, 2), {}),
            # 中支座负弯矩（N*m）
            ("hogging_e10_M2", "force", 10, "M2", ("elem", 10, 6), {}),
            ("hogging_e11_M1", "force", 11, "M1", ("elem", 11, 3), {}),
            # 跨内最大正弯矩（N*m）
            ("sagging_e4_M2", "force", 4, "M2", ("elem", 4, 6), {}),
            ("sagging_e17_M1", "force", 17, "M1", ("elem", 17, 3), {}),
            # 端弯矩应为 0（ANSYS 输出 ~1e-9 舍入噪声）
            ("zero_end_moment_e1_M1", "force", 1, "M1", ("elem", 1, 3),
             {"force_ref_zero": True, "atol": 1e-6}),
            ("zero_end_moment_e20_M2", "force", 20, "M2", ("elem", 20, 6),
             {"force_ref_zero": True, "atol": 1e-6}),
            # 端剪力：仅比幅值（两边符号约定相反）
            ("end_shear_e1_V1", "force", 1, "V1", ("elem", 1, 2),
             {"force_ref_abs": True}),
            ("end_shear_e10_V2", "force", 10, "V2", ("elem", 10, 5),
             {"force_ref_abs": True}),
        ],
    },
    "cb02_3span": {
        "id": "cb02_3span",
        "dat": "cb02_3span.dat",
        "inp": "cb02_3span.inp",
        "exe": "example_mixed.exe",
        "description": (
            "三跨等截面连续梁 3x30 m（总长 90 m），箱形截面 h=2.0 m b=8.0 m，"
            "A=0.478 m2 I=0.4148 m4，E=2.06e11 Pa，30 个 2D BEAM 单元 / 31 节点，"
            "全桥均布荷载 q=20 kN/m 向下，自重关闭。支座节点 1/11/21/31。"
        ),
        "ansys": {"release": "18.0", "element": "BEAM3", "mode": "batch -b"},
        "checks": [
            ("span1_mid_uy", "nodal_disp", 6, 2, ("disp", 6, 2), {}),
            ("span1_max_uy", "nodal_disp", 5, 2, ("disp", 5, 2), {}),
            ("span2_mid_uy", "nodal_disp", 16, 2, ("disp", 16, 2), {}),
            ("span3_mid_uy", "nodal_disp", 26, 2, ("disp", 26, 2), {}),
            ("end_support_uy", "nodal_disp", 1, 2, ("disp", 1, 2), {"atol": 1e-9}),
            ("support2_uy", "nodal_disp", 11, 2, ("disp", 11, 2), {"atol": 1e-9}),
            ("support3_uy", "nodal_disp", 21, 2, ("disp", 21, 2), {"atol": 1e-9}),
            ("far_support_uy", "nodal_disp", 31, 2, ("disp", 31, 2), {"atol": 1e-9}),
            ("reaction_R1", "reaction", 1, 2, ("reac", 1, 2), {}),
            ("reaction_R2", "reaction", 11, 2, ("reac", 11, 2), {}),
            ("reaction_R3", "reaction", 21, 2, ("reac", 21, 2), {}),
            ("reaction_R4", "reaction", 31, 2, ("reac", 31, 2), {}),
            ("hogging_e10_M2", "force", 10, "M2", ("elem", 10, 6), {}),
            ("hogging_e11_M1", "force", 11, "M1", ("elem", 11, 3), {}),
            ("hogging_e20_M2", "force", 20, "M2", ("elem", 20, 6), {}),
            ("hogging_e21_M1", "force", 21, "M1", ("elem", 21, 3), {}),
            ("sagging_e4_M2", "force", 4, "M2", ("elem", 4, 6), {}),
            ("sagging_e5_M1", "force", 5, "M1", ("elem", 5, 3), {}),
            ("zero_end_moment_e1_M1", "force", 1, "M1", ("elem", 1, 3),
             {"force_ref_zero": True, "atol": 1e-6}),
            ("zero_end_moment_e30_M2", "force", 30, "M2", ("elem", 30, 6),
             {"force_ref_zero": True, "atol": 1e-6}),
            ("end_shear_e1_V1", "force", 1, "V1", ("elem", 1, 2),
             {"force_ref_abs": True}),
            ("end_shear_e10_V2", "force", 10, "V2", ("elem", 10, 5),
             {"force_ref_abs": True}),
        ],
    },
    "cb03_selfweight": {
        "id": "cb03_selfweight",
        "dat": "cb03_selfweight.dat",
        "inp": "cb03_selfweight.inp",
        "exe": "example_mixed.exe",
        "description": (
            "两跨连续梁 2x20 m（与 CB01 同几何/截面/均布荷载），另开启材料自重："
            "DENS=7850 kg/m3、A=0.72 m2、g=9.81 -> 自重线荷载 55446.12 N/m，"
            "与均布 20 kN/m 叠加后 q=75446.12 N/m，总荷载 3.01784e6 N。"
        ),
        "load_case": ("UDL 20 kN/m + self weight; FEM *GRAVITY 0 -9.81 <-> "
                      "ANSYS MP,DENS + ACEL,0,9.81,0 (opposite sign by definition)"),
        "ansys": {"release": "18.0", "element": "BEAM3", "mode": "batch -b"},
        "checks": [
            ("left_span_quarter_uy", "nodal_disp", 6, 2, ("disp", 6, 2), {}),
            ("left_span_max_uy", "nodal_disp", 5, 2, ("disp", 5, 2), {}),
            ("right_span_quarter_uy", "nodal_disp", 16, 2, ("disp", 16, 2), {}),
            ("end_support_uy", "nodal_disp", 1, 2, ("disp", 1, 2), {"atol": 1e-9}),
            ("mid_support_uy", "nodal_disp", 11, 2, ("disp", 11, 2), {"atol": 1e-9}),
            ("far_support_uy", "nodal_disp", 21, 2, ("disp", 21, 2), {"atol": 1e-9}),
            ("reaction_R1", "reaction", 1, 2, ("reac", 1, 2), {}),
            ("reaction_R2", "reaction", 11, 2, ("reac", 11, 2), {}),
            ("reaction_R3", "reaction", 21, 2, ("reac", 21, 2), {}),
            ("hogging_e10_M2", "force", 10, "M2", ("elem", 10, 6), {}),
            ("hogging_e11_M1", "force", 11, "M1", ("elem", 11, 3), {}),
            ("sagging_e4_M2", "force", 4, "M2", ("elem", 4, 6), {}),
            ("sagging_e17_M1", "force", 17, "M1", ("elem", 17, 3), {}),
            ("zero_end_moment_e1_M1", "force", 1, "M1", ("elem", 1, 3),
             {"force_ref_zero": True, "atol": 1e-6}),
            ("end_shear_e1_V1", "force", 1, "V1", ("elem", 1, 2),
             {"force_ref_abs": True}),
            ("end_shear_e10_V2", "force", 10, "V2", ("elem", 10, 5),
             {"force_ref_abs": True}),
        ],
    },
    "cb04_point_loads": {
        "id": "cb04_point_loads",
        "dat": "cb04_point_loads.dat",
        "inp": "cb04_point_loads.inp",
        "exe": "example_mixed.exe",
        "description": (
            "两跨连续梁 2x20 m（CB01 几何），无均布荷载、无自重；"
            "左右两跨跨中节点 6 / 16 各施加 200 kN 向下集中力，总荷载 400 kN。"
        ),
        "load_case": ("two -200 kN point loads at nodes 6/16; "
                      "FEM *CLOAD <-> ANSYS F,<node>,FY,-200000"),
        "ansys": {"release": "18.0", "element": "BEAM3", "mode": "batch -b"},
        "checks": [
            ("point_load_node6_uy", "nodal_disp", 6, 2, ("disp", 6, 2), {}),
            ("point_load_node16_uy", "nodal_disp", 16, 2, ("disp", 16, 2), {}),
            ("span_quarter_uy", "nodal_disp", 3, 2, ("disp", 3, 2), {}),
            ("end_support_uy", "nodal_disp", 1, 2, ("disp", 1, 2), {"atol": 1e-9}),
            ("mid_support_uy", "nodal_disp", 11, 2, ("disp", 11, 2), {"atol": 1e-9}),
            ("far_support_uy", "nodal_disp", 21, 2, ("disp", 21, 2), {"atol": 1e-9}),
            ("reaction_R1", "reaction", 1, 2, ("reac", 1, 2), {}),
            ("reaction_R2", "reaction", 11, 2, ("reac", 11, 2), {}),
            ("reaction_R3", "reaction", 21, 2, ("reac", 21, 2), {}),
            ("hogging_e10_M2", "force", 10, "M2", ("elem", 10, 6), {}),
            ("hogging_e11_M1", "force", 11, "M1", ("elem", 11, 3), {}),
            ("sagging_at_load_e6_M1", "force", 6, "M1", ("elem", 6, 3), {}),
            ("zero_end_moment_e1_M1", "force", 1, "M1", ("elem", 1, 3),
             {"force_ref_zero": True, "atol": 1e-6}),
            ("end_shear_e1_V1", "force", 1, "V1", ("elem", 1, 2),
             {"force_ref_abs": True}),
            ("shear_e7_V1", "force", 7, "V1", ("elem", 7, 2),
             {"force_ref_abs": True}),
        ],
    },
    "cb05_elastic_support": {
        "id": "cb05_elastic_support",
        "dat": "cb05_elastic_support.dat",
        "inp": "cb05_elastic_support.inp",
        "exe": "example_mixed.exe",
        "description": (
            "两跨连续梁 2x20 m（CB01 几何），均布 20 kN/m；中支座节点 11 "
            "由刚性改为竖向接地弹簧 K=2e6 N/m。中支座下沉约 32.6 mm，"
            "弹簧反力约 65.15 kN，两端刚性反力各约 367.42 kN。"
        ),
        "load_case": ("elastic mid support K=2e6 N/m; FEM *BOUNDARY 11 2 0 K=2e6 "
                      "<-> ANSYS COMBIN14 + KEYOPT(2)=2 + R,2,2e6 + grounded node"),
        "ansys": {"release": "18.0", "element": "BEAM3 + COMBIN14", "mode": "batch -b"},
        "checks": [
            ("spring_node_uy", "nodal_disp", 11, 2, ("disp", 11, 2), {}),
            ("left_span_quarter_uy", "nodal_disp", 6, 2, ("disp", 6, 2), {}),
            ("left_span_max_uy", "nodal_disp", 5, 2, ("disp", 5, 2), {}),
            ("end_support_uy", "nodal_disp", 1, 2, ("disp", 1, 2), {"atol": 1e-9}),
            ("far_support_uy", "nodal_disp", 21, 2, ("disp", 21, 2), {"atol": 1e-9}),
            ("reaction_R1", "reaction", 1, 2, ("reac", 1, 2), {}),
            ("reaction_R4", "reaction", 21, 2, ("reac", 21, 2), {}),
            # FEM-Eigen 不把弹簧自由度列入反力表，故由 -K*UY 反算后与
            # ANSYS COMBIN14 SMISC1 的幅值比较。
            ("spring_force_node11", "spring_force", 11, 11, ("spring", 11, 1),
             {"force_ref_abs": True, "k": 2.0e6, "spring_value": 0.0}),
            ("hogging_e10_M2", "force", 10, "M2", ("elem", 10, 6), {}),
            ("hogging_e11_M1", "force", 11, "M1", ("elem", 11, 3), {}),
            ("sagging_e4_M2", "force", 4, "M2", ("elem", 4, 6), {}),
            ("zero_end_moment_e1_M1", "force", 1, "M1", ("elem", 1, 3),
             {"force_ref_zero": True, "atol": 1e-6}),
            ("end_shear_e1_V1", "force", 1, "V1", ("elem", 1, 2),
             {"force_ref_abs": True}),
        ],
    },
}

CSV_NAMES = {
    "disp": "{stem}_ansys_disp.csv",
    "reac": "{stem}_ansys_reac.csv",
    "elem": "{stem}_ansys_elem.csv",
    "spring": "{stem}_ansys_spring.csv",   # 仅弹性支承算例有
}
CSV_LABEL = {"disp": "UX/UY/ROTZ", "reac": "FX/FY",
             "elem": "NI/VI/MI/NJ/VJ/MJ", "spring": "node/SMISC1"}


def read_csv_table(path):
    """读 ANSYS *VWRITE 导出的 CSV，返回 {首列 int: 整行 float 列表}。

    保留首列，使下面 src 的 col_idx 就是整行 0-based 索引：
      disp [node,UX,UY,ROTZ] / reac [node,FX,FY] / elem [elem,NI,VI,MI,NJ,VJ,MJ]
    """
    rows = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        toks = raw.replace(",", " ").split()
        if len(toks) < 2:
            continue
        try:
            key = int(round(float(toks[0])))
        except ValueError:
            continue
        try:
            rows[key] = [float(t) for t in toks]
        except ValueError:
            continue
    return rows


def load_tables(cases_dir, stem):
    tables = {}
    for kind, pattern in CSV_NAMES.items():
        path = cases_dir / pattern.format(stem=stem)
        if not path.exists():
            if kind == "spring":
                tables[kind] = {}      # 只有弹性支承算例才有弹簧 CSV
                continue
            raise FileNotFoundError("缺少 ANSYS 结果文件: %s" % path)
        tables[kind] = read_csv_table(path)
    return tables


def build_case(spec, cases_dir):
    stem = spec["id"]
    tables = load_tables(cases_dir, stem)
    checks = []
    for name, kind, key, fem_sel, src, opts in spec["checks"]:
        csv_kind, row_key, col_idx = src
        row = tables[csv_kind].get(row_key)
        if row is None or col_idx >= len(row):
            raise ValueError("%s: 在 %s 中找不到 %s 的第 %d 列" %
                             (name, csv_kind, row_key, col_idx))
        value = row[col_idx]
        note = None
        if opts.get("force_ref_abs"):
            value = abs(value)
            note = ("ANSYS COMBIN14 SMISC1 与弹簧对结构的反力反号，此处比幅值"
                    if kind == "spring_force" else
                    "FEM-Eigen 与 ANSYS SMISC 剪力符号约定相反，此处比幅值")
        if opts.get("force_ref_zero"):
            note = "ANSYS 实测 %g（舍入噪声），参考值取精确 0" % value
            value = 0.0
        chk = {
            "name": name,
            "kind": kind,
            "ref": value,
            "rtol": 0.001,
            "atol": opts.get("atol", 0.0),
            "source": "%s:%s%d.col%d" %
                      (CSV_NAMES[csv_kind].format(stem=stem),
                       "node" if csv_kind != "elem" else "elem", row_key, col_idx),
        }
        if kind == "nodal_disp":
            chk["node"] = key
            chk["dof"] = fem_sel
        elif kind == "reaction":
            chk["node"] = key
            chk["dof"] = fem_sel
        elif kind == "force":
            chk["elem"] = key
            chk["col"] = fem_sel
        elif kind == "spring_force":
            chk["node"] = key
            chk["k"] = float(opts["k"])
            chk["value"] = float(opts.get("spring_value", 0.0))
        if opts.get("force_ref_abs"):
            chk["abs"] = True
        if note:
            chk["note"] = note
        # 字典键顺序：把 name/kind/定位字段放前面，ref/容差放后面
        order = ["name", "kind", "node", "dof", "elem", "col", "abs", "k", "value",
                 "ref", "rtol", "atol", "source", "note"]
        checks.append({k: chk[k] for k in order if k in chk})

    return {
        "id": spec["id"],
        "dat": spec["dat"],
        "inp": spec["inp"],
        "exe": spec["exe"],
        "description": spec["description"],
        "load_case": spec.get("load_case", "equal-span continuous beam under full UDL"),
        "reference_source": (
            "ANSYS %s / %s 批处理实测，冻结 CSV 见 cases/%s_ansys_*.csv；"
            "由 build_reference.py 从这些 CSV 抽取，非解析解、非 FEM-Eigen 自测值。"
        ) % (spec["ansys"]["release"], spec["ansys"]["element"], stem),
        "ansys": dict(spec["ansys"], csv_columns={k: v for k, v in CSV_LABEL.items()}),
        "same_source_apdl": (
            "cases/%s.inp 与 cases/%s.dat 由同一个 BridgeSpec 对象生成"
            "（BridgeStudio/ansys_bridge.py），节点/单元/材料/截面/约束/荷载完全一致。"
        ) % (stem, stem),
        "tolerances": {"default_rtol": 0.001, "comment": "|got-ref| <= atol + rtol*|ref|"},
        "checks": checks,
    }


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default=str(here / "cases"))
    ap.add_argument("--reference", default=str(here / "reference"))
    ap.add_argument("--check", action="store_true",
                    help="只校验现有 JSON 与 CSV 是否一致，不写文件")
    args = ap.parse_args()

    cases_dir, ref_dir = Path(args.cases), Path(args.reference)
    ref_dir.mkdir(parents=True, exist_ok=True)
    failed = 0
    for stem, spec in CASES.items():
        expected = build_case(spec, cases_dir)
        out = ref_dir / (stem + "_ansys.json")
        if args.check:
            if not out.exists():
                print("[FAIL] %s 参考 JSON 缺失" % stem)
                failed += 1
                continue
            actual = json.loads(out.read_text(encoding="utf-8"))
            if actual == expected:
                print("[PASS] %s 与冻结 CSV 一致（%d 检查项）" %
                      (stem, len(expected["checks"])))
            else:
                print("[FAIL] %s 与冻结 CSV 不一致，请重跑 build_reference.py" % stem)
                failed += 1
            continue
        out.write_text(json.dumps(expected, ensure_ascii=False, indent=2) + "\n",
                       encoding="utf-8")
        print("[WRITE] %s  %d 个检查项" % (out, len(expected["checks"])))
        for chk in expected["checks"]:
            print("    %-28s ref=%-16.8g  %s" %
                  (chk["name"], chk["ref"], chk["source"]))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
