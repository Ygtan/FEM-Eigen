"""CB3D29 动力 snap-through / Newmark-β benchmark 对比脚本。

验证五层：
  1. 冻结输入 / 输出 sha256 + 当前 dyn_probe.exe 重跑逐字节一致；
  2. 独立闭式锚点（von Mises 静力极限 P_cr、一致/集中等效质量）；
  3. 独立 SDOF RK4 时程锚点（阶跃荷载，三种质量/荷载幅值；
     位移与速度逐采样点比对）；
  4. 力学行为（低于静力极限的动力跃越、荷载时程差异、阻尼耗散、
     时间步加密收敛、能量非负/外功一致）；
  5. 负例与零污染（DYN=0 / 零步长 / 零质量 / 非法 β；无 T 行、
     无附加 S 键）。
不依赖 ANSYS；口径为 FEM 自洽（ANSYS 瞬态动力学不同源，见 README）。
"""
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_reference as refmod  # noqa: E402

CHECKS = 0
FAILED = []


def check(name, ok, detail=""):
    global CHECKS
    CHECKS += 1
    if ok:
        print(f"  [PASS] {name}" + (f"  {detail}" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  [FAIL] {name}" + (f"  {detail}" if detail else ""))


def close(a, b, rel=1e-6, abs_=1e-10):
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse(text):
    T, S = {}, {}
    other = []
    for line in text.splitlines():
        t = line.split(",")
        if not t:
            continue
        if t[0] == "T" and len(t) >= 13:
            T.setdefault(t[1], []).append({
                "step": int(t[2]), "t": float(t[3]), "lam": float(t[4]),
                "v": float(t[5]), "vel": float(t[6]), "acc": float(t[7]),
                "res": float(t[8]), "iter": int(t[9]),
                "kin": float(t[10]), "str": float(t[11]), "work": float(t[12]),
            })
        elif t[0] == "S" and len(t) >= 4:
            S.setdefault(t[1], {})[t[2]] = float(t[3])
        else:
            other.append(line)
    return T, S, other


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="目录，含 dyn_probe.exe")
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    a = ap.parse_args()

    bin_dir = Path(a.bin)
    cases = Path(a.cases)
    ref_path = Path(a.reference) / "cb3d29_dynamic_snapthrough.json"
    inp = cases / "cb3d29_dynamic_snapthrough.dyn"
    frozen = cases / "cb3d29_dynamic_snapthrough_fem.csv"
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    an = ref["analytic"]

    exe = bin_dir / "dyn_probe.exe"
    check("dyn_probe.exe 存在", exe.is_file(), str(exe))
    if not exe.is_file():
        return 1
    check("输入 .dyn 在档", inp.is_file())
    check("冻结 CSV 在档", frozen.is_file())
    check("参考 JSON 在档", ref_path.is_file())

    check("输入 sha256 与冻结一致", sha256(inp) == ref["sha256"]["input"])
    check("冻结 CSV sha256 与参考一致", sha256(frozen) == ref["sha256"]["fem_csv"])

    frozen_text = frozen.read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "fresh.csv"
        r = subprocess.run([str(exe), str(inp), str(out)],
                           capture_output=True, text=True, timeout=600)
        check("当前 exe 重跑退出码 0", r.returncode == 0, r.stderr.strip()[:200])
        fresh_text = out.read_text(encoding="utf-8") if out.is_file() else ""
    check("当前 exe 重跑与冻结 CSV 逐字节一致",
          fresh_text == frozen_text,
          f"frozen={len(frozen_text)} fresh={len(fresh_text)}")

    T, S, other = parse(fresh_text)
    check("无异常记录行（只含 T/S）", not other, str(other[:2]))
    all_ids = {"dyn_snap", "dyn_lump", "dyn_low", "dyn_ramp", "dyn_sine",
               "dyn_damp", "dyn_fine", "dyn_off", "dyn_baddt",
               "dyn_badmass", "dyn_badbeta"}
    check("全部 11 个 case 均有 S 汇总", all_ids <= set(S),
          f"missing={sorted(all_ids - set(S))}")

    # ---- 2. 独立闭式锚点 ----
    vcr = refmod.limit_v()
    pcr = refmod.p_of_v(vcr)
    check("闭式 v_cr 独立重算一致", close(vcr, an["v_cr"], rel=1e-9),
          f"{vcr:.10f}")
    check("闭式 P_cr 独立重算一致", close(pcr, an["P_cr"], rel=1e-9),
          f"{pcr:.6f}")
    # PREF 写成 7 位有效数字，比例允许 1e-6 舍入差
    check("P_snap/P_cr ≈ 0.85", abs(an["ratio_snap"] - 0.85) < 1e-6,
          f"{an['ratio_snap']:.9f}")
    check("P_low/P_cr ≈ 0.60", abs(an["ratio_low"] - 0.60) < 1e-6,
          f"{an['ratio_low']:.9f}")
    L0 = refmod.math.hypot(refmod.B, refmod.H)
    rhoA = refmod.DENS * refmod.A
    m_c = 2.0 * rhoA * L0 / 3.0
    m_l = rhoA * L0
    check("一致质量闭式 m=2·rhoA·L0/3 一致", close(m_c, an["m_consistent"], rel=1e-12))
    check("集中质量闭式 m=rhoA·L0 一致", close(m_l, an["m_lumped"], rel=1e-12))
    check("m_lumped/m_consistent = 1.5",
          close(m_l / m_c, 1.5, rel=1e-12))
    with_mass = [c for c in all_ids if c != "dyn_badmass"]  # 负例组 DENS=0
    check("S,m_consistent 与闭式一致",
          all(close(S[c]["m_consistent"], m_c, rel=1e-10) for c in with_mass))
    check("S,m_lumped 与闭式一致",
          all(close(S[c]["m_lumped"], m_l, rel=1e-10) for c in with_mass))
    # C++/Python 黄金分割在极值附近 v 的条件数差，放宽容差到 1e-8
    check("S,analytic_v_cr 与闭式一致",
          all(close(S[c]["analytic_v_cr"], vcr, rel=1e-8) for c in with_mass))
    check("S,analytic_P_cr 与闭式一致",
          all(close(S[c]["analytic_P_cr"], pcr, rel=1e-10) for c in with_mass))
    check("一致质量矩阵 mass_total = 2m_c",
          all(close(S[c]["mass_total"], 2.0 * m_c, rel=1e-10)
              for c in ("dyn_snap", "dyn_low", "dyn_ramp", "dyn_sine",
                        "dyn_damp", "dyn_fine")))
    check("集中质量矩阵 mass_total = 2m_l",
          close(S["dyn_lump"]["mass_total"], 2.0 * m_l, rel=1e-10))

    # ---- 3. 独立 RK4 时程锚点 ----
    nsteps = an["steps"]
    rk = ref["rk4"]
    check("参考 RK4 每个 case 冻结 21 个采样点",
          all(len(rk[c]) == 21 for c in ("dyn_snap", "dyn_lump", "dyn_low")))
    recomputed = {
        "dyn_snap": refmod.rk4_step_load(m_c, refmod.PREF_SNAP, nsteps),
        "dyn_lump": refmod.rk4_step_load(m_l, refmod.PREF_SNAP, nsteps),
        "dyn_low": refmod.rk4_step_load(m_c, refmod.PREF_LOW, nsteps),
    }
    check("参考 RK4 采样点可由独立脚本重现",
          all(max(abs(row[1] - recomputed[c][i * refmod.SAMPLE][0]),
                  abs(row[2] + recomputed[c][i * refmod.SAMPLE][1]))
              < 1e-9
              for c in recomputed
              for i, row in enumerate(rk[c])
              if (i * refmod.SAMPLE) in recomputed[c]))

    for cid, key in (("dyn_snap", "dyn_snap"),
                     ("dyn_lump", "dyn_lump"),
                     ("dyn_low", "dyn_low")):
        rows = T.get(cid, [])
        check(f"{cid}: T 行数 = steps+1",
              len(rows) == nsteps + 1, f"{len(rows)}")
        check(f"{cid}: 全部步收敛",
              S[cid].get("converged") == 1.0 and S[cid].get("error_empty") == 1.0
              and S[cid].get("exception") == 0.0)
        check(f"{cid}: 阶跃荷载 λ=1（步1起）",
              all(close(r["lam"], 1.0, abs_=1e-12) for r in rows[1:]))
        check(f"{cid}: v(0)=0",
              close(rows[0]["v"], 0.0, abs_=1e-12) if rows else False)
        check(f"{cid}: 时间步均匀且严格递增",
              all(close(rows[i]["t"] - rows[i - 1]["t"], an["dt"], rel=1e-8)
                  for i in range(1, len(rows))))
        worst_v = max(abs(rows[i * refmod.SAMPLE]["v"] - row[1])
                      for i, row in enumerate(rk[key])
                      if i * refmod.SAMPLE < len(rows))
        worst_vel = max(abs(rows[i * refmod.SAMPLE]["vel"] - row[2])
                        for i, row in enumerate(rk[key])
                        if i * refmod.SAMPLE < len(rows))
        check(f"{cid}: FEM 位移逐点对齐独立 RK4（<2e-3 m）",
              worst_v < 2.0e-3, f"worst={worst_v:.3e} m")
        check(f"{cid}: FEM 速度逐点对齐独立 RK4（<5 m/s）",
              worst_vel < 5.0, f"worst={worst_vel:.3f} m/s")
        check(f"{cid}: Newton 迭代 1~maxIter（步1起）",
              all(1 <= r["iter"] <= S[cid]["MAXITER"] for r in rows[1:]))
        check(f"{cid}: 自由残差 < tol·|Fref|",
              all(r["res"] <= S[cid]["TOL"] * refmod.PREF_SNAP * 1.01
                  for r in rows))
        check(f"{cid}: 动能/应变能非负",
              all(r["kin"] >= 0.0 and r["str"] >= 0.0 for r in rows))

    # ---- 4. 动力行为 ----
    check("阶跃荷载幅值低于静力极限",
          refmod.PREF_SNAP < pcr, f"{refmod.PREF_SNAP:.0f} < {pcr:.0f}")
    check("动力 snap-through 发生（v_max 远超 v_cr）",
          S["dyn_snap"]["v_max"] > 0.6,
          f"v_max={S['dyn_snap']['v_max']:.4f} v_cr={vcr:.4f}")
    check("snap 时刻落在 (0, t_end) 内",
          0.0 < S["dyn_snap"]["snap_time"] < an["t_end"],
          f"t_snap={S['dyn_snap']['snap_time']:.5f}")
    # 离散时间步首次越阈会过冲，允许 0.3~0.31
    check("snap 阈值位移落在 0.3~0.31（离散过冲）",
          0.3 <= abs(S["dyn_snap"]["snap_disp"]) <= 0.31,
          f"{abs(S['dyn_snap']['snap_disp']):.6f}")
    check("S,snap_time 与 T 行首次越阈一致",
          any(r["t"] == S["dyn_snap"]["snap_time"] and abs(r["v"]) >= 0.3
              for r in T["dyn_snap"]))
    check("低荷载不跃越（v_max < 0.25）",
          S["dyn_low"]["v_max"] < 0.25,
          f"v_max={S['dyn_low']['v_max']:.4f}")
    check("低荷载无 snap_time", S["dyn_low"]["snap_time"] < 0.0)
    check("跃越幅值远大于低荷载幅值",
          S["dyn_snap"]["v_max"] > S["dyn_low"]["v_max"] + 1.0)
    check("一致/集中质量 snap 时刻不同（质量敏感性）",
          abs(S["dyn_snap"]["snap_time"] - S["dyn_lump"]["snap_time"]) > 5.0e-4,
          f"{S['dyn_snap']['snap_time']:.5f} vs {S['dyn_lump']['snap_time']:.5f}")
    check("集中质量仍发生 snap-through",
          S["dyn_lump"]["snap_time"] > 0.0 and S["dyn_lump"]["v_max"] > 0.6)
    check("爬升荷载不跃越", S["dyn_ramp"]["snap_time"] < 0.0
          and S["dyn_ramp"]["v_max"] < 0.2)
    check("半正弦荷载不跃越", S["dyn_sine"]["snap_time"] < 0.0
          and S["dyn_sine"]["v_max"] < 0.2)
    check("阶跃幅值远大于爬升幅值",
          S["dyn_snap"]["v_max"] > S["dyn_ramp"]["v_max"] + 1.0)
    check("阻尼工况收敛", S["dyn_damp"]["converged"] == 1.0)
    check("阻尼工况仍跃越", S["dyn_damp"]["snap_time"] > 0.0)
    check("阻尼降低峰值动能",
          S["dyn_damp"]["kinetic_max"] < S["dyn_snap"]["kinetic_max"])
    check("阻尼降低峰值位移",
          S["dyn_damp"]["v_max"] < S["dyn_snap"]["v_max"])
    check("加密时间步 npoints=801",
          len(T["dyn_fine"]) == 2 * nsteps + 1)
    check("加密时间步收敛", S["dyn_fine"]["converged"] == 1.0)
    check("dt 加密 v_max 收敛（<2e-3 m）",
          abs(S["dyn_fine"]["v_max"] - S["dyn_snap"]["v_max"]) < 2.0e-3,
          f"Δ={abs(S['dyn_fine']['v_max']-S['dyn_snap']['v_max']):.3e}")
    check("dt 加密 v_final 收敛（<2e-3 m）",
          abs(S["dyn_fine"]["v_final"] - S["dyn_snap"]["v_final"]) < 2.0e-3)
    check("阶跃外功 ≈ P_ref·v_final（相对 2%）",
          all(close(S[c]["work_final"], S[c]["PREF"] * abs(S[c]["v_final"]),
                    rel=2e-2)
              for c in ("dyn_snap", "dyn_lump", "dyn_low", "dyn_damp")))
    check("初加速度 a(0) = −P_ref/m（阶跃）",
          close(T["dyn_snap"][0]["acc"], -refmod.PREF_SNAP / m_c, rel=1e-8),
          f"{T['dyn_snap'][0]['acc']:.3f}")

    # ---- 5. 负例 / 零污染 ----
    for c in ("dyn_off", "dyn_baddt", "dyn_badmass", "dyn_badbeta"):
        check(f"{c}: 无 T 行", not T.get(c))
    check("dyn_off: DYN=0 标记", S["dyn_off"]["dyn"] == 0.0)
    check("dyn_off: 无 npoints 等运行时键", "npoints" not in S["dyn_off"])
    check("dyn_baddt: error 非空、converged=0",
          S["dyn_baddt"]["error_empty"] == 0.0
          and S["dyn_baddt"]["converged"] == 0.0)
    check("dyn_badmass: error 非空、converged=0",
          S["dyn_badmass"]["error_empty"] == 0.0
          and S["dyn_badmass"]["converged"] == 0.0)
    check("dyn_badbeta: error 非空、converged=0",
          S["dyn_badbeta"]["error_empty"] == 0.0
          and S["dyn_badbeta"]["converged"] == 0.0)
    anchor_keys = {"dyn", "E", "A", "B", "H", "L0", "PREF", "DENS",
                   "ADDMASS", "analytic_v_cr", "analytic_P_cr",
                   "analytic_lambda_cr", "rhoA", "m_consistent", "m_lumped"}
    check("dyn_off 只输出锚点 S 键（无 DT/TEND 等运行时键）",
          set(S["dyn_off"]) == anchor_keys,
          f"extra={sorted(set(S['dyn_off'])-anchor_keys)}")

    print(f"\nALL CHECKS PASSED ({CHECKS}/{CHECKS})" if not FAILED
          else f"\nFAILED {len(FAILED)}/{CHECKS}: {FAILED}")
    return 0 if not FAILED else 1


if __name__ == "__main__":
    sys.exit(main())
