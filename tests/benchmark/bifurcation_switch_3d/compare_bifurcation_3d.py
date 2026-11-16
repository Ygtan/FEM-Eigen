"""CB3D28 分岔检测 / 分支切换 benchmark 对比脚本。

验证四层：
  1. 冻结输入 / 输出 sha256 + 当前 exe 重跑逐字节一致；
  2. 独立闭式锚点（反对称分岔点 v_bif/λ_bif、对称极值点 v_sym/λ_sym）；
  3. 分支切换行为（反对称进入分支、承载力下降、幅值敏感性）；
  4. 负例 / 零污染（对称模态、零幅值、越界模态、BIFURC=0 全部不漂移）。
不依赖 ANSYS；口径为 FEM 自洽（与 ANSYS 的 RIKS 类算法不同源，见 README）。
"""
import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

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
    P, M, S = {}, {}, {}
    for line in text.splitlines():
        t = line.split(",")
        if not t:
            continue
        if t[0] == "P" and len(t) >= 9:
            P.setdefault(t[1], []).append({
                "step": int(t[2]), "lam": float(t[3]), "v": float(t[4]),
                "res": float(t[5]), "iter": int(t[6]), "ds": float(t[7]),
                "piv": float(t[8]), "raw": line,
            })
        elif t[0] == "M" and len(t) >= 4:
            M.setdefault(t[1], []).append((int(t[2]), float(t[3])))
        elif t[0] == "S" and len(t) >= 4:
            S.setdefault(t[1], {})[t[2]] = float(t[3])
    return P, M, S


def p_lines(text, cid):
    """只取 step 之后的数值字段，去掉 case id 前缀再比较路径数值。"""
    return [",".join(ln.split(",")[2:])
            for ln in text.splitlines() if ln.startswith(f"P,{cid},")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="目录，含 snap_probe.exe")
    ap.add_argument("--cases", required=True)
    ap.add_argument("--reference", required=True)
    a = ap.parse_args()

    bin_dir = Path(a.bin)
    cases = Path(a.cases)
    ref_path = Path(a.reference) / "cb3d28_bifurcation_switch.json"
    inp = cases / "cb3d28_bifurcation_switch.snap"
    frozen = cases / "cb3d28_bifurcation_switch_fem.csv"
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    cf = ref["closed_form"]

    exe = bin_dir / "snap_probe.exe"
    check("snap_probe.exe 存在", exe.is_file(), str(exe))
    if not exe.is_file():
        return 1

    check("输入 sha256 与冻结一致", sha256(inp) == ref["sha256"]["input"])
    check("冻结 CSV sha256 与参考一致", sha256(frozen) == ref["sha256"]["frozen_csv"])

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

    P, M, S = parse(fresh_text)
    cids = ["bif_off", "bif_anti", "bif_small",
            "bif_sym", "bif_zero", "bif_badmode"]
    for cid in cids:
        check(f"{cid} 产出路径点", len(P.get(cid, [])) > 100,
              f"npoints={len(P.get(cid, []))}")

    # ---- 闭式锚点 ----
    check("反对称分岔闭式有两个正根", len(cf["bifurcation_roots_h"]) >= 2,
          str([round(x, 6) for x in cf["bifurcation_roots_h"]]))
    check("反对称分岔先于对称极限点", cf["bif_before_limit_point"],
          f"v_bif={cf['v_bif']:.6f} < v_sym={cf['v_sym']:.6f}")

    v_bif, lam_bif = cf["v_bif"], cf["lambda_bif"]
    v_sym, lam_sym = cf["v_sym"], cf["lambda_sym"]

    # 主路径（BIFURC=0）应越过对称极值点，λ_max 对齐闭式 λ_sym。
    off_lmax = S["bif_off"].get("lambda_max_fit", 0.0)
    off_vmax = S["bif_off"].get("v_at_lambda_max_fit", 0.0)
    check("主路径 λ_max 对齐闭式对称极值 λ_sym",
          close(off_lmax, lam_sym, rel=2e-3),
          f"FEM={off_lmax:.6f} closed={lam_sym:.6f}")
    check("主路径 v(λ_max) 对齐闭式 v_sym",
          close(off_vmax, v_sym, abs_=5e-3),
          f"FEM={off_vmax:.6f} closed={v_sym:.6f}")
    check("主路径 P(v) 逐点满足对称解析",
          S["bif_off"].get("worst_dP_rel", 9.9) < 1e-9,
          f"worst_rel={S['bif_off'].get('worst_dP_rel')}")

    # ---- 分岔检测（反对称切换）----
    an = S["bif_anti"]
    check("bif_anti 检测到分岔", an.get("bif_detected") == 1.0)
    check("bif_anti 实际切换分支", an.get("bif_switched") == 1.0)
    check("bif_anti 切换步号为正", an.get("bif_step", 0) > 0,
          f"step={an.get('bif_step')}")
    check("bif_anti 一致切线特征值确认临界点",
          0.0 < an.get("bif_eig_ratio", 9.9) < 1e-3,
          f"eig_ratio={an.get('bif_eig_ratio'):.6e}")
    check("bif_anti 切换点 pivotRatio 远小于 1",
          0.0 <= an.get("bif_pivot", 9.9) < 1e-5,
          f"pivot={an.get('bif_pivot'):.6e}")

    v_sw = -an.get("bif_ctrl_disp", 0.0)   # snap_probe: v = -ctrlDisp
    lam_sw = an.get("bif_lambda", 0.0)
    check("FEM 切换点 v 对齐闭式 v_bif",
          close(v_sw, v_bif, abs_=5e-3),
          f"FEM={v_sw:.6f} closed={v_bif:.6f}")
    check("FEM 切换点 λ 对齐闭式 λ_bif",
          close(lam_sw, lam_bif, rel=5e-3),
          f"FEM={lam_sw:.6f} closed={lam_bif:.6f}")

    # ---- 分支进入行为 ----
    amp_before = an.get("bif_max_amp_before", 0.0)
    amp_at = an.get("bif_amp_at_switch", 0.0)
    amp_after = an.get("bif_max_amp_after", 0.0)
    amp_final = abs(an.get("bif_amp_final", 0.0))
    check("切换后模态幅值显著放大",
          amp_after > 3.0 * max(amp_before, 1e-8),
          f"before={amp_before:.5f} after={amp_after:.5f}")
    # amp_at 已含切换当步的路径自然推进（同一推进量在两个 case 中相同），
    # 故用「两 case 的 Δamp 之差」精确核对施加幅值之差，避免把推进量算进来。
    check("切换瞬间含路径推进，增量须大于施加幅值",
          amp_at - amp_before > 0.02,
          f"Δamp={amp_at - amp_before:.6f} > 0.02")
    check("分支终点保持显著破缺", amp_final > 1.0,
          f"|amp_final|={amp_final:.5f}")

    an_lmax = an.get("lambda_max_fit", 0.0)
    check("分岔后分支承载力显著低于主路径",
          an_lmax < 0.5 * off_lmax,
          f"branch={an_lmax:.4f} main={off_lmax:.4f}")
    check("分岔后分支 v 方向反转（离开对称下降段）",
          S["bif_anti"].get("v_final", 0.0) < 0.0,
          f"v_final={S['bif_anti'].get('v_final'):.5f}")
    check("分支不再满足对称路径 P(v) 解析",
          an.get("worst_dP_rel", 0.0) > 1e-3,
          f"worst_rel={an.get('worst_dP_rel'):.6e}")

    # ---- 幅值敏感性（切换瞬间）----
    sm = S["bif_small"]
    check("bif_small 检测到并切换",
          sm.get("bif_detected") == 1.0 and sm.get("bif_switched") == 1.0)
    d_small = sm.get("bif_amp_at_switch", 0.0) - sm.get("bif_max_amp_before", 0.0)
    check("小扰动切换瞬间增量含同一路径推进",
          d_small > 0.005, f"Δamp={d_small:.6f} > 0.005")
    check("两 case 的增量差精确等于施加幅值差 0.015",
          close((amp_at - amp_before) - d_small, 0.015, abs_=1e-4),
          f"ΔΔamp={(amp_at - amp_before) - d_small:.6f} expect=0.015")

    # ---- 负例 ----
    sy = S["bif_sym"]
    check("bif_sym 对称模态不产生分支破缺",
          abs(sy.get("bif_amp_final", 9.9)) < 1e-6,
          f"amp_final={sy.get('bif_amp_final'):.3e}")
    # 对称模态扰动只改变采样位置（同样落在对称路径上），故用解析闭合
    # 残差 + 宽松承载力窗口判定，不要求与 off 采样点逐位一致。
    check("bif_sym 仍在同一对称平衡路径上",
          sy.get("worst_dP_rel", 9.9) < 1e-9 and
          close(sy.get("lambda_max_fit", 0.0), off_lmax, rel=2e-3),
          f"worst_rel={sy.get('worst_dP_rel'):.3e} lmax={sy.get('lambda_max_fit'):.5f}")

    ze = S["bif_zero"]
    check("bif_zero 检测到但不切换",
          ze.get("bif_detected") == 1.0 and ze.get("bif_switched") == 0.0)
    bm = S["bif_badmode"]
    check("bif_badmode 越界模态不切换",
          bm.get("bif_detected") == 1.0 and bm.get("bif_switched") == 0.0)

    for cid in ("bif_zero", "bif_badmode"):
        check(f"{cid} 路径 P 行与主路径逐位一致",
              p_lines(fresh_text, cid) == p_lines(fresh_text, "bif_off"))

    # ---- 零污染 ----
    check("bif_off 不输出任何 M 行", "bif_off" not in M)
    check("bif_off 无分岔诊断 S 行",
          not any(k.startswith("bif") for k in S["bif_off"]))
    check("bif_off 与既有冻结输出不漂移（P 行数稳定）",
          len(P["bif_off"]) == 801, f"n={len(P['bif_off'])}")
    check("所有平衡点残差在工程容差内",
          all(p["res"] < 0.2 for cid in cids for p in P[cid]))

    print("\n" + "=" * 78)
    print(f"  {CHECKS - len(FAILED)}/{CHECKS} checks passed")
    if FAILED:
        print("  FAILED:", ", ".join(FAILED))
    print("  [PASS] CB3D28 bifurcation detection & branch switching benchmark"
          if not FAILED else "  [FAIL] CB3D28 benchmark")
    print("=" * 78)
    return 0 if not FAILED else 1


if __name__ == "__main__":
    raise SystemExit(main())
