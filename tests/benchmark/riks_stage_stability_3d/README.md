# 第 27 项 CTest：`riks_stage_stability_3d`

**第二十次（2026-11-16）** —— 在 ③b 施工阶段中间态起点上，新增
**Wempner–Riks 正交轨迹法**作为 opt-in 第二路径约束，并用独立闭式解 +
同问题双方法交叉验证闭环。BS 规格 §2.5 里「位移控制 / Riks」挂账项的
**Riks 部分**由此收口；位移控制、分支屈曲、初始缺陷、材料非线性、
动力 snap-through 仍未覆盖，故 **B3 仍 `P` 不 `Y`**。

## 这一项验证什么

1. **闭式独立锚点**：von Mises 两杆浅桁架的解析 `P(v)` 逐点核验 RIKS
   平衡路径（不读求解器内部量）；
2. **同问题双方法**：同一 STAGE 起点、同一模型，球面弧长与正交弧长
   各自独立跟踪 `R(u,λ)=0`，按 `|ΔUY|` 插值交叉比较；
3. **跨基准锚点**：CB3D23 的 RIKS 路径终点与第 26 项冻结的 SPHERE
   终点一致；
4. **零污染**：第 20 / 26 项原有全部冻结输出逐字节不变；
5. **负例不静默**：STAGE 起点给未知 `PATHMETHOD` 直接
   `ERROR=unsupported_path_method`；REFERENCE / COMPLETION 入口仍
   **拒绝** RIKS，原有负例保持原样。

## Riks 口径（必须写清，避免误读）

本项目 `PATHMETHOD=RIKS` 指 **Wempner–Riks 正交轨迹法**：每轮校正增量
满足

```text
duR·prevDu + dλ·(duF·prevDu + ψ²·prevDlam) = 0
```

`dλ` 由**线性方程**唯一确定，无二次方程选根、无夹角启发式。

它**不是** ANSYS `ARCLEN` / ABAQUS `*STATIC,RIKS` 的球面型弧长，
两者不可当作同一算法对拍；后者对应本项目的 `PATHMETHOD=SPHERE`
（Crisfield 球面约束 `duᵀdu + ψ²dλ² = ds²`）。第一步无前序增量时
RIKS 一律用球面约束启动，因此球面默认路径逐位不变。

残差装配、一致切线、LU 分解、约束施加、步长自适应全部复用；唯一
改动是 `stability_path_3d.hpp` 校正步里 `dλ` 的约束方程。

> 实现注记：正交约束必须落在**校正增量**上。若错误地落在「本步总
> 增量」上，预测点沿当前切线走（恰是正交面的法向），校正会把整步
> 增量投影掉，路径原地不动 —— 实测 2000 步 `CTRL` 恒为
> `-0.00076374043`。现实现已避开该陷阱。

## 算例构型

| 文件 | 说明 |
| --- | --- |
| `cases/vm_stage_riks.dat` | von Mises 两杆，STAGE-1 预载 `P0=1e6 N`，`PATHMETHOD=RIKS` |
| `cases/vm_stage_sphere.dat` | 同模型同起点，`PATHMETHOD=SPHERE`（双方法对照） |
| `cases/cb3d23_riks.dat` | 第 26 项 CB3D23 三杆浅拱改 `PATHMETHOD=RIKS` |
| `cases/st_bad_method.dat` | `PATHMETHOD=FOO` → `ERROR=unsupported_path_method` |
| `cases/*_stage_stability.txt` | 冻结输出（mingw64 工具链产出） |
| `reference/cb3d24_riks_stage.json` | 3109 checks + 20 strings + 20 criteria + 8 文件清单 |

von Mises 两杆：`B=2 m, H=0.5 m, E=2.1e11 Pa, A=1e-3 m²`，
节点 1/3 全固、顶点 2 约束 `UZ`，`STRAIN=1`（对数应变）。
`*STABILITY_LOAD` 给顶点 `FY=−2e6 N`（`PREF`），`CTRL=(2,2)`、
`CTRLTARGET=0.8`、`DS0=0.005` / `DSMAX=0.01`、`PSI=1`。

闭式解析（STRAIN=1）：

```text
L(v) = hypot(B, H−v)
P(v) = 2·EA·ln(L0/L)·(H−v)/L          L0 = hypot(B, H)
λ(v) = [P(v) − P0] / PREF             v = V0 − CTRL
```

`V0 = 0.12443409386935` 来自 `STATE CTRL0`；起点闭合残差
`RES0 = P(V0) − P0 = 6506.4665684 N` 也由同一闭式给出，与输出逐位吻合
（这不是求解器误差，而是阶段 1 平衡态的解析荷载差）。

## 冻结结果

```text
vm_stage_riks     KIND=wempner_riks_orthogonal BASIS=stage_state STAGE=1 NELEM=2 NDOF=9
                  184 步 / 185 点，λ_final=-0.85723836，CTRL_end=-0.80111
                  闭式逐点 worst|Δλ| = 4.57e-14（184 点）
vm_stage_sphere   KIND=spherical_arc_length，184 步 / 185 点
                  双方法按 CTRL 插值 worst|Δλ| = 1.47e-4（≥150 点）
cb3d23_riks       STAGE=2 NELEM=3 NDOF=12，133 步 / 134 点
                  λ_final=2.5757758309，与第 26 项冻结 SPHERE 终点一致
负例              PATHMETHOD=FOO -> ERROR=unsupported_path_method，NPOINT=0
```

容差：闭式 `|λ−λ_closed| ≤ 1e-8`、双方法 `≤ 3e-4`、跨基准 `≤ 1e-5`、
重放要求输出 **逐字节一致**（不使用 `replay_tolerance.py`，工具链
变更须重新冻结，见 `tests/benchmark/README.md`）。

## 实测结果

```
[PASS] 63/63 checks: closed_form_worst=4.508e-14,
       cross_method_worst=1.473e-04, cb_lambda_end=2.5757758309
build_reference_riks_stage.py --check: 3109 checks + 20 strings + 20 criteria
第 20 / 26 项 CTest：均 Passed（SPHERE 冻结输出逐字节不变）
```

## 复现

```bash
# 生成 / 校验黄金值
python build_reference_riks_stage.py
python build_reference_riks_stage.py --check

# 离线守卫（--bin 指向含 cable_solver.exe 的构建目录）
python compare_riks_stage_3d.py --bin <build_dir> --cases cases --reference reference

# 第 27 项 CTest
ctest --test-dir <build_dir> -R riks_stage_stability_3d --output-on-failure
```

## 如实边界

- RIKS 只在 `STABASIS=STAGE` 入口 opt-in；默认 `PATHMETHOD=SPHERE`
  与 `STABILITY=0` 的全部既有输出逐位不变。
- 验证构型是 von Mises 两杆 / 三杆浅拱**研究构型**，不是真实桥梁
  施工阶段；不覆盖 CR 梁 / 索真实施工阶段、分支屈曲与分支切换、
  初始缺陷、材料非线性、动力 snap-through。
- `λ` 是增量荷载因子，不是规范稳定安全系数；JTG 3362 无 `λ`
  合格限值，本项不做「通过 / 不通过」判定。
- 与 ANSYS 的关系：本项**不跑 ANSYS**。ANSYS `ARCLEN` 是球面型，
  与本项目 RIKS 正交型不同源；如要做 ANSYS 对拍，应对
  `PATHMETHOD=SPHERE`，证据在第 17 / 26 项。
