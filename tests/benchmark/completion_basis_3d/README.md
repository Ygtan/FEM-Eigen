# 第 24 项 CTest：`completion_basis_3d`

**第十九次（2026-11-15，a 阶段 ③ / `③a`）** —— 收口 ② 阶段登记的
`MODE=COMPLETION` 起点不同源缺陷，并核对施工阶段链路与成桥态稳定
主链路的共存性。本项**不跑 ANSYS**，完全离线。

> **编号口径**：本项覆盖 `CB3D10` / `CB3D13` 的既有冻结输出，不是新的
> BridgeStudio 算例；黄金值 `case=CB3D10_BASIS_FIX`、文件名
> `completion_basis_3d.json`。`CB3D22` 专属于 BS 侧
> `cases/CB3D22_limit_point_vs_ansys/` 与 FE 侧第 25 项
> `bridge_limit_point_3d`（真实桥梁极值点 / 下降段），两者不复用编号。

## 缺陷与修复

`fem_3d_nl.hpp::solveCompletionState3D(Model3D m, ...)` 与
`iterateCableUnstrainedLengths(Model3D m, ...)` 都**按值接收模型**：

- 前者把索 `L0` 改成下料长度，并由 `iterateErectionGeometry()` 把参考
  构形改成制造线形 `X0 + ΔX`（`[ERECTED]`），返回的 `u` 是相对该构形的
  成桥位移；
- 后者把 `L0` 迭代到命中目标索力 `T*`；
- `cable_solver.cpp` 手里的 `m` 却仍是**设计线形 + `.dat` 原始 `L0`**。

`STABASIS=COMPLETION` 用这份外部模型装配恒定荷载 `F0` 并计算
`Fint(u0)`，起点不闭合。CB3D10 实测：

| 量 | 修复前 | 修复后 |
|---|---|---|
| `RES0` | `3.1492e6 N` | `7.4001e-3 N` |
| `RES0 / F0_NORM`（`F0_NORM=2.9572e6 N`） | `1.06` | `2.50e-9` |
| `[NEWTON] RESLAST` | `6.1573e-3 N` | `6.1573e-3 N` |
| 首步 `λ` | `-5.693e-3`（变负） | `+9.984e-3` |
| `LAMBDA_FINAL` | `2.6833` | `4.5760` |

修复：给 `CableLengthResult3D` / `ErectionGeometryResult3D` /
`CompletionStateResult3D` 增加 `Model3D basisModel`，回传与 `u` 同源的
最终模型；`cable_solver.cpp` 只在 `STABASIS=COMPLETION`（opt-in）时把
`basisModel` 交给 `formatCompletionStabilitySection()`。
`MODE=STATIC` 与悬索桥找形的主链路 `m` 本身同源，不替换；
`STABILITY=0` / `STABASIS=REFERENCE` 路径零改动。
`formatCompletionStabilitySection()` 另加自由度结构一致性检查，
不一致写 `ERROR=invalid_initial_state`，不静默兜底。

## 四个 case 与独立锚点

| case | MODE | 基准模型 | 用途 |
|---|---|---|---|
| `cb3d10_completion_basis` | `COMPLETION` | 制造线形 + 下料 `L0` | 修复主体 |
| `cb3d10_static_basis` | `STATIC` | 设计线形 + 冻结 `[CABLE]` 显式 `L0` | **独立物理锚点** |
| `cb3d10_default_basis` | `COMPLETION` | `STABILITY=0` | 默认零污染 |
| `cb3d13_erection_basis` | `COMPLETION`（悬索桥） | 悬索桥成桥态闭环 | **施工阶段链路共存核对** |

`MODE=COMPLETION` 与 `MODE=STATIC` 是两条**独立成桥态闭环**：参考构形
不同（制造线形 vs 设计线形），`CTRL0` 分别是 `-0.0016328 m` 与
`-0.0044246 m`，本来就不该相等；可比较的是路径增量
`ΔCTRL = CTRL − CTRL0` 与 `λ`。实测逐点对拍：

```
NPOINT                    129 / 129
worst |Δλ|                1.081e-6
worst |Δ(ΔCTRL)|          3.131e-6 m
LAMBDA_FINAL              4.5760060833765 / 4.5760050023174（相对差 2.36e-7）
F0_NORM                   2957212.2214688（两者完全相同）
RES0/|F0|                 2.502e-9 / 2.128e-9
```

这不是「把黄金值抄一遍」：`MODE=STATIC` 的 `L0` 从历史冻结 `[CABLE]`
段独立取出写进 `.dat`，与 `COMPLETION` 内部的制造线形迭代没有任何
共享状态；两条路径能重合到 `1e-6` 量级，说明 `u0 / F0 / Finc /
Fint(u0) / 切线刚度` 五者确实同源。

其它判据：

- `RES0/|F0| < 1e-8` 且与 `[NEWTON] RESLAST` 同量级（内层 Newton
  独立残差）；`CTRL0 ≡ [U](1016,3)`（`< 1e-12`）；首步 `λ > 0`；
- `FINC_NORM ≡ |*STABILITY_LOAD|`：从 `.dat` 独立读段算范数，不读求解器；
- `λ` 路径单调、每点 `RES < 100 N`、`CONVERGED=1`、`ERROR=none`；
- 闭环段齐备：`[CABLE] / [U] / [ERECTED] / [LENGTH] / [NEWTON] /
  [BUCKLING] / [STABILITY]`，`[BUCKLING] ERROR=nonlinear_state_not_supported`
  原样保留；`COMPLETION` 的 `[LENGTH] KIND=erection_geometry STATUS=ok`，
  `STATIC` 写 `KIND=none STATUS=not_applicable`；
- 默认零污染：无 `[STABILITY]`，离线重跑逐行一致，与历史 CB3D10
  冻结基线的段集合 / `[U](1016,3)` / `[NEWTON] RESLAST` 一致。
  跨编译器的 `1e-9` 量级浮点噪声按 `1e-5` 相对容差吸收，
  **不做**跨二进制逐字符比对。

## 施工阶段链路接入核对（CB3D13）

`cable_solver.cpp` 的施工阶段链路（`*ERECTION` / `SYSTEM=SELF`）是**解析
索多边形找形**：`findErectionState()` 二分主缆水平力 `H` 命中固定下料
长度，逐级给出 `[ERECTION] / [ERECTION_SHAPE] / [ERECTION_CHAIN] /
[ERECTION_LOAD] / [ERECTION_FORCE]` 五段；它**不走**
`solveModel3DNonlinear()`，因此每个施工阶段中间态没有 FEM 收敛位移
`u0`、也没有一致切线刚度。

在 CB3D13（自锚式悬索桥、5 个施工阶段）上同时开
`STABASIS=COMPLETION`，实测：

- 段集合 = 历史基线 + `[STABILITY]`（`[NEWTON]` 之后、`[END]` 之前），
  5 个 `[ERECTION*]` 段与冻结基线**逐字符一致**；
- `[STABILITY] BASIS=completion_state` 起点闭合
  `RES0=2.011e-3 N`、`RES0/|F0| = 1.112e-9`，
  `CTRL0 ≡ [U](1014,3)`，首步 `λ>0`，26 个平衡点、
  `LAMBDA_FINAL=0.456467`、`CONVERGED=1`；
- `[NEWTON] RESLAST` 是迭代停止残差，跨工具链有 `8e-3` 相对浮点
  噪声，只校核同量级；物理口径由 5 个施工阶段段逐位一致与
  「同工具链重跑逐行一致」共同保证。

**核对结论（如实边界）**：施工阶段链路与成桥态稳定主链路是**并列输出、
互不污染**，已具备「在含施工阶段的算例上做成桥态全过程稳定」的能力；
但 `[STABILITY]` 的起点仍是**成桥态**，不是某个施工阶段中间态。
「每个施工阶段都具备全过程稳定分析」需要为各阶段另建 FEM 平衡态
（阶段激活单元 / 阶段荷载 / 阶段 `L0` / 阶段 Newton 收敛），属实质新增，
**本项未覆盖，不得据此把 B3 涂成 `Y`**。

合计 **110 项检查**，四个 `.dat` 离线重跑（归一化 `elapsed_ms` 后）
逐行一致。

> **工具链口径（③a 全量回归时补记）**：本项最初用 `ucrt64/bin/c++.exe`
> 的构建生成冻结输出，与仓库历史冻结证据的 `mingw64/bin/g++.exe`
> 不同；两套工具链同为 GCC 16.2 / `-O3 -DNDEBUG`，但牛顿残差、
> 模态向量末位浮点不同，导致异工具链重放时**非逐点行**也漂移
> （约 1100 行、量级 `1e-19`~`1e-3`）。已用规范 mingw64 构建
> `_d_build/cable_solver.exe` **重新生成本项四个 case 与黄金值**，
> 恢复严格逐字节重放守卫，不做跨二进制容差兜底。跨工具链的结构化
> 比较只用于第 25 项 `bridge_limit_point_3d`（见
> `tests/benchmark/replay_tolerance.py`）。

## 生成器

```
python build_reference_basis.py <cable_solver.exe> <cases_dir> [reference_dir]
```

`cases/cb3d10_baseline_cable.txt` 与 `cases/cb3d13_baseline_cable.txt`
分别是从 `nonlinear_ansys_3d/cases/cb3d10_completion_cable.txt`、
`suspension_erection_3d/cases/cb3d13_balanced_erection_cable.txt`
复制的历史冻结基线。正常回归不需要重跑生成器（会覆写冻结输出）。

CB3D13 的稳定性键必须落在 `*CABLE_BRIDGE` 段内（`*SUSP_BRIDGE` 之前）；
插到文件末尾 `*ERECTION` 段后的 `*CABLE` 前会被解析成施工阶段参数，
`STABILITY=1` 根本不生效。

## 能力边界

本项只修复**成桥态稳定路径的起点同源性**，并核对施工阶段链路的
**共存性**（见上节）。它不覆盖分支屈曲 / 分支切换、初始缺陷、
材料非线性、**施工阶段中间态的全过程稳定分析**，也不改变
`MODE=COMPLETION` 其余输出段的历史口径（`[U]` / `[ERECTED]` /
`[REACTION]` 等仍按主链路 `m` 写出）。B3 是否推 `Y` 见
`BridgeStudio/B类后续项接入规格.md` §2.5 / §2.6，本项**不**代替该结论。
