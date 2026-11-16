# 第 26 项 CTest：`stage_stability_3d`

**第二十次（2026-11-15，a 阶段 ③b）** —— **施工阶段中间态起点**的全过程
稳定 FEM ↔ ANSYS 对标。BS 侧正式算例为
`BridgeStudio/cases/CB3D23_stage_state_vs_ansys/`，本项是它的**纯离线
回归守卫**：15 份冻结证据同步到 `cases/`，跑 CTest 时既不重跑 ANSYS，
也不重新建模，只用当前 `cable_solver.exe` 重放两个 `.dat`。

## 为什么需要这一项

③a 的三处收口（起点同源修复、CB3D13 五个 `[ERECTION*]` 段与
`[STABILITY]` 并列输出、CB3D22 真实斜拉桥越过极值点）**起点全部是
成桥态**。`STABASIS=COMPLETION` 拿到的 `u0` 是找形 / 索长迭代收敛后的
全桥平衡位移，`F0` 是成桥态恒定荷载；即便模型里有施工阶段段，它们
也只是「在同一份结果文件里共存」，并没有把弧长路径的起点放到某个
施工阶段的中间平衡态上。③b 补的正是这个缺口：

- `fem_3d_nl.hpp::solveModel3DStagesNonlinear()` 逐阶段给出**真实共转
  列式全量平衡态** `u_stage`（含已激活单元、累计阶段荷载、当前 `L0`），
  但 `NLStageResult3D` 原先只回传 `u` / `dm` / `cableL0` / 激活集合，
  **不回传同源模型**。③b 增加 `Model3D model` 并在每阶段求解后写
  `row.model = mc`，让 `[STABILITY]` 能用与 `u_stage` 同一份模型装配
  `F0` 与 `Fint(u0)`。
- `cable_solver.cpp` 新增 `STABASIS=STAGE` + `STAGEID=<id>` 独立 opt-in：
  先跑完施工阶段链路，按 `STAGEID` 选中指定阶段的 `(model, u)`，
  再以 `u0 = u_stage` 走 Crisfield 球面弧长路径。未知阶段写
  `ERROR=stage_not_found`，**不静默降级**成成桥态或参考构形。
- `[STABILITY] BASIS=stage_state` 与 `BASIS=completion_state` /
  `BASIS=reference_configuration` 三者分列；新增 `STAGE=<id> NELEM=<n>
  NDOF=<n>` 行与 `[STAGES]` 逐阶段摘要，证明选中的确实是**中间态**。

## 构型（如实边界）

4 节点 / 3 根 `TRUSS` 的对称浅拱，`E = 2.1e11 Pa`、`A = 0.001 m²`，
节点 1/4 在 `(±2, 0, 0)`、节点 2/3 在 `(±0.7, 0.35, 0)`：

| 阶段 | 激活单元 | 新增阶段荷载 |
|---|---|---|
| 1 | 1, 3（两根边杆） | 节点 2 / 3 各 `FY = −10000 N` |
| 2 | 2（中间杆） | 无（累计恒载仍是阶段 1 的两条） |
| 3 | 无新增 | 节点 2 / 3 各 `FY = −5000 N` |

`STAGEID=2` 选**阶段 2**：此时 3 根杆全部激活、`NELEM=3`、
`NDOF=12`，但阶段 3 的荷载尚未施加——这是几何与荷载都处于中间态的
平衡点，不是成桥态。`*STABILITY_LOAD` 给节点 2 / 3 各
`FY = −20000 N`（`PINC`），`CTRL=(2, 2)`、`CTRLTARGET=0.005`、
`PATHMETHOD=SPHERE`、`PSI=1`、`DS0=0.01` / `DSMAX=0.02`。

约束：节点 1 / 4 全固；节点 2 / 3 的 `UZ` 与 `UX` 为 0。后者是 FEM
与 ANSYS **同源**加的侧向导向约束——3 杆浅拱在零轴力下 4 个平动
自由度只受 3 个杆长约束，是一阶运动机构，ANSYS 稀疏求解器会报
`Sparse solver minimum pivot ... at node 3 UY` 的 small pivot。
曾试过加对角三角化杆 `E,1,3`，结果 FEM 阶段求解器
`非线性求解未收敛：达到最大迭代次数`，故回退 3 杆 + 导向约束。

**这是为验证「起点是施工阶段中间态」构造的研究构型，不是实际桥梁
尺寸，也不是实际施工工序。** 它不覆盖 CR 梁 / 索的真实施工阶段、
分支屈曲与分支切换、初始缺陷、材料非线性或动力 snap-through。

## 求解口径

| | FEM | ANSYS 18.0 |
|---|---|---|
| 起点 | `STABASIS=STAGE STAGEID=2`：`solveModel3DStagesNonlinear()` 的阶段 2 平衡态 `u0` 与 `NLStageResult3D::model` | `step1 TIME=1` 求同一阶段平衡态（只建阶段 2 已激活的 3 个单元、只加累计恒载） |
| 恒定荷载 | `F0` = 阶段 2 累计节点荷载 + 弹性支承 + CABLE 自重 | `F,2,FY,-10000` / `F,3,FY,-10000` |
| 增量 | `R(u,λ) = F0 + λ·Finc − Fint(u)`，`Finc` 只来自 `*STABILITY_LOAD` | `step2` `F_end = F_orig + LAMMAX·PINC`（`−64000 N`），`TIME = 1 + LAMMAX = 3.7`，故 `λ = TIME − 1` |
| 跟踪 | Crisfield 球面弧长 `duᵀdu + ψ²dλ² = ds²` | `st_load` 荷载控制（`NSUBST,400`）；`st_arc` `ARCLEN,ON` 交叉验证 |
| 终止 | `CTRLTARGET=0.005`（`|ΔUY|`） | `LAMMAX=2.7` |

**不用反力求和**：`λ` 直接由 `TIME` 换算，避免完整结构上反力求和
口径不可靠。逐点对拍以递增自变量 `|ΔUY| = |UY − UY_step1|` 插值 `λ`
（本例路径单调，`λ` 未越过极值点），而不是按 `λ` 反查。

**应变档**：阶段荷载 / `PINC` / `CTRLTARGET` 取到轴向应变约 `2e-4`。
初版放大 10 倍时轴向应变约 `2e-3`，ANSYS `LINK180` 默认 Almansi
应变与 FEM 共转列式工程应变出现系统性差异
（`u0` 差 `2.89e-5 m`、`worst|Δλ|≈0.038`），缩小 10 倍后
`worst|Δλ|≈0.005`，属两套应变度量的残差而非力学回归。

## 冻结结果

```
[STAGES] [1, 2, 3]  NELEM 2 -> 3 -> 3；STAGE=2 NELEM=3 NDOF=12
FEM    BASIS=stage_state  STATUS=ok  NPOINT=134
       RES0/|F0| = 6.34e-13   CTRL0 = -0.000952173872（非零初态）
       lambda_end = 2.57577583   |dUY|_end = 0.00502513021（单调）
负例   STAGEID=99 -> ERROR=stage_not_found，NPOINT=0，不降级
ANSYS  st_load  step1 UY=-0.000951910681；step2 27 点，0 error
       st_arc   step1 UY=-0.000951910681；step2 25 点，0 error
       两 job 共同区间 |dUY|=[1.28988e-05, 0.00502513]，覆盖 FEM 路径 99.7%
       st_load worst|dlambda|=0.005124 mean=0.001891
       st_arc  worst|dlambda|=0.005109 mean=0.001881
BS 报告 72 checks / 0 failed；本项离线守卫 214 checks / 0 failed
```

## 独立锚点（`compare_stage_state_3d.py`，214 项检查）

1. **冻结证据完整性**：15 个文件逐个 `sha256` / 字节数核对，集合不得增删；
2. **黄金值自洽**：从冻结 `.dat` / `.txt` / `.inp` / `.out` / `.csv`
   用 `build_reference_stage_state.py` 的同一套解析函数重新派生，与
   `reference/cb3d23_stage_state.json` 逐字段相等（含 134 个 FEM 平衡点、
   2 条 ANSYS 曲线、`[STAGES]` 3 阶段摘要、负例与报告结论）；
3. **物理判据**：`BASIS=stage_state`、`STAGE=2 NELEM=3 NDOF=12`、
   阶段 1→2 激活单元递进、`RES0/|F0| < 1e-9`、`CTRL0 != 0`、
   FEM 路径单调且收敛、负例 `ERROR=stage_not_found` 且无平衡点；
4. **同源锚点 / ANSYS 证据**：APDL `N,` / `E,` / `F,` / `D,` 与
   `.dat` / `.cfg.json` 逐项一致（坐标最大偏差 `0`、4 节点 / 3 单元 /
   10 条约束；`loads_step1` 与 `loads_step2` 分列，不互相覆盖）；
   两个 job 各 2 次 `SOLVE`、`λ = TIME − 1`、`0 error`；
   按 `|ΔUY|` 插值 `λ` 逐点比较，逐点容差 `6e-3`、均值容差 `2.5e-3`；
5. **离线重放**：用当前二进制重跑 `cb3d23_stage_state.dat` 与
   `st_bad_stage.dat`，要求 `<stem>_stage_stability.txt` **逐字节一致**。
   本项冻结证据与历史黄金值同由 `mingw64/bin/g++.exe` 产出，
   因此不放宽为 `replay_tolerance.py` 的结构化比较；换 CRT 工具链
   必须重新冻结，见 `tests/benchmark/README.md`。

## 生成器

```
python build_reference_stage_state.py [--cases <dir>] [--out <reference dir>]
```

**不重跑 FEM / ANSYS**：`cases/` 是唯一冻结源，脚本只做哈希、派生指标
与同源锚点计算。新增 / 更新证据时先在 BS 侧跑
`cases/run_cases_stage_state.py --check-only --sync-to <cases>` 同步
15 个文件，再重新生成黄金值。正常回归不需要跑生成器。

生成器 `apdl_anchor()` 的 `F,<node>,FY,<val>` 需按**列表**收集：两段
求解会在同一节点写两次 `F,`，初版只留最后一个值导致 `step1` 恒载
被 `step2` 的 `−64000` 覆盖；现输出 `loads_step1` / `loads_step2` 两列。

## 能力边界

本项证明「`cable_solver.exe` 的球面弧长路径可以把起点放到**指定施工
阶段的共转列式全量平衡态**，且该路径与 ANSYS 两段求解定量一致」。
它**不**改变以下事实：`STAGE` 仍是独立 opt-in（默认 `STABASIS=REFERENCE`、
`STABILITY=0` 的输出逐位不变）；算例是 3 杆 `TRUSS` 浅拱研究构型，不是
真实桥梁施工阶段；位移控制 / Riks、分支屈曲与分支切换、初始缺陷、
材料非线性、动力 snap-through 均未闭环。因此 **B3 仍 `P` 不 `Y`**，
依据见 `BridgeStudio/B类后续项接入规格.md` §2.5 / §2.6；本项不代替
该结论，也不得靠改文档提前打勾。
