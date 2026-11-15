# 第 25 项 CTest：`bridge_limit_point_3d`

**第十九次（2026-11-15，a 阶段 ③ / `③a`）** —— 真实斜拉桥**越过极值点
并跟踪下降段**的 FEM ↔ ANSYS 对标。BS 侧正式算例为
`BridgeStudio/cases/CB3D22_limit_point_vs_ansys/`，本项是它的**离线
回归守卫**：冻结证据复制到 `cases/`（15 份），跑 CTest 时既不重跑
ANSYS，也不重新建模，只用当前 `cable_solver.exe` 重放两个 `.dat`。

> **编号口径**：`CB3D22` 专属于本项与 BS 侧极限点算例。第 24 项
> `completion_basis_3d` 覆盖的是既有 `CB3D10` / `CB3D13` 冻结输出，
> 其黄金值已改为 `case=CB3D10_BASIS_FIX` /
> `completion_basis_3d.json`，不再复用 `CB3D22`。

## 为什么需要这一项

② 阶段打通了「稳定路径从成桥态起步」的两段 ANSYS 求解，但路径始终停在
**上升段**，`λ` 单调增加；球面弧长法真正的能力是越过 `dλ/ds = 0` 的
极限点继续跟踪下降段（`λ` 减小而位移继续增大）。此前只有
`snap_through_ansys_3d` / `frame_postbuckling_ansys_3d` 这类抽象构件
探针，**没有真实桥梁算例**证明主链路 `cable_solver` 能做到。

## 构型（如实边界）

在 ② 的 `CB3D21` 同源双塔三跨斜拉桥（50 节点 / 63 单元 / 16 索）上：

- **不改**几何、材料、恒载、索面积与下料长度；
- 只把桥塔实常数 `RealConst3D(id=2)` 的 `A / IY / IZ / J / AW` 按
  `TOWER_SCALE=0.012` 缩放，让塔在水平力下出现明显 P-Δ 软化；
- 增量荷载由竖向改为塔顶 `2009` 节点 `+X`，`PINC = 1.0e6 N`；
- 控制量 `CTRL = (2009, 1)`，`CTRLTARGET=1.5`，`PSI=0.5`，
  `PATHMETHOD=SPHERE`，`MAXPATHSTEP=2500`。

初版曾用 `tower_b=tower_h=3*sqrt(scale)` **重新建模**，塔自重缩小 83 倍
（塔顶恒载从 `-344882.8125` 变成 `-4138.59375`），成桥态不再同源，
FEM 路径不收敛、ANSYS 也不越极限点。现版本保留 `COMPLETION_CFG` 原模型，
只缩放实常数，恒载 / 索 / 几何逐位一致。

**这是为观测极值点构造的研究构型，不是实际设计尺寸。** 它不覆盖分支
屈曲 / 分支切换、初始缺陷、材料非线性，也不代表施工阶段中间态稳定。

## 求解口径

| | FEM | ANSYS 18.0 |
|---|---|---|
| 成桥态 | `MODE=STATIC` + 显式 16 根 `L0`，`STABASIS=COMPLETION` | `step1` 默认 `TIME=1`，`INISTATE` 赋 `EPEL,X` |
| 增量 | `R(u,λ)=F0 + λ·Finc − Fint(u)`，`Finc` 只来自 `*STABILITY_LOAD` | `step2` `F,2009,FX`，`λ = TIME − 1` |
| 跟踪 | Crisfield 球面弧长 `duᵀdu + ψ²dλ² = ds²` | `lp_arc` `ARCLEN,ON,0.005,1.0E-8`；`lp_load` 荷载控制到 `λ=1.5` |
| 终止 | `CTRLTARGET=1.5` / `MAXPATHSTEP=2500` | `ARCTRM,U,1.6,2009,UX` |

`λ` 在极值点后**非单调**，因此逐点对拍不能以 `λ` 为插值自变量，而是按
单调量 `ΔUX = UX − UX_step1` 插值 `λ`（`ARCLEN` 预测步偶有 `ΔUX`
微小回退，插值前先取严格递增子序列）。

## 冻结结果

```
FEM    NPOINT=553  lambda_max=1.62494051 @ dUX=0.200410507
       descending=207  lambda_end=1.24293771
       RES0/|F0| = 1.264e-9
ANSYS  lp_load  step2=31 pts, lambda_end=1.5（上升段，0 error / 1 warning）
       lp_arc   step2=2039 pts, lambda_max=1.62274808 @ dUX=0.201182045
                lambda_end=1.24868184（0 error / 3 warnings）
极值点相对差  d-lambda=1.349e-3, d-UX=3.850e-3
逐点 |Δλ|     上升段 125 pts worst=5.493e-3 mean=1.749e-3
               下降段 427 pts worst=2.054e-2 mean=3.273e-3
               lp_load 上升段 worst=9.438e-3 mean=1.526e-3
```

`lp_arc_ansys.out` 出现
`One indefinite tangent matrix is detected for TIME = 1.001`：这是极值点
附近切线刚度失去正定性的**物理标志**，不是求解失败；ANSYS job 是否成功
以 `NUMBER OF ERROR MESSAGES ENCOUNTERED=0` 判定。`lp_load` 在极值点前
终止，不应出现该消息。

## 独立锚点（`compare_bridge_limit_point_3d.py`，140 项检查）

1. **冻结证据完整性**：15 个文件逐个 `sha256` / 字节数核对，集合不得增删；
2. **黄金值自洽**：从冻结 `.dat`/`.txt`/`.inp`/`.out`/`.csv` 重算全部
   派生量，与 `reference/cb3d22_limit_point.json` 逐字段相等（含 553 个
   FEM 平衡点、2039 个 ANSYS 平衡点整条曲线）；
3. **同源锚点**：`.dat` 16 个 `L0` ↔ `.cfg.json` 的 `l0`；
   ANSYS `N,`/`E,` 弦长 ↔ `L0·(1+ε)`；`INISTATE eps ≡ 0.0008`。
   实测 `eps` 最大偏差 `1.83e-16`、弦长最大偏差 `2.43e-11 m`，
   66 节点 / 63 单元 / 16 索全部匹配。ANSYS 单元号按 `E,` 行创建顺序
   （索单元 48–63），不是该行第一个节点号；
4. **口径锚点**：`FINC_NORM ≡ |*STABILITY_LOAD|`、
   `CTRL0 ≡ [U](2009,1)`（`< 1e-10`，`[U]` 段按 10 位有效数字写出）、
   `[BUCKLING] ERROR` 段保留、`STABILITY=0` 的 `lp_default` 不泄漏
   `[STABILITY]` 且 `[U](2009,1)` 与主算例一致；
5. **极值点物理判据**：`RES0/|F0| < 1e-8`、`RES0` 与 `[NEWTON] RESLAST`
   同量级、极值点在路径中段、下降段 ≥ 100 步、`λ` 降幅 > 10%、
   `LAMBDA_FINAL` 与末点一致、每点残差 < 100 N；
6. **ANSYS 证据**：两个 job 各 2 次 `SOLVE`、1 条 `F,2009,FX`、
   16 条 `INISTATE`、0 error、`step1` 成桥态 `UX` 一致；
   `lp_arc` 越过极值点且出现不定切线刚度，`lp_load` 只走上升段；
7. **FEM ↔ ANSYS 逐点对拍**：按 `ΔUX` 插值 `λ`，上升 / 下降 /
   荷载控制三段分别限制最大偏差与均值（容差来自 BS 侧
   `run_cases_limit_point.py --check-only`，实测值见上节）；
8. **离线重放**：用当前二进制重跑 `cb3d22_limit_point.dat` 与
   `lp_default.dat`。本项 FEM 证据由 `ucrt64/bin/c++.exe` 产出并与 BS
   算例目录共享，而仓库历史证据用 `mingw64/bin/g++.exe`；两套工具链在
   `[STABILITY]` 逐点行的末位浮点不同，故走
   `tests/benchmark/replay_tolerance.py` 的**结构化比较**：非逐点行
   逐位一致（仅掩 `elapsed_ms`）、`STEP/ITER/CONVERGED` 逐位一致、
   `LAMBDA/CTRL/DS/PIVOT` 相对差 `< 1e-9`、`RES` 绝对差 `< 0.1 N`，
   且点数必须相同。冻结文件字节本身仍由第 1 项 sha256 与第 2 项
   黄金值自洽锁死。实测 mingw64 重放 `worst rel=9.546e-14`、
   `worst |dRES|=4.063e-3 N`。

## 生成器

```
python build_reference_limit_point.py <cases_dir> [reference_dir]
```

**不重跑 FEM / ANSYS**：`cases/` 是唯一冻结源，脚本只做哈希、派生指标
与同源锚点计算。新增 / 更新证据时先从 BS 算例目录同步 15 个文件，再重新
生成黄金值。正常回归不需要跑生成器。

早期生成器在 `apdl_anchor()` 用
`^N,(\d+),([^,]+),([^,]+),([^,]+)\s*$` 解析节点，`[^,]` 会跨行吞掉
`\n` 与 `! ---- Orientation nodes K ----` 注释，把说明文字当成 Z 坐标
（`ValueError: could not convert string to float`）。现改为逐行
`split(',')` + 严格 `float()` 的数值化解析，忽略注释行。

## 能力边界

本项证明「真实斜拉桥构型上，主链路球面弧长法能越过极限点并跟踪下降段，
且与 ANSYS `ARCLEN` 定量一致」。它**不**覆盖分支屈曲 / 分支切换、
初始缺陷、材料非线性，也不改变 `[STABILITY]` 起点仍是**成桥态**这一
事实（不是施工阶段中间态）。B3 是否推 `Y` 见
`BridgeStudio/B类后续项接入规格.md` §2.5 / §2.6，本项不代替该结论。
