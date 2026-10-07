# 更新记录

本文件记录 FEM-Eigen 项目的新增功能和主要修改内容。

## 2026-09-18

### 新增：CABLE（只受拉索）与 STRUT（只受压杆）

- 新增 `*ET` 单元类型 `CABLE`、`STRUT`，与 `TRUSS` 一样使用每节点 `u/v` 两个自由度。
- 新增实常数 `ISTRN`（初应变），用于输入预应力/索力/压杆预紧。
- 轴向机械轴力统一为：

  ```text
  N = E·A/L·ΔL + E·A·ISTRN
  ```

- `CABLE` 在 `N ≤ 0` 时退出；`STRUT` 在 `N ≥ 0` 时退出；退出单元的刚度、等效荷载、内力均置零。
- 求解采用**方案 A**：线性系统 + cable/strut 开关迭代，默认最多 100 次。
- 新增 `*ANALYSIS` 命令：`LINEAR` 为方案 A（默认）；`NONLINEAR` 为方案 C（几何非线性）。
- 校验数据：
  - `verify_cable_tension.dat`：预张紧索，轴力应为 `EA·ISTRN`。
  - `verify_cable_slack.dat`：初始松弛索，应退出工作。
  - `verify_strut.dat`：预压缩压杆，轴力应为 `-EA·ISTRN`。

### 新增：方案 C（Newton-Raphson 几何非线性）

- `*ANALYSIS TYPE=NONLINEAR` 现已可用，不再只是占位。
- 采用 Newton-Raphson 迭代，内层继续使用可切换的线性求解器（LDLT/CHOLMOD/AMGCL 等）。
- 几何非线性目前覆盖 `TRUSS`/`CABLE`/`STRUT`；`BEAM` 暂保持线性欧拉梁切线。
- 轴向单元使用当前构形的工程应变与一致切线：

  ```text
  N = EA·((L - L0)/L0 + ISTRN)
  K_t = (EA/L0 - N/L)·n·nᵀ + (N/L)·I
  ```

- `*ANALYSIS` 支持可选参数 `TOL`（缺省 `1e-6`）和 `MAXITER`（缺省 `50`）。
- 校验数据：`verify_cable_nonlinear.dat`（水平预张紧索跨中竖向荷载）。

### 新增：专用 `*PRESTRESS` 命令

- `*PRESTRESS` 可以自动完成混凝土梁预应力加载：

  ```text
  *PRESTRESS
  1 ELEM=1-10 P=1500e3 F=0.20 E1=0.0 E2=0.0 LOSS=0.85 EPS=195e9 APS=0.001 METHOD=AUTO
  ```

- `METHOD=EQUIV`：自动生成等效横向均布荷载、端部弯矩和轴力。
- `METHOD=TENDON`：自动生成 `CABLE` 预应力筋单元和 `RIGID_PIN` 内部销接刚臂，接入 `*ANALYSIS TYPE=NONLINEAR`。
- `METHOD=AUTO`：线性分析使用 `EQUIV`，非线性分析使用 `TENDON`。
- 端节点/内部节点由 `ELEM` 顺序自动识别，输入文件无需特别注明。
- 新增内部单元类型 `RIGID_PIN`，只传平动不向筋节点传弯矩，不参与应力结果输出。
- 校验数据：
  - `verify_prestress_equiv.dat`
  - `verify_prestress_tendon.dat`
  - `verify_prestress_ecc_tendon.dat`（端节点偏心弯矩）

### 新增：影响线与车道荷载包络

- 新增 `*INFLUENCE` 命令，可输出：
  - 节点位移影响线；
  - 支座反力影响线；
  - 单元内力影响线；
  - 单元应力影响线。
- 新增 `*LANE_LOAD` 命令，按方案 A/B 计算车道荷载内力、应力包络：

  ```text
  *LANE_LOAD
  1 NODES=1-21 Q=-10.5 P=-360 DIR=GY
  ```

- 新增示例程序 `example_influence`：

  ```powershell
  .\build\example_influence.exe .\verify_influence.dat
  ```

- 输出：
  - `<模型>_influence_lane<N>.txt`
  - `<模型>_envelope_lane<N>.txt`

- 校验数据：`verify_influence.dat`。

### 新增：约束方程 `*CE`

- 采用罚函数法实现 APDL 风格约束方程：

  ```text
  *CE_PARAMETER PENALTY=1e12
  *CE
  1 0.0 2 1 1.0 3 1 -1.0
  ```

- 形式：

  ```text
  Σ C_i·U_i = CONST
  ```

- 罚函数项：

  ```text
  K_pen = α·c·cᵀ
  F_pen = α·CONST·c
  ```

- 兼容 `SimplicialLDLT`、`CHOLMOD`、`AMGCL`。
- 校验数据：`verify_ce.dat`。

### 新增：混凝土徐变与收缩 AEMM

- 第一阶段实现 `*CREEP` 和 `*SHRINKAGE`：

  ```text
  *CREEP
  ALL PHI=2.0 CHI=0.8

  *SHRINKAGE
  ALL ECS=3.0e-4
  ```

- 徐变采用龄期调整有效模量法：

  ```text
  E_eff = E / (1 + χ·φ)
  ```

- 收缩应变 `ECS` 按等效初应变节点荷载施加，正值为缩短。
- 校验数据：`verify_creep_shrinkage.dat`。

### 新增：规范竖向温度梯度 `*TEMP_PROFILE`

- 输入格式：

  ```text
  *TEMP_PROFILE
  1-10 T1=14.0 T2=5.5 A=0.10 B=0.80 H=0.80
  ```

- 对矩形截面近似积分：

  ```text
  N_T = ∫ E·α·T(z)·b dz
  M_T = ∫ E·α·T(z)·(H/2 − z)·b dz
  ```

- 自动转成梁单元轴向力和端弯矩等效节点荷载。
- 校验数据：`verify_temp_profile.dat`。

### 新增：80 m 综合系杆拱桥算例

- 新增综合算例 `bridge_tied_arch_80m.dat`，覆盖：

  - `BEAM` 主梁与拱肋
  - `CABLE` 吊杆
  - `STRUT` 压杆
  - `*PRESTRESS` 主梁预应力
  - `*CE` 拱脚刚接
  - `*CREEP`、`*SHRINKAGE`
  - `*TEMP_PROFILE`
  - `*GRAVITY`
  - `*ANALYSIS NONLINEAR`
  - `*INFLUENCE`
  - `*LANE_LOAD`

- 可分别用 `example_mixed`、`example_influence`、`example_viz` 运行。
- 新增 APDL 基本模型：`bridge_tied_arch_80m_apdl.inp`。

### 新增：车道荷载非线性直接扫描

- `*LANE_LOAD` 增加 `METHOD` 选项：

  ```text
  *LANE_LOAD
  1 NODES=1-9 Q=-10500.0 P=-360000.0 DIR=GY METHOD=LINEAR
  ```

  或：

  ```text
  METHOD=NONLINEAR
  ```

- `LINEAR`：保留单位荷载影响线叠加方法。
- `NONLINEAR`：对每个集中荷载位置施加实际 `Q`、`P`，调用非线性求解器直接扫描包络。
- 校验数据：`verify_lane_nonlinear.dat`。

### 新增：施工阶段正装分析

- 新增 `*STAGE`、`*STAGE_LOAD`、`*STAGE_ELEMENT_LOAD`。
- 新增 `*STAGE_PRESTRESS`、`*STAGE_BOUNDARY`、`*STAGE_UNBOUNDARY`。
- 新增 `*STAGE_CREEP`、`*STAGE_SHRINKAGE`。
- 新增 `*STAGE_TEMP_PROFILE`。
- 新增 `*CREEP_CURVE`、`*SHRINKAGE_CURVE`，阶段 `TIME` 按龄期插值徐变/收缩。
- 收缩时变曲线改为增量收缩应变，逐步累计 `Δε_cs`。
- 徐变时变曲线增加轴向徐变逐步积分，根据上一阶段轴力施加增量徐变荷载。
- 梁单元增加弯曲徐变逐步积分，根据上一阶段弯矩施加增量徐变弯矩荷载。
- 新增 `*CREEP_MODEL`、`*SHRINKAGE_MODEL`，内置简化 CEB-FIP 1990 徐变收缩公式。
- 三维桁架基础：`fem_3d.hpp`、`example_3d_truss`、`verify_3d_truss.dat`。
- 三维施工阶段基础：`Model3D::Stage3D`、`solveModel3DStages()`。
- 三维 `*STAGE` / `*STAGE_LOAD` 输入解析和 `example_3d_stage`。
- 三维施工阶段包络和预拱度基础：`compute3DStageEnvelope`、`compute3DCamber`、`example_3d_stage_envelope`。
- `solveModel3D` 支持混合自由度：桁架 3 DOF/节点，梁 6 DOF/节点。
- 三维梁基础验证：`verify_3d_beam.dat`、`example_3d_beam`。
- 三维梁施工阶段和预拱度支持 6 DOF/节点。
- 三维梁施工阶段位移包络：`compute3DDisplacementEnvelope`、`example_3d_beam_envelope`。
- 三维桁架/梁轴向预应力初应变 `ISTRN`。
- 三维 AEMM 徐变/收缩基础：`*CREEP PHI/CHI`、`*SHRINKAGE ECS`。
- 三维影响线/车道荷载基础：`*INFLUENCE_3D`、`*LANE_3D`、`example_3d_influence`。
- 三维多车道基础：`*LANE_3D` 支持多条车道并合并包络。
- 三维影响面基础：`*LANE_SURFACE_3D`、`compute3DInfluenceSurface`、`example_3d_influence_surface`。
- 三维影响面包络基础：`envelope3DFromSurface`。
- 三维多车道折减系数 `REDUCTION`。
- 三维梁单元内力包络基础：`compute3DBeamEnvelope`。
- 三维施工阶段+活载组合基础：`compute3DStageLiveEnvelope`、`example_3d_stage_live`。
- 三维梁活载轴力影响基础：`*ELEMENT_INFLUENCE_3D`、`compute3DBeamAxialInfluence`、`example_3d_beam_live`。
- 三维横向车道自动折减基础：`LANE_COUNT`。
- 三维活载冲击系数基础：`IMPACT`。
- 三维影响面 CSV 输出。
- 三维影响面 SVG 可视化。
- 三维梁活载弯矩影响基础：`compute3DBeamMomentInfluence`。
- 三维梁活载剪力影响基础：`compute3DBeamShearInfluence`。
- 三维梁活载应力影响基础：`compute3DBeamStressInfluence`，支持 `WT/WB`。
- 三维梁活载包络报告 CSV 输出。
- 三维梁格综合算例：`bridge_3d_grillage.dat`、`example_3d_grillage`。
- 新增 `*STAGE_UNLOAD`，支持移除节点荷载和单元均布荷载。
- 新增倒拆分析：`example_stage --backward`，输出倒拆阶段结果和正装-倒拆状态残差。
- 倒拆分析改为线性正装结果逆序，正装-倒拆残差为零。
- 非线性倒拆初步版本：按非线性正装结果逆序输出。
- 非线性倒拆升级为逆序移除荷载、反向激活/钝化单元并继承非线性状态。
- 新增正装-倒拆迭代：`*STAGE_TARGET`、`*STAGE_CONTROL`，`example_stage --iterate` 调整阶段预应力。
- 正装-倒拆迭代扩展为多变量/多目标，采用有限差分 Newton 法。
- 正装-倒拆迭代升级为 Levenberg-Marquardt 最小二乘，并支持控制变量 `LOWER`、`UPPER`。
- 迭代目标增加支座反力 `REACTION` 类型。
- 迭代目标增加单元内力 `FORCE` 和单元应力 `STRESS` 类型。
- 迭代目标支持权重 `WEIGHT`，实现加权最小二乘。
- 新增 `*STAGE_CONSTRAINT`，支持位移、反力、内力、应力上下界不等式约束。
- 新增 `*CAMBER`，计算节点预拱度并输出 `_camber.txt`。
- 预拱度结果增加原始坐标和安装标高。
- 新增非线性预拱度迭代：`example_stage --camber-iterate`。
- `*CAMBER` 支持目标位移，预拱度 = 目标位移 − 最终位移。
- 新增施工反馈修正：`*FEEDBACK` + `example_stage --feedback`。
- 新增 `example_stage_viz`，输出每阶段模型和位移 SVG。
- `example_stage_viz` 补齐弯矩、剪力、轴力、上下翼缘应力图。
- 新增 `example_stage_envelope`，输出施工阶段内力/应力包络。
- `example_stage_envelope` 增加弯矩、上翼缘应力、下翼缘应力包络 SVG。
- `example_stage_envelope` 增加轴力、剪力包络 SVG。
- 实现单元激活/钝化和分阶段求解。
- 阶段节点荷载和单元均布荷载按正装顺序累积。
- 施工阶段状态继承：弹性模量、预应力等效荷载、累积荷载、支座变化在阶段间继承。
- 非线性施工阶段继承上一阶段位移作为初值。
- 非线性施工阶段继承上一阶段 cable/strut 激活状态。
- 新增示例程序 `example_stage`：

  ```powershell
  .\build\example_stage.exe .\verify_stage.dat
  ```

- 输出 `<模型>_stages.txt`。
- 校验数据：`verify_stage.dat`。
- 阶段预应力/支座变化校验：`verify_stage_advanced.dat`。

### 修改：输出与文档

- 轴向单元（`TRUSS`/`CABLE`/`STRUT`）在控制台和结果文件中输出真实单元类型，不再统一显示为 `TRUSS`。
- `input_format.md`、`README.md` 补充 `CABLE`/`STRUT`、`ISTRN`、`*ANALYSIS`、方案 C 的使用说明。

## 2026-09-17

### 新增：求解器运行时选择、计时与对比

`example_mixed` 现在可以在命令行选择求解器，并统计耗时：

```powershell
.\build\example_mixed.exe .\model.dat                      # 默认 SimplicialLDLT
.\build\example_mixed.exe .\model.dat --solver CHOLMOD     # 指定 CHOLMOD
.\build\example_mixed.exe .\model.dat --solver AMGCL       # 指定 AMGCL
.\build\example_mixed.exe .\model.dat --compare            # 对比 LDLT/CHOLMOD/AMGCL
```

- 支持的后端：`LDLT`、`LLT`、`LU`、`QR`、`CG`、`BICGSTAB`、`AMGCL`、`CHOLMOD`。
- `--compare` 会依次求解并输出各自耗时、相对 LDLT 的倍数和最大位移，用于判断“多少单元适用哪个求解器”。
- 计时为“组装 + 求解”总时间。

### 新增：单元应力输出

- 桁架：轴向应力 `σ = N/A`，输出两端 `σ1`、`σ2`。
- 梁：上下翼缘正应力 `σ_top = N/A − M/WT`、`σ_bottom = N/A + M/WB`，输出两端共 4 个值。
- 梁需要在 `*ELEMENT` 上输入上下翼缘抗弯模量 `WT`、`WB`。
- 应力符号：受拉为正；上翼缘 = 局部 +y 侧，下翼缘 = 局部 −y 侧。
- 结果文件 `*.result` 中新增“单元应力”段；`example_mixed` 控制台同步输出。

### 新增：可视化程序输出 7 张独立 SVG

`example_viz` 现在为每个模型生成 7 张独立图，均带节点编号、自动居中、四边留边距：

```text
<名称>_model.svg          模型图
<名称>_disp.svg           位移图
<名称>_moment.svg         弯矩图
<名称>_shear.svg          剪力图
<名称>_axial.svg          轴力图
<名称>_stress_top.svg     上翼缘应力图
<名称>_stress_bottom.svg  下翼缘应力图
```

### 修复：单元内力符号（重要）

`elementForces2D` 之前把“杆端力”直接当作“截面内力”，导致：

- 弯矩、剪力、上下翼缘应力在相邻单元共用节点处**不连续**；
- 弯矩图/剪力图看起来“反了”。

修正：把梁单元右端剪力、左端弯矩的符号换算成真正的截面内力：

```cpp
f.V1 = fphys(1);    // 左端剪力
f.V2 = -fphys(4);   // 右端剪力
f.M1 = -fphys(2);   // 左端弯矩
f.M2 = fphys(5);    // 右端弯矩
```

修正后，弯矩图、剪力图、上下翼缘应力图均连续，且悬臂梁端部应力符号恢复物理正确（向下荷载时上翼缘受拉）。

### 修改：CMake 链接与运行库

- MinGW 下静态链接 C++/GCC 运行库：

  ```cmake
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
  ```

  这样示例 exe 不再依赖 `libstdc++-6.dll`、`libgcc_s_seh-1.dll`。

- CHOLMOD 只链接到实际使用它的目标（`example_cholmod`、`example_amgcl`、`example_solvers`、`example_mixed`），避免 `example_viz` 等产生不必要的 DLL 依赖。

## 2026-09-16

### 新增：四种荷载类型

- `*ELEMENT_LOAD`：单元均布/梯形荷载，方向支持 `GX/GY/LX/LY` 及 `G<θ>`/`L<θ>`（角度支持小数）。
- `*GRAVITY` + `*SELF_WEIGHT`：自重，规则 `RHO`（面积×密度×g）或 `MASS`（线质量×g）。
- `*TEMPERATURE`：体系温度升降（全结构统一）。
- `*ELEMENT_TEMPERATURE`：指定单元温度，`UNIFORM`（整体升温）或 `GRADIENT`（梁温度梯度，需 `H`）。

### 新增：单元内力输出

- 桁架：两端轴力 `N1`、`N2`。
- 梁：两端轴力 `N1/N2`、剪力 `V1/V2`、弯矩 `M1/M2`。
- 结果文件 `*.result` 输出节点位移、单元内力、支座反力。

### 新增：二维欧拉梁单元与输入文件

- `EulerBeam2D`：二维欧拉-伯努利梁（每节点 3 自由度 u/v/θ）。
- `fem_io.hpp`：`*.dat` 关键字格式读取、`Model2D`、`solveModel2D`、`elementForces2D`、`writeResult`。
- 修复 `*ELEMENT` 关键字下多行数据只保留最后一个单元的 bug。

## 更早版本

- 可切换求解后端：`SimplicialLDLT`（默认）、`SimplicialLLT`、`SparseLU`、`SparseQR`、`ConjugateGradient`、`BiCGSTAB`、`AMGCL`、`CHOLMOD`。
- `CholmodSolver` 分阶段求解：`analyzePattern` / `factorize` / 单右端项 / 多右端项。
- 单元：一维/二维/三维桁架、二维平面应力三角形（CST）、三维四面体（T4）。
- 内置 Eigen 3.3.9 与 AMGCL 1.5.0 头文件。
