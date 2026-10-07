# FEM-Eigen 修改记录

本文件记录最近一次集中开发中新增和修改的功能。

## 2026-09-20

### 新增：约束方程 `*CE`

- 采用罚函数法实现 APDL 风格约束方程。
- 输入示例：

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
- 校验文件：`verify_ce.dat`

### 新增：混凝土徐变与收缩 AEMM

- 新增 `*CREEP`、`*SHRINKAGE`。
- 第一阶段采用龄期调整有效模量法：

  ```text
  E_eff = E / (1 + χ·φ)
  ```

- 收缩应变 `ECS` 按等效初应变节点荷载施加，正值为缩短。
- 校验文件：`verify_creep_shrinkage.dat`

### 新增：规范竖向温度梯度 `*TEMP_PROFILE`

- 输入示例：

  ```text
  *TEMP_PROFILE
  1-10 T1=14.0 T2=5.5 A=0.10 B=0.80 H=0.80
  ```

- 按矩形截面近似积分：

  ```text
  N_T = ∫ E·α·T(z)·b dz
  M_T = ∫ E·α·T(z)·(H/2 − z)·b dz
  ```

- 自动转成梁单元轴向力和端弯矩等效节点荷载。
- 校验文件：`verify_temp_profile.dat`

### 新增：80 m 综合系杆拱桥算例

- 新增 `bridge_tied_arch_80m.dat`。
- 覆盖功能：

  - `BEAM` 主梁与拱肋
  - `CABLE` 吊杆
  - `STRUT` 压杆
  - `*PRESTRESS`
  - `*CE`
  - `*CREEP`
  - `*SHRINKAGE`
  - `*TEMP_PROFILE`
  - `*GRAVITY`
  - `*ANALYSIS NONLINEAR`
  - `*INFLUENCE`
  - `*LANE_LOAD`

### 新增：APDL 基本模型

- 新增 `bridge_tied_arch_80m_apdl.inp`。
- 仅包含 ANSYS APDL 可直接覆盖的基本模型、约束和基本荷载。
- 未展开影响线、车道包络、徐变收缩、规范温度梯度和自动预应力筋生成。

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
- 校验文件：`verify_lane_nonlinear.dat`

### 修复：影响线基准响应扣除

- `computeInfluenceCurve()` 现在先求解无单位移动荷载的基准响应，再计算单位荷载增量。
- 避免自重、预应力、温度等基准荷载污染影响线和车道荷载包络。

### 修复：吊杆节点连接

- 修正 `bridge_tied_arch_80m.dat` 吊杆连接：

  ```text
  吊杆 201：节点 2 → 102
  吊杆 202：节点 3 → 103
  ...
  吊杆 207：节点 8 → 108
  ```

### 修复：车道荷载量纲

- 统一采用 SI 单位制。
- 原车道荷载：

  ```text
  Q=-10.5 kN/m
  P=-360 kN
  ```

- 已改为：

  ```text
  Q=-10500 N/m
  P=-360000 N
  ```

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
- 新增倒拆分析：`example_stage --backward`。
- 倒拆分析改为线性正装结果逆序，正装-倒拆残差为零。
- 非线性倒拆初步版本：按非线性正装结果逆序输出。
- 非线性倒拆升级为逆序移除荷载、反向激活/钝化单元并继承非线性状态。
- 新增正装-倒拆迭代：`*STAGE_TARGET`、`*STAGE_CONTROL`，`example_stage --iterate`。
- 正装-倒拆迭代扩展为多变量/多目标，采用有限差分 Newton 法。
- 正装-倒拆迭代升级为 Levenberg-Marquardt 最小二乘，支持变量上下界。
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
- 实现单元激活/钝化、分阶段荷载和状态输出。
- 阶段节点荷载和单元均布荷载按正装顺序累积。
- 施工阶段状态继承：弹性模量、预应力等效荷载、累积荷载、支座变化在阶段间继承。
- 非线性施工阶段继承上一阶段位移作为初值。
- 非线性施工阶段继承上一阶段 cable/strut 激活状态。
- 新增示例程序 `example_stage`。
- 校验文件：`verify_stage.dat`、`verify_stage_advanced.dat`

## 当前可用编译目录

由于原 `build` 目录部分目标文件被系统占用，当前使用：

```text
build_cable
```

锁释放后仍可执行：

```powershell
cmake --build "D:\AI-work\FEM-Eigen\build"
```

## 文档

- `input_format.md`
- `README.md`
- `CHANGELOG.md`
