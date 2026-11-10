# FEM-Eigen 使用说明书

## 1. 概述

FEM-Eigen 是一个轻量、无第三方运行依赖的 C++ 有限元程序，面向桥梁结构静力分析，包含两套内核：

- **二维内核**（`fem.hpp` + `fem_2d.cpp`）：平面杆系结构（梁、桁架、索、压杆、刚臂），覆盖施工阶段、影响线、车道荷载、预应力、徐变收缩、温度、几何非线性。
- **三维内核**（`fem_3d.hpp`）：空间梁格结构（桁架、梁、薄壁梁），覆盖施工阶段包络、影响线/影响面、多车道荷载、模态/屈曲。

### 单位约定

| 物理量 | 单位 |
|---|---|
| 长度 | m |
| 力 | N |
| 应力 / 弹性模量 | Pa |
| 密度 | kg/m³ |
| 线质量（附加质量） | kg/m |
| 面积 | m² |
| 惯性矩 | m⁴ |
| 翘曲常数（扇形惯性矩）`IW` | m⁶ |
| 转角 / 翘曲率 | rad / rad/m |

输入文件、输出文件、命令行均使用上述单位，**程序不做单位换算**。

## 2. 快速开始

### 2.1 编译

环境要求：MSYS2 的 **mingw64** 工具链（编译器与 SuiteSparse/CHOLMOD 必须同属 mingw64，混用 ucrt64 会导致运行库不一致）。CMake 4.x 为独立安装。

```powershell
$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
cmake -S . -B build -G "MinGW Makefiles" `
  -DCMAKE_C_COMPILER=C:/msys64/mingw64/bin/gcc.exe `
  -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe
cmake --build build -j8
```

### 2.2 运行

二维模型（`*NODE` 两坐标）：

```powershell
.\build\example_mixed.exe model.dat
```

三维模型（`*NODE` 三坐标）：

```powershell
.\build\example_3d_beam.exe model.dat      # 单个静力工况
.\build\example_3d_grillage.exe bridge.dat # 梁格 + 施工阶段包络
```

### 2.3 回归测试

```powershell
.\run_tests.ps1      # 全部测试（解析解 + 逐字节回归 + ANSYS 对标）
ctest --test-dir build  # 等价入口
```

`dat_byte_regression` 把 37 个 `verify_*.dat` 算例的输出与
`tests/regression/baselines/` 下的黄金基线逐字节比对。**任何改动
`example_*.cpp` 输出行、`defaultEnvelopeSpecs` 规格、
`envelope3DFromInfluence` 的乘子或积分权重的提交，都必须同步刷新基线**：

```powershell
python tests/regression/compare_regression.py --freeze --bin build --src . --baseline tests/regression/baselines
```

刷新前先跑一次不带 `--freeze` 的比对，确认差异全部能解释；刷新时**只提交
真正有内容变化的文件**——`--freeze` 会把 37 个文件全量 `copy2`，CRLF/LF
行尾差异会把无变化文件也标成 `M`。每条基线变动都要在 `CHANGELOG.md` 写明
数值依据。2026-11-09 曾因 `6a2d439`（二维包络新增位移行）与 `09c56bc`
（补线荷载梯形积分权重）两处有意变更后未同步基线，导致该项长期 3/37
失败；已刷新并恢复 6/6 全绿，详见 CHANGELOG 当日条目。

## 3. 两套输入格式

程序根据节点坐标维度自动区分二维/三维模型，二者关键字风格一致但有差异：

| 项目 | 二维 | 三维 |
|---|---|---|
| 节点 | `*NODE 1 x y` | `*NODE 1 x y z` |
| 单元类型 | `*ET`：`BEAM`/`TRUSS`/`CABLE`/`STRUT`/`RIGID`/`RIGID_PIN` | `*ET`：`TRUSS`/`BEAM`/`BEAM7` |
| 材料 | `*MP`（`E`/`ALPX`/`DENS`） | `*MP`（`E`/`DENS`） |
| 实常数 | `*R`（`A`/`I`/`WT`/`WB`/`H` 等） | `*R`（`A`/`IY`/`IZ`/`J`/`G`/`WT`/`WB`/`WY`/`WTORS`/`IW`/`ADDMASS` 等） |
| 自由度 | 每节点 3 个：u、v、θ | 每节点 6 个：u、v、w、θx、θy、θz（`BEAM7` 为 7 个，多一个翘曲率 φ） |

预处理宏（二维）：`*FILL a b`（等分生成节点）、`*DO … *ENDDO`（循环）。

完整语法见 `input_format.md`（二维）与 `input_format_3d.md`（三维），本节只讲“怎么选”与新增功能。

## 4. 单元类型怎么选

| 需求 | 选择 | 说明 |
|---|---|---|
| 常规梁 / 刚架 | `BEAM` | 默认，二维 3 自由度、三维 6 自由度 |
| 薄壁开口截面（工字梁、槽钢）受约束扭转 | `BEAM7`（三维） | 7 自由度，需在 `*R` 给 `IW` |
| 只受拉构件（吊杆、斜拉索） | `CABLE` | 受压自动退出 |
| 只受压构件（压杆） | `STRUT` | 受拉自动退出 |
| 可拉可压杆件 | `TRUSS` | 双向 |
| 刚体连接 | `RIGID` / `RIGID_PIN` | `RIGID_PIN` 不传弯矩 |

`CABLE`/`STRUT` 的初应变 `ISTRN`：`CABLE` 配 `ISTRN>0` 表示预张紧；`STRUT` 配 `ISTRN<0` 表示预压缩。

## 5. 分析类型怎么选

`*ANALYSIS TYPE=LINEAR`（默认）或 `*ANALYSIS TYPE=NONLINEAR TOL=1e-6 MAXITER=50`。

- `LINEAR`（方案 A）：线弹性 + `CABLE`/`STRUT` 单向开关迭代。绝大多数常规结构用它。
- `NONLINEAR`（方案 C）：Newton-Raphson 几何非线性（覆盖 `TRUSS`/`CABLE`/`STRUT`），适合大位移、索的几何刚度效应、非线性施工阶段；`BEAM` 用线性欧拉梁切线 + 一阶 P-Δ 几何刚度，**端部释放（`*RELEASE`）同样生效**（局部切线凝聚后装配，见 `input_format.md` §7.14）。

判断依据：位移很小、索始终张紧时两者结果几乎一致，用 `LINEAR` 更快；存在大位移或索可能松弛时用 `NONLINEAR`。

## 6. 求解器怎么选（二维 `example_mixed`）

命令行：`--solver 名称`。可用值：`LDLT`（默认）、`LLT`、`LU`、`QR`、`CG`、`BICGSTAB`、`AMGCL`、`CHOLMOD`。

| 情形 | 建议 |
|---|---|
| 常规对称正定（默认） | `LDLT` |
| 需要纯 Cholesky | `LLT` |
| 切线矩阵非正定 / 病态（几何非线性） | `LU` 或 `BICGSTAB` |
| 最小二乘 / 欠定 | `QR` |
| 超大规模稀疏（高性能） | `CHOLMOD`（需 SuiteSparse）或 `AMGCL`（需 AMGCL 头文件） |

求解器对比：`example_mixed.exe --compare model.dat`。

三维内核内部使用 `SimplicialLDLT`（失败自动回退 `SparseLU`），暂不通过命令行切换。

## 7. 预应力方法怎么选

`*PRESTRESS` 的 `METHOD`：`AUTO`（默认）/`EQUIV`/`TENDON`。

- `AUTO`：`LINEAR` 分析自动用 `EQUIV`，`NONLINEAR` 分析自动用 `TENDON`，推荐。
- `EQUIV`：把有效预应力换算成等效荷载（横向均布 + 端部弯矩 + 轴力），只对 `BEAM`。梁组只有一端 `dof1` 自由时在该端施加 `Pe`；两端都自由时施加自平衡力对 `+Pe/-Pe`（连续梁中间跨逐跨张拉走这条）；两端都固定时报错，改用 `TENDON`。不做筋-混凝土协同，也不含分项损失。
- `TENDON`：自动生成 `CABLE` 筋单元 + `RIGID_PIN` 刚臂，保留真实筋几何，适合非线性。

损失用 `LOSS` 系数一次性折减；若需按规范逐项计算（锚固、摩擦、弹性压缩、松弛、收缩徐变），用 `prestress_loss.hpp` 的 `computePrestressLosses`（见 §9.7）。

## 8. 约束与荷载

### 8.1 支座约束 `*BOUNDARY`

刚性约束：`<节点|区间> <dof|区间> <值>`。值非零表示支座沉降/强制位移。

**弹性支承（新增）**：第四列 `K=<刚度>`，该自由度变为接地弹簧。

```text
*BOUNDARY
1 1-3 0.0
10 2 0.0 K=2.0e6   ! 节点10 的 v 由竖向弹簧支承
```

三维同样支持：`<节点> <dof 或 a-b> [值] [K=刚度]`。

### 8.2 常见荷载

| 关键字 | 含义 |
|---|---|
| `*CLOAD` | 节点集中力 / 力矩 |
| `*ELEMENT_LOAD` | 单元分布荷载（二维 `GY` 等方向；梁上梯形/均布） |
| `*TEMPERATURE` / `*ELEMENT_TEMPERATURE` | 均匀温度 / 单元温度 |
| `*TEMP_PROFILE` | 规范竖向温度梯度 |
| `*GRAVITY` | 重力加速度（配合 `DENS` 自动自重） |
| `*SELF_WEIGHT` | 显式自重（旧格式兼容） |
| `*CE` | 约束方程（主从自由度关系） |
| `*CREEP` / `*SHRINKAGE` | 徐变 / 收缩（AEMM） |
| `*RELEASE` | 梁端释放（铰 `M`、剪力 `V`、轴力 `N`） |

### 8.3 影响线与车道荷载

二维：`*INFLUENCE` 指定响应，`*LANE_LOAD` 给车道 q/p。

三维：

```text
*LANE_3D
NODES=1,2,3 Q=-10.0 P=-100.0 DOF=2 IMPACT=1.15 REDUCTION=0.9

*LANE_SURFACE_3D
X_NODES=1,2,3 Y_NODES=1,2 Q=-10.0 P=-100.0 DOF=2 LANE_COUNT=2 IMPACT=1.15
```

**自动车道荷载（新增）**：`CODE=JTG_D60` 按车道跨径自动生成 `Q=10.5 kN/m`、`P=Pk(L)`（L≤5m→270kN，L≥50m→360kN，中间线性内插），影响面再按 `LANE_COUNT` 自动横向折减。

```text
*LANE_3D
NODES=1,2,3 DOF=2 CODE=JTG_D60
```

冲击系数可先用 `example_modal` / `example_3d_modal` 求基频，再由 JTG D60 公式（`impactCoefficientFromFrequency`）换算后填 `IMPACT`。

## 9. 新增功能专项说明

### 9.1 弹性支承

`*BOUNDARY` 第四列 `K=`（见 §8.1）。`K>0` 时刚度加到对角线上、基准位移产生等效荷载 `K·value`；不再作为刚性约束。适用桥墩弹性支座、桩基弹簧、地基反力等。

### 9.2 七自由度薄壁梁 `BEAM7`（约束扭转/翘曲）

三维专用。每节点 7 个自由度（多一个翘曲率 `φ=dθx/dx`），用于薄壁**开口**截面（工字梁、槽钢）的约束扭转。

```text
*ET
1 BEAM
2 BEAM7
*R
1 A=0.01 IY=1e-8 IZ=1e-8 J=1e-8 G=80e9
2 A=0.01 IY=1e-8 IZ=1e-8 J=1e-8 G=80e9 IW=1e-10
```

`IW` 为扇形惯性矩（m⁶）。`IW=0` 时退化为圣维南扭转，与 `BEAM` 等价。`example_3d_beam_live` 会为 `BEAM7` 额外输出双力矩 `B` 包络。

### 9.3 变截面

- 二维：`section.hpp` 的 `VariableSection`（控制点 + 线性插值），配合纤维截面/组合梁使用，采用“单元中点截面”近似。
- 三维：`beamElementFromSection` 从 `VariableSection` 取中点截面生成 `BeamElement`；等价做法是给每根梁单元直接赋插值后的 `A/IY/IZ/J`，网格加密即收敛。

### 9.4 纤维截面与材料非线性

`section.hpp`：把截面离散为纤维（y、z、面积、材料），统一计算面积/惯性矩/形心/上下缘，支持组合截面（各纤维不同 `E`）、温度按纤维积分、材料非线性（`STEEL` 理想弹塑性屈服、`CONCRETE` 受拉开裂）。

### 9.5 自振频率 / 屈曲 / 冲击系数

- 二维：`example_modal.exe model.dat`（`computeNaturalFrequencies2D`、`computeBucklingFactors2D`、`impactCoefficientFromFrequency`）。
- 三维：`example_3d_modal.exe model.dat`（`computeNaturalFrequencies3D`、`computeBucklingFactors3D`）。质量矩阵需在 `*MP` 给 `DENS`（或 `*R` 给 `ADDMASS`）。
- 自振频率默认走稀疏 shift-invert 特征求解、屈曲走稀疏 Lanczos（`A=L⁻¹KgL⁻ᵀ` 的最大 |μ|，μ<0 对应屈曲模态，λ=-1/μ），均基于 `external/spectra/` 头文件库，大模型从稠密的 O(n³) 降为近似 O(n·k)，几千自由度可秒级完成；若稀疏求解未收敛或分解失败会自动回退稠密求解，结果一致。

### 9.6 规范荷载库

`code_loads.hpp`（JTG D60-2015）：`jtgD60LaneLoadQk()`（10.5 kN/m）、`jtgD60LaneLoadPk(L)`（270~360 kN）、`jtgD60PedestrianLoad()`（3.0 kN/m²）、`jtgD60LaneReduction(n)`（1/1/0.78/0.67/0.60/0.55/0.52/0.50）。输入侧通过 `CODE=JTG_D60` 直接调用（§8.3）。

### 9.7 预应力损失分项

`prestress_loss.hpp` 的 `computePrestressLosses`（JTG 3362-2018 第 6.2 条）分项输出：σl1 摩擦（6.2.2）、σl2 锚具变形/回缩/接缝压缩（6.2.3）、σl3 温差（6.2.4）、σl4 弹性压缩（6.2.5）、σl5 松弛（6.2.6）、σl6 收缩徐变（6.2.7），以及有效预应力 σpe。

- σl5 分两类：钢丝/钢绞线按 6.2.6-1（Ψ·ζ·(0.52σpe/fpk−0.26)σpe）；预应力螺纹钢筋按 6.2.6-2/3（一次张拉 0.05σcon、超张拉 0.035σcon），用 `*PRESTRESS_LOSS` 的 `REBAR=1` 选择。
- 后张法 σpe=σcon−σl1−σl2−σl4，先张法 σpe=σcon−σl2，用 `POST=0/1` 选择。
- 已对照正式条文核对 σl1~σl6 系数。仍为简化式的是：σl2 曲线筋按附录 G 的反摩擦影响未实现（仅直线筋 ΣΔl/l·Ep）；σl6 的 εcs、φ 需由用户按附录 C 计算后输入，程序不自动算收缩/徐变随龄期发展。
- **2D / 3D 同源**：`parsePrestressLossLine()` 是两边共用的唯一键名映射，`computePrestressLosses()` 是唯一公式。2D 用 `*PRESTRESS ... LOSSREF=<编号>` 引用；3D 用 `*STAGE ... PRESTRESS=<elems>:P=...:E=...:LOSSREF=<编号>` 引用（只在关键字模式可用，与 `LOSS=` 互斥）。3D 侧 `resolvePrestressLosses3D()` 在整份文件读完后把 `σpe/σcon` 写回 `loss`，所以表写在 `*STAGE` 前后都行；没解析过就直接求解会报错，不静默用 1.0。格式与边界见 `input_format_3d.md`「阶段预应力分项损失」。
- 只有引用了表的束才写 `stage prestress loss` 结果段（σl1~σl6、σpe、RATIO，12 位有效数字）；只给一个总损失比的束不产生这段。σl3/σl5/σl6 属后期损失，不进 σpe，也不进等效荷载的 `Pe = P·(σpe/σcon)`。

### 9.8 索长迭代与制造线形几何闭合（`cable_solver`）

成桥态索承体系有两条**互不等价**的外层迭代，由 `.dat` 的 `MODE=` 选择，
都基于共转列式几何非线性内核（`fem_3d_nl.hpp`）：

| MODE | 未知量 | 控制误差 / 判据 | 方法 | 产出 |
| --- | --- | --- | --- | --- |
| `TENSION` | 各索无应力长度 L0 | max 相对索力偏差（无量纲，`tolTension`） | 影响矩阵 `J = ∂T_i/∂L0_j = M·D` 的 Newton，一次分解 + nc 次回代，冻结复用、残差停滞才重建 | 命中目标索力的下料 L0；**不保证线形闭合** |
| `COMPLETION` | 节点坐标偏移 ΔX（预抬高） | `max|ΔX + u|`（m，`tolClosure`） | 低秩 Newton / Woodbury + 回溯线搜索 + 不动点兜底 | 制造线形（预抬高）与成桥态坐标 |
| `STATIC` | 无（给定 L0 一次求解） | 无外层迭代 | 一次几何非线性求解 | 成桥态；`[LENGTH]` 写 ERROR 段 |

要点：

- **索力命中 ≠ 线形闭合**。`TENSION` 的唯一判据是索力相对偏差，成桥弦长
  相对设计弦长的偏差 `geomErr` 只是报告量。实测斜拉桥算例索力收敛到
  ~3e-8，而 `geomErr` 稳定在 0.15 m 左右；需要制造线形必须另跑
  `MODE=COMPLETION`。
- 两种迭代的控制误差**单位不同**（无量纲 / m），不能画在同一根轴上。
- 悬索桥主缆找形是「索多边形找形 → 平衡态弦长闭式反算 L0 → 一次求解」，
  同样没有外层迭代。
- 逐轮过程（误差、弦长偏差、步长、线搜索次数、影响矩阵重建/限幅/兜底
  事件、每索终值、完整 J）写入结果文件 `[LENGTH]` 段，格式见
  `cable_length_output_3d.hpp`；过程量只说明数值如何收敛，**不是设计合格
  判定**。

**这是外层迭代，不是内层。** 每一次 `[LENGTH]` 轮次里都要重跑一次完整的
`solveModel3DNonlinear`（共转列式 Newton + 荷载步）求平衡。内层过程单独写
`[NEWTON]` 段（见 §10）：成桥态只暴露收敛那一轮，非线性施工阶段每阶段
各暴露一段。两层迭代的控制量、未知量与收敛判据都不同，不要混读。

**`P3-3b`：完整形状敏感度（2026-11-09 已实现，默认仍走低秩）**

低秩敏感度 `G = -W·Kc·Bᵀ` 只计入「移动节点 → 索弦长变 → 索力变 →
位移变」，漏掉「梁参考轴方向、无应力长度、局部刚度随参考构形 X 变」以及
「缆索自重外荷载随几何变」。它影响 **Newton 收敛率**，不影响结果正确性：
线搜索按实际闭合误差做 Armijo 判据，不动点步 `ΔX ← -u` 是收敛保证。

实现（都在 `fem_3d_nl.hpp`，不新建力学求解）：

```
平衡      R(X,u) = F_ext(X) - F_int(X,u) = 0
全微分    ∂R/∂X·dX - K_t·du = 0      （K_t = ∂F_int/∂u）
⇒         du/dX = K_t⁻¹·∂R/∂X = -K_t⁻¹·K_shape
其中      K_shape = ∂(F_int - F_ext)/∂X |_{u fixed} = -∂R/∂X
Newton    (I + G_full)·δ = -(ΔX + u)，G_full = du/dX
```

- `crShapeResidual3D()` 固定 `u` 装配 `R`，与 `solveModel3DNonlinear` 的
  `Fext` / `assemble` 逐项同源（节点荷载、弹簧 `k·value`、缆索自重
  `w·Lc/2`、弹簧内力 `k·u`；初应变已在单元内力里）。
- `computeShapeSensitivity3D()` 对每个 `isT` 平动自由度做中心差分
  （默认 h=1e-5 m），只重建受扰动节点相连的 `CRElemData`，是
  `2·nt` 次**装配**而不是 `2·nt` 次非线性求解；再用收敛构形的
  `NLTangent3D` 一次分解 + `nt` 次回代给出 `G_full`。
- **行与列要分开处理**（本轮最重要的口径）：未知量只有节点坐标，所以
  **列**只取 `isT`；但**行必须保留全部自由自由度（含转角、弹性支承）**，
  改变梁参考轴方向会改变端弯矩。只有硬约束行清零（与 `applyDirichlet`
  一致）。初版误把「非 isT 行」一并清零，小模型差分对拍能过，实桥上
  Newton 却退化成线性收敛、比低秩还慢；修正后才恢复二次收敛。
- `K_shape` 不再是 nc 秩 ⇒ 不能用 Woodbury，改为在 `isT` 子空间做
  `nt×nt` 稠密 LU；奇异时显式退到不动点步，由线搜索按实际误差裁决，
  不静默。

`ErectionGeometryOptions3D::sensitivity` / `CompletionStateOptions3D::sensitivity`
三选一，三条路径**只换 G**，线搜索与兜底完全相同：`LowRank`（默认）、
`Full`、`FixedPoint`（G=0）。`cable_solver` 用
`*CABLE_BRIDGE` 的 `SENSITIVITY=LOWRANK|FULL|FIXEDPOINT` 选择，缺省
`LOWRANK`；BridgeStudio 前端不写该关键字，生产结果不受影响。

**验收（`test_shape_sensitivity_3d.cpp`，44/44）**

1. 收敛态自由自由度 `|R|/2P = 3.1e-11`（残差口径同源）；非索锚点列
   `max|K_shape| = 2.14e6 N/m`、自由转角行 `1.81e6 N·m/m`，证明梁参考
   形状项确实存在；`||G_full - G_low||/||G_full||` 实测 **0.17%**
   （源码旧注释「约 20%」是错的）；`σ(I+G_full)∈[0.058,1.17]`，cond≈20。
2. `G_full·v` 与「扰动 X + 完整重求平衡」的 `du/dX·v` 对拍：
   h=1e-2/1e-3/1e-4 三档 L2 相对差均 `3.6e-5`、相关系数 `1.00000000`；
   单自由度方向 `1.7e-8`。
3. 完整 / 低秩 / 纯不动点三条路径收敛到同一解：制造线形差
   `1.6e-13 m`（Full vs LowRank）、`1.4e-11 m`（FixedPoint），索力差
   `5.5e-13` / `7.3e-11`；三条闭合误差都 ≤ `tolClosure`。
4. 默认 LowRank 与「显式 LowRank」的 camber、轮数、闭合误差逐位相同；
   线搜索与不动点兜底未改。

**几何闭合的刚体规范自由度**（本轮发现，见该测试的 `testGaugeFreedom`）：
若整段梁的节点都没有竖向 / 面内转动约束，把参考构形作刚体平移（或绕
固结点转动）δ、位移取 `u = -δ`，当前构形、CR 梁内力、索力全不变，
于是 `(I+G)` 有零空间。此时**成桥构形与索力唯一，但制造线形 ΔX 是一族
解**，不同 Newton 路径会停在同一族的不同点上；实测两路径 camber 差
均值 `-0.6246 m`、去掉该均匀平移后偏差仅 `1.1e-4 m`、成桥构形差
`6.3e-11 m`。实桥的塔墩固结天然打断它；测试模型必须显式固结，否则
「三路径制造线形一致」这条 oracle 本身不成立。

**真实斜拉桥 `MODE=COMPLETION`（48 索、nt=258、TOL=1e-8）实测**：

| 路径 | 轮数 | 逐轮闭合误差 (m) | COMPLETION 耗时 |
| --- | --- | --- | --- |
| `LOWRANK`（默认） | 6 | `0.158 → 3.03e-3 → 8.73e-5 → 2.20e-6 → 5.55e-8 → 1.42e-9` | 317 ms |
| `FULL` | 4 | `0.158 → 1.06e-3 → 5.15e-7 → 9.53e-12` | 720 ms |
| `FIXEDPOINT` | 400 未收敛 | 停滞 `1.98e-3`，每轮 FALLBACK | 63.7 s |

`FULL` 的逐轮比值 `0.0067 → 4.9e-4 → 1.9e-5` 是**二次收敛**，闭合误差
比默认再小 150 倍；但每轮要 `2·nt` 次装配 + `nt×nt` 稠密 LU，
实桥上总耗时仍是默认的约 2.3 倍，且大自由度模型会按 `O(nt³)` 恶化。
`FIXEDPOINT` 停滞说明 Newton 步不可省，兜底只能作安全网。
**结论：四条验收 oracle 全部达标，但默认路径保持 `LOWRANK` 不变**——
它已稳定 6 轮收敛、零 FALLBACK，且代价随索数而非总自由度增长。
需要更高闭合精度、或某算例低秩路径反复兜底时，再显式用
`SENSITIVITY=FULL`。

### 9.9 带约束成桥目标索力优化（`OPTIMIZE=1`）

`METHOD=RIGID` 与 `METHOD=MINENERGY` 是「按一个固定准则直接算一组索力」；
本节是「给定可行域，求一组更优的索力」。默认 `OPTIMIZE=0` **完全不进入
这条路径**，既有算例的索力口径逐位不变。

力学全部复用内核已有积木，不重跑非线性求解、不在 Python 侧另算一套
索力分配：去掉待定索后的裸结构线性刚度 `K` → 单位索力对 `Ld` →
`Factorization3D::solveMulti` 给 `U = K⁻¹Ld` 与恒载位移 `u0 = K⁻¹Fdl`
→ 弯曲应变能柔度 `Fm = LdᵀU`、交叉项 `g = Ldᵀu0`，以及
`beam3DLocalEndForces` + `beamEndForcesAt` 给出单位索力下的主梁弯矩
影响矩阵 `A` 与恒载弯矩 `M0`。

- **目标**：`MINENERGY`（`min ½TᵀFmT + gᵀT`）或 `MINMAXM`（`min max|M主梁|`，
  IRLS 权重 `w ∝ |M|^(p−2)`，`p = 2,4,8,16` 逐轮压平峰值）。
- **约束**：索力箱约束 `T ∈ [TMIN, TMAX]`；逐塔 `ΣT·cos = 0` 作为
  **软约束罚项**并入二次型。不给 `TMAX` 时程序取
  `1.5 × RIGID 起点最大索力`（`TMAXSRC=DEFAULT_1.5_RIGID`，属程序口径，
  **不是规范限值**）——实测完全不限上界会把 48 索从 0.68~1.40 MN 推到
  0.14~6.67 MN，工程上不可用。
- **起点**：`*CABLE` 段把目标索力全给齐且都为正时从这组索力出发
  （`START=GIVEN`），否则用刚性支承连续梁法解析解（`START=RIGID`）。

**量纲必须归一**：弯矩项 ~(N·m)²、平衡项 ~N²、变量 ~1e6 N，直接拼进
同一个二次型会差十几个数量级。实现用 `x = T/Tscale`、
`Ã = (Tscale/Mscale)·A`、`m̃0 = M0/Mscale` 把目标与约束都缩到 O(1)，
箱约束、步长、收敛判据全在 `x` 空间，输出时还原成 N。早期版本把平衡罚
乘上弯矩项的尺度，罚项比弯矩项大几百倍，实测「二次型收敛了但主梁弯矩
反而变大」；修正后 `balanceScale=0.1` 才是真正的相对权重。

**两组弯矩峰值不是同一个量**：`MPEAK0`/`MPEAK` 是裸结构线性化影响矩阵
口径的优化目标，`MPEAK_FINAL` 是本次成桥态 `[GIRDEREND]` 实际 `max|My|`
（由 `cable_solver` 回填）。实测 `MINENERGY` 下二者方向一致、`MINMAXM`
开塔顶平衡时可能反向。**不得互相相减当作改善量**；未优化基线需另跑
`OPTIMIZE=0` 取 `[GIRDEREND] max|My|` 对拍。

**数学守卫**：若线性化弯矩峰值不降反升，写
`STATUS=objective_not_improved`、`APPLIED=0`，**不把 T\* 写回 `*CABLE`**，
保留原索力口径；`T` 行只是过程量。未收敛也如实写
`STATUS=not_converged`。塔顶平衡与切比雪夫弯矩目标存在实测权衡：
48 索算例 `MINMAXM + OPTBALANCE=1` 会被守卫拦下。

**体系守卫（`SYSGUARD=1`，默认开启）**：目标函数里只有弯矩、约束只有
索力箱，只盯目标会得到「数学收敛、体系不可用」的 T\*——48 索算例
`MINENERGY + OPTBALANCE=0` 把边跨索压到 `TMIN`，索承担桥面恒载从
99.31 % 掉到 89.04 %、同跨相邻索力比从 1.08 涨到 31.98，计算书工程
校核表出现两处不通过。内核在**写回 T\* 之前**用与计算书同一口径
（`100·Σ(T·ez)/DECKLOAD` 与按支座断面分跨的相邻索力比）预估三条判据：
`DECKSHARE > 90 %`、`DECKSHARE ≤ 100 %`（支座不上拔）、
`ADJRATIO ≤ 1.5`。任一不满足就写
`STATUS=system_check_rejected`、`APPLIED=0`、**不写回 T\***，
`REJECT` 列出 `system_deck_share_low` / `system_support_uplift` /
`system_adjacent_jump`；`SYSCHECK_PRED` 是 T\* 预估（守卫依据），
`SYSCHECK_FINAL` 是成桥态实测回填，`APPLIED=0` 时后者标
`ORIGINAL_TENSIONS`（原方案，不是过程量）。
`SYSGUARD=0` 仅供诊断，只记录不拦截（`SYSCHECK=OFF`、`PASS=NA`）；
`DECKLOAD` 缺失或支座断面不足两个时写 `NOT_AVAILABLE`，不假判。
限值可用 `SYSDECKMIN` / `SYSDECKMAX` / `SYSADJMAX` 覆盖。
这三条是体系合理性判据，**不是规范条文**，也不是设计合格判定。
实测 `MINMAXM + MINT=0.30` 三条全过（97.92 % / 748 kN / 1.39），
守卫正常放行，不会误杀可用方案。

优化轮数、收敛率、目标下降幅度都不构成规范意义上的设计合格判定。
段格式、错误码与实测数值见 `input_format_3d.md` 的 `[OPT]` 专节。

## 10. 输出与可视化

- 二维 `example_mixed`：写 `<模型名>.result`（节点位移、单元内力、单元应力、支座反力）。
- 二维 `example_viz`：额外生成 `*_model.svg`、`*_disp.svg`、`*_moment.svg`、`*_shear.svg`、`*_stress_top/bottom.svg` 等矢量图。
- 三维 `example_3d_grillage`：写 `*_summary.txt`（位移包络、梁内力包络、成桥态恒载/活载内力、施工阶段索力与预应力、冲击系数，以及 `modal` 段：前 8 阶自振频率/周期/参与质量/节点振型；质量全为零时写 `ERROR=no_mass`）。`LOSSREF` 路径另写 `stage prestress loss` 段（JTG 3362-2018 §6.2 的 σl1~σl6 / σpe / 有效比），只给 `LOSS=` 的束不写。成桥态索承体系 `cable_solver` 的对应段是 `[MODAL]`（12 阶，成桥态一致切线刚度）。段格式见 `modal_output_3d.hpp`。
- 线性屈曲：`example_3d_grillage` 另写 `buckling` 段（前 6 阶 `L <i> <λ> <Pcr_kN>`、`REF` 行的参考荷载与几何刚度诊断、`STATUS`）。参考荷载 = 成桥态构型 + 成桥态恒载合力（阶段增量累加），`Pcr = λ·Σ|Fi|`；`STATUS=no_compression` 表示参考荷载下无受压杆件、本就不存在压屈模态，不是求解失败。`cable_solver` 写 `[BUCKLING] ERROR=nonlinear_state_not_supported`：共转列式成桥态的切线刚度已含几何刚度，线弹性 λ 不适用，需另做非线性稳定全过程分析。段格式与口径见 `buckling_output_3d.hpp`。
- **λ 不是规范稳定安全系数**：JTG 3362 没有 λ 的合格限值；线性屈曲不含初始缺陷、材料非线性与几何非线性，是真实稳定承载力的上限。构件稳定按 §5.3 的稳定系数 φ 另行验算。
- 带约束成桥目标索力优化另写 `[OPT]` 段（`OPTIMIZE=1` 时）：
  `OBJECTIVE/START/NCABLES/NMOMENTS/NBALANCE/IRLS/NITER/CONVERGED/CONVBY/
  APPLIED/TMIN/TMAX/TMAXSRC` + `I`（逐轮）与 `T`（逐索）。`MPEAK_KIND=
  LINEAR_BARE` 与 `MPEAK_FINAL_KIND=GIRDEREND_FINAL` 分列两种弯矩口径；
  `APPLIED=0` 表示 T\* 未写回。默认 `OPTIMIZE=0` 写
  `ERROR=no_optimization_requested`。段格式见 `cable_opt_3d.hpp` 与
  `input_format_3d.md`。
- 成桥态索承体系 `cable_solver` 另写 `[LENGTH]` 段：`KIND=cable_length`
  （MODE=TENSION 的索长迭代，`ERRKIND=tension_rel`，含完整 `NJ` + `J` 影响
  矩阵）或 `KIND=erection_geometry`（MODE=COMPLETION 的几何闭合，
  `ERRKIND=closure_m`，不产出 J）；每轮 `I` 行、每索 `E` 行与
  `SENS` 诊断见 §9.8。`MODE=STATIC` / 悬索桥找形写
  `ERROR=no_iteration_static_mode` / `no_iteration_suspension_formfinding`：
  本来就无外层迭代，不是求解失败。段格式见 `cable_length_output_3d.hpp`。
- 内层共转列式 Newton：`cable_solver` 写 `[NEWTON]`（最终成桥态那一次
  平衡），`example_3d_stage_nl` 在 `<stem>.dat.nlstage.txt` 的每个
  `STAGE` 块里各写一段并带 `STAGE=<id>`。段内是
  `NSTEP/NITER/MAXITER/TOLF/TOLD/CONVERGED/REFACT/RESFIRST/RESLAST`
  + 逐迭代 `S` 行（荷载步、步内迭代、λ、|R|、|du|、|Ftarget|、
  重建分解、本步收敛）。收敛判据是残差**或**位移增量二选一，每个
  荷载步都要满足；换步时残差允许跳升。段格式见
  `newton_output_3d.hpp`。线性内核链路不产出此段，下游应写
  `no_newton_section` 说明，不摆空段。
- 三维 `example_3d_beam_live`：写 `*_live_envelope.csv`（轴力/弯矩/剪力/应力/双力矩包络）。

## 11. 常见问题

- **报“零刚度自由度”**：存在未约束、且未被任何单元连接的自由度；检查 `*BOUNDARY` 是否约束了全部刚体模式（二维梁纯桁架节点的 θ，三维纯轴向节点的转动分量会自动处理，但刚体平动/转动仍需手动约束）。
- **`CABLE` 全部退出导致奇异**：索的初应变或荷载方向有误，检查 `ISTRN` 与荷载符号。
- **非线性不收敛**：放宽 `TOL`、增大 `MAXITER`，或改用 `--solver LU`。
- **三维斜杆结果异常**：确认使用的是修复后的构建（早期版本三维梁坐标变换矩阵被转置，斜杆刚度算错）。
- **CHOLMOD 不可用**：确认 CHOLMOD 与编译器同属 mingw64，且 `PATH` 含 `C:\msys64\mingw64\bin`。
