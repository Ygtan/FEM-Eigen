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
- `NONLINEAR`（方案 C）：Newton-Raphson 几何非线性（覆盖 `TRUSS`/`CABLE`/`STRUT`），适合大位移、索的几何刚度效应、非线性施工阶段。

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
- `EQUIV`：把有效预应力换算成等效荷载（横向均布 + 端部弯矩 + 轴力），只对 `BEAM`，要求梁组至少一端轴向可动。
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

## 10. 输出与可视化

- 二维 `example_mixed`：写 `<模型名>.result`（节点位移、单元内力、单元应力、支座反力）。
- 二维 `example_viz`：额外生成 `*_model.svg`、`*_disp.svg`、`*_moment.svg`、`*_shear.svg`、`*_stress_top/bottom.svg` 等矢量图。
- 三维 `example_3d_grillage`：写 `*_summary.txt`（位移包络、梁内力包络）。
- 三维 `example_3d_beam_live`：写 `*_live_envelope.csv`（轴力/弯矩/剪力/应力/双力矩包络）。

## 11. 常见问题

- **报“零刚度自由度”**：存在未约束、且未被任何单元连接的自由度；检查 `*BOUNDARY` 是否约束了全部刚体模式（二维梁纯桁架节点的 θ，三维纯轴向节点的转动分量会自动处理，但刚体平动/转动仍需手动约束）。
- **`CABLE` 全部退出导致奇异**：索的初应变或荷载方向有误，检查 `ISTRN` 与荷载符号。
- **非线性不收敛**：放宽 `TOL`、增大 `MAXITER`，或改用 `--solver LU`。
- **三维斜杆结果异常**：确认使用的是修复后的构建（早期版本三维梁坐标变换矩阵被转置，斜杆刚度算错）。
- **CHOLMOD 不可用**：确认 CHOLMOD 与编译器同属 mingw64，且 `PATH` 含 `C:\msys64\mingw64\bin`。
