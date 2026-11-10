# FEM-Eigen 开发说明书

本文档面向维护者/二次开发者，说明架构、数据模型、算法流程、测试体系与扩展方式。面向使用者的说明见 `USER_MANUAL.md`。

## 1. 项目结构与模块

```
.
├── fem.hpp            # 共用基础（2D/3D 共享）：求解器、DofMap、BeamElement、单元刚度
├── fem_3d.hpp         # 三维内核：Model3D、装配/求解、包络/影响线、模态/屈曲、BEAM7、读输入
├── fem_2d.hpp/.cpp    # 二维装配/求解/内力、输入解析、施工阶段、模态/屈曲（二维）
├── section.hpp        # 纤维截面：Fiber/Section、组合截面、变截面、材料非线性
├── code_loads.hpp     # JTG D60 规范荷载库
├── prestress_loss.hpp # 预应力损失分项（JTG 3362-2018）
├── example_*.cpp      # 可执行示例（每个对应一种分析能力）
├── test_2d.cpp        # 解析解断言测试（2D，129 项）
├── test_3d.cpp        # 解析解断言测试（3D + 统一内核，115 项）
├── tests/regression/  # 逐字节回归（37 个黄金基线 + 比对脚本）
├── tests/benchmark/ansys/   # 100m 简支梁 ANSYS 对标
├── tests/benchmark/bridge/  # 80m 系杆拱 ANSYS 对标
├── external/          # 自带 Eigen3、AMGCL、Spectra 头文件
└── run_tests.ps1      # 全量测试入口
```

依赖：Eigen3（自带）、AMGCL（自带，可选）、Spectra（自带，稀疏特征求解）、SuiteSparse CHOLMOD（系统可选）。`fem_2d.cpp` 编译进 `fem_2d` 静态库；`fem_3d.hpp` 与共享核心 `fem.hpp` 为纯头文件，3D 层以 INTERFACE 目标 `fem_3d` 表达。

## 2. 数据模型

### 2.1 二维（`fem_2d.hpp`）

- `Model2D`：`nodes`（map 节点号→二维坐标）、`elements`、`constraints`、`loads`、`elemLoads`、`prestresses`、`laneLoads`、`stages` 等。
- `Model2D::Element`：`type`（`BEAM`/`TRUSS`/`CABLE`/`STRUT`/`RIGID`/`RIGID_PIN`）、`E/A/I/alpha/h/Wt/Wb/dens/addmass/istrn/nu/K`、`conn`、`releases`。
- `Model2D::Constraint`：`node/dof/value` + `k`（弹性支承刚度，`>0` 即弹簧）。
- `ModelResult`：`U`（位移）、`R`（反力）、`nodeIds/coords`、非线性诊断。

### 2.2 三维（`fem_3d.hpp`）

- `Model3D`：`nodes`（三维坐标）、`elements`、`constraints`、`loads`、`stages`、`influences`、`elementInfluences`、`lanes`、`laneSurface`。
- `Element3D`：`type`（`TRUSS`/`BEAM`/`BEAM7`）、`E/A/Iy/Iz/J/G`、`Iw`（翘曲常数）、`dens/addmass`、`istrn`、`Wtop/Wbot/Wy/Wt`、`n1/n2`。
- `Constraint3D`：`node/dof/value` + `k`。
- `Lane3D` / `LaneSurface3D`：含 `autoCode`（`CODE=JTG_D60` 自动车道荷载）。

### 2.3 共用结构

- `DofMap`（`fem.hpp`）：统一自由度映射。`index`（节点号→序号）、`nodeIds`、`dofStart`、`dofCount`（每节点 3/6/7）、`nDofs`；`dof(node,local)`、`safeDof(...)`。二维 `buildDofMap2D` 给每节点 3 自由度；三维 `buildDofMap3D` 给纯轴向节点 3、`BEAM` 节点 6、`BEAM7` 节点 7。
- `BeamElement`（`fem.hpp`）：统一的 12 自由度三维梁单元（二维梁是其平面内特例）。局部自由度顺序 `[u1,v1,w1,θx1,θy1,θz1, u2,v2,w2,θx2,θy2,θz2]`。提供 `localStiffness()`（含 Timoshenko `PhiY/PhiZ`）、`transformation()`、`stiffness()`、`consistentMassMatrix(rhoA)`、`geometricStiffness(N)`、`localEndForces(ug)`。
- `SectionProperties`（`section.hpp`）：`A/cy/cz/Iy/Iz/Iyz/J/yTop/yBot/zLeft/zRight`，以及截面模量访问器。
- `BeamEndForces`（`fem_3d.hpp`）：统一梁端内力 `N/Vy/Vz/T/My/Mz/B`，由 `beamEndForcesAt(fl, dofsPerElem, end)` 从 12/14 自由度局部杆端力提取。

## 3. 装配与求解流程

### 3.1 二维线性（`solveModel2D` → 线性分支）

1. 建立节点序号映射，`nDofs = 3 × 节点数`。
2. `assembleLinear2D`：按单元类型组装 `K`（梁用 `EulerBeam2D`，含剪切参数 `Phi`；轴向单元用 `TrussElement`；端部释放用 `condenseBeamGlobal` 静态凝聚），组装荷载向量 `F`（节点力、分布荷载等效、CE 罚函数、弹性支承 `K(d,d)+=k`、`F(d)+=k·value`）。
3. `autoConstrainFreeRotations` 自动约束零刚度转角；`applyDirichlet` 施加硬约束（弹性支承不进 `cdofs`）。
4. `solve(K, F, SolverType)` 求解；`ReusableSolver` 支持结构不变时复用符号分解。

### 3.2 二维非线性（`solveModel2DNonlinear`）

Newton-Raphson：`assembleState` 每步组装 `Fint/Fext` 与一致切线（轴向单元用当前构形工程应变切线，梁叠加一阶 P-Δ 几何刚度），弹性支承在切线中加 `k`、在 `Fint/Fext` 中加 `k·u`/`k·value`；内层按残差收敛，外层按荷载分级。

### 3.3 三维线性（`assembleModel3D` → `solveModel3D`）

1. `buildDofMap3D` 建立自由度映射；`elementGdofs3D` 给出单元全局自由度（BEAM 12、BEAM7 14、轴向 6）。
2. 组装 `K`：`BEAM` 用 `beam3DLocalStiffness` + `beam3DTransformation`；`BEAM7` 用 `beam7LocalStiffness` + `beam7Transformation`；轴向单元用 `truss3DStiffness`。
3. 组装 `F`（节点力、初应变/收缩等效荷载、弹性支承）。
4. 硬约束进 `cdofs/cvals`，零刚度诊断，`applyDirichlet`；`solveSparse3D`（LDLT 失败回退 SparseLU）。

## 4. 单元实现与约定

### 4.1 坐标变换约定

三维变换 `λ = R^T`（全局→局部），其中 `R=[ex,ey,ez]`（列向量为局部轴在全局坐标的分量）。`Ke = T^T Kl T`，`fl = Kl·(T·ug)`。注意：早期版本误用 `R`（转置），导致非轴线对齐梁刚度算错，已在合并阶段修复；`BeamElement::transformation` 与 `beam3DTransformation`/`beam7Transformation` 保持一致。

### 4.2 梁刚度符号约定

- 绕 z 轴弯曲（v,θz）：与二维 `EulerBeam2D` 的 (v,θ) 一致，`θz=+dv/dx`。
- 绕 y 轴弯曲（w,θy）：耦合项符号与 z 轴相反，`θy=-dw/dx`。
- 扭转：圣维南 `GJ`（6 自由度）或 `GJ + E·Iw`（7 自由度，`φ=dθx/dx` 为翘曲率）。

### 4.3 七自由度梁 `BEAM7`（14 自由度局部刚度）

局部自由度 `[u1,v1,w1,θx1,θy1,θz1,φ1, u2,v2,w2,θx2,θy2,θz2,φ2]`。扭转块 `K_t = E·Iw·K_beam + G·J·K_geo`，其中 `K_beam` 为梁弯曲刚度形（12/6L/…），`K_geo` 为 `∫N'^T N'` 形（6/5、L/10、…）。`φ` 不参与旋转变换（标量）。

### 4.4 内力

- 二维：`elementForces2D`（`Kl·ul`，扣除分布荷载固端力，处理端部释放、索/压杆松弛、符号）。
- 三维：`beam3DLocalEndForces`（12 自由度）、`beam7LocalEndForces`（14 自由度，含双力矩 `B`）。

## 5. 2D/3D 内核合并的边界

合并分三阶段，已完成部分见 `kernel_merge_design.md`：

- A 统一单元接口：`BeamElement` 统一刚度/质量/几何刚度/变换/内力，2D 为其平面内特例。
- B 统一自由度表：`DofMap` 共用，`buildDofMap2D`/`buildDofMap3D` 分别产生编号。
- C 统一内力：`BeamEndForces` + `beamEndForcesAt` 统一 12/14 自由度内力提取。

保留边界：二维求解器 `solveModel2D` 额外承载端部释放、分布荷载固端力、索/压杆松弛、施工阶段、预应力、徐变收缩等三维不具备的特性，整体重写风险大于收益，故只统一到“单元/自由度/内力”层，求解器各自保留。

## 6. 测试体系

### 6.1 解析解断言（`test_2d.cpp` / `test_3d.cpp`）

用材料力学闭式解校验数值结果，当前 **244 项**（`test_2d` 129 项 + `test_3d` 115 项），
拆成两个可执行，避免一个维度崩溃吞掉另一个维度的结果；另有 `test_3d_nl` 79 项、
`test_shape_sensitivity_3d` 44 项专项断言（计数以 `CHANGELOG.md` 为准）。
公共断言辅助在 `test_common.hpp`。新增测试：写 `testXxx()` 并在对应 `main()` 中调用。

### 6.2 逐字节回归（`tests/regression/`）

`compare_regression.py` 对 37 个 `verify_*.dat` 用当前构建重算，与黄金基线逐字节/数值比对。刷新基线：`python compare_regression.py --freeze ...`（先确认当前结果正确再刷新）。

### 6.3 ANSYS 对标（`tests/benchmark/`）

- `ansys/`：100m 简支梁，参考值为 ANSYS BEAM4 的欧拉-伯努利闭式解。
- `bridge/`：80m 系杆拱，参考值来自用户提供的 ANSYS 结果（`bridge_tied_arch_80m_apdl_disp.txt`、`..._elemforce.txt`），`compare_bridge.py` 比对节点位移与单元内力。

三者均接入 CTest（`analytical_exact_solutions_2d`、`analytical_exact_solutions_3d`、`dat_byte_regression`、`ansys_benchmark`、`bridge_benchmark`）与 `run_tests.ps1`。

## 7. 构建与工具链

- 必须使用统一的 mingw64 工具链。编译器与 CHOLMOD/OpenBLAS 必须同属 `C:\msys64\mingw64`；混用 ucrt64 会导致 C 运行库不一致（历史上曾出现“PATH 陷阱”）。CHOLMOD 仅安装在 mingw64。
- CMake 目标：`fem_2d`（静态库）+ `fem_3d`（INTERFACE，纯头文件）+ 各 `example_*` + `test_2d` / `test_3d`。可选依赖：CHOLMOD（`FEM_USE_CHOLMOD`）、AMGCL。
- 生成器：`MinGW Makefiles`；CMake 4.x 独立安装，`mingw32-make` 可用 CodeBlocks 自带或 MSYS2 `make`。

## 8. 提交规范

遵循 `CONTRIBUTING.md` 的中文提交风格：修复类 `P0-3/P1-5：…`、功能类 `新增…：…`、`完善…：…`、`docs: …`。`.git` 有 ACL 时由用户手动 `git add` + `git commit`。

## 9. 如何扩展

### 9.1 新增单元类型

1. 在 `Model2D::Element.type` / `Element3D.type` 增加取值，解析器 `*ET` 读入。
2. 实现局部刚度 `localStiffness()`（与既有 `BeamElement`/`TrussElement` 保持同一局部自由度顺序与符号约定）、`transformation()`（三维）、`consistentMassMatrix`/`geometricStiffness`（如需）。
3. 在装配分支（`assembleLinear2D` / `assembleModel3D`）按类型分发；三维同时改 `buildDofMap3D` 的自由度数（若每节点自由度数不同）与 `elementGdofs3D`。
4. 实现内力提取（`elementForces2D` / `beamXxxLocalEndForces`）。
5. 加解析解/一致性测试，跑 `run_tests.ps1` 全绿。

### 9.2 新增约束/荷载

1. 在模型结构体加字段（如 `Constraint.k`、`Lane.autoCode` 的做法）。
2. 解析器读关键字；装配处施加贡献（如弹性支承在 `assembleLinear2D`/`assembleModel3D` 中加刚度与等效荷载）。
3. 若该约束影响自由度数或收敛（如非线性），同步改 `solveModel2DNonlinear`。
4. 加测试并补 `input_format*.md`、`USER_MANUAL.md`。

### 9.3 新增输出/可视化

在示例 `example_*` 中复用 `ModelResult`/`BeamEnvelope3D` 等结构输出；结构体新增字段（如 `BeamEnvelope3D.Bmin/Bmax`）后，同步更新所有消费该结构的输出点。

## 10. 关键设计决策备忘

- 先验后改：任何重构前先建黄金基线（解析解 + 逐字节回归 + ANSYS 对标），改动后三者必须全绿。
- 输入向后兼容：新增关键字/字段一律带缺省值，旧 `*.dat` 不修改仍可运行。
- 6 自由度与 7 自由度并存：不强行把 `BEAM` 扩成 7 自由度，新增 `BEAM7` 类型由用户选择，老路径零风险。
- 求解器边界：二维非线性/施工阶段逻辑复杂，保留独立求解器，只在单元/自由度/内力层统一。
