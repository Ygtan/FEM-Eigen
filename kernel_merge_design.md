# 2D/3D 内核合并设计

## 目标与原则

把 `fem.hpp`（2D 内核）与 `fem_3d.hpp`（3D 内核）中重复的梁单元、自由度映射、
装配、求解、内力后处理收敛到同一套实现，让“二维梁 = 三维梁位于 x-y 平面内、
仅绕 z 轴弯曲”的特例。

原则：**先验后改**。每一步都要保证
- `test_2d` / `test_3d`（解析解断言）全过；
- `dat_byte_regression`（31 个黄金基线逐字节回归）全过；
- `ansys_benchmark`（ANSYS BEAM4 对标）全过。

## 约定

统一 12 自由度梁单元局部自由度顺序：

    [u1, v1, w1, θx1, θy1, θz1, u2, v2, w2, θx2, θy2, θz2]

- 绕局部 z 轴弯曲（v, θz）对应 2D 梁的 (v, θ)，截面参数 `Iz` 对应 2D 的 `I`；
- 绕局部 y 轴弯曲（w, θy）对应 `Iy`，耦合项符号与 z 轴相反（右手系约定）；
- 扭转（θx）对应 `J`、`G`；
- 剪切参数 `PhiZ/PhiY = 12EI/(G·As·L²)`，0 时退化为欧拉-伯努利梁。

## 分阶段

### A. 统一单元接口（已完成）

新增 `BeamElement`（`fem.hpp`）：`localStiffness()`、`transformation()`、`stiffness()`、
`consistentMassMatrix(rhoA)`、`geometricStiffness(N)`。

等价性由 `test_3d.cpp::testBeamElementUnification()` 逐元素校验：
- 3D：与 `beam3DLocalStiffness` / `beam3DTransformation` 一致；
- 2D：平面内 (u,v,θz) 子空间与 `EulerBeam2D`（含 Timoshenko）、一致质量、
  几何刚度一致。

### B. 统一自由度表（已完成）

把 `DofMap3D` 提升为 `fem.hpp` 中共用的 `DofMap`（`DofMap3D` 保留为别名），
并为 2D 增加 `buildDofMap2D`（每节点 3 自由度 [u,v,θ]，编号与 `solveModel2D`
内部一致）。等价性由 `testDofMapUnification()` 校验：
“2D 装配刚度 == 3D 刚度限制到平面内 [u,v,θz] 的子矩阵”。

### C. 统一装配/求解/内力（进行中）

已落地：`BeamElement::localEndForces` 统一单元级内力 `fl = Kl·(T·ug)`，并由
`testUnifiedEndForces()` 校验与 `beam3DLocalEndForces` 逐元素一致。

剩余（需单独评估，见下）：把 `solveModel2D` 的装配/求解改用 `BeamElement` + `DofMap`。
注意 2D 内核额外承载了端部释放、分布荷载固端力、索/压杆松弛、施工阶段、预应力、
徐变收缩等三维内核不具备的特性，整体重写风险高、收益有限，建议保持“单元级与
自由度级统一、求解器各自保留”的边界，避免破坏已稳定的回归体系。

## 当前状态

- 阶段 A、B 已落地并通过全部回归；
- 阶段 C 单元级内力已统一；求解器级装配重写保留原边界，不强行合并。
