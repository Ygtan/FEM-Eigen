# 更新记录

本文件记录 FEM-Eigen 项目的新增功能和主要修改内容。

## 2026-11-09

### 实现：P3-3b 完整形状敏感度 `K_shape`（默认路径仍为低秩）

承接下方「取证」条目，本轮按 `USER_MANUAL.md §9.8` 的规格实现并通过四条
验收 oracle，**但没有改默认数值路径**。

- `fem_3d_nl.hpp` 新增 `ShapeSensitivityOptions3D` /
  `ShapeSensitivityResult3D`、`crShapeResidual3D()`（固定 u 装配
  `R = F_ext - F_int`）、`computeShapeSensitivity3D()`（中心差分装配
  `K_shape`，收敛构形切线一次分解 + nt 次回代得到 `G_full`）。
  只格式化 / 装配内核已有量，不在别处另算一遍力学。
- `SensitivityMode3D {LowRank, Full, FixedPoint}` 加到
  `ErectionGeometryOptions3D` / `CompletionStateOptions3D`：三条路径统一
  `δ = -(I+G)⁻¹(ΔX+u)`，只换 G，线搜索与不动点兜底完全相同。
  默认 `LowRank`，与既有算例逐位一致（`dat_byte_regression` 37/37 不变）。
- `cable_solver.cpp` 的 `*CABLE_BRIDGE` 新增可选关键字
  `SENSITIVITY=LOWRANK|FULL|FIXEDPOINT`（缺省 LOWRANK）；未知值明确报错，
  不静默回退。BridgeStudio 前端不写该关键字。
- **修正一个隐蔽口径错**：`K_shape` 的**列**只取 isT 平动自由度，但**行**
  必须保留全部自由自由度（含转角、弹性支承），只有硬约束行清零。初版把
  「非 isT 行」一起清零，等于假设转角不随 X 变；小模型差分对拍能过，
  实桥上 Full 却只线性收敛、6 轮且比低秩慢。修正后实桥 4 轮、
  `0.158 → 1.06e-3 → 5.15e-7 → 9.5e-12`，逐轮比值二次收敛。
  `verbose` 现打印 `nt` / `|K_shape|max` / `|G_full|max` /
  `σ(I+G)` 与 `cond`，便于诊断。
- 新增 `test_shape_sensitivity_3d.cpp`（CTest 名 `shape_sensitivity_3d`，
  44/44）：残差口径同源、非索锚点列与自由转角行非零、
  `G_full·v` 与差分 `du/dX·v` 三档步长对拍、三条路径逐轮误差与终值一致、
  默认路径逐位不变回归保护。数值见 §9.8。
- **发现几何闭合的刚体规范自由度**：无端部固结时，参考构形作刚体平移
  δ、位移取 -δ 不改变当前构形与索力，`(I+G)` 有零空间 ⇒ 制造线形是一族
  解而物理量唯一。新增 `testGaugeFreedom` 把它固化成回归；实桥塔墩固结
  天然打断，测试模型必须显式固结，否则「三路径 camber 一致」不成立。
- **决策：默认保持 `LOWRANK`**。四条 oracle 全达标，但实桥
  `FULL` 720 ms / 默认 317 ms（每轮 `2·nt` 装配 + `nt×nt` 稠密 LU，
  大自由度按 `O(nt³)` 恶化），而默认已 6 轮零 FALLBACK 稳定收敛；
  `FIXEDPOINT` 400 轮停滞在 1.98e-3 m，说明 Newton 步不可省。
  `FULL` 作为 opt-in 高精度 / 诊断路径（闭合误差再小 150 倍）。

验证：全量 `mingw32-make -j2` rc=0、`ctest` 7/7（新增该测试）、
`dat_byte_regression` 37/37；真实斜拉桥三路径结果见 §9.8 表。

### 新增：内层共转列式 Newton 过程量输出（`[NEWTON]`，不改算法）

以前 `solveModel3DNonlinear` 的荷载步 / 逐次 Newton 残差只在 `verbose`
控制台里打印，结果文件完全没有；`[LENGTH]` 只记录**外层**迭代（改 L0 或
制造线形 ΔX）。这次把**内层**过程量补上，两层迭代在计算书里分开读：

- `fem_3d_nl.hpp` 新增 `NewtonStep3D` / `NewtonTrace3D`，在 Newton 循环里
  逐次记录荷载步号 `k`、步内迭代号、λ=k/NSTEP、|R|、|du|、|Ftarget|、
  是否重建符号分解、本步是否收敛。`NonlinearResult3D.trace` 带出；
  `CableLengthResult3D` / `ErectionGeometryResult3D` /
  `CompletionStateResult3D` 用 `lastNewton` 透传**收敛那一轮**，
  `NLStageResult3D.newton` 每阶段各带一段。
- 新增 `newton_output_3d.hpp`（只格式化，不重算力学）：`[NEWTON]` 段含
  `KIND/STATUS/STAGE/NSTEP/NITER/MAXITER/TOLF/TOLD/CONVERGED/REFACT/
  RESFIRST/RESLAST` + 逐迭代 `S` 行；无内层过程时写 `ERROR=` 段，不静默。
- `cable_solver.cpp` 在 `[LENGTH]` 后写最终成桥态那一段；
  `example_3d_stage_nl.cpp` 在每阶段 `[ENDSTAGE]` 前各写一段（带 `STAGE=`）。
  线性梁格 / 桁架 / 拱桥走线性内核，本来没有这一段。
- **口径**：收敛判据是 `|R| ≤ tolForce·max(|Ftarget|,1e-12)` **或**
  `|du| ≤ tolDisp·max(1,|u|)`，每个荷载步都要满足；换荷载步时 |Ftarget|
  变大、残差允许跳升，整条曲线不是单调的。`REFACT` 只统计符号分解重建，
  数值 `refactorize` 不计。规范没有 Newton 迭代轮数的合格限值，过程量
  不构成设计合格判定。

留痕字段一律**不参与**判据、步长或收敛计算，与 `IterationTrace3D` 同一
条纪律；改动前后 Newton 收敛判据逐位一致。

### 核销：P3-2b「影响矩阵复用同一次分解」早已实现，注释过时

`solveModel3DNonlinear` 与索长迭代上方两处注释仍写「影响矩阵靠 nc 次扰动
重跑非线性求解」「见 TODO(P3-2b)」，与代码不符。实际 `cableSensitivity3D`
接收收敛构形的 `NLTangent3D`，闭式给出 `J = M·D`（一次分解 + nc 次回代），
冻结复用直到残差停滞才重建，重建次数已写进 `[LENGTH]` 的 `SENS JREBUILD`。
本次只改注释并写明「P3-2b 已实现」，不动算法。

### 取证：P3-3b（梁参考形状未计入敏感度 G）暂不改数值路径

> 历史记录：本节的「本版未实现」已由本文件顶部「实现：P3-3b 完整形状
> 敏感度」条目承接；下方「约 20%」是源码旧注释的估计，实测
> `||G_full - G_low||/||G_full|| = 0.17%`，详见 `USER_MANUAL.md §9.8`。

`fem_3d_nl.hpp` 的 `TODO(P3-3b)` 说低秩敏感度 G 只计索项、漏了梁参考
形状，误差约 20%，会拖累 Newton 收敛率。用真实斜拉桥 `MODE=COMPLETION`
算例复核当前 `[LENGTH]` 段：6 轮闭合误差
`0.158 → 0.0030 → 8.7e-5 → 2.2e-6 → 5.6e-8 → 1.4e-9 m`，逐轮收敛比
约 0.025，**每轮都是 Newton 步被线搜索接受，没有一次触发 FALLBACK 不动点
兜底**。即该算例里 G 的近似误差没有导致 Newton 失效，注释里的「20%」不能
直接推导出「必须改算法」。

完整计入需构造 `∂F_int/∂X` 的 nd×nd 算子并放弃现有 Woodbury 低秩结构，
属于架构级改动；现有不动点步 `ΔX ← -u` 仍是收敛保证。本次不盲改已验证
的成桥态数值路径，改为把实现规格与验证要求写进 `USER_MANUAL.md` §9.8，
下一轮先做差分 G 对拍再决定是否动手。

### 修复：刷新 3 处陈旧黄金基线，ctest 恢复 6/6 全绿（不改算法）

`dat_byte_regression` 曾报 3/37 失败。逐行取证确认三处都是**口径变更后基线
未同步**，不是算法回归，因此刷新基线而不改内核：

- `verify_influence_envelope_lane1.txt`（+5 行）、
  `verify_lane_nonlinear_envelope_lane1.txt`（+3 行）：`6a2d439` 给二维车道
  包络新增竖向位移规格 `D*_2`，原有 43 / 11 行内力数值一字未变，失败原因是
  「数字个数不一致 243→273 / 51→69」。
- `verify_3d_beam_live_live_envelope.csv`（`1,N,1,-0,-120` → `-115`）：
  `09c56bc` 修正 `envelope3DFromInfluence` 线荷载漏积分权重（补梯形
  `laneWeights3D`），`Q=-10` 项由 `10×2=20` 变 `10×1.5=15`，
  `20+100=120` → `15+100=115`，`max_rel=0.0417` 正好 `5/120`。
  **与冲击系数无关**：该算例 `lane.impact=1.0`；`8145a1f` 只改
  `compute3DStageLiveEnvelope`（位移包络），本算例走
  `envelope3DFromInfluence` 直接算轴力，从未乘 `1+μ`。内力包络按 JTG D60
  计冲击的口径没有变化。此前把它猜成「差一个 μ、需要人工拍板是否该含
  1+μ」是误判，记录在此，避免以后照着错误结论去改内核。

刷新前 34 个 PASS 全部是 `identical`（字节相同），故只重写这 3 个文件，
其余 34 个基线一字未动。`--freeze` 会把 37 个文件全量 `copy2`，CRLF/LF
行尾差异会把 34 个无变化文件也标成 `M`——**不要直接全量提交**，按实际
差异逐个刷新并保持原行尾。`ctest` 全套 6/6 通过（含
`dat_byte_regression` 37/37）。

**防呆**：任何改动 `example_*.cpp` 输出行、`defaultEnvelopeSpecs` 规格、
`envelope3DFromInfluence` 乘子或积分权重的提交，都要同步
`python tests/regression/compare_regression.py --freeze --bin <build>
--src <source> --baseline <dir>`，并在 CHANGELOG 写明每条基线变动的数值
依据；否则字节回归会长期挂红灯，掩盖真正的算法回归。

### 新增：索长迭代 / 制造线形几何闭合的逐轮过程量输出（不改算法、不改输入格式）

`fem_3d_nl.hpp` 里的两类外层迭代（影响矩阵 Newton 求下料索长 L0、
Woodbury 低秩 Newton + 回溯线搜索求制造线形预抬高 ΔX）一直把每轮误差
打到控制台，但 `cable_solver` 的结果文件只给 `converged / iterations /
closure_err / tension_err` 四个汇总标量，BridgeStudio 看不到收敛过程，
也无法区分「索力命中」与「线形闭合」。本轮把**内核已经算出来的**过程写成
`[LENGTH]` 段，数值口径一个字没改。

- 新增 `cable_length_output_3d.hpp`：`formatCableLengthSection3D(tr, ...,
  header)` 与 `cableLengthErrorSection3D(header, code)`；STATUS 由
  `cableLengthStatus3D` 单点判定（`no_iteration` / `not_converged` / `ok`），
  下游只读不算。段内 `KIND/STATUS/ERRKIND/NCABLES/NITER/MAXITER/TOL/
  CONVERGED`、`SENS JAVAIL/JREBUILD/JDIAGMIN/JDIAGMAX`、每轮
  `I <it> <err> <geomErr> <stepMax> <stepRelMax> <alpha> <lineSearch>
  <flags>`、每索 `E <elem> <Ttarget_N> <Tfinal_N> <L0_m> <chord_m>
  <Ltar_m>`，以及 `NJ=<n>` + `J <i> <j> <value>` 完整影响矩阵。
  flags 为 `REBUILD/CAPPED/FALLBACK/CONVERGED` 逗号串。
- `fem_3d_nl.hpp`：新增 `IterationStep3D` / `IterationTrace3D`，由
  `iterateCableUnstrainedLengths`（`CableLengthResult3D.trace`）与
  `iterateErectionGeometry`（`ErectionGeometryResult3D.trace`）在循环里顺手填。
  **留痕不参与判据、步长或收敛计算**；`solveCompletionState3D` 把同一份
  trace 透传出来。补记了几何闭合的「设计弦长」（只供报告量 geomErr）。
- `cable_solver.cpp`：在 `[END]` 前写 `[LENGTH]`。`MODE=TENSION` 走
  `KIND=cable_length`（`ERRKIND=tension_rel`，控制误差是最大相对索力偏差，
  输出 J = M·D）；`MODE=COMPLETION` 走 `KIND=erection_geometry`
  （`ERRKIND=closure_m`，控制误差是 max|ΔX + u|，单位 m，**不产出 J**）；
  `MODE=STATIC` 与悬索桥找形是闭式反算 L0 后一次解完，本来就没有外层
  迭代，分别写 `ERROR=no_iteration_static_mode` /
  `ERROR=no_iteration_suspension_formfinding`，不摆空表冒充「0 轮收敛」。
- 口径要点（必须随段进报告）：两种迭代的控制误差**单位不同、判据不同**，
  不能画同一根轴；`geomErr` 是成桥变形量级、只是结果量不参与判据——
  实测 TENSION 下索力收敛到 3e-8 而弦长偏差稳定在 0.15 m 左右，说明
  「索力对了但线形没闭合」；过程量只说明数值如何收敛，**不构成规范意义上
  的设计合格判定**，未收敛如实标 `STATUS=not_converged`。

### 新增：线性屈曲荷载因子全量输出（不改算法、不改输入格式）

`computeBucklingFactors3D` 一直在解 `(K+λ·Kg)·φ=0`，但只在
`example_3d_modal.cpp` / `bench/modal_bench.cpp` 里把 λ 打到控制台，
`example_3d_grillage.exe` 的结果文件里没有，BridgeStudio 拿不到。本轮
把**同一次特征求解**的结果与诊断写进结果文件：

- 新增 `buckling_output_3d.hpp`：`formatBucklingSection3D(ref, dm, result,
  header)` 写 `buckling`（`<stem>_summary.txt`）/ `[BUCKLING]`
  （`<stem>_cable.txt`）段：`NMODES/NDOF/BASIS`、`REF NLOAD/PREF/SUMFZ/
  NBEAM/NAXIAL/KGNZ/NFREE/COMPRESSION/CHOLESKY/STATUS`、
  每阶 `L <i> <λ> <Pcr_kN>`。`P_ref = Σ|Fi|`（N，不限方向），
  `Pcr = λ·P_ref/1000`（kN）；`SUMFZ` 只作诊断，竖向是 dof 3。
  `bucklingErrorSection3D` 在真正做不了时写 `ERROR=<code>`，不静默。
- `fem_3d.hpp`：
  - 抽出 `assembleGeometricStiffness3D(m, dm, U)`，除 BEAM 的一致几何
    刚度外补齐 **TRUSS/CABLE** 的 `Kg=(N/L)·[[P,-P],[-P,P]]`
    （`P=I-n⊗n` 横向投影）。原实现只遍历 BEAM，纯桁架桥 Kg 恒为零、
    永远报"没有屈曲模态"；`BEAM7` 的 14 自由度翘曲项无对应式子，
    显式跳过而不是塞 12×12 错位。
  - `computeBucklingResult3D` 在 λ 之外返回 `kgNonZeros / nFree /
    choleskyOk / compression`，`STATUS` 由 `bucklingStatus3D` 单点判定：
    `no_compression`（参考荷载下无受压软化）是物理预期，不是求解失败；
    `Kg` 全零时早退，避免稀疏无解后掉进 nf³ 稠密回退。旧接口
    `computeBucklingFactors3D` 保持不变。
  - 新增 `finalStageBucklingModel3D(m)`：成桥态构型 + 成桥态恒载合力。
    施工阶段是**增量法**，按 `(node,dof)` 累加 `st.loads`；有阶段时
    不读 `m.loads`。**坑**：`finalStageModel3D` 无阶段时 `return m`
    是连 `loads` 一起拷的，有阶段时只新建构型；不清空 `ref.loads`
    会把 `*CLOAD` 算两遍，实测悬臂柱 λ 正好减半（10.17 vs Euler
    20.33）而 Pcr 不变。
- `example_3d_grillage.cpp`：末尾写 `buckling` 段，前 6 阶。无参考荷载 /
  特征求解异常分别写 `no_reference_load` / `solver_exception`。
- `cable_solver.cpp`：写 `[BUCKLING] ERROR=nonlinear_state_not_supported`。
  成桥态是共转列式几何非线性平衡态，一致切线刚度已含成桥态几何刚度，
  基于未变形构型的一次特征值 λ 不适用；真实需求是非线性稳定全过程
  （荷载–位移极值点）分析，本版本未覆盖，不拿线弹性 λ 冒充。
- 结果段向后兼容：新段放在既有段之后，段内行不匹配任何旧解析正则；
  两个结果文件的其它字段一字未改。

**口径**：λ 是参考荷载整体放大倍数，不是规范稳定安全系数；JTG 3362
没有 λ 的合格限值。线性屈曲不含初始缺陷、材料非线性与几何非线性，
是稳定承载力的上限。

**验证**：20 段悬臂柱（E=206 GPa、Iy=0.4、L=100 m、轴压 1 MN）
`λ₁ = 20.33138614`，Euler 闭式解 `π²EI/(4L²)/P = 20.33138507`，
误差 5×10⁻⁸；`example_3d_modal.exe` 旧接口输出一致。空间桁架桥
（109 根 TRUSS）Kg 非零元 1030、求出 5 阶 λ。
构建：`mingw32-make -C build example_3d_grillage cable_solver`。

## 2026-11-08

### 新增：模态结果全量输出（不改算法、不改输入格式、不改 CMake 目标）

内核一直在解 Kφ=λMφ（`computeNaturalFrequencies3D`、`cable_solver.cpp` 的
`modalFromTangent`），但只为 JTG D60 汽车冲击系数留下一个竖向基频
`f1_vertical`，各阶频率与振型只存在于局部变量里。本轮把**同一次特征求解**
的结果写进结果文件，供 BridgeStudio 解析：

- 新增 `modal_output_3d.hpp`：`formatModalSection3D(dm, M, modal, header)`
  输出 `NMODES/NDOF`、`MASS X/Y/Z`（各方向平动总质量 rᵀMr，kg）、每阶
  `F <i> <f_Hz> <T_s> <mx> <my> <mz> <mnorm>`（有效参与质量 Γ_d²/(φᵀMφ)
  与 φᵀMφ 自检值）、每节点 `S <i> <node> <ux> <uy> <uz>`；
  `modalErrorSection3D` 在算不出时写 `ERROR=<code>`，不静默跳过。
  自由度顺序 [u,v,w,θx,θy,θz]，故 Z = 竖向 w。
- `fem_3d.hpp`：把 `compute3DLiveBeamEnvelope` 里的「取各 *STAGE
  activate/deactivate 后的成桥态构型 + `constrainDisconnected3D`」抽成
  `finalStageModel3D(m)`，两个调用方共用，保证模态与 `resolveLaneImpact3D`
  的 f1 使用**同一个构型**。无施工阶段时返回全模型。
- `example_3d_grillage.cpp`：末尾写 `modal` 段，成桥态构型 + 线弹性刚度
  前 8 阶。无质量 / 无自由自由度 / 特征求解异常分别写 `no_mass` /
  `no_free_dof` / `solver_exception`。
- `cable_solver.cpp`：在算 `f1Vertical` 的同一次 `modalFromTangent` 里
  就地格式化 12 阶（成桥态一致切线刚度，含索的 Ernst 效应与几何刚度），
  输出 `[MODAL]` 段；未收敛 / 无车道 / 无质量写 ERROR 段。
- 结果段向后兼容：新段放在既有段之后，段内行不匹配任何旧解析正则；
  `<stem>_summary.txt` 与 `<stem>_cable.txt` 的其它字段一字未改。

**验证**：10 m 简支梁 `verify_3d_modal.dat` 前 8 阶频率
119.954 / 128.199 / 146.913 / 387.766 / 479.865 / 587.712 Hz 与黄金基线
`tests/regression/baselines/verify_3d_modal.modal.txt` 逐项吻合，
f₁ 与 Euler-Bernoulli 闭式解 119.9534 Hz 误差 0.0007%。
构建：`mingw32-make -C build example_3d_grillage cable_solver`。

## 2026-11-07（续）

### 修正：活载位移包络不计冲击

- 二维车道包络（`example_influence.cpp`）与三维活载位移包络（`fem_3d.hpp`
  `compute3DStageLiveEnvelope`、`fem_3d_nl.hpp` `fillTangentGirderLive`）不再对
  **位移**乘 1+μ。
- 依据：JTG D60-2015 表 4.3.2 的挠度限值针对汽车荷载「不计冲击」本身的挠度；
  冲击系数只放大**内力/应力**（承载力组合用）。此前位移与内力同乘 1+μ，验算时
  再除掉，报告展示的「含冲击挠度」会误导。
- 内力/应力包络仍乘 1+μ，`impact_factor`/`f1_vertical`/`impact_mu` 仍正常输出。

## 2026-11-07

### 新增：二维活载位移、规范徐变与规范收缩

- 二维车道包络增加竖向位移。线性影响线的内力乘自动 1+μ、位移不计冲击（见上一条修正）。未写冲击系数时按竖向基频计算；写了不等于 1 的冲击系数时仍按手写值。
- `*CREEP` 同时给出 `RH`、`H`、`FCM` 时，名义徐变系数和 βH 按 JTG 3362 附录 C.2.1 计算，收缩按附录 C.1。没写这三项时仍是 φ0=2、βH=700。手写 `PHI=` 或 `ECS=` 优先，不与规范值相加。
- 非线性阶段里，混凝土梁的规范收缩只计该梁参与受力之后新发展的部分。`TS=` 是干燥开始龄期，缺省与加载龄期相同。索不收缩。线性阶段仍不读 `TIME`。
- 没写统一的 `H` 时，各混凝土截面用自己的 `HTH`（理论厚度，mm）。拉索和钢杆不折减。手写 `H` 仍覆盖全桥。

## 2026-11-06

### 新增：非线性施工阶段、切线活载、多龄期徐变

- 新程序 `example_3d_stage_nl` 调用 `solveModel3DStagesNonlinear`。每个阶段是当前构形上的全量平衡，位移不能按线性阶段累加。`TENSION=` 用无应力索长命中目标索力。
- 没有拉索时，阶段文件末尾的活载包络与线性梁格相同，并自动乘 1+μ。有拉索时，末阶段用收敛构形的切线刚度：竖直索与梁并联，斜索含方向余弦和几何刚度，垂度索用 Ernst 切线而不是 EA/L0。
- 未写冲击系数时，基频在这套刚度上按 JTG D60 计算。弹性支承留在模态里，不当成固结。
- 收缩是初应变，正值缩短，只缩短一次，自由缩短不留下轴力。徐变增量作用在已有应力上并留在后续平衡里。阶段写 `TIME=`、模型写 `PHI0=` 时，Δφ = φ(t)−φ(上一阶段)；本阶段新加的荷载要到下一龄期才开始徐变。
- 曲线仍是 φ=φ0[(t−t0)/(βH+t−t0)]^0.3（φ0=2，βH=700，t0 默认 7 天），不是按湿度、理论厚度取值的规范徐变。

## 2026-11-02

### 新增：悬索桥施工阶段分级找形（空缆 → 对称吊装 → 合拢 → 铺装）

- 悬索桥此前只做成桥态。本次 `cable_solver` 新增施工阶段模拟：主缆下料长度 L0
  一旦固定，从空缆态（只受主缆自重）起逐级加载，主缆逐级下沉、矢高从 f0 增大到
  成桥矢高 f。施工顺序：空缆 → 主跨加劲梁自跨中向两侧对称吊装一期恒载 Q1 →
  合拢（全桥 Q1）→ 二期铺装 Q2（成桥态）。
- 核心力学：成桥找形是「给定设计矢高、二分 H 命中矢高」的**设计态**；施工阶段
  找形是「给定 L0、二分 H 命中 L0」的**结果态**——矢高是结果不是输入。新增
  `ErectionState` / `findErectionState(m, cfg, L0target)`，`totalL0(H)` 沿三跨
  索多边形积分 + 闭式 `cableUnstrainedLengthForTension` 反算无应力长度，二分命中
  L0（L0(H) 单调递减）。
- 成桥 L0 反算用**纯找形口径**（`fit.zNew` 弦长 + `fit.T` 索力），不复用
  `solveModel3DNonlinear` 迭代出的 L0，否则铺装态不能精确重现成桥矢高。
- 输出新增 `[ERECTION]`（每阶段一行 stage/名称/H/矢高/跨中高程/吊索力范围）与
  `[ERECTION_SHAPE]`（每阶段主跨主缆高程序列，画下沉曲线用）。
- 关键修复：合拢/铺装阶段此前只用主跨吊点（`pts`）激活荷载，**边跨加劲梁恒载被
  清零**，铺装态 H 偏低 2.2%、主跨矢高偏大 2.3%（102.3 m vs 100 m）。改为用全部
  吊索梁端节点（`deckSet`）激活后，铺装态 H=226.95 MN、矢高 100.00 m 与成桥态
  精确一致——这是施工阶段模拟正确性的自洽判据。

## 2026-11-01

### 新增：汽车冲击系数（JTG D60，竖向一阶弯曲基频）

- 索承体系活载包络此前 `impact=1.0`（未计冲击力）。本次 `cable_solver` 在成桥态
  一致切线刚度上求结构自振频率，识别**竖向一阶弯曲基频**，按 JTG D60 公式算冲击
  系数 μ，活载乘 `1+μ`。
- 新增两个辅助函数（`cable_solver.cpp`）：`modalFromTangent`（在成桥态切线刚度
  `NLTangent3D.K` + 一致质量矩阵上做广义特征值 `Kφ=λMφ`，复用 `solveSparseModal` +
  稠密回退）、`verticalFlexuralFrequency`（识别主梁节点竖向位移 w 主导的最低阶，
  w² 占主梁平动 u²+v²+w² 比例 > 0.5 判定，识别不出回退第一阶 f1）。
- 索承体系基频**必须**用成桥态切线刚度：张紧索切线刚度远大于弹性刚度，线弹性装配
  刚度会把索算软、基频算低、μ 失真。质量矩阵用材料 `DENS` + `ADDMASS` 的一致质量。
- 冲击系数公式复用 `impactCoefficientFromFrequency`（JTG D60）：f<1.5→0.05，
  1.5~14→0.1767·ln f−0.0157，>14→0.45。注意 `impact` 字段是**放大倍数**（`1+μ`），
  不是纯 μ。默认 `IMPACT=1.0` 时用自动 `1+μ`，显式非 1.0 则尊重手动值。
- meta 输出 `f1_vertical`/`impact_mu`/`impact_factor` 三键，前端报告据此标注。
  默认算例（斜拉 200m、悬索 300m）竖向基频均 <1.5 Hz，μ 取下限 0.05、活载放大 1.05。

### 新增：横向分布系数（活载多车道横向布载）

- 索承体系活载梁内力包络此前是「纵向单车道 + 多车道折减近似」，横向布载只靠一个
  折减系数近似、不做横向影响面。本次在 `Lane3D` 增两个字段 `laneCount`（横向设计
  车道数）、`eccentricity`（偏载增大系数），`envelope3DFromInfluence` 的乘子从
  `impact·reduction` 扩为 `impact·reduction·laneCount·eccentricity`——单梁模型横向
  只有一个节点，严格横向影响面需要梁格/板壳，规范上用横向分布系数把多车道并排
  布载 + 偏载放大折进纵向影响线包络。
- 默认值（1, 1.0）等价于旧的单车道近似，向后兼容；`*LANE_3D` 新增 `LANE_COUNT=` /
  `ECCENTRICITY=` 关键字，`REDUCTION=` 语义明确为多车道横向折减系数（JTG D60-2015
  表 4.3.1-5）。索承体系建造器默认 3×0.78×1.15=2.691。
- 能力边界如实修正（写入 `input_format_3d.md`）：横向影响面**严格版**仍不做（需梁格），
  但横向布载已按系数法计入，偏载增大系数是经验值、写入假设清单。乘子为线性，内核侧
  由 BridgeStudio `verify_cable_3d.py` 的放大对拍锁定（lane_count=3 对 1 的 MyMax
  比值 = 2.6910，相对偏差 <1e-3）。

### 新增：索承体系活载梁内力包络（`[GIRDER_LIVE]`）

- 斜拉桥/悬索桥的 `cable_solver` 此前**彻底不含活载**，梁内力只有成桥态恒载单一
  状态。本次补齐活载梁内力包络（N/My/Vz），并在 `[GIRDEREND]` 之后新增
  `[GIRDER_LIVE] elem end NMin NMax MyMin MyMax VzMin VzMax`（注意 N 分量键统一为
  `NMin`/`NMax` 大写，与 My/Vz 命名一致——早期曾混用小写 `Nmin`，导致前端拼键取不到值）。
- 活载影响线必须在**成桥态切线刚度**上算，不能用线弹性装配。`solveModel3DNonlinear`
  返回的 `NLTangent3D.K`（张紧索切线刚度远大于弹性刚度）直接分解回代；三条新函数
  （`fem_3d.hpp`）：`scanTangentLane3D`（对非约束自由度加单位荷载，约束行必须跳过，
  否则回代成位移 1）、`tangentBeamResponseInfluence`（切线位移后用
  `beam3DLocalEndForces` 提取梁端内力增量）、`autoIdentifyGirderLane`（按
  `|Δx| ≥ max(|Δy|,|Δz|)` 自动识别水平纵向主梁，排除塔柱/横梁，缺省车道时兜底）。
- 包络公式复用梁格内核：`envelope3DFromInfluence` + `laneWeights3D` 梯形积分权重 +
  `jtgD60LaneLoadQk/Pk/Reduction`（JTG D60 车道荷载、多车道折减）。前端显式
  `*LANE_3D` 优先，空则自动识别主梁（`autoCode=true`）。
- 能力边界（如实标注，写入 `input_format_3d.md`）：①纵向单车道 + 多车道折减近似，
  不做横向影响面；②梁 P-Δ 一阶略去（恒载由非线性精确算入，活载增量用切线刚度求
  位移、再用线性弹性刚度提取梁端内力）；③组合应力验算（1.2S_G+1.4S_Q）未接入。
- `test_3d_nl.cpp` 第 11 节 `testCableLiveInfluence`：简支梁（纯梁无索）切线刚度必退化
  线性弹性，断言切线影响线命中手算 `x(L−x)/L`、与线性 `compute3DBeamMomentInfluence`
  逐位一致、两端支座 z 自由度被约束跳过。`test_3d_nl` 79/79。

### 新增：施工阶段曲线筋（`PRESTRESS=` 关键字模式，二次抛物线）

- 直线筋（`s=0`）只有轴压 + 端弯矩，**不产生上拱**；混凝土梁上拱来自曲线筋的
  跨内反向均布力。此前 3D 内核只做直线筋，本次补齐曲线筋。
- `Prestress3D` 增 `s`（矢高 = 跨中偏心 − 端部偏心，向上为正；`s<0` → 上拱）。
- 语法：直线筋保持位置模式 `<P>:<e>[:<loss>]`（向后兼容），曲线筋走关键字模式
  `P=...:E=...:S=...[:LOSS=...]`——矢高必须命名，否则 `P:e:0.85` 里 0.85 是 loss
  还是矢高会歧义。两者由值段是否含 `=` 自动判定。
- 等效荷载在直线筋三样之上叠加两样，必须同时施加才静力自平衡（竖向合力 = 0）：
  两端锚固切线竖向力 `4Pe·S/L`（两端同向）+ 跨内竖向均布 `q = −8Pe·S/L²` 节点化
  （首尾各半跨、中间节点全跨）。曲线筋要求链共线（抛物线定义在首尾弦上），不共线
  报错；上拱依赖跨内单元细分（单元内弯矩仍线性）。
- `test_3d.cpp` [26]：跨中上拱闭式解 `−Pe·S·L²/(12EI)`、与手工 CLOAD 逐位等价、
  `s=0` 退化为直线筋、曲线筋链不共线报错（直线筋链不共线仍按旧近似）、关键字/
  位置两种写法解析。`test_3d` 78 → 88 项。

## 2026-10-31

### 新增：施工阶段预应力（`*STAGE` 行内 `PRESTRESS=`，直线筋等效荷载法）

- 3D 内核此前**没有预应力**：`*PRESTRESS` 只存在于 2D 内核，3D 施工阶段只能
  算自重与调索，混凝土连续梁的纵向预应力无处施加。
- 现在 `Model3D::Stage3D` 增 `prestresses`，`*STAGE` 行内支持
  `PRESTRESS=<单元列表>:<P>:<e>[:<loss>]`：一条直线筋（常偏心），翻译成三组
  等效节点荷载——两端一对轴压 `Pe`、两端一对偏心弯矩 `Pe·e`（`Pe=P·loss`）。
- **不产生上拱**：直线筋只在锚固端给轴压与端弯矩，跨内无反弯；上拱来自曲线筋
  的跨内反向均布力（本版不做）。这一点写进 `input_format_3d.md` 与测试注释。
- 硬校验：`P<=0`、`loss` 越界、非梁单元、未激活、单元不构成开链，全部报错。
- 轴向取锚固端连线单位向量（要求直线链）；偏心弯矩绕横轴，方向由
  `Pe·e·(ẑ×c)` 给出——符号在 `test_3d.cpp` [25] 用悬臂梁闭式解锁死。
- `StageResult3D` 增 `prestressForces`，`example_3d_grillage` 增
  `stage prestress force` 段（`STAGE<n> P<束号> Pe=<有效预加力>`）。


### 新增：`[GIRDEREND]` 梁端内力，斜拉桥也有了

- 成桥态内核原先只在**悬索桥**、且只在 `DECKE0..DECKE1` 范围内算加劲梁内力，
  斜拉桥一根梁也没有——于是"索力很准、但答不出主梁扛不扛得住"。
- 现在两种桥型都算，新增 `[GIRDEREND]` 段：每端一行 `elem end N My Vz`
  （`end=1` i 端、`2` j 端），覆盖全部 `BEAM`/`BEAM7`。
- `[GIRDER]` 保持原格式与原语义（悬索桥摘要，各分量跨端取绝对值最大者），
  但它的 `(N, My)` 可能来自不同端，**不能用于应力验算**——这一点写进了文件头
  与 `input_format_3d.md`。
- 顺带修掉一个潜伏缺陷：原代码调用 `beamEndForcesAt(fl, 12, end)` 写死 12
  个自由度，`BEAM7`（14 自由度）的杆端力会被取错位置。现在按单元类型取。
- 梁内力对着**制造线形**（`[ERECTED]`）算，不是设计构形：`u` 是相对制造线形
  的位移，用设计构形会把预拱度当成变形。

### 新增：三维施工阶段调索（`*STAGE` 的 `TENSION=` / `L0=`）

- `Stage3D` 增加 `tensions`（单元 → 目标索力 N）与 `unstrained`（单元 → 无应力索长 m）。
  `*STAGE` 行内写法为逗号列表：`TENSION=12:3.2e6,13:3.1e6`、`L0=15:10.234`。
- `TENSION=` 语义是「本阶段末必须命中这个力」，内核在阶段内迭代无应力索长
  （默认 6 轮、相对容差 1e-4）；`L0=` 语义是「直接给下料长度」，索力是结果。
  同一单元两者互斥，另有四条硬校验（非索单元 / 未激活 / ≤0 / 互斥）直接报错。
- **为什么线性内核也要迭代**：`assembleStiffness3D` 把 `CABLE` 当桁架（只按 `EA/L`
  组装），根本不读 `e.L0`，改索长对结果毫无影响。目标索力必须翻译成初应变
  （`istrn += ΔP/(Eeff·A)`）才能进入方程。
- **张拉量按索力状态算**：`ΔP = N目标 − N当前`。早期写成
  `f(L_geom, L0新) − f(L_geom, L0旧)`，在分批张拉中把先张的索带到 +89%（结构已变形
  0.5 m 时两者差近一倍），已废弃。
- 跨阶段继承三样东西：索力状态 `Tcur`、无应力索长 `L0state`、**按节点累积的位移**
  （各阶段自由度编号不同，但「节点 n 沿 x 的位移」连续；弦长必须用累积变形后的位置，
  否则逐阶段调索会越调越大）。
- `StageResult3D` 增加 `cables`（`CableStageState3D`：`elem / L0 / chord / T /
  target / dev / Tdes / devDes`）。`Tdes` 是**跨阶段继承的设计目标**，`devDes` 用来
  量化「张好之后被后来的索带跑多少」——分批张拉这是物理事实，不是缺陷。
- 输出段 `stage cable tension`（`example_3d_grillage.exe`）。
- 索力是**状态量、不累加**，与梁内力（成桥态 = Σ 各阶段增量）性质相反；接口上没有
  「索力包络」这种东西。
- 校验：`test_3d` 新增 `test3DStageCableTension()`（竖向索 + 弹簧手算模型，13 条判据，
  含 4 条负例）。

## 2026-10-30

### 修复：成桥态内力/位移改为逐施工阶段累加

- `compute3DFinalStageBeamForces` / `compute3DFinalStageDisplacement` 原先只取
  `stages.back()`，等于丢掉前面各阶段的荷载贡献；荷载全在前一阶段时成桥态直接为 0。
- 改为 `Σ_k（第 k 阶段荷载在第 k 阶段构型上的增量）`，每阶段用各自的 `DofMap3D`。
- 单阶段模型输出不变（`dat_byte_regression` 逐字节回归通过）；多阶段模型结果修正。
- 详见 `MODIFICATION_LOG.md` 的「判决性反例」两算例。

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
