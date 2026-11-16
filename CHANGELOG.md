# 更新记录

本文件记录 FEM-Eigen 项目的新增功能和主要修改内容。

## 2026-11-15

### 新增：③b 施工阶段中间态起点的全过程稳定（CB3D23，第二十次）

用户 `按建议做`：继续 ③ 的 §2.5 边界收口，先做**施工阶段中间态
全过程稳定**，并补齐 Python 依赖。B3 结论**仍 `P` 不 `Y`**——本轮
只交付一个 `STABASIS=STAGE` opt-in 研究构型算例，§2.5 的位移控制 /
Riks、分支屈曲与分支切换、初始缺陷、材料非线性、动力 snap-through
仍未闭环，且算例是 3 杆 `TRUSS` 浅拱而非真实桥梁施工阶段。

- **同源阶段模型回传**：`fem_3d_nl.hpp::NLStageResult3D` 增加
  `Model3D model`，`solveModel3DStagesNonlinear()` 每阶段求解后写
  `row.model = mc`。此前只回传 `u` / `dm` / `cableL0` 与激活集合，
  主链路无法用与 `u_stage` 同源的模型装配 `F0` / `Fint(u0)`。
  这是 ③a `basisModel` 设计决策在施工阶段链路上的对应补齐。
- **主链路独立 opt-in**：`cable_solver.cpp` 新增
  `STABASIS=STAGE` + `STAGEID=<id>`（第 15 个 `*CABLE_BRIDGE` 键）与
  `runStageStabilityMode()`：先逐阶段跑共转列式全量平衡，按
  `STAGEID` 选中 `(pick->model, pick->u)`，再以 `u0 = u_stage` 走
  Crisfield 球面弧长路径；`F0` 取该阶段累计节点荷载 + 弹性支承 +
  CABLE 自重，`Finc` 只来自 `*STABILITY_LOAD`。未知阶段写
  `ERROR=stage_not_found`，**不静默降级**成成桥态或参考构形。
  输出文件 `<stem>_stage_stability.txt`，`[STAGES]` 逐阶段摘要 +
  `[STABILITY] BASIS=stage_state` + `STAGE=<id> NELEM=<n> NDOF=<n>`；
  `stability_output_3d.hpp::StabilityStateHeader3D` 增加
  `stageId` / `stageElems` / `stageDofs`，只在 `stageId >= 0` 时写
  `STAGE=` 行，其余两种基准输出逐位不变。
  三种 `BASIS`（`reference_configuration` / `completion_state` /
  `stage_state`）继续分列。
- **BS 侧同源导出与正式算例**：`ansys_stage_stability.py`（唯一
  常量表 `CASE_ID=cb3d23_stage_state`）与
  `cases/run_cases_stage_state.py`（`--generate-only` /
  `--check-only` / `--sync-to` / `--sync-only`）生成 4 节点 / 3 根
  `TRUSS` 浅拱的 FEM `.dat` 与两个 ANSYS `.inp`。3 个施工阶段
  （`STAGE_ACTIVE={1:[1,3], 2:[2], 3:[]}`），`STAGEID=2` 取
  3 杆全激活、阶段 3 荷载未施加的**中间态**。两侧同源给节点 2/3
  加 `UX` 导向约束，消除 4 平动自由度 / 3 杆长约束的一阶运动机构
  （否则 ANSYS 报 `Sparse solver minimum pivot`）。阶段荷载 /
  `PINC` / `CTRLTARGET` 取到轴向应变约 `2e-4`，避免 10 倍档下
  ANSYS `LINK180` 默认 Almansi 应变与 FEM 共转工程应变的系统性差异。
  `step1 TIME=1` 求同一阶段平衡态，`step2` 用
  `F_end = F_orig + LAMMAX·PINC`、`TIME = 1 + LAMMAX`，故
  `λ = TIME − 1`，不靠反力求和；`st_load` 荷载控制、
  `st_arc` `ARCLEN` 交叉验证，按递增自变量 `|ΔUY|` 插值 `λ` 对拍。
- **正式结果**：`[STAGES] [1,2,3]`、`STAGE=2 NELEM=3 NDOF=12`、
  起点 `RES0/|F0|=6.34e-13`、`CTRL0=-0.000952173872`（非零初态）、
  FEM `NPOINT=134` / `λ_end=2.57577583` / `|ΔUY|_end=0.00502513021`；
  两个 ANSYS job 均 0 error，`st_load` 27 点
  `worst|Δλ|=0.005124 mean=0.001891`、`st_arc` 25 点
  `worst|Δλ|=0.005109 mean=0.001881`，共同区间覆盖 FEM 路径
  99.7%。负例 `STAGEID=99` 写 `ERROR=stage_not_found` 且无平衡点。
  BS 报告 `72 checks, 0 failed`。
- **第 26 项 CTest `stage_stability_3d_benchmark`**（214 项检查）：
  `tests/benchmark/stage_stability_3d/` 冻结 15 份证据 + 黄金值
  `reference/cb3d23_stage_state.json`（由
  `build_reference_stage_state.py` 从冻结证据派生，不重跑 FEM/ANSYS）。
  `compare_stage_state_3d.py` 五层判据：证据 `sha256` / 字节数、
  黄金值自洽、物理判据、同源 APDL 锚点 / 两段 ANSYS 证据 /
  FEM↔ANSYS 逐点对拍、当前二进制重放两个 `.dat` **逐字节一致**
  （证据与历史黄金值同为 mingw64 工具链，故不放宽为
  `replay_tolerance.py`）。`CMakeLists.txt` 与 `run_tests.ps1`
  已注册，`tests/benchmark/README.md` 更新为合计 26 个 CTest。
- **依赖补齐**：`anastruct 1.7.1` / `python-docx 1.2.0` /
  `pypdf 7.1.0`（清华镜像；默认 PyPI 走代理报
  `ProxyError 502 Bad Gateway` 与 wheel 文件占用）。补齐前 53 个
  BS 守卫有 25 项因缺依赖 / 会话态失败，补齐后 51/53。余下两项是
  **与 ③b 力学无关的历史遗留 / 环境适配**，均按等价口径修正并在
  `B类后续项接入规格.md` §6 第 19 条登记：
  `verify_app.py` 的 `AppTest.session_state` 在 Streamlit 1.45.1
  改为 `SafeSessionState`（只支持 `in` / `[]`，不再有 `dict.get()`），
  新增 `_ssget()` 用 `in` + `[]` 严格等价替换 12 处旧写法，
  **功能断言一字未放宽**；`verify_opt_2d.py` 的「`## 2026-11-12`
  必须在 CHANGELOG 前 12 行」是 ③a 起过时的**位置**断言
  （当日小节持续追加到顶部），改为锁「该日期小节仍完整在档」，
  其余判据不放宽。修正后 `verify_opt_2d.py` 286/286，另 4 个 ③b
  关键守卫全绿（`verify_completion_stability_3d.py` 100/100、
  `verify_b_class_spec.py` 185/185、`verify_capability_matrix.py`
  282/282、`verify_stability_3d.py` 30/30）。
- **UI 兼容修复（与力学无关，但影响真实界面）**：`_ssget()` 让
  `verify_app.py` 从「第 202 行直接崩溃」恢复为可跑完，也**暴露出
  原本被崩溃掩盖的 7 条界面失败**——Streamlit 1.45.1 明确禁止
  `Expanders may not be nested inside other expanders`，
  `ui_3d.py` 的「与 MIDAS 互导」外层 `st.expander` 内嵌
  `_interop_block()` 的 `c1/c2.expander`。最小修复是把外层改为
  普通 `st.container()` + 一行 `caption`（保留内层 4 个折叠块、
  全部按钮 / 下载 / 试算功能与守卫断言，仅去掉「互导整块可折叠」
  这一层包装）。修复后 `verify_app.py` **109/109 通过**
  （此前 script run 中断，后续断言根本没跑到，63/70 只是已执行部分）。
  这是业务 UI 与新版 Streamlit 的兼容性问题，不是 ③b 引入，
  不涉及任何力学判据或能力矩阵；如实登记，不放宽断言。

### 新增：③a 起点同源修复 + 施工阶段共存核对 + 真实桥梁越过极值点

用户选 `3a`：先做实质技术收口，再评估 B3 是否具备推 `Y` 条件。
结论 **B3 仍 `P` 不 `Y`**（§2.5 的分支屈曲 / 初始缺陷 /
材料非线性 / 施工阶段中间态稳定未收口）。

- **修复 ② 登记的 `MODE=COMPLETION` 起点不同源缺陷**：
  `fem_3d_nl.hpp` 的 `CableLengthResult3D` /
  `ErectionGeometryResult3D` / `CompletionStateResult3D` 增加
  `Model3D basisModel`，把与 `u0` 同源的最终模型（制造线形 +
  迭代收敛后的 `L0`）回传给主链路；`cable_solver.cpp` 只在
  `STABASIS=COMPLETION`（opt-in）且 `MODE=TENSION/COMPLETION` 时
  用它装配 `F0` / `Fint(u0)`，`MODE=STATIC` 与悬索桥找形仍用主
  链路 `m`；`formatCompletionStabilitySection()` 增加自由度结构
  一致性检查，不一致写 `ERROR=invalid_initial_state`。
  CB3D10 `RES0` 由 `3.1492e6 N` 降到 `7.4001e-3 N`
  （`RES0/|F0|` `1.06 → 2.502e-9`），首步 λ 由负转正。
- 连带改进：第 22 项 CTest 的 `vshape_tension`（`MODE=TENSION`）
  起点 `RES0` 从 `0.0151 N` 降到 `2.84e-8 N`，λ 路径只有末位变化；
  已重新冻结 `cases/vshape_tension_cable.txt` 与
  `reference/cb3d20_completion_stability.json`（真实改进，不是浮点噪声）。
- **第 24 项 CTest `completion_basis_3d`**（110 项检查）：4 个 case，
  跨 `MODE` 独立锚点（`MODE=COMPLETION` 制造线形基准 vs
  `MODE=STATIC` 设计线形 + 冻结 `L0` 基准逐点
  `worst|dλ|=1.081e-6`、`worst|dΔCTRL|=3.131e-6 m`）、
  默认零污染、以及 CB3D13 悬索桥 5 个 `[ERECTION*]` 段与
  `[STABILITY]` **并列输出互不污染**（`RES0/|F0|=1.112e-9`、
  26 个平衡点）。核对结论是「共存」，不是「每个施工阶段中间态
  已具备全过程稳定分析」。
- **第 25 项 CTest `bridge_limit_point_3d`**（140 项检查）：
  BS 侧 `ansys_bridge_limit_point.py` /
  `cases/run_cases_limit_point.py` / `CB3D22_limit_point_vs_ansys/`
  在 CB3D21 同源桥上只缩塔截面（`TOWER_SCALE=0.012`）、增量荷载
  改为塔顶 2009 节点 `+X`，构造稳定敏感研究构型。FEM 球面弧长
  `NPOINT=553`、`λ_max=1.624941 @ dUX=0.200411`、下降段 207 步、
  `λ_end=1.242938`；ANSYS `lp_arc` 两段 `ARCLEN`
  `λ_max=1.622748 @ dUX=0.201182`，极值点 `dλ=1.349e-3`、
  `dUX=3.850e-3`，两个 job 均 0 error。λ 在极值点后下降，
  按单调量 `ΔUX` 插值 λ 逐点对拍（上升段 worst `5.493e-3`、
  下降段 worst `2.054e-2`）。同源锚点：16 根索
  `L0·(1+ε)` 与 ANSYS 弦长偏差 `2.43e-11 m`、`ε` 偏差 `1.83e-16`。
- **构建口径**：历史冻结证据由 `mingw64/bin/g++.exe` 产出，
  ③a 初版证据由 `ucrt64/bin/c++.exe` 产出，两套 CRT 的末位浮点
  不同。第 24 项冻结证据已用规范 mingw64 构建重新生成并保持
  严格逐字节重放；第 25 项（与 BS 算例目录共享 FEM 证据）改用
  `tests/benchmark/replay_tolerance.py` 的结构化跨工具链比较
  （`λ/CTRL/DS/PIVOT` 相对 `<1e-9`、`RES` 绝对 `<0.1 N`，
  非逐点行与 STEP/ITER/CONVERGED 仍逐位一致）。
  `run_tests.ps1` 新增 `-BuildDir` 参数（默认 `build`），
  干净 Ninja 构建用 `_d_build`。`.gitignore` 排除 `_c_build/` /
  `_d_build/`。
- **全量回归**：`_d_build`（显式 `-DCMAKE_CXX_COMPILER=
  C:/msys64/mingw64/bin/g++.exe`）全量编译后
  `ctest` **25/25**、`run_tests.ps1 -BuildDir _d_build`
  `ALL TESTS PASSED`。
- 编号口径：`completion_basis_3d` 覆盖既有 CB3D10 / CB3D13，
  黄金值改为 `case=CB3D10_BASIS_FIX` / `completion_basis_3d.json`；
  `CB3D22` 专属于 BS 侧极限点算例与 FE 侧第 25 项。
- **文档同步（③a 收尾）**：`USER_MANUAL.md` 补 §2.3 的 25 个 CTest /
  `-BuildDir` / mingw64 工具链口径与 §9.10 `STABASIS=COMPLETION`
  专项说明；`DEVELOPMENT_GUIDE.md` 更新 §6.3 测试索引、§7 工具链
  浮点漂移教训、§10 `basisModel` 同源模型设计决策。BS 侧
  `verify_completion_stability_3d.py::stage3_contract()` 新增
  5 条断言锁住这四份 FE 文档（68 → 73 项），
  `B类后续项接入规格.md` §6 第 18 条同步登记。

### 新增：真实斜拉桥成桥态两段对标（CB3D21，第十九次 / a 阶段 ②）

- 承接 ① 的 `STABASIS=COMPLETION`：把同一条主链路 opt-in 放到**真实双塔
  三跨斜拉桥**（50 节点 / 63 单元 / 16 索，与 CB3D10 同几何 / 材料 /
  恒载），与 ANSYS 18.0「step1 恒载成桥态 + step2 增量荷载」两段求解
  端到端对拍。② 完成后 B3 能力矩阵**仍 `P` 不标 `Y`**（③ 才推 `Y`）。
- BridgeStudio 侧新增 `ansys_bridge_completion.py`（唯一数据源
  `COMPLETION_CFG + EPS0=8e-4 + PINC=1e6 + LAMMAX=4.55`）与
  `cases/run_cases_completion.py`（`--generate-only` / `--sync-to` /
  `--sync-only`）：同一份 `ModelSpec3D(mode=STATIC, 设计线形)` 同时生成
  FEM `.dat`（16 根索显式 `L0 = L_design/(1+EPS0)`、
  `STABASIS=COMPLETION` 键、`*STABILITY_LOAD`）与 ANSYS 两段 `.inp`
  （`INISTATE,DEFINE,...,EPEL,X,eps ≡ EPS0`，step1 `TIME=1`，
  step2 `F_orig → F_orig − LAMMAX·PINC`）。
- 物理定标（四个探针已归档到 `completion_stability_3d/probes/`）：
  **不用反力求和**——完整斜拉桥 50 个 `F` / 34 个约束节点上
  `ΣRF,FZ` 随加载反而下降；荷载控制主 job `cs_load` 用
  `λ=(TIME−1)·LAMMAX` 与 `ΔUZ=UZ−UZ_step1` 对 FEM `[STABILITY]`
  逐子步插值对拍，`ARCLEN` job `cs_arc` 的 `TIME` 不是 λ、只比终点。
  两个 job 都恰有 2 个 `SOLVE`、2 条 `F,1016,FZ`、16 条 `INISTATE`。
- **已知缺陷（独立后续项，② 未修）**：`MODE=COMPLETION` 的索长迭代
  会改写 `m.elements[*].L0`，与收敛位移 `u0` 不同源（CB3D10 实测
  `RES0=3.15 MN`）；② 用 `MODE=STATIC` + 显式 `L0` 绕开，
  `RES0=7.97e-3 N` 与 `[NEWTON] RESLAST=5.96e-3 N` 同量级。
- 新增 `tests/benchmark/completion_bridge_3d/`：第 **23** 项 CTest
  `completion_bridge_3d_benchmark`，**2392 checks + 12 strings +
  35 criteria + 4 cases**；离线重跑 `cable_solver.exe` 复现 2 个 `.dat`
  （除 `elapsed_ms` 外逐行一致），同源独立锚点（16 根索
  `L0*(1+EPS0) ≡ 设计弦长`、`INISTATE eps` 偏差 `1.83e-16`、
  `F_end` 偏差 `3.31e-7 N`）、默认零污染（既有 CB3D10 / von Mises
  REFERENCE 逐位不变）、报告 `ALL CHECKS PASSED`、两 job 0 error。
  实测 28 子步 `worst|ΔUZ| = 1.85e-5 m`（容差 `1e-4`）、
  初始斜率相对差 `1.37e-4`、`ARCLEN` 与荷载控制终点
  `|d| = 6.9e-8 m`；`--live-ansys` 临时目录复跑两 job
  `max|d| = 0 / 1e-11`，不覆写冻结 CSV。
- `CMakeLists.txt` / `run_tests.ps1` 注册第 23 项；本机 `ctest` 23/23、
  `run_tests.ps1` `ALL TESTS PASSED`、干净 Ninja 自动生成
  23 条 `add_test` / 40 个 exe。同步 `tests/benchmark/README.md`
  （22→23 CTest）与 BridgeStudio 六份文档 / 守卫。
- 如实边界：CB3D21 路径**单调硬化、无极限点**，ANSYS 荷载控制即可
  覆盖全程；越过极限点 / 下降段能力由 L / N 阶段独立基准题与
  FEM 球面弧长内核证明。施工阶段链路、分支屈曲 / 初始缺陷 /
  材料非线性仍未覆盖。

### 新增：成桥态起点的全过程稳定分析（`STABASIS=COMPLETION`，第十九次 / a 阶段 ①）

- 用户按「先 b 后 a」完成 b（CB3D19）后选 **a1 三阶段**，本轮只交付 ①：
  规格/守卫先行 + FE 内核与主链路 opt-in + 第 22 项 CTest。
  ②（真实桥梁 ANSYS 两段对标）与 ③（B3 推 `Y`）**尚未做**，能力矩阵
  B3 该行**仍 `P` 不标 `Y`**。
- `stability_path_3d.hpp` 新增 `ArcLengthState3D`（`u0` / `F0` / `Finc`），
  给 `solveArcLengthPath3D` 加**带默认值的第三参数**。平衡方程从
  `R = λ·Fref − Fint(u)` 推广为 `R = F0 + λ·Finc − Fint(u)`、`u(0)=u0`；
  Crisfield 球面弧长约束仍对 `du`，一致切线与 `crElementInternal` 完全复用。
  段内 `CTRL` 改成**相对起点的控制位移增量**、`CTRLTARGET` 也按增量判定；
  参考构形基准起点为 0，**第 19/20/21 项 CTest 与全部冻结输出逐位不变**。
  `ArcLengthResult3D::res0` 记录起点闭合残差 `|F0 − Fint(u0)|`，不补 0。
- `stability_output_3d.hpp` 给 `formatStabilitySection3D` /
  `stabilityErrorSection3D` 增加带默认值的 `basis` 参数，新增
  `StabilityStateHeader3D`。`completion_state` 段多一行
  `STATE CTRL0= / F0_NORM= / FINC_NORM= / RES0=`，`CTRL0` 可与 `[U]` 自洽对拍。
  `BASIS=reference_configuration` 与 `BASIS=completion_state` 的 λ
  **不是同一个量**（前者放大总参考荷载、后者只放大增量荷载），段里分列。
- `cable_solver.cpp` 新增 `STABASIS`（默认 `REFERENCE`）与独立段
  `*STABILITY_LOAD`（`node dof value`，**不进 `m.loads`**，成桥态恒载 `F0`
  不被增量荷载污染）。`STABASIS=COMPLETION` 时**不提前 return**，先走完
  找形 / 索长迭代 / 内层 Newton，再在主链路 `[END]` 前写
  `[STABILITY] BASIS=completion_state`；`[U]` / `[NEWTON]` / `[LENGTH]` /
  `[BUCKLING] ERROR=nonlinear_state_not_supported` 全部保留。
  `F0` 与 `solveModel3DNonlinear` 的 `Fext` 逐项同源。新增 4 类 `ERROR=`：
  `zero_incremental_load` / `unsupported_stability_basis` /
  `completion_not_converged` / `invalid_initial_state`。
  `fem_3d_nl.hpp` 仍是荷载控制 Newton、`NonlinearOptions3D` 默认值零改动。
- 新增 `tests/benchmark/completion_stability_3d/`：V 形双索单自由度体系
  （`B=2.0`、`H=0.6`、参考构形带预拉 `L0=0.995·hypot(B,H)`），
  `MODE=STATIC`（显式 `L0`）与 `MODE=TENSION`（索长迭代命中 `T*`）
  两条成桥态闭环各 241 个平衡点。独立闭式
  `λ(v)=[P(v)−F0]/PINC` 逐点对拍 `worst|dλ|=1.505e-13` / `3.515e-08`；
  第二套独立锚点 `K_yy(V0)` 与中心差分 `dP/dv` 相对差 `3.964e-10`，
  路径首点割线斜率 `51.156` 对解析 `51.150`。
  冻结 **2963 checks + 18 strings + 18 criteria + 6 cases**，另含
  「无增量荷载 / 未知基准 / 成桥态未收敛 / 默认零污染」四类负例。
- 默认零污染实测：既有 CB3D10 结果文件 1257 行除 `elapsed_ms` 逐位一致、
  M 阶段 von Mises `REFERENCE` 基准 `_stability.txt` 逐字节一致、
  `*STABILITY_LOAD` 存在时默认态 `[U](2,2)` 仍是 `−0.1`（未混进 `m.loads`）。
- `CMakeLists.txt` / `run_tests.ps1` 注册第 22 项 CTest（离线、不跑 ANSYS）；
  `input_format_3d.md` 新增「成桥态起点基准」小节（第 14 个键、
  `*STABILITY_LOAD` 段、`STATE` 行、4 类错误码、两种 λ 分列口径）；
  `tests/benchmark/README.md` 同步 21→**22 个 CTest**（exe 仍 **40**，
  本轮不新增 exe）。BridgeStudio 侧新增
  `verify_completion_stability_3d.py`（32/32 纯离线守卫），并更新
  `stability_3d.py` 唯一解析器（ERROR 段 `BASIS` 通配、`STATE` 行解析、
  4 类中文错误码、`BOUNDARY_NOTE` 分列两种基准）。

## 2026-11-14

### 新增：共转 CR 梁刚接框架后屈曲独立对标（`frame_probe.exe`，第十八次 / N 阶段）

- 把 L 阶段（CB3D17）的 `stability_path_3d.hpp::solveArcLengthPath3D` 与
  K 阶段的共转 CR 梁单元组合起来，做 **Williams toggle 型刚接框架**的
  后屈曲全过程：两根各 24 个 `BEAM` 单元、支座端铰支释放 `ROTZ`、
  顶点与中间节点刚接、全节点约束 `UZ/ROTX/ROTY` 保持 XY 平面，
  顶点施加 `-Y` 参考荷载。新增 `frame_probe.exe` / `frame_probe.cpp`
  读 `*.frame` 纯 ASCII 参数文件，输出 `P` 行（逐平衡点
  λ / 控制位移 / 残差 / 迭代数 / 弧长 / LU 量级比）、`S` 行派生量与
  `E` 行错误码；**不修改** `fem_3d_nl.hpp` / `stability_path_3d.hpp`
  / `cable_solver.cpp` 任何一行，力学完全复用内核。
- 与 ANSYS `BEAM188 + NLGEOM + ARCLEN` 端到端对标（2 个 job）：
  `KEYOPT(3)=2` + `MP,GXY = G×1000` 退化到 Euler-Bernoulli 口径，
  `IY=IZ` 对称截面消除与 BEAM188 缺省局部轴约定的差异；
  两边同源取 `NELEM=24`（8→12→16→24 离散收敛诊断归档在
  `tests/benchmark/frame_postbuckling_3d/probes/`）。
  正式判据限定 `0 ≤ v ≤ VMAX=0.95`，避开 `v→1` 强硬化区
  `dλ/dv` 过大导致的插值/采样密度放大误差。
- **极限点取第一个 λ 局部极大**：刚接框架翻越后弯曲刚度强烈再硬化，
  全局 `max(λ)` 通常落在路径终点而不是极限点；`lambda_global_max`
  只作对照量输出，不作判据。
- 独立锚点：梁框架没有桁架 `P(v)` 闭式解，改用**初始切线刚度 `k0`**，
  三套实现互相对拍 —— `build_reference_frame.py` 自带纯 Python 2D
  线性刚架刚度矩阵、`ansys_bridge_frame.py::linear_frame_stiffness`
  第二套纯 Python 实现、`frame_probe.cpp::linearSlope` 在参考构形上
  组装内核 `K0` 的解。实测 `k0_analytic=1.211946322189e7 N/m`，
  FEM `k0_fem` 相对差 `1.907e-13`，`nodes=49 / apex=25`。
- 第 21 项 CTest `frame_postbuckling_3d_benchmark`：
  `tests/benchmark/frame_postbuckling_3d/` 冻结 **12513 项黄金值**、
  4 个字符串、37 条 criteria、6 个 case（含 `bad_noload` 负例与
  ANSYS 日志 `errors=0` 断言），离线重跑 `frame_probe.exe` 复现
  13958 个数值 `worst|d|=0`；`--live-ansys` 在临时目录复核、
  不覆写冻结 CSV。`run_tests.ps1` 同步追加该段。
- **如实边界**：这是**独立 `frame_probe.exe` 探针**，不是
  `cable_solver.exe` 成桥态闭环 / 施工阶段链路；能力矩阵 B3 仍
  **标 `P` 不标 `Y`**，把全过程稳定分析接进成桥态闭环属于后续阶段。
  既有 20 项 CTest 与全部冻结算例逐位不变，`STABILITY=0` 默认
  零污染口径不受影响。

### 新增：B3 非线性稳定全过程 `[STABILITY]` 段（`STABILITY=1` opt-in，第十七次）

- 把 L 阶段（CB3D17）已用 ANSYS `ARCLEN` 对标过的
  `stability_path_3d.hpp` 球面弧长法，以 `STABILITY=1` 接入
  `cable_solver.exe` 主程序：`*CABLE_BRIDGE` 增 13 个 opt-in 键，
  在 `readModel3D` 参考构形上跟踪完整荷载–位移平衡路径，
  输出 `<stem>_stability.txt` 的 `[STABILITY]` 段（逐平衡点
  λ / 控制位移 / 残差 / 切线奇异指标 / 收敛标志 + `ERROR=` 双函数）。
- **默认全关**，既有成桥态闭环 / `[NEWTON]` / `[LENGTH]` / `[OPT]` /
  默认参数与全部既有结果文件逐位不变；
  `fem_3d_nl.hpp::solveModel3DNonlinear` 仍是荷载控制 Newton，
  `NonlinearOptions3D` 默认值一字未动。
  `[BUCKLING] ERROR=nonlinear_state_not_supported` **原样保留**：
  线性屈曲对共转成桥态不适用是独立的物理结论，不与全过程极值点混淆。
- 上面 §2026-11-09 附近「真实需求是…本版本未覆盖，不拿线弹性 λ 冒充」
  是当时的历史事实记录，一字不改仍予保留；**自第十七次起该结论已被
  `STABILITY=1` opt-in 取代**，但能力矩阵仍标 `P` 不标 `Y`：
  只有球面弧长法，位移控制 / Riks / 分支屈曲 / 初始缺陷 / 材料非线性 /
  CR 梁后屈曲 / 真实桥梁全过程对标均未覆盖，且不等于成桥态闭环已具备
  全过程稳定分析。
- `input_format_3d.md` 新章登记 13 键、`[STABILITY]` 段格式与 7 类
  `ERROR=` 码；第 20 项 CTest 冻结 von Mises 两杆浅桁架解析锚点、
  默认零污染断言与 ERROR 负例。

## 2026-11-12

### 新增：2D 平面索力优化 [OPT2D] 段（`--opt-detail` / `--opt-detail-only`，opt-in）

- 第十四次（B2「2D 平面索力优化 [OPT2D] 一层」）是 2D **首次新造第二层
  力学**：第十三次已新造 `[LENGTH2D]` 平面索长闭合一层，本轮在其无垂度
  线性杆本构之上再新造 `[OPT2D]` 平面索力优化一层。第十一次「本轮不新造
  力学」、第十三次「`[OPT]` 在 2D 仍不做、仍标 `2D=N`、不造占位」都是
  当时的历史事实记录，一字不改、仍予保留；**自第十四次起该结论已被
  `[OPT2D]` 取代**：能力矩阵该格 `2D=N` → `2D=P`，但仍 `P ≠ Y`。
- 用户拍板的力学口径（逐条如实落到内核与文档）：
  1. **目标**：只做 `MINMAXM`，最小化**用户显式登记的 `BEAM` 单元端弯矩
     集合**的 `max|M|`；内核不自动识别「主梁」、不补零；**不做
     `MINENERGY`**（2D 无裸结构弯曲应变能柔度 `Fm` / 影响矩阵）。
  2. **约束**：只做**显式索力箱** `TMIN` / `TMAX`（必填、
     `0<TMIN<TMAX`、`BND_SRC=USER`），缺一写 `ERROR=missing_tension_bounds`；
     **绝不猜 3D 的 `1.5×RIGID`**；2D 无 `DECKLOAD` / 塔顶 `ΣT·cos` /
     `SYSGUARD` 三判据，**不造占位**。
  3. **算法**：变量是无量纲名义索力 `x=T/Tscale`，经
     `Lu(x)=Lgeom/(1+Tscale·x/EA−istrn)` 落到第十三次无垂度线性杆本构；
     **中心有限差分 Jacobian + IRLS(p=2/4/8/16) + 投影梯度**。段内明写
     「**不是 3D 的 `J=M·D` 闭式影响矩阵**」，也**不等于 3D
     `[OPT]`→`[LENGTH]`→`[NEWTON]` 三层迭代**。
  4. **写回边界**：`APPLIED=1` 只表示建议 `T*` 被采纳，命中 `T*` 属
     `[LENGTH2D]` 职责；`[OPT2D]` 只定目标、**不写回 `.dat` / `.result`**；
     `T` 行只是过程量不能当设计索力；峰值不降写
     `objective_not_improved` / `APPLIED=0`。`Nfinal ≠ T*`：结构柔性使
     `Lu(T)` 闭式只是近似。
- 新增输入段 `*CABLE_OPT2D`：参数行
  `OPT2D=1 TMIN=... TMAX=... TOL=... MAXITER=... IRLS=... ROUNDTOL=...`
  + 逐行 `<CABLE单元号> <T0_N>` 显式登记待优化索 + 逐行
  `M <BEAM单元号> <end=1|2>` 显式登记纳入 `max|M|` 的梁端弯矩。
  `fem_2d.cpp::iterateCableOpt2D()` 求建议索力 `T*`，不做命中。
- 新增 `cable_opt_output_2d.hpp` 写 `[OPT2D]` 段，字段按 2D 事实取舍：
  `I` 行（逐 IRLS 轮 `p` / `max|M|` / 投影梯度步）、`T` 行（每索
  `T` / `TMIN` / `TMAX` / `BND_SRC`）、`M` 行（每个登记梁端 `M`），
  及 `STATUS` / `CONVERGED` / `CONVBY` / `APPLIED` / `REJECT` 摘要字段；
  无契约时写 `ERROR=no_opt_request`，**不摆空表**。
  `example_mixed.exe <dat> --opt-detail` 把该段写进独立
  `<stem>.opt.txt`；`--opt-detail-only` 求解后直接返回、
  **不写 `.result`**、不打印结果明细。
- 验证：`mingw32-make -C build -j4` rc=0；`ctest` 7/7（含
  `dat_byte_regression` 37/37）。两个新算例：`verify_opt_2d.dat`
  正例 `T*=TMIN=200000`、`MPEAK0=1304347.826`、
  `MPEAK=260869.5652`、`APPLIED=1`、`CONVERGED=1`；
  `verify_opt_2d_negative.dat` 负例 `objective_not_improved` /
  `APPLIED=0`。**三条零污染证据通过**：默认 `.dat` / 默认 `.result` /
  默认 HTML 逐位不变，12 格 `DEFAULT_DAT_SHA` 不刷新，37 个黄金基线
  不刷新；双层 opt-in 默认全关。
- 口径：`[OPT2D]` 是 2D opt-in 平面索力优化一层，**不等于 3D `[OPT]`**
  （2D 仍缺 `SYSGUARD` / `MINENERGY` / 影响矩阵 / 垂度 / Ernst /
  共转成桥态 / 三维制造线形），能力矩阵 `2D=P` 不是 `Y`。
  BridgeStudio 侧 `opt_2d.py` 只解析、不重算、不补零，
  `verify_opt_2d.py` 守卫由用户另行编写。

## 2026-11-11

### 新增：2D 平面索长闭合 [LENGTH2D] 段（`--length-detail` / `--length-detail-only`，opt-in）

- 第十三次（B2「2D 平面索长闭合 [LENGTH2D] 一层」）是**首次真正新造
  2D 力学**：此前 2D 连索长输入契约都没有，不是第九 / 十 / 十一次那种
  「内核已算、只补留痕 / 补段」。第十一次「本轮不新造力学」是历史条目，
  一字不改、仍为当时的事实记录。用户拍板只做 `[LENGTH2D]` **一层**，
  `[OPT]` 带约束目标索力优化在 2D **仍不做、仍标 `2D=N`、不造占位**。
- 新增输入段 `*CABLE_LENGTH2D`：参数行
  `LENGTH2D=1 TOLT=... MAXLENITER=...` + 逐行
  `<CABLE单元号> <目标索力T*_N>`（T*>0、单元必须为 CABLE）；
  `Model2D` 增 `Element::Lu`（无应力下料长度 m，独立于几何弦长，与
  `ISTRN` 不是同一个量）、`LengthTarget` / `lengthTargets` /
  `lengthEnabled` / `lengthTolT` / `lengthMaxIter` /
  `LengthStep2D` / `LengthCable2D` / `LengthIterationResult2D`。
  **`Lu<=0`（未指定）时本构逐位不变**：仍是原线性杆
  `N = E·A/L·ΔL + E·A·ISTRN`；只有给定 `Lu>0` 才切到
  `N = EA·((L-Lu)/Lu + istrn)` 的无垂度线性杆本构。
- `fem_2d.cpp::iterateCableLength2D()`：闭式初值 `Lu0 = L*/(1+T*/EA)`，
  **中心有限差分 Jacobian + LM 阻尼最小二乘**求根（不是 3D 的闭式影响
  矩阵 `J = M·D`）；唯一判据 `max|N/T*-1| <= TOLT`，
  `ERRKIND=tension_rel`（无量纲），不掺几何闭合误差 / 线形判据。
  无垂度、无 Ernst 等效模量、无共转成桥态，不等于 3D 的
  `[OPT]` → `[LENGTH]` → `[NEWTON]` 三层迭代。
- 新增 `cable_length_output_2d.hpp::formatCableLengthSection2D()` 写
  `[LENGTH2D]` 段，字段按 2D 事实取舍：**不写 sag / Ernst / SENS /
  NJ / J / geomErr / STAGE 占位**；无契约时写
  `KIND=none STATUS=not_applicable` + `ERROR=no_length_request`，
  **不摆空表**。`example_mixed.exe model.dat --length-detail` 把该段
  写进独立 `<stem>.length.txt`；`--length-detail-only` 求解后直接返回，
  **不写 `.result`**、不打印结果明细。默认用例不带开关，`.result` /
  控制台 / 37 个黄金基线**逐位不变**；
  `tests/regression/compare_regression.py` 的 CASES 不加该算例、基线不刷新。
- 验证：`mingw32-make -C build -j4` rc=0；`ctest` 7/7（含
  `dat_byte_regression` 37/37）。新算例 `verify_length_2d.dat`
  （两端固定单根 CABLE + 非零 `ISTRN`，使闭式初值不直接命中）以
  `Lu = 10/(1+T*/EA−istrn)` 为**闭式解 oracle**：首轮残差
  `0.04 = EA·istrn/T*`，真实迭代 3 轮收敛到 `3.83e-13`，
  `Lu = 9.952229299 m`、`N = T* = 1e6 N`。
- 口径：过程量只说明数值如何收敛、不是规范意义上的设计合格判定；
  2D `[LENGTH2D]` 不等于 3D 三层迭代已移植，`[OPT]` 在 2D 仍没有
  对应力学、不造占位。BridgeStudio 侧 `length_2d.py` 只解析、不重算、
  不补零，下料长度由内核迭代求出。

### 新增：2D 内层 Newton 逐迭代过程量段（`--newton-detail` / `--newton-detail-only`，opt-in）

- 第十一次（B2「三层迭代过程量」）**只闭环内层 Newton 一层**。只读调研
  确认：2D 内核有 `*ANALYSIS TYPE=NONLINEAR` 的方案C Newton-Raphson
  （`fem_2d.cpp::solveModel2DNonlinear`），但逐迭代的 `dU` / `alpha` /
  `lamTarget` / `Fint` / `Fext` 全是局部变量用完即丢，`writeResult()` 只写
  5 个纯文本段；`[LENGTH]` 索长 / 几何闭合迭代与 `[OPT]` 带约束成桥索力
  优化在 2D **连输入契约都没有、内核根本没算**（无 `OPTIMIZE=` / `L0=`
  迭代 / 影响矩阵 / 体系守卫）。本轮**不新造力学、不写不适用段占位**：
  两格能力边界继续如实标 `2D=N`，Python 侧也不得从位移反推过程量。
- `fem_2d.hpp` 新增 `NewtonStep2D` / `NewtonTrace2D`，`ModelResult` 挂
  `newtonTrace` 字段；`fem_2d.cpp` 的 Newton 循环里只做 `push_back` 留痕，
  **判据 `resRef = tol*max(1,|Fext|)`、`duOk`、`alpha` 线搜索、收敛与荷载
  步细分计算一字未改**，三个求解入口签名都不动。残差已达标、本步未求解
  增量时 `|du|=0` 是事实值，不是补零冒充。
- 新增 `newton_output_2d.hpp::formatNewtonSection2D()` 写独立 `[NEWTON2D]`
  段：`BASIS=2d_global`、`KIND/STATUS`、头部 `NSTEP NITER MAXITER TOL
  CONVERGED RESFIRST RESLAST LINECUTS STEPREFINE LOADRATIO`、逐迭代
  `S <loadStep> <iter> <lambda> <|R|> <|du|> <|Fext|> <alpha> <lineCuts>
  <converged>`，`setprecision(10)` 与 3D `[NEWTON]` 同源。字段按 2D 事实
  取舍：2D 用直接求解、**不写 `REFACT`**（无符号分解重建）；只有一个 tol
  同时管残差与位移增量、**写 `TOL` 不拆 `TOLF`/`TOLD`**；主算例没有阶段
  维度、**不设阶段占位、不造阶段号**。无 Newton 迭代时写
  `ERROR=no_newton_solution`，不摆空表、不输出零行冒充收敛。
- `example_mixed.exe model.dat --newton-detail` 把该段写进独立
  `<stem>.newton.txt`，`.result` / 控制台 / 默认黄金基线**逐位不变**；
  `--newton-detail-only` 求解后直接返回，不写 `.result`、不打印结果明细，
  供前处理在 `example_stage.exe` 主路径之外单独调用（避免覆盖主结果）。
  过程量必须求解后才有，与 `--loss-detail-only` 的「不求解只取明细」不同。
  `tests/regression/compare_regression.py` 的默认用例不带新开关、基线不刷新。
- 口径：**过程量只说明数值如何收敛，不是规范意义上的设计合格判定**
  （规范没有 Newton 迭代轮数的合格限值）；本段只覆盖内层 Newton，
  **不等于 3D 的 `[OPT]` → `[LENGTH]` → `[NEWTON]` 三层迭代已移植到 2D**。
  BridgeStudio 侧 `newton_2d.py` 只解析、不重算 `K·ΔU=R`。
- 验证：`mingw32-make -C build -j4` rc=0；`ctest` 7/7（含
  `dat_byte_regression`）；`verify_cable_nonlinear.dat` 默认 / opt-in
  `.result` 逐位一致，`[NEWTON2D]` 出 3 条 `S` 行（|R| 2000 →
  0.2497499531 → 1.180023901e-08，末行收敛），`--newton-detail-only`
  不写 `.result`；线性算例 `verify_prestress_loss.dat` 如实写
  `ERROR=no_newton_solution`。

### 新增：2D 分项预应力损失明细段（`--loss-detail` / `--loss-detail-only`，opt-in）

- 第十次（B2）：`prestress_loss.hpp::computePrestressLosses()` 的六项公式本来
  就被 2D 调用，但 `fem_2d.cpp::processPrestress` 只取
  `effectivePrestressRatio()` 标量，`PrestressLossResult` 当场丢弃，2D 状态层
  没有可写的明细容器。本轮在 `Model2D` 增加 `prestressLossDetails`（每次
  `processPrestress` 先清空、再按 `LOSSREF` 束重建），
  `loss = effectivePrestressRatio(...)` 一行一字不动——明细是**另取一份
  同源结果**，不改任何数值路径。
- 新增 `prestress_loss_output_2d.hpp::formatPrestressLossSection2D()` 写
  `[PRESTRESSLOSS2D]` 段：`NPRESTRESS=<n> BASIS=2d_global`、`NOTE=`、
  `P<id> SIGMA_CON= SL1=..SL6= SPE= RATIO=`，字段 / 顺序 / 单位（Pa，
  RATIO 无量纲）/ 12 位有效数字与 3D `stage prestress loss` 段同源。
  **2D 没有阶段维度**：不设 `STAGE` 占位、不凭空造阶段号；只收录
  `LOSSREF` 束，只给 `LOSS=` 的束不输出一堆 0 冒充「各项损失为 0」；
  没有引用时写 `ERROR=no_lossref_prestress`，不摆空表。
- `example_mixed.exe model.dat --loss-detail` 把该段写进 `<stem>.loss.txt`；
  `--loss-detail-only` 在求解前 return，**不求解、不写 `.result`**，供前处理
  在 `example_stage.exe` 主路径之外单独取明细。两个开关都不改 `.result` /
  控制台输出，`tests/regression/compare_regression.py` 的
  `verify_prestress_loss.dat` 默认用例不带开关、黄金基线不刷新；默认输出
  **逐位不变**。函数只读 `m.prestressLossDetails`，不重算力学量。
- 口径：σpe 后张 = σcon−σl1−σl2−σl4、先张 = σcon−σl2；σl3/σl5/σl6 是
  后期损失，不进 σpe 也不进等效荷载 `Pe = P·(σpe/σcon)`；σl2 只按直线筋，
  σl6 的 εcs / φ 需按附录 C 查表输入。Python 侧只解析、不现算、不补零
  （见 BridgeStudio `prestress_loss_2d.py`）。
- 验证：`cmake --build build --target example_mixed` rc=0；`ctest` 7/7
  （含 `dat_byte_regression` 37/37）；默认 `.result` 与黄金基线规范化后
  逐位一致；opt-in 段与 3D 官方示例数值逐位一致
  （`SL1=87790888.5885 … RATIO=0.820782595059`）；`--loss-detail-only`
  不写 `.result` 且与 `--loss-detail` 的明细逐位一致。

### 新增：2D 模态明细段（`--modal-detail`，opt-in）

- `modal_output_2d.hpp::formatModalSection2D(model, modal)` 把内核**已经算出**的
  `ModalResult2D.modes` 与 `assembleMassMatrix2D()` 重建的同一份 M 写成
  `[MODAL2D]` 段：`NMODES/NDOF/BASIS=xy`、`NOTE`、`MASS X/Y`、
  `F <i> <f_Hz> <T_s> <mx> <my> <mnorm>`、`S <i> <node> <ux> <uy> <theta>`。
  公式 `Γ_d = φᵀMr_d`、`eff_d = Γ_d²/(φᵀMφ)`、`mtot_d = r_dᵀMr_d` 与
  `modal_output_3d.hpp::formatModalSection3D` 逐项同构；**不求解新特征值**，
  不重算 K/M，不改 `.dat` 解析，不动静力 / 活载数值路径。
- 方向口径按用户拍板 `BASIS=xy`：2D 每节点 3 dof（1=u 水平、2=v 竖向、
  3=θz 转角），刚体平动只有 X/Y 两向；2D 无出平面自由度，**不设 Z 向
  参与质量（不存在，不是 0，不写 mz、不补零）**；θz 只作振型分量展示，
  不设参与质量。周期 `T=1/f` 由内核写出，BridgeStudio 侧只解析、不得现算。
- `example_modal.exe model.dat --modal-detail` 才在 `<stem>.modal.txt` 末尾
  追加该段；**不带开关时文件与控制台输出逐字节不变**，
  `tests/regression/baselines/verify_modal.modal.txt` 黄金基线不刷新。
  段头以 `[` 开头、置于既有两个 `==` 段之后，既有 2D 解析器零污染。
- 不静默：`nm==0` / 缺振型向量 / 振型维度与 `NDOF` 不符，写 `ERROR=`
  而不摆空表或补零。
- 验证：`cmake --build build --target example_modal` rc=0；
  `ctest --output-on-failure` 7/7（含 `dat_byte_regression`）；
  opt-in 段冒烟实测 6 阶 / 11 节点 / 66 行振型、`T=1/f`、`mnorm≈1`、
  Σmx/MASS X≈94.9%、Σmy/MASS Y≈81.1%（6 阶截断判据，如实写出）。

## 2026-11-10

### 实现：3D 阶段预应力分项损失（JTG 3362-2018 §6.2）

- **2D / 3D 同源**：把 `*PRESTRESS_LOSS` 的键名映射抽成
  `prestress_loss.hpp::parsePrestressLossLine()`，`fem_2d.cpp` 与 `fem_3d.hpp`
  共用；公式仍是 `computePrestressLosses()` 唯一来源。同一组参数在两个内核
  必得同一个 `σpe/σcon`，不再由两份重复的 `KEY=VALUE` 映射各写一套。
- `fem_3d.hpp`：新增 `Model3D::prestressLosses` 表、`Prestress3D::lossRef /
  lossDetail / hasLossDetail`，以及 `resolvePrestressLosses3D()`（整份文件读完
  后把 `σpe/σcon` 写回 `loss`，幂等）；`*STAGE` 的 `PRESTRESS=` 关键字模式
  新增 `LOSSREF=<编号>`，与 `LOSS=` 互斥，位置模式不支持。引用不存在的编号、
  `σpe/σcon` 不在 `(0,1]`、或没解析就直接求解，全部明确抛错，不静默退回
  `loss=1.0`。等效荷载与结果段仍只认 `loss` 一个量。
- `example_3d_grillage.cpp`：新增 `stage prestress loss` 段（σl1~σl6、σpe、
  RATIO，12 位有效数字，写完立即恢复原精度），只在 LOSSREF 路径出现；
  `stage prestress force` 段及其余段逐位不变。只给 `LOSS=` 的束不产生空段，
  不拿一堆 0 冒充「各项损失为 0」。段头需登记进 `fem_eigen_3d.py::sections`，
  否则 STAGE/P 行会被上一个 section 的正则收走。
- `test_3d.cpp` 新增 `[27]` **27 项**断言：六项损失与独立闭式解逐项对拍、
  表写在 `*STAGE` 前后两种顺序、默认 `LOSS=` 路径逐位不变、四类硬报错、
  未 resolve 直接求解报错 + 解析后 `Pe=P·ratio`。`test_3d` 88 → 115/115。
- 口径：`σpe` 后张 = `σcon−σl1−σl2−σl4`、先张 = `σcon−σl2`；
  **σl3/σl5/σl6 是后期损失，不进 σpe，也不进 `Pe = P·(σpe/σcon)`**。
  σl2 仍只按直线筋，εcs/φ 需查表输入，不含筋-混凝土协同与次内力。
- 验证：全量编译 rc=0、`ctest` 7/7（含 `dat_byte_regression` 37/37；
  基线 `.dat` 均不含 `*PRESTRESS_LOSS`，无需刷新）。

### 修复：`EQUIV` 预应力支持梁组两端轴向自由（逐跨张拉）

- `fem_2d.cpp::processPrestress` 的 `EQUIV` 原先只允许「梁组恰好一端 `dof1` 固定」，其余情形一律抛错；连续梁**中间跨**两端 `dof1` 都自由，逐跨阶段张拉因此无法建模。现在该分支施加一对自平衡轴压力 `+Pe`（起点）/ `-Pe`（终点）：合力为零，只有本梁组受压，外围结构与桥台不承担附加轴力，正是「该阶段只张拉本跨」的等效荷载口径。两端都固定仍报错（`TENDON` 才能模拟弹性缩短）。
- `test_2d.cpp` 新增 `[20c]` 5 项 oracle：等效力对 `+Pe/-Pe`、合力为零、本跨 `N=-Pe`、外围跨 `N=0`、桥台轴向反力为 0，全部精确命中；`test_2d` 124 → 129/129。
- 边界：`EQUIV` 只模拟轴向缩短与端部偏心/曲线束等效荷载；`LOSS` 是总损失系数，不含摩阻/锚固/松弛/收缩徐变分项损失，也不做筋-混凝土协同。
- 验证：全量编译 rc=0、`ctest` 7/7（含 `dat_byte_regression` 37 基线；既有 `*PRESTRESS` 算例均不触发新分支，逐位不变）。文档 `input_format.md` §6.3 / 施工阶段段、`USER_MANUAL.md` §7 同步。

### 修复：2D 端部释放 + NONLINEAR（方案 C），以及凝聚散射自由度错位

- `fem_2d.cpp::solveModel2DNonlinear` 去掉「端部释放暂不支持 NONLINEAR」的显式抛错。BEAM/RIGID 带 `*RELEASE` 时，在局部坐标系先算 `Kl + Kgg(N,L)`，再用 `releaseBeam2D()` 静态凝聚，按 `Tc`（nf×6）投影回**完整 6 个全局自由度**：`Fint = Tcᵀ(rel.Kl·Tc·u)`、`Fext = Tcᵀrel.Fq`、`Kt = Tcᵀrel.Kl·Tc`。释放端内力强制为 0，切线与内力共轭；几何刚度用的轴力仍按凝聚前物理轴力估计（一阶近似）。**无释放单元的装配代码原样保留，默认路径数值逐位不变**。
- 顺带修同源缺陷：`condenseBeamGlobal()`（LINEAR 路径）原先只散射 nf 个「局部自由自由度」编号，而 `Ke` 是 6×6。释放 j 端 θ 时两者编号恰好重合，但释放 i 端自由度会错位（斜梁尤其明显）。现改为按 `elementGdofs2D()` 的完整 6 个自由度散射；被释放行列为零，既有水平梁结果不变。37 个 `dat_byte_regression` 基线 `.dat` 均不含 `*RELEASE`，基线无需刷新。
- `test_2d.cpp` 新增 `[20b]` 16 项断言：水平/斜 30° 梁的旋转不变性、固定端 `|M1|=qL²/8`、铰端 `M2=0`、i 端释放 `M1=0`、LINEAR 与方案 C 一致性、Newton 收敛。验证：`test_2d` 124/124、`ctest` 7/7（含 `dat_byte_regression`）、全量编译 rc=0。
- 文档：`input_format.md` §7.14、`USER_MANUAL.md` §6.1 同步去掉「暂不支持」的过期说明。

### 实现：`[OPT]` 体系守卫（`SYSGUARD`，默认开启，OPTIMIZE=0 逐位不变）

承接计算书工程校核表的两处不通过：`MINENERGY + OPTBALANCE=0` 数学上收敛
（`NITER=196`、成桥态弯矩降 26.5 %），但把边跨索压到 `TMIN` 后索承担
桥面恒载 99.31 % → 89.04 %、同跨相邻索力比 1.08 → 31.98。根因是优化
目标只有弯矩、约束只有索力箱，采纳守卫只查「线性化目标不降反升」。
本节把三条体系判据接进**写回 T\* 之前**的采纳守卫。

- `cable_opt_3d.hpp` 新增 `CableOptSystemCheck3D` 与
  `evaluateCableOptSystem3D()`：按与 `deck_share_pct` /
  `_same_span_smoothness` 完全相同的口径算
  `100·Σ(T·ez)/DECKLOAD`（> 90 %、≤ 100 % 不上拔）与按支座断面
  （`node<2000` 且约束 dof=2/3）分跨的相邻索力比（≤ 1.5）。
  支座约束在解析时已展开为逐 dof，无需另造分跨规则。
- `cable_solver.cpp` 新增 `SYSGUARD/SYSDECKMIN/SYSDECKMAX/SYSADJMAX`；
  守卫不合格写 `STATUS=system_check_rejected`、`APPLIED=0`、**不写回
  T\***、`REJECT=<原因码>`，结果仍按原索力求解。成桥态求解后用实测
  `T` 回填 `SYSCHECK_FINAL`（`ORIGINAL_TENSIONS` /
  `OPTIMIZED_TENSIONS` 区分口径）。`SYSGUARD=0` 只记录不拦截；
  `DECKLOAD` 缺失或支座断面不足写 `NOT_AVAILABLE`，不假判。
- 数学守卫 `objective_not_improved` 优先于体系守卫，但 `REJECT` 行
  始终如实列出体系判据原因码，避免 `PASS=0 / REJECT=none` 的自相矛盾。
- 实测：默认守卫拦下原示例（原索力 99.31 % / 1.08，工程校核表两条
  体系判据回到通过）；`SYSGUARD=0` 诊断模式仍如实标出两条；
  `MINMAXM + MINT=0.30` 三条全过（97.92 % / 748 kN / 1.39）正常放行；
  `MINT=0.40` 触发 `system_support_uplift`。内核全量编译 rc=0、
  `ctest` 7/7、`dat_byte_regression` 37/37 未变。
- 前端 `verify_cable_opt_3d.py` 离线 72 / live 119；新增 `[OPT]`
  体系守卫表与 `SYSCHECK` 解析、原因码中文说明、8 类故障注入。

## 2026-11-09

### 实现：带约束成桥目标索力优化 `[OPT]`（默认 OPTIMIZE=0，数值口径不变）

承接「阶段二 D」。此前 `METHOD=RIGID|MINENERGY` 只是「按固定准则算一组
索力」，`cable_opt.py` 只改索面积；本节补上真正的带约束优化。

- 新增 `cable_opt_3d.hpp`：`optimizeCableTensions3D()` 复用内核裸结构
  线性刚度、`solveMulti` 柔度与 `beamEndForcesAt` 弯矩影响矩阵，解带箱约束
  + 逐塔 `ΣT·cos=0` 软约束的二次型；`MINMAXM`（IRLS 逼近 min max|M|）与
  `MINENERGY` 两种目标，`START=GIVEN|RIGID` 两种起点。优化结果 T\* 写回
  `*CABLE`，随后仍走既有 TENSION / COMPLETION 命中路径，不另造闭环。
- `cable_solver.cpp` 新增 `OPTIMIZE/OPTOBJ/OPTBALANCE/TMAX` 关键字与 `[OPT]`
  段输出。默认全关，既有算例逐位不变（`dat_byte_regression` 37/37 通过）。
- **修量纲而非调参**：弯矩项 ~(N·m)²、平衡项 ~N²、变量 ~1e6 N，早期把
  平衡罚乘弯矩尺度导致罚项大几百倍，实测「收敛了但弯矩反而变大」。
  改为 `x=T/Tscale` 全量无量纲化 + 相对配平后回归正常。
- **守卫不静默**：线性化目标不降反升时写 `STATUS=objective_not_improved`、
  `APPLIED=0`，不写回 T\*、保留原索力；未收敛写 `not_converged`。
  `MPEAK_KIND=LINEAR_BARE` 与 `MPEAK_FINAL_KIND=GIRDEREND_FINAL` 分列两种
  弯矩口径，禁止互相相减冒充改善量。
- 48 索真实斜拉桥实测（相对 `OPTIMIZE=0` 基线成桥态 `max|My|` 39.64 MN·m）：
  `MINENERGY+OPTBALANCE=0` 收敛到 29.12 MN·m（降 26.5%），
  `MINMAXM+OPTBALANCE=0` 到 21.67 MN·m（降 45.3%）；
  `MINMAXM+OPTBALANCE=1` 因塔顶平衡与切比雪夫目标冲突被守卫拦下，
  如实标注未采纳。塔顶平衡是软约束、与弯矩目标存在权衡，不是「越严越好」。
- 段格式与错误码见 `input_format_3d.md` 的 `[OPT]` 专节，实现口径见
  `USER_MANUAL.md §9.9`；前端唯一解析器 `BridgeStudio/cable_opt_3d.py`，
  专项回归 `verify_cable_opt_3d.py` 82/82。

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
