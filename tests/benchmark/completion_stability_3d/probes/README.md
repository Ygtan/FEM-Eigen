# 诊断探针（`completion_stability_3d/probes/`）

只入库可复现的 `.py` / `.md`，不进 CTest、不产出黄金值。

| 文件 | 用途 | 运行方式 |
|---|---|---|
| `make_cases.py` | **冻结 cases 的唯一生成器**：从同一份几何/材料常量写出 6 个 `.dat`，调 `cable_solver.exe` 产出 `cases/*_cable.txt` / `*_stability.txt`。内含独立闭式 `N(v)` / `P(v)` / `F0` / `T*` 与成桥态起点自洽检查。正常回归**不需要**重跑（会覆写冻结输出） | `python make_cases.py <cable_solver.exe> <out_dir>` |
| `two_step_vshape_probe.py` | ② 阶段：V 形双索上验证 ANSYS 两段求解（`step1` 成桥态静力 + `step2` `ARCLEN`），用反力反推总荷载、与闭式 `λ(v)` 对拍，定标 `λ_FEM = TIME-1` 语义与 LINK180 应变度量差 | `python two_step_vshape_probe.py`（需 ANSYS 18.0） |
| `two_step_fem_probe.py` | ② 阶段：在 CB3D10 完整双塔三跨斜拉桥上试 `STABASIS=COMPLETION`，扫描 `STABTOL` / `DS0` / `PSI` 参数组 | `python two_step_fem_probe.py` |
| `two_step_static_probe.py` | ② 阶段：把 CB3D10 改成 `MODE=STATIC` + 冻结 `[CABLE]` 段显式 `L0`，验证起点是否闭合 | `python two_step_static_probe.py` |
| `two_step_bridge_probe.py` | ② 阶段：完整斜拉桥上的 ANSYS 两段求解方向/口径定标（荷载控制 vs `ARCLEN`、重施全部 `F`、加大目标荷载到 `4.6×PINC`），并做 FEM↔ANSYS 逐点对拍 | `python two_step_bridge_probe.py`（需 ANSYS 18.0） |

## 已固化的可行性结论（早期探针不入库）

1. **参考构形必须带预拉**。最早一版把低点放在 `(0,0)`、支座放在
   `(±2, 0.5)`，并用 `L0 = L(V0)/(1+T*/EA)`：参考弦长
   `hypot(2,0.5)=2.0616 < L0=2.0768`，`CABLE` tension-only 在 `u=0`
   松弛，`solveModel3DNonlinear` 的零刚度诊断直接报
   `节点 2 的第 2 个自由度未约束且未被任何单元连接`。
   改成 `H=0.6`、`L0=0.995·hypot(B,H)` 后参考态
   `N_ref=1.055e6 N > 0`，闭环正常收敛。
2. **`MODE=STATIC` 与 `MODE=TENSION` 都与闭式吻合**：
   `v_fem = 0.1`（误差 0）与 `0.09999999941`（误差 `5.9e-10`，
   索长迭代命中 `T*` 的残差）。
3. **纯索无梁体系 `MODE=COMPLETION` 不收敛**：制造线形几何闭合需要梁，
   实测 `converged=0`、`closure_err=1.7e-2`。这一事实被用作
   `vshape_unconverged` 负例（验证 `ERROR=completion_not_converged`），
   **不是**能力声明；真实斜拉 / 悬索桥的 `COMPLETION` 闭环见
   CB3D10~CB3D13。
4. **不归档的临时脚本**：早期 `probe.py` / `probe2.py` 含本机绝对路径
   `D:\AI-work\FEM-Eigen\build\cable_solver.exe` 与两版参数的中间试算，
   口径已被 `make_cases.py` 统一取代，按项目约定不入库。

## 与 ANSYS 的分工

本目录**没有** ANSYS `.inp` / `.out` / CSV：① 阶段只做 FEM 侧自洽与
独立闭式锚点。真实桥梁算例的两段对标属阶段 ②，届时另建 benchmark
目录，不回填本项的冻结黄金值。

## ② 阶段物理定标结论（2026-11-15，探针已归档、尚未写产品代码）

1. **ANSYS 两段求解可行**：`step1` = 参考构形 + `INISTATE` 初应变 +
   恒载 + `NLGEOM` 静力（`TIME=1`）；`step2` = 在 `step1` 荷载基础上
   加增量荷载，`TIME=2`。V 形双索上 `ARCLEN` 的自动步长确实生效
   （`TIME` 增量非均匀），且 `F_total = F0 + (TIME-1)·PINC`，即
   **`λ_FEM = TIME − 1`**。完整斜拉桥上 `0 error`。
2. **应变度量**：FEM `CABLE` 用工程应变 `N=EA/L0·(L−L0)`，
   ANSYS LINK180 `KEYOPT(2)=0` 是 Almansi、`=1` 是对数应变。
   V 形双索（成桥态 ε≈2%）实测 `KEYOPT(2)=0` 起点位移差
   `2.67e-3 m`、`=1` 差 `1.01e-3 m`；CB3D10 冻结算例 ε≈1e-3，
   既有对标相对差 `7.7e-4`。② 的正式算例必须把索应变控制在
   `1e-3` 量级，或按同一份 `L0` 同源生成 `INISTATE` eps。
3. **反力求和不可靠**：V 形双索（单一荷载）用 `RF1+RF3` 反推总荷载
   成立，但完整斜拉桥（50 个 `F`、34 个约束节点）在 `step2` 中
   `ΣRF,FZ` 随加载反而下降（`19.84 MN → 15.33 MN`），荷载控制与
   `ARCLEN` 两个 job 均如此。故 ② 的对拍**不用反力**，改用
   `λ = (TIME−1)·LAMMAX` 与 `ΔUZ = UZ − UZ_step1`。
4. **FEM 路径单调硬化**：`MODE=STATIC` + 显式 `L0` 时
   `RES0=6.3e-3 N`（与 `[NEWTON] RESLAST=7.0e-3 N` 同量级，起点
   严格闭合），`λ` 从 0 单调升到 `4.553`、`CTRL` 到 `-0.05001`，
   **无极限点 / 无下降段**。因此 ANSYS `step2` 用荷载控制 Newton
   即可覆盖全程：目标荷载放大到 `4.6×PINC` 后，27 个子步与 FEM
   弧长路径逐点对拍 `worst|ΔUZ| = 1.76e-5 m`；`ARCLEN` job 与荷载
   控制 job 的终点 `UZ` 仅差 `8.6e-8 m`（可作交叉验证，但其 `TIME`
   不是 λ，不能直接换算）。
5. **发现 ① 阶段遗留缺陷（如实登记）**：`MODE=COMPLETION` 的索长
   迭代会改写 `m.elements[*].L0`，而传给
   `formatCompletionStabilitySection()` 的 `u0` 与最终 `m` 不同源，
   导致 CB3D10 上 `RES0=3.15 MN`（≈`F0_NORM` 量级）、路径第一步
   `λ` 变负。① 的 V 形双索（`MODE=STATIC`/`TENSION`）不受影响。
   ② 的正式算例因此采用 `MODE=STATIC` + 冻结 `L0`；
   「`COMPLETION` 起点同步修复」列为独立后续项，不得靠改文档略过。
6. **同源要求**：CB3D10 冻结 APDL 用 `[ERECTED]` 制造线形作参考
   构形、eps 按 `(L_ERECTED−L0)/L0`，而 `.dat` 用设计线形。
   ② 必须让 FEM 与 ANSYS 用**同一参考构形**（建议设计线形），
   从同一份 `L0` 同时生成 `*EN ... L0=` 与 `INISTATE ... EPEL,X`。
