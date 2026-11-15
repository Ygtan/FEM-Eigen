# 诊断探针（`completion_stability_3d/probes/`）

只入库可复现的 `.py` / `.md`，不进 CTest、不产出黄金值。

| 文件 | 用途 | 运行方式 |
|---|---|---|
| `make_cases.py` | **冻结 cases 的唯一生成器**：从同一份几何/材料常量写出 6 个 `.dat`，调 `cable_solver.exe` 产出 `cases/*_cable.txt` / `*_stability.txt`。内含独立闭式 `N(v)` / `P(v)` / `F0` / `T*` 与成桥态起点自洽检查。正常回归**不需要**重跑（会覆写冻结输出） | `python make_cases.py <cable_solver.exe> <out_dir>` |

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
独立闭式锚点。真实桥梁算例的「成桥态等效初应变 + 两段 `ARCLEN`」
对标属阶段 ②，届时另建 benchmark 目录，不回填本项的冻结黄金值。
