# cable_arch_3d — 3D 斜拉桥 / 拱桥混合单元 benchmark

本类目检验 **「同一个 `ModelSpec3D` 里同时存在 BEAM 与轴向杆单元」** 的
建模 → 求解 → 商业软件对标链路，以及 `cable_solver.exe` 非线性成桥态链路。

与既有三个 benchmark 正交（互不污染黄金参考值）：

| 类目 | 内核 / 单元 | 参考值来源 |
|---|---|---|
| `ansys` / `bridge` | 2D `example_mixed.exe` | 闭式解 / 离线冻结结果 |
| `continuous_beam` | 2D BEAM3，5 个工况 | 真实 ANSYS 18.0 CSV |
| `continuous_beam_3d` | 3D **纯 BEAM**（BEAM4 / BEAM188） | 真实 ANSYS 18.0 CSV |
| **`cable_arch_3d`** | 3D **BEAM + CABLE/TRUSS**（BEAM4 + LINK180） | 真实 ANSYS 18.0 CSV |
| `cable_arch_3d/nonlinear_smoke` | `cable_solver.exe` / `fem_3d_nl.hpp` | 冻结的 `*_cable.txt` |

## 目录结构

```
cases/
  cb3d03_cable_stayed.{dat,inp,model.json,report.txt}
  cb3d03_cable_stayed_ansys_{disp,reac,beam,link}.csv   # 冻结的 ANSYS 实测
  cb3d03_cable_stayed_ansys.out / _compare_report.txt / _manifest.txt
  cb3d04_arch_bridge.*                                  # 同上
  cb3d09_real_tower.*                                   # 真实竖直塔柱算例
  cb3d09_real_tower_ansys_{disp,reac,beam,link}.csv
reference/
  cb3d03_cable_stayed_ansys.json   162 个检查项
  cb3d04_arch_bridge_ansys.json    299 个检查项
  cb3d09_real_tower_ansys.json     234 个检查项
nonlinear_smoke/
  cb3d05_completion.{dat,model.json}
  cb3d05_completion_cable.txt      # cable_solver 原始输出（冻结）
  cb3d05_completion_{smoke_report,manifest}.txt
build_reference_cable.py           # 从冻结 CSV 生成 reference JSON
compare_cable_arch_3d.py           # CTest 入口：线性对标
compare_cable_solver_smoke.py      # CTest 入口：非线性成桥态
README.md
```

## 算例

| 编号 | 模型 | 关键结果（FEM-Eigen 3D 与 ANSYS 一致到 ~1e-8） |
|---|---|---|
| CB3D03 | 独塔双跨斜拉桥，L=60 m，12 个 BEAM 主梁单元；塔顶 z=18 m 刚性锚点；6 根 CABLE 斜拉索（A=4000 mm²）；q=20 kN/m 节点化 | 主梁 UZ 峰值 −5.8339 mm；索力 94.5~213.6 kN（峰值在 elem 15）；塔顶竖向反力 857.2 kN；ΣFZ 平衡残差 0 |
| CB3D04 | 下承式抛物线拱桥，L=60 m / 矢高 12 m；12 主梁 BEAM + 12 拱肋 BEAM + 11 竖向吊杆 TRUSS；q=20 kN/m 节点化 | 主梁 UZ 峰值 −8.7997 mm；拱肋轴力峰值 840.0 kN；吊杆力 45.6~98.9 kN（关于跨中对称）；端反力 600 / 600 kN |
| CB3D05 | 非线性成桥态 smoke：主跨 200 m + 两侧各 80 m，塔高 60 m，4×12=48 根索，`MODE=COMPLETION METHOD=DECK` | 6 轮收敛；`closure_err=1.4218e-9 m`；`tension_err=5.2875e-9`；索长闭合 max \|dL\|=3.89e-10 m；索力 682.5~1402.8 kN |
| CB3D09 | **真实竖直塔柱斜拉桥**（F 阶段）：L=60 m，12 主梁 BEAM + 6 段竖直塔柱 BEAM + 6 根 CABLE；塔底固结、塔梁分离 | 20 节点 / 24 单元；位移、18 根梁 i 端内力、6 根索轴力共 234 项与 ANSYS 逐项一致；方向节点 K 解决纯竖直 BEAM 局部坐标分叉 |

**参考值不手工填数**：`build_reference_cable.py` 逐项从冻结 ANSYS CSV 抽取，
每项带 `source = 文件名:行键:列号`；`--check` 可随时验证 JSON 与 CSV 脱节与否。
检查项覆盖全部节点 6 自由度位移、全部 BEAM 单元 i 端 6 分量内力（幅值）、
全部轴向杆轴力（幅值）。容差沿用项目惯例 `rtol=1e-3`，位移 `atol=1e-9`、
力 `atol=1e-6`。支座反力不进离线比对（FEM-Eigen 3D `.report.txt` 无反力段），
由 `--live-ansys` 逐行复核冻结 CSV 并检查支座 ΣFZ。

## 实测确认的语义（本机 ANSYS v180）

- **FEM-Eigen 线性 3D 内核**（`fem_3d.hpp`）按 `BEAM / BEAM7 / 其它` 分派，
  `CABLE` 与 `TRUSS` **都走 `truss3DStiffness` = `(EA/L)·c⊗c`**，不区分只拉不压、
  不读 `L0`、无垂度。ANSYS **LINK180** 同为 3D 轴向杆（可拉可压），口径严格对应。
  实测最小混合模型：位移、梁端力、杆轴力全部逐位一致。
- 只受拉 / 松弛退出 / 弹性悬链线 / Ernst 等效弹模 / 索力迭代 / 成桥态几何闭合
  属于 **`fem_3d_nl.hpp` + `cable_solver.exe`** 非线性链，即 CB3D05，不与 ANSYS 对标。
- **BEAM4** 与 FEM-Eigen 3D 同为 12×12 Euler-Bernoulli 空间梁
  （`R` + `RMORE`，`SHEARZ/SHEARY=0`）；BEAM188 是 Timoshenko，差异已在
  `continuous_beam_3d` 单独量化，本类目不重复。
- `LINK180` 轴力 = **SMISC 1**；`R,<id>,A,ISTRN` 即可，无需 `SECTYPE`。
  `ETABLE` 必须用 `ESEL,S,TYPE` 分组，否则对杆单元请求梁的 SMISC 分量会产生 warning。
- **局部坐标陷阱**：`beam3DTransformation` 用 `ez_ref=global Z`（单元轴平行 Z 时
  改用 global Y）、`ey = ez_ref × ex`。ANSYS BEAM4 默认规则与之一致，故 **X-Z 平面
  斜梁逐位一致**（拱肋折线实测通过）；但**纯竖直 BEAM 两边局部 y 轴差 90°**，
  IY/IZ 互换，会导致结果不一致。CB3D03/04 的竖向构件用轴向杆规避；
  **F 阶段已在 `ansys_bridge_cable.py::_vertical_beam_orientations()` 解决**：
  对纯竖直 BEAM 自动在 I 点全局 +Y 生成方向节点 K，用 `E,I,J,K` 对齐 FEM
  局部轴，K 节点不需要约束（`../nonlinear_ansys_3d/probes/vertical_beam_probe.py` 实测 0 warning / 0 error）。
  该逻辑只在检测到竖直梁时插入，CB3D03/04 的 APDL 逐字节不变；CB3D09 是首个
  用到它的算例。
- `example_3d_report.exe` 对「仅由轴向杆连接的节点」只输出 3 个平动位移；
  APDL 导出时对这些节点的转动分量补 0，保证两张表可逐项比对。
- ANSYS `*GET,NODE,RF,FZ` 在**无约束但有外荷载**的节点返回节点荷载残差，
  力平衡只能统计支座节点。
- `*VWRITE,ARR(1,1),...` 一次写出整个数组，外面不能再套 `*DO`；APDL 必须纯 ASCII。

## 已知边界（不是缺陷，但不声称已验证）

- CB3D03 的**塔按刚性锚点简化**（塔顶节点全约束），不建塔柱梁单元；本轮校核对象是
  「混合 BEAM + 轴向杆」链路，不是塔梁相互作用。F 阶段的 CB3D09 补上了真实竖直
  塔柱梁单元，但仍不建塔梁连接 / 塔梁相互作用（塔梁分离）。
- 线性对标链**不验证**只受拉、索垂度、Ernst 模量、索力优化；这些由 CB3D05 的
  非线性链覆盖，CB3D05 本身只做内核自洽校核。G 阶段的
  **`tests/benchmark/nonlinear_ansys_3d`（CB3D10）** 才首次为该非线性链提供
  真实 ANSYS 18.0 `NLGEOM + LINK180 + INISTATE` 数值基准，两者不是同一条目。
- CB3D05 头部 `balance_pct=-52.79%` 是「索力竖向分量 / 总恒载」的比值偏差，
  差额由桥塔 / 支座承担，属斜拉桥成桥态正常受力分配，**不是不平衡残差**。
  `support_v_sum_kN` 与 `applied_v_kN` 相差 0.265%，源于 `cable_solver.cpp:570-585`
  的 `appliedV` 只累加 `*CLOAD` 显式荷载与 CABLE 自重等效节点力，未计入
  `DECKLOAD` 分配的桥面荷载（与 `support_share_kN=250 kN` 量级吻合）。
  真正的收敛判据是 `closure_err` / `tension_err` / 索长闭合 `|dL|`。
- FEM-Eigen 3D 线性内核不认识 `*ELEMENT_LOAD` / `*GRAVITY` / `*SELF_WEIGHT`
  （静默忽略），故均布荷载节点化为 `*CLOAD`，ANSYS 用同一组 `F,node,FZ`。

## 复现

```powershell
# 1) 重跑两个算例的完整「同源导出 → 双求解器 → 对比报告」链（需 ANSYS 许可）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_cable.py --case all

# 1b) F 阶段：重跑真实竖直塔柱算例 CB3D09（需 ANSYS 许可）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_tower.py

# 2) 重跑非线性成桥态 smoke（只需 cable_solver.exe + msys64，不需要 ANSYS）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cable_solver_smoke.py

# 2b) G 阶段：非线性成桥态 vs ANSYS NLGEOM 对标（需 ANSYS 许可）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_nl.py

# 3) CTest：离线回归，不依赖 ANSYS
ctest --test-dir D:\AI-work\FEM-Eigen\build -R "cable_arch_3d_benchmark|cable_solver_smoke|nonlinear_ansys_3d_benchmark" --output-on-failure
ctest --test-dir D:\AI-work\FEM-Eigen\build -L benchmark --output-on-failure

# 4) 可选：实时 ANSYS 复核归档的 .inp（需许可，约 1 分钟）
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\cable_arch_3d\compare_cable_arch_3d.py `
  --bin D:\AI-work\FEM-Eigen\build `
  --cases D:\AI-work\FEM-Eigen\tests\benchmark\cable_arch_3d\cases `
  --reference D:\AI-work\FEM-Eigen\tests\benchmark\cable_arch_3d\reference --live-ansys

# 5) 校验黄金参考值与冻结 CSV 一致
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\cable_arch_3d\build_reference_cable.py --check

# 6) 全量（33 个 exe + 10 组基准/回归）
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1
```

CTest 里 `cable_arch_3d_benchmark` / `cable_solver_smoke` /
`nonlinear_ansys_3d_benchmark` 登记在
`CMakeLists.txt` 与 `run_tests.ps1`；`build/CTestTestfile.cmake` 已同步
（本机 `CMAKE_MAKE_PROGRAM` 指向损坏的 `mingw32-make.exe/mingw32-make.exe`，
无法直接重新 configure，改用临时目录 + `MSYS Makefiles` 验证语法后同步条目；
环境修复后重新 configure 会自动生成同样内容）。

脚本运行输出统一 ASCII，PowerShell / CTest 日志 / 重定向都不会中文乱码。

## 验证记录（2026-11-13）

| 项目 | 结果 |
|---|---|
| `ctest -L benchmark` | 9/9 Passed（F/G 收尾后的全量 CTest 为 14/14 Passed） |
| `run_tests.ps1` 全量 | `ALL TESTS PASSED` |
| 线性对标离线 | **695/695**（162 + 299 + CB3D09 的 234） |
| `--live-ansys` | CB3D03/04 共 246 / 443 个值与冻结 CSV maxdiff 5.8e-10；0 error / 0 warning |
| `build_reference_cable.py --check` | 三份 JSON 与冻结 CSV 一致 |
| `cable_solver_smoke` | 15/15，重跑与冻结 `_cable.txt` 关键量逐位一致 |
| CB3D09 回归兼容性 | CB3D03 / CB3D04 / CB3D09 的 APDL 逐字节不变（方向节点与 `*CLOAD` 合并只在需要时生效） |
| 既有链路 | 2D 84/84、3D 纯梁 113/113 均照常通过，`.dat/.inp` 未改动 |
| G 阶段对标 | 见 `nonlinear_ansys_3d/README.md`：325/325 离线、616 个 CSV 值实时复核零差异 |
