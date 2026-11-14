# nonlinear_ansys_3d — 非线性成桥态 vs ANSYS NLGEOM + LINK180 对标基准

本类目验证 **FEM-Eigen 几何非线性成桥态**（`fem_3d_nl.hpp` + `cable_solver.exe`）
与 **ANSYS 18.0 大变形静力**（`NLGEOM,ON` + `LINK180` + `INISTATE`）的数值一致性。

这是全项目第一次把「非线性共转列式索 + 索长迭代 + 制造线形 + 成桥目标态」
这条链拿到商业软件上做逐项对标；此前 `cable_arch_3d/nonlinear_smoke`（CB3D05）
只做内核自洽校核，不涉及 ANSYS。

## 1. 与其它 benchmark 的边界

| 类目 | 内核 / 单元 | 参考值来源 |
|---|---|---|
| `cable_arch_3d` | 3D **线性** `example_3d_report.exe`，BEAM + CABLE/TRUSS 按 `EA/L·c⊗c` 装配 | 真实 ANSYS 18.0 BEAM4 + LINK180 CSV |
| `cable_arch_3d/nonlinear_smoke` | `cable_solver.exe` 共转列式非线性 | 冻结的 `*_cable.txt`，**无**商业软件基准 |
| **`nonlinear_ansys_3d`** | `cable_solver.exe` 共转列式非线性，`MODE=COMPLETION METHOD=DECK` | 真实 ANSYS 18.0 **NLGEOM + LINK180 + INISTATE** CSV |

`nonlinear_ansys_3d` 与 `cable_solver_smoke` 是两个独立 CTest 条目：前者对外部
商业软件，后者检查 48 索大模型的收敛 / 几何闭合 / 索长闭合，互不替代。

## 2. 数据链路（单一数据源）

```
BridgeStudio/cable_stayed_3d.build_cable_stayed_3d(mode=COMPLETION, cable_dens=0)
    -> ModelSpec3D.to_dat()                         -> cb3d10_completion.dat
    -> cable_solver.exe                             -> cb3d10_completion_cable.txt
         [ERECTED]   制造线形绝对坐标  ──────────────┐
         [CABLE]     逐索下料长度 L0 ─────────────┐  │
         [U]         相对制造线形的成桥位移       │  │
         [REACTION]  支座反力                     │  │
                                                 │  │
    -> ansys_bridge_nl.spec_nl_to_apdl()  <────────┴──┘
         BEAM4（节点取 [ERECTED] 几何）
         LINK180 KEYOPT(3)=1 tension-only
         INISTATE,SET,DTYP,EPEL
         INISTATE,DEFINE,<elem>,,EPEL,X, eps=(L_ERECTED - L0)/L0
         NLGEOM,ON + NSUBST,20,100,10 + NEQIT,50
                                                   -> cb3d10_completion_nl.inp
    -> ANSYS180.exe -b -j ... -i ... -o ...        -> _disp/_reac/_link.csv
    -> ansys_bridge_nl.compare_nl()                -> compare_report.txt
    -> build_reference_nl.py                        -> reference/cb3d10_completion_ansys.json
    -> compare_nonlinear_ansys_3d.py（CTest 离线重跑 FEM）
```

FEM 与 ANSYS 只共享同一个 `ModelSpec3D`，ANSYS 侧的初始几何与初应变来自
FEM 已收敛的成桥态结果；不存在手工重敲节点坐标或索力。

## 3. 本机 ANSYS v180 实测语义（G 阶段探测结论）

探测脚本已归档在仓库内 `probes/`（A~G 共 16 个，见 `probes/README.md`；
脚本含本机绝对路径常量，跨设备运行前需按该 README 调整）：

| 结论 | 说明 |
|---|---|
| `LINK180` 的 `R,id,A,ISTRN` 第 3 参数**静默无效** | 初应变必须用 `INISTATE,SET,DTYP,EPEL` + `INISTATE,DEFINE,<e>,,EPEL,X,eps`；零位移下实测 `N0=E·A·eps` |
| `KEYOPT(3)=1` = tension-only，`2` = compression-only | 压缩时轴力归零；ANSYS 18.0 会给出一条 `KEYOPT(3) = 1 is an undocumented option` **warning**，属已知基线，不隐藏、不作失败 |
| `KEYOPT(2)=0` 用 Almansi 应变，`1` 用对数应变 | FEM 共转列式索为 Lagrangian 工程应变 `N=EA/L0·(Ln−L0)`；三者一阶一致，差异 `O(eps²)` |
| 竖直 / 近竖直 `BEAM4` 需方向节点 `K(+Y)` | F 阶段 `probes/vertical_beam_probe.py` 实测：`E,I,J,K` 且 K 取 I 点全局 +Y，可对齐 FEM `ey=+X / ez=+Y`；K 节点无需约束 |
| `*CLOAD` 累加 vs ANSYS `F` 覆盖 | 塔柱自重按「每段两端各半」写两遍，同一 `(node,dof)` 会出现多条；必须先合并再输出 `F`。未合并时两塔共丢 4.83 MN，反力总和从 18.03 MN 掉到 13.21 MN |
| `*CFOPEN` + `*VWRITE` | 一次写出整个数组，外层不能套 `*DO`；APDL 命令流保持纯 ASCII |

## 4. 算例 CB3D10

| 项 | 值 |
|---|---|
| 体系 | 双塔三跨斜拉桥，主跨 80 m + 两侧各 30 m，塔高 25 m，每塔每侧 4 索，共 16 索 |
| 离散 | 50 节点 / 63 单元 / 300 自由度 |
| 荷载 | `q_deck=50 kN/m`，`cable_dens=0`（关闭垂度 / Ernst） |
| 求解 | `MODE=COMPLETION METHOD=DECK`，4 轮收敛，`closure_err=2.006e-10 m`，`tension_err=3.530e-09` |
| 索力 | 423.7 ~ 731.1 kN；初应变 5.58e-4 ~ 1.49e-3 |
| 方向节点 | 16 个（塔柱近竖直 BEAM4） |

对标结果（`cb3d10_completion_nl_compare_report.txt`）：

| 分组 | 项数 | 最大绝对差 | 分组量级 | 相对 | 判定 |
|---|---|---|---|---|---|
| 位移（6 DOF） | 300 | 9.466e-06 | 1.224e-02 | 7.731e-04 | PASS |
| 索力 | 16 | 1.877e+02 N | 7.312e+05 N | 2.567e-04 | PASS |
| 反力 | 9 | 2.996e+03 N | 8.972e+06 N | 3.339e-04 | PASS |

## 5. 冻结参考值与容差口径

`build_reference_nl.py` 从冻结 ANSYS CSV 逐项抽取 `ref`，每项带
`source = 文件名:行键:列号`；检查项清单由 FEM `[U] / [CABLE] / [REACTION]`
决定，共 **325 项**（300 位移 + 16 索力 + 9 反力）。

容差与 `ansys_bridge_nl.py::compare_nl()` **完全一致**：

```
对每组取 scale = max(|FEM|, |ANSYS|)
逐项判定  |got - ref| <= atol + rtol * scale      rtol = 2e-3
```

不能退化成 `rtol·|ref|`：转动自由度只有 1e-5 rad 量级，与 ~1e-2 m 的平动
共用一个分组量级；若按逐项相对容差，ANSYS/FEM 的 `O(eps²)` 一阶差异会被
误判为 FAIL。JSON 中每项都显式写入 `scale` 与 `group_scale`，便于审计。

`--check` 可随时验证 JSON 与冻结 CSV 是否脱节。CTest 入口除 325 项对标外还会：

- 校验归档文件清单、冻结 CSV 行数（disp 50 / reac 50 / link 16）、FEM 结果段行数；
- 校验每个 `ref` 仍等于冻结 CSV 对应单元格（防止冻结表被改动而 JSON 未更新）；
- 重跑 `cable_solver.exe`，校验 `converged=1`、`closure_err/tension_err<1e-6`；
- 用 `1e-6` 相对容差把重跑的 325 个 FEM 值与冻结 `_cable.txt` 逐项比对
  （实测 `worst rel = 0`，完全复现）。

## 6. 目录结构

```
cases/
  cb3d10_completion.dat                    COMPLETION 前处理输出
  cb3d10_completion.model.json             模型快照
  cb3d10_completion_cable.txt               cable_solver 原始结果（冻结）
  cb3d10_completion_nl.inp                  ANSYS NLGEOM APDL
  cb3d10_completion_nl_ansys.out           ANSYS 作业日志
  cb3d10_completion_nl_ansys_{disp,reac,link}.csv   冻结的 ANSYS 实测
  cb3d10_completion_nl_compare_report.txt  FEM vs ANSYS 对标报告
  cb3d10_completion_nl_{manifest,run_log}.txt
reference/
  cb3d10_completion_ansys.json              325 个检查项
build_reference_nl.py                      从冻结 CSV 生成 reference JSON
compare_nonlinear_ansys_3d.py               CTest 入口（离线 / --live-ansys）
probes/                                     A~G 全部 16 个 ANSYS 语义探测脚本 + README
README.md
```

`probes/` 是原始探测证据，不进 CTest；其 README 记录了每条结论的实测数字、
脚本清单，以及跨设备运行前需要调整的本机绝对路径常量。

归档 CSV 命名把 ANSYS 现场生成的 `<apdl_stem>_{disp,reac,link}.csv`
改名为 `<apdl_stem>_ansys_{...}.csv`，`--live-ansys` 内部做了名称映射。

## 7. 运行方法

```powershell
# 1) BridgeStudio 完整复现：同源导出 -> cable_solver -> ANSYS -> 对标报告（需 ANSYS 许可）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_nl.py

# 2) 重新生成 / 校验冻结参考值
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\build_reference_nl.py
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\build_reference_nl.py --check

# 3) CTest 离线回归（默认，不需要 ANSYS 许可）
ctest --test-dir D:\AI-work\FEM-Eigen\build -R nonlinear_ansys_3d_benchmark --output-on-failure

# 4) 可选：现场重跑归档 APDL，逐格复核冻结 CSV（需许可）
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\compare_nonlinear_ansys_3d.py `
  --bin D:\AI-work\FEM-Eigen\build `
  --cases D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\cases `
  --reference D:\AI-work\FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\reference --live-ansys

# 5) 全量回归（33 个 exe + 10 组 benchmark / regression）
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1
```

`--live-ansys` 会核对 ANSYS 日志计数：期望 **0 error / 1 warning**，且警告文本
必须含 `KEYOPT(3) = 1 is an undocumented`；任何多余 error / warning 都判失败。

`nonlinear_ansys_3d_benchmark` 已登记在 `CMakeLists.txt`、`run_tests.ps1`
与 `build/CTestTestfile.cmake`（本机 CMake cache 的 `CMAKE_MAKE_PROGRAM`
指向损坏的 `mingw32-make.exe/mingw32-make.exe`，只能手工同步 CTest 条目；
环境修复后重新 configure 会自动生成同样内容）。所有脚本运行输出统一 ASCII。

## 8. 已知边界（不是缺陷，但不声称已验证）

- **索自重 W=0**：ANSYS LINK180 无悬链线 / 垂度模型，FEM 侧关闭垂度与
  Ernst 等效弹模，否则两边口径不同。带垂度的成桥态尚未与商业软件对标。
- **应变度量差异**：FEM 用 Lagrangian 工程应变，ANSYS NLGEOM 用 Almansi
  应变；本算例把索应变控制在 ~1e-3，实测差异 ~8e-4，属 `O(eps²)`。
  大应变 / 大转动工况不在此基准覆盖范围。
- **BEAM4 为小变形线性梁**：主梁成桥位移为 mm 量级，ANSYS 侧不参与大转动
  更新的影响为二阶；未验证几何非线性梁柱效应显著的工况。
- **`KEYOPT(3)=1` 是 ANSYS 18.0 的 undocumented 选项**：实测压缩时轴力归零
  确为 tension-only，但依赖该版本行为，升级 ANSYS 后需重跑探测脚本。
- **ANSYS 只对归档 `.inp` 提供 `--live-ansys` 复核**；默认 CTest 完全离线，
  不依赖许可证，也不重新生成 APDL。
- 对标只覆盖位移 / 索力 / 支座反力；主梁单元内力、动力特性、稳定性不在本
  基准内。`cable_solver` 头部 `balance_pct` 同样是「索力竖向分量 / 总恒载」
  的分配比值，不是力平衡残差，详见 `cable_arch_3d/README.md`。

## 9. 验证记录（2026-11-13）

| 项目 | 结果 |
|---|---|
| `build_reference_nl.py --check` | 325 项与冻结 CSV 一致，分组量级 disp=0.0122444 / force=731245 / reac=8.9715e6 |
| 离线 CTest | `325/325 checks passed`，`worst d/tol=0.387`；FEM 重跑 325 值 `worst rel=0` |
| `--live-ansys` | 0 error / 1 已知 warning；616 个 CSV 数值 `max abs diff = 0` |
| `ctest --output-on-failure` | 14/14 Passed（含新增 `nonlinear_ansys_3d_benchmark` 0.40 s） |
| `run_tests.ps1` 全量 | `ALL TESTS PASSED` |
| 既有基线无回归 | 2D 84/84、3D 纯梁 113/113、cable/arch 695/695、element loads 1098/1098、kernel 67/67、smoke 15/15 均照常通过 |

## 10. 相关文件

| 文件 | 作用 |
|---|---|
| `BridgeStudio/ansys_bridge_nl.py` | `spec_nl_to_apdl()` / `initial_cable_strain()` / `compare_nl()` |
| `BridgeStudio/cases/run_cases_nl.py` | CB3D10 一键复现脚本 |
| `BridgeStudio/ansys_bridge_cable.py` | `spec_mixed_to_apdl()`；F 阶段加入竖直 BEAM 方向节点与 `*CLOAD` 合并 |
| `BridgeStudio/ansys_bridge_tower.py` | F 阶段 CB3D09 真实竖直塔柱斜拉桥构造器 |
| `FEM-Eigen/tests/benchmark/cable_arch_3d/` | CB3D03/04/09 线性混合单元 + CB3D05 非线性 smoke |
| `probes/README.md` 及目录下 16 个 `*_probe.py` / `trial.py` / `diag_reaction.py` | A~G 全部 ANSYS 语义探测脚本与实测结论归档 |
