# 悬索桥找形 vs ANSYS NLGEOM 基准（CB3D11 / CB3D12）

H 阶段：把「BridgeStudio 前处理 → FEM-Eigen 计算 → ANSYS 校核」链路扩展到
**悬索桥**，并给出**地锚式与自锚式各一个算例**。两个算例的几何、材料、
荷载参数完全相同，只改 `*SUSP_BRIDGE SYSTEM`，因此两者差异可以直接归因
于结构体系，而不是参数凑出来的巧合。

---

## 1. 数据链路（单一数据源）

```
suspension_3d.build_suspension_3d(system=EARTH|SELF, **CFG)   -> ModelSpec3D
  |
  +-- to_dat()  ----------------------------> <stem>.dat
  |                                            含 *SUSP_BRIDGE / *MAINCABLE /
  |                                            *HANGER / *CABLE
  +-- cable_solver.exe ---------------------> <stem>_cable.txt
  |      [CABLE]     T / L0 / chord（平衡态弦长）
  |      [MAINCABLE] / [HANGER] 主缆分段与吊索
  |      [U] / [REACTION] / [CAMBER] / [ERECTED]
  |
  +-- ansys_bridge_susp.design_geometry(res) -> 设计几何 D = [ERECTED] + [U]
  +-- ansys_bridge_nl.spec_nl_to_apdl(
  |        coords=D, l0=[CABLE].L0)  ------> <stem>_nl.inp
  |        BEAM4 + LINK180 + NLGEOM + INISTATE
  |        eps = (L_D - L0) / L0
  |
  +-- ANSYS180.exe -b ----------------------> *_nl_ansys_{disp,reac,link}.csv
  +-- compare_susp() -----------------------> *_nl_compare_report.txt
```

所有文件都由同一个 `ModelSpec3D` 派生，不存在两边手工建模的偏差。

---

## 2. 两个算例

| | **CB3D11 地锚式（EARTH）** | **CB3D12 自锚式（SELF）** |
|---|---|---|
| 主缆锚固 | 桥台外侧锚碇（全约束） | 加劲梁端节点 |
| 水平力 H | 锚碇传给地基，梁不承受 | 以轴压形式压进加劲梁 |
| 纵向约束 | 加劲梁端 | **跨中**（对称荷载下支反力恒为零） |
| 节点 / 单元 | 50 / 65 | 48 / 65 |
| 主缆段 / 吊索 | 20 / 17（共 37 根 CABLE） | 20 / 17（共 37 根 CABLE） |
| 位移 DOF | 243 | 237 |
| 反力分量 | 12 | 8 |
| FEM `H` | 3000 kN | 3000 kN |
| 主跨矢高 | 12 m（目标 12 m） | 12 m（目标 12 m） |
| `closure_err` | 6.073e-04 m | 1.841e-03 m |
| `tension_err` | 1.028e-03 | 2.201e-03 |
| 索力范围 | 199.8 ~ 3315.7 kN | 199.6 ~ 3383.3 kN |
| **位移** | 243 项，相对差 **6.145e-04** PASS | 237 项，相对差 **1.013e-03** PASS |
| **索力** | 37 项，相对差 **7.129e-06** PASS | 37 项，相对差 **3.681e-05** PASS |
| **反力** | 12 项，相对差 **4.496e-06** PASS | 8 项，相对差 **1.942e-05** PASS |

共用参数（`BridgeStudio/cases/run_cases_susp.py::CFG`）：主跨 120 m、
边跨 40 m、矢跨比 1/10、吊索间距 10 m、塔高分 4 段、`q_deck=20 kN/m`、
主缆 30000 mm²、吊索 5000 mm²、`cable_dens=0`（W=0）、
`main/tower` 材料按 `suspension_3d.build_suspension_3d` 默认。
生成后清空 `model.lanes`：本批只做静力成桥态对标，不引入车道活载 /
影响线 / 模态，ANSYS 侧没有对应荷载。

---

## 3. 关键口径：ANSYS 初始几何不是 `[ERECTED]`

这是本批最重要的一处，也是与 G 阶段斜拉桥（CB3D10）**不同**的地方。

G 阶段的斜拉桥走 `solveCompletionState3D`：制造线形是闭环迭代的未知量，
L0 相对 `[ERECTED]` 反算，所以 ANSYS 必须从 `[ERECTED]` 出发。

悬索桥走 `cable_solver.cpp` 的 **FormFinding 分支**（源码注释明确）：
主缆节点坐标不是可制造量，可制造的是**无应力索长**；把节点坐标也当
未知量会让方程数少于未知量、迭代漂到索力只剩几十 kN 的松弛构型。
因此该分支是「找形给出设计几何与索力 → 由**变形后弦长**闭式反算 L0 →
解一次几何非线性」。于是：

- `[U]` 是**从设计几何 D 出发**的一次 Newton 位移（毫米级弹性压缩量），
  不是「从制造线形出发的成桥位移」；
- `[CAMBER] = -[U]`、`[ERECTED] = D - [U]` 只是**预拱度建议值**，
  不参与 L0 反算，二者不构成闭环；
- `closure_err` 定义为 `max|[U]|`，是位移量级，**不是迭代残差**；
  `tension_err` 是索力达标偏差。

所以 ANSYS 必须以设计几何 D 为初始构形、使用同一批绝对 L0。
D 由 `[ERECTED] + [U]` 逐节点重建，不依赖 `.dat` 里那条只作初值用的
抛物线，也不依赖内核是否把找形线形写回 `m.nodes`。

**实测对照**（同样参数、同样 L0，只改 ANSYS 初始几何）：

| ANSYS 初始几何 | 位移相对差 | 结论 |
|---|---|---|
| `[ERECTED]`（照搬斜拉桥做法） | 0.248（EARTH）/ 0.385（SELF） | FAIL |
| 设计几何 D = `[ERECTED] + [U]` | 6.1e-04（EARTH）/ 1.0e-03（SELF） | PASS |

自锚式对二阶差异更敏感：H 压进加劲梁，梁的几何刚度差异被放大。
把 `q_deck` 与索面积配到初应变 2e-4 ~ 6e-4 量级后，两个体系都在
统一判据 `rtol=2e-3` 下通过，没有为某个算例单独放宽容差。

---

## 4. ANSYS 语义（v180 实测）

| 项目 | 结论 |
|---|---|
| LINK180 初应变 | `R,id,A,ISTRN` 的 `ISTRN` 静默无效；必须用 `INISTATE,SET,DTYP,EPEL` + `INISTATE,DEFINE,<e>,,EPEL,X,eps` |
| tension-only | `KEYOPT(3)=1`（压缩时轴力归零）；ANSYS 报 `undocumented option` 警告，两个 CABLE 类型各一条，共 2 条，如实保留 |
| 应变度量 | `KEYOPT(2)=0` 用 Almansi 应变；FEM 共转列式索用 Lagrangian 工程应变 `N=EA/L0*(Ln-L0)`，一阶一致、差异 `O(eps²)` |
| 梁单元 | BEAM4（12×12 Euler-Bernoulli）；竖直塔柱用方向节点 `K(+Y)` 对齐 FEM 的 `ey=+X/ez=+Y`，本批每例 8 个方向节点 |
| 求解 | `NLGEOM,ON` + `NSUBST,20,100,10` + `NEQIT,50`；索轴力取 SMISC 1 |
| 荷载 | `*CLOAD` 为累加语义、ANSYS `F` 为覆盖语义，导出时按 `(node,dof)` 合并 |
| 垂度 | `cable_dens=0 → W=0`，关闭垂度 / Ernst；ANSYS LINK180 无垂度模型 |

探测脚本与实测结论归档在 `../nonlinear_ansys_3d/probes/`。

---

## 5. 黄金参考值与容差

`build_reference_susp.py` 逐项从冻结 ANSYS CSV 取值，检查项清单由 FEM
`_cable.txt` 决定，每项带 `source=<文件>:<行键>:<列号>`：

- `cb3d11_susp_earth_ansys.json` — **292** 项（243 位移 + 37 索力 + 12 反力）
- `cb3d12_susp_self_ansys.json` — **282** 项（237 位移 + 37 索力 + 8 反力）

容差与 `compare_susp()` 完全一致：每组取 `scale = max(|FEM|, |ANSYS|)`，
逐项判 `|got-ref| <= atol + 2e-3 × scale`。**不能**退化成逐项
`rtol×|ref|`：转动自由度量级只有 1e-5，和平动共用分组量级，逐项相对
容差会把 ANSYS/FEM 的一阶差异误判成 FAIL。

`compare_suspension_ansys_3d.py` 的收敛检查用 JSON 里的
`closure_tol=1e-2` / `tension_tol=1e-2`（FormFinding 的这两个量不是迭代
残差），但另有 `frozen_reproduction` 对 `[U]` / `[CABLE]` / `[REACTION]`
**逐位对拍**重跑结果与冻结文件，宽松阈值不会掩盖任何漂移。

---

## 6. 目录

```
cases/
  cb3d11_susp_earth.dat / .model.json / _cable.txt / _nl.inp
  cb3d11_susp_earth_nl_ansys.out / _disp.csv / _reac.csv / _link.csv
  cb3d11_susp_earth_nl_compare_report.txt / _manifest.txt / _run_log.txt
  cb3d12_susp_self.*                                    （同上 11 个文件）
reference/
  cb3d11_susp_earth_ansys.json    292 项
  cb3d12_susp_self_ansys.json     282 项
build_reference_susp.py           从冻结 CSV 生成 / --check 校验
compare_suspension_ansys_3d.py    CTest 入口（离线 / --live-ansys）
README.md
```

---

## 7. 复现

```powershell
# 1) 重新生成两个算例（需要 ANSYS 许可 + cable_solver.exe）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_susp.py

# 2) 重建 / 校验黄金参考值
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\suspension_ansys_3d\build_reference_susp.py --check

# 3) CTest 离线回归（不需要 ANSYS）
ctest --test-dir D:\AI-work\FEM-Eigen\build -R suspension_ansys_3d_benchmark --output-on-failure
ctest --test-dir D:\AI-work\FEM-Eigen\build -L benchmark --output-on-failure

# 4) 现场重跑归档 APDL，逐格复核冻结 CSV（需要 ANSYS 许可）
C:\Users\user\anaconda3\python.exe D:\AI-work\FEM-Eigen\tests\benchmark\suspension_ansys_3d\compare_suspension_ansys_3d.py --bin D:\AI-work\FEM-Eigen\build --cases D:\AI-work\FEM-Eigen\tests\benchmark\suspension_ansys_3d\cases --reference D:\AI-work\FEM-Eigen\tests\benchmark\suspension_ansys_3d\reference --live-ansys

# 5) 全量
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1
```

CTest 里 `suspension_ansys_3d_benchmark` 登记在 `CMakeLists.txt`，
`run_tests.ps1` 有一键入口；本机 `build/CTestTestfile.cmake` 因 CMake
cache 的 make 路径损坏而手工同步，另一台机器重新 configure 会自动生成。

---

## 8. 实测回归结果

| 项目 | 结果 |
|---|---|
| `build_reference_susp.py --check` | 两份 JSON 与冻结 CSV 一致 |
| 离线回归 | **574/574**（292 + 282），`worst d/tol` 0.307 / 0.506 |
| FEM 重跑复现 | 574 个值 `worst rel = 0` |
| `--live-ansys` | 两算例 637 + 613 = 1250 个 CSV 数值，max abs diff 1.0e-06 / 3.7e-09；**0 error / 2 已知 warning** |
| 既有链路 | 2D 84/84、3D 纯梁 113/113、混合单元 695/695、分布荷载 1098/1098、CB3D10 325/325 均未改动 |

---

## 9. 已知边界（未声称验证）

- 索垂度 / Ernst 等效模量：`W=0` 刻意关闭，ANSYS LINK180 没有对应模型。
- 大转动下的 BEAM4 梁柱效应：成桥位移为毫米量级，影响为二阶；
  自锚式梁承受 H 轴压，对这一项更敏感（位移差比地锚式大约一个量级）。
- 活载 / 影响线 / 模态 / 施工阶段（`*ERECTION`）：本批算例清空 `lanes`、
  未启用 `erection`，只校核静力成桥态。
- 主缆线形找形本身（`[CABLESHAPE]` 的索多边形解）是 FEM 内核能力，
  ANSYS 侧没有独立复算；本批验证的是**给定找形结果与下料长度后，
  两套软件对同一平衡态的力学响应一致**。
- ANSYS 版本升级后需重跑 `probes/link180_nl_probe.py` 复核
  `KEYOPT(3)=1` 行为，并同步 `ansys_log_expect` 的 warning 计数。
