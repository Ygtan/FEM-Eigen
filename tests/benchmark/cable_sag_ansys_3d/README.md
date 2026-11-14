# `cable_sag_ansys_3d` — 索垂度 / Ernst 等效弹模 + 带索自重空缆态 vs ANSYS 18.0

J 阶段。**CB3D14** 是单根 CABLE 单元本构的端到端对标，**CB3D15** 是
`W>0` 悬索桥的空缆态验证。此前 A~I 阶段所有 ANSYS 对标算例都用
`cable_dens=0 → W=0` 主动关闭垂度，因此「Ernst 等效弹模 / 索垂度」虽在
内核中已实现，却从未被外部软件独立验证；空缆 stage 0 也因无横向荷载
退化成直线态。本 benchmark 把这两个缺口一起补上。

```
BridgeStudio/ansys_bridge_sag.py::SAG_CFG   (唯一数据源)
  -> cb3d14_cable_sag.sag   -> cable_sag_probe.exe -> *_fem.csv
       fem_3d.hpp::cableTensionFromChord
       fem_3d_nl.hpp::cableErnstModulus
  -> sag_01..05 / rod_00 / sag_half .inp -> ANSYS 18.0 -> *_disp.csv / *_link.csv
  -> reference/cb3d14_cable_sag_ansys.json   573 项黄金值

BridgeStudio/ansys_bridge_sag.py::BRIDGE_CFG
  -> cb3d15_cable_selfweight.dat -> cable_solver.exe -> *_cable.txt
       [ERECTION] stage 0  (空缆态)
```

## 为什么必须用「分段索链」

FEM-Eigen 把**一根** CABLE 建成浅垂度弹性悬链线（`fem_3d.hpp:160-222`）：

```
Ln = L0·(1 + T/EA) − w²·Lh³/(24·T²)          s(T) = w²Lh³/(24T²)
dT/dLn = 1 / (L0/EA + w²Lh³/(12T³)) = Eeq·A/L0
Eeq = E / (1 + w²·Lh³·E·A/(12·T³·L0))         ← Ernst 等效弹模
```

而 ANSYS 18.0 的 LINK180 是**无垂度轴向杆**，单根单元不可能自行产生
Ernst 软化；直接给它 `Eeq` 又只是把 FEM 公式抄进 ANSYS，属循环论证。
唯一能独立验证的口径是把一根索离散成 N 段 LINK180、按**水平均布**把
自重加到节点上，让 `NLGEOM,ON` 自己收敛出悬链线——此时索链的等效轴向
柔度天然包含垂度效应。这与 I 阶段 `ansys_bridge_erection.stage_apdl()`
用的是同一套技术，但对象是**单根本构单元**而非全桥主缆。

## 算例参数（CB3D14）

| 项 | 值 |
|---|---|
| 弹性模量 / 截面 | E = 195 GPa，A = 0.006 m²，EA = 1.17e9 N |
| 索重度 | W = 462.051 N/m = 7850 kg/m³ × A × 9.81 |
| 水平跨径 / 分段 | Lh = 100 m（两端等高，弦长 Ln = Lh），NSEG = 80 |
| 基准目标张力 | T = 1.0 MN → L0 = 99.923490771 m，f = 0.577564 m，f/L = 0.0058 |
| 初始几何 | 直线 + 0.20 m 正弦扰动（两端不动），**不从答案出发** |
| 初应变 | 每段 `eps = L_init/L0_seg − 1`，`INISTATE,SET,DTYP,EPEL` 逐段注入 |
| 荷载 | 内部节点 `−W·dx`，两端 `−W·dx/2`（水平均布，对应抛物线） |
| ANSYS | `LINK180`, `KEYOPT(2)=0`（Almansi）、`KEYOPT(3)=0`，`NLGEOM,ON`，`NSUBST,100,500,20 / NEQIT,200` |

7 个 ANSYS job：`sag_01..05` 为 L0 扫描（base ±0.002 / ±0.004 m），
`rod_00` 为 `W=0` 直杆基线，`sag_half` 为 `W/2` 低重度工况。

## 实测结果

| 工况 | W (N/m) | L0 (m) | Eeq/E | T_FEM (N) | H_ANSYS (N) | 相对差 |
|---|---|---|---|---|---|---|
| sag_01 | 462.051 | 99.919491 | 0.843427 | 1039172.830 | 1038952.066 | 2.124e-4 |
| sag_02 | 462.051 | 99.921491 | 0.835705 | 1019493.327 | 1019268.963 | 2.201e-4 |
| sag_03 | 462.051 | 99.923491 | **0.827600** | 1000000.000 | 999772.054 | 2.279e-4 |
| sag_04 | 462.051 | 99.925491 | 0.819102 | 980701.933 | 980470.412 | 2.361e-4 |
| sag_05 | 462.051 | 99.927491 | 0.810200 | 961608.453 | 961373.392 | 2.444e-4 |
| rod_00 | **0** | 99.923491 | **1.000000** | 895843.386 | 895843.375 | **1.268e-8** |
| sag_half | 231.026 | 99.916825 | 0.950496 | 1000000.000 | 999950.664 | 4.934e-5 |

- **Ernst 切线**：用 `sag_02/sag_04` 的 ANSYS 水平力做中心差分
  `dT/dL0 = −9.699638e6 N/m`，内核解析值 `−(1+T/EA)·Eeq·A/L0 = −9.698616e6 N/m`，
  相对差 **1.053e-4**。这是本 benchmark 最敏感的一项。
- **垂度**：6 个 `W>0` 工况 `f` 相对差 3.824e-5 ~ 1.966e-4。
- **W=0 直杆基线**相对差 **1.268e-8**（几乎逐位一致），`Eeq/E` 精确为 1。
- 7 个 job 全部 **0 error / 0 warning**；`--live-ansys` 复核 7 job /
  **2261 个 CSV 值**，max abs diff **1.0e-10**。

## 已知差异的定量解释（不是错误）

内核闭式是**浅垂度一阶近似**（用单一 T 代替沿索变化的张力），ANSYS
分段索链是更精细的离散解。差异随 `(V/H)²` 缩放，实测：

| V/H | (V/H)² | dT/T | df/f |
|---|---|---|---|
| 0.0116 | 1.33e-4 | 4.934e-5 | 3.824e-5 |
| 0.0231 | 5.34e-4 | 2.279e-4 | 1.837e-4 |
| 0.0462 | 2.13e-3 | 7.061e-4 | 5.292e-4 |

分段数 N=20/40/80/120 时误差收敛到 ~2.2e-4 平台；把 EA 提高 10 倍使
`eps=T/EA` 降到 8.5e-5 后偏差仍为 1.7e-4 —— 说明**不是** Almansi /
Lagrangian 应变度量差，也不是离散误差，而是浅垂度近似本身的二阶项。
`W=0` 时两套求解器逐位一致（1.27e-8）是这一结论最直接的证据。

Ernst 修正是 **17.2% ~ 67.5%** 的一阶效应，与该二阶偏差量级分离约
3 个数量级，因此统一取 `rtol=2e-3`（与 H/I 阶段相同），逐项
`|got−ref| <= 1e-9 + 2e-3·max(|FEM|,|ANSYS|)`。实测
`573/573 checks passed, worst d/tol = 0.113`。

## CB3D15：带索自重的空缆态

与 CB3D13 **完全相同**的自锚式悬索桥参数（主跨 120 m、自动平衡边跨
71.414 m、矢跨比 1/10、`q_deck=20 kN/m`、主缆 30000 mm²），只把
`cable_dens` 从 0 改成 7850（`W>0`）。

| 工况 | 空缆 stage 0 H | 空缆矢高 | 跨中高程 | 成桥态 H |
|---|---|---|---|---|
| `W=0`（CB3D13 口径） | **0.1 kN**（命中二分下界 `lo=1e2`） | 0 m | 35 m（直线） | 3000 kN |
| `W>0`（CB3D15） | **354.451 kN** | **11.732 m** | 5.268 m | 3346.54 kN（**+11.55%**） |

索自重提供了横向荷载，空缆态自然成为真实悬链线，不再退化到
`cable_solver.cpp:881` 的 100 N 数值下界。因此 I 阶段「空缆不参与对标」
的边界在 `W>0` 时被解除；本 benchmark 用 FEM 重跑冻结 `.dat` 的
6644 个数值逐位复现（`worst rel=0`）并断言 `H>1e3 N && sag>1e-6 m`。

> 注意：CB3D15 是 **FEM 侧自洽回归**，没有 ANSYS 全桥对标。ANSYS
> LINK180 无垂度模型，全桥主缆若要复现需把每一段再细分并施加节点
> 自重，属于另一个量级的算例，未纳入本阶段。

## 归档文件

```
cases/
  cb3d14_cable_sag.sag / .model.json / _fem.csv
  cb3d14_cable_sag_compare_report.txt / _manifest.txt / _run_log.txt
  sag_01..05, rod_00, sag_half : .inp / _ansys.out / _disp.csv / _link.csv
  cb3d15_cable_selfweight.dat / .model.json / _cable.txt / _report.txt
reference/
  cb3d14_cable_sag_ansys.json       573 项黄金值
build_reference_sag.py             从冻结 CSV 生成 / --check 校验
compare_cable_sag_ansys_3d.py       CTest 入口（离线 / --live-ansys）
README.md
```

573 = 7×H + 6×f + 7×80 段索力 = 7 + 6 + 560。每项 `ref` 来自冻结
ANSYS CSV；H / f / T 都是派生量，派生口径写在 JSON 的 `derived_rules`，
由回归脚本从冻结 CSV 重算校验，防止黄金值漂移。

## 运行

```powershell
# 1) 重生成算例（需要 ANSYS 许可 + BridgeStudio）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_sag.py `
  --sync-to D:\AI-work\FEM-Eigen\tests\benchmark\cable_sag_ansys_3d\cases

# 2) 生成 / 校验黄金值
C:\Users\user\anaconda3\python.exe `
  D:\AI-work\FEM-Eigen\tests\benchmark\cable_sag_ansys_3d\build_reference_sag.py --check

# 3) CTest 离线回归（只需 cable_sag_probe.exe + cable_solver.exe）
ctest --test-dir D:\AI-work\FEM-Eigen\build -R cable_sag_ansys_3d_benchmark --output-on-failure

# 4) 实时 ANSYS 复核 7 个 job
C:\Users\user\anaconda3\python.exe `
  D:\AI-work\FEM-Eigen\tests\benchmark\cable_sag_ansys_3d\compare_cable_sag_ansys_3d.py `
  --bin D:\AI-work\FEM-Eigen\build --live-ansys
```

## 实测验证结果

| 项目 | 结果 |
|---|---|
| 离线回归 | **573/573 checks passed**，`worst d/tol=0.113` |
| FEM 内核 CLI 复现 | 70 个值，`worst rel=0` |
| 冻结 CSV 重算 | 573/573 与黄金值一致 |
| Ernst 切线 | `rel=1.053e-4` PASS |
| CB3D15 FEM 复现 | 6644 个值 `worst rel=0`，空缆态判据 PASS |
| `--live-ansys` | 7 job / 2261 CSV 值，max abs diff `1.0e-10`，0 error / 0 warning |
| 既有链路 | 16 项 CTest 全部未改动，A~I 冻结 `.dat` 逐字节不变 |

## 已知边界

- **大转动梁柱效应未纳入**：CR 内核（`fem_3d_nl.hpp` 共转 + 解析几何刚度）
  已具备该能力，但需要 BEAM188 + NLGEOM 的独立算例，属 K 阶段。
- ANSYS 侧 `KEYOPT(2)=0` Almansi 应变；换用 `KEYOPT(2)=1` 实测可把
  `W>0` 的系统偏差从 2.28e-4 降到 1.63e-4，但不改变浅垂度近似这一主因。
- CB3D15 只做 FEM 自洽回归，无 ANSYS 全桥对标（见上文说明）。
- 索链用 `KEYOPT(3)=0`（各段始终受拉）；低张力/松弛工况需另行处理。
