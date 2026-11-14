# `suspension_erection_3d` — 自锚式自动平衡边跨 + 施工阶段找形 vs ANSYS 18.0

CB3D13 一个算例，三层验证。唯一数据源仍是 `BridgeStudio` 的
`ModelSpec3D`；所有 `.dat` / `.inp` 都由它同源导出。

```
suspension_3d.build_suspension_3d(system="SELF", side_span=0,
                                  q2_deck=5e3, erection=True)
  -> ModelSpec3D.to_dat()  ->  cb3d13_balanced_erection.dat
  -> cable_solver.exe      ->  *_cable.txt
       [ERECTED]/[U]/[CABLE]/[REACTION]        [A] 全桥成桥态
       [ERECTION] / [ERECTION_SHAPE]
       [ERECTION_CHAIN] / [ERECTION_LOAD]
       [ERECTION_FORCE]                       [C] 施工阶段（I 阶段新增）
  -> [A] ansys_bridge_susp + ansys_bridge_nl  -> *_nl.inp -> ANSYS CSV
  -> [B] ansys_bridge_erection.balance_report()  纯 FEM 判据
  -> [C] ansys_bridge_erection.stage_apdl()  -> *_stage1..5.inp -> ANSYS CSV
```

## 算例参数

| 项 | 值 |
|---|---|
| 体系 | `SELF` 自锚式，主缆锚于加劲梁端，H 压进梁 |
| 主跨 / 边跨 | 120 m / **自动平衡 71.414284 m**（边主跨比 0.5951） |
| 矢跨比 / 塔顶高程 | 1/10（矢高 12 m）/ z=17 m |
| 恒载 | `q_deck=20 kN/m`，一期 Q1=15、二期 Q2=5 kN/m |
| 主缆 / 吊索 | 30000 / 5000 mm²，E=195 GPa，`cable_dens=0` |
| 规模 | 60 节点 / 83 单元 / 49 根 CABLE（26 主缆段 + 23 吊索） |
| 成桥水平力 | H = 3000 kN，主跨矢高 12.000 m |
| 施工阶段 | 6 个：空缆、吊装×3、合拢、铺装；**5 个有荷载阶段参与对标** |

## [A] 全桥成桥态对标

与 `suspension_ansys_3d` 同口径：ANSYS 初始几何 = 设计几何
`D = [ERECTED] + [U]`（**不是** `[ERECTED]`），逐索绝对 L0 由
`INISTATE,SET,DTYP,EPEL` + `INISTATE,DEFINE,<e>,,EPEL,X,eps` 注入，
`eps=(L_D-L0)/L0`；`BEAM4 + LINK180(KEYOPT(3)=1) + NLGEOM,ON`。
竖直塔柱 BEAM4 用方向节点 K(+Y) 对齐 FEM 局部轴。

| 组 | 项数 | 最大绝对差 | 量级 | 相对 |
|---|---|---|---|---|
| 位移 | 291 | 1.700e-6 | 2.421e-3 | 7.023e-4 |
| 索力 | 49 | 118.0 N | 3.277e6 | 3.601e-5 |
| 反力 | 8 | 84.5 N | 4.682e6 | 1.805e-5 |

初应变范围 2.106e-4 .. 5.880e-4；ANSYS 0 error / 2 warning
（两个 CABLE 单元类型各一条 `KEYOPT(3)=1 is an undocumented`）。

## [B] 自动平衡边跨

自锚式主缆锚于加劲梁端，梁端支座能承压不能承拉。令主缆在锚点的
竖向分力 `V=0`：

```
V = H·(Δz/L) − qL/2 = 0   =>   L_bal = sqrt(2·H·Δz/q)
H0 = q·L_main²/(8f) = 3000 kN, Δz = 17 m  =>  L_bal = 71.414284 m
```

生成器在 `side_span=0` 时自动取 `L_bal`。量化证据（支座反力）：

| 边跨 | 梁端 Rz | `max|Rz| / (q·L_main)` | 含义 |
|---|---|---|---|
| **自动 71.414 m** | **+0.975 kN** | **4.064e-4** | 无需抗拔支座 |
| 手动 40 m（CB3D12） | −872.621 kN | 0.3636 | 梁端被上拔 873 kN |

塔顶纵向残差 Rx = ±6.934 kN（相对 H=3000 kN 约 2.3e-3），说明塔两侧
主缆水平力基本自平衡。回归判据 `max|Rz|/(q·L_main) <= 1e-3`。

## [C] 施工阶段逐阶段对标

### 内核新增输出段

I 阶段在 `cable_solver.cpp` 里给 `*ERECTION` 增加 4 段**仅在启用施工
阶段时输出**的同源数据，使 ANSYS 侧无需重演内核内部的吊索吊装规则：

| 段 | 内容 |
|---|---|
| `[ERECTION]` | 每阶段 `H / sag / midz / hmin / hmax` |
| `[ERECTION_SHAPE]` | 每阶段主跨主缆 13 个节点高程 z |
| `[ERECTION_CHAIN]` | 主跨主缆节点 `id / x`（各阶段相同） |
| `[ERECTION_LOAD]` | 每阶段吊索搬运到主缆节点的竖向荷载 q（N） |
| `[ERECTION_FORCE]` | 每阶段主跨 12 段主缆轴力 T（N） |

重编译 `cable_solver.exe` 后，对 17 个 A~H 阶段冻结 `.dat` 逐字节对拍：
**数值/结构零差异**（唯一差别是 `elapsed_ms` 运行耗时，benchmark 本来就跳过）。
`build/` 被 `.gitignore` 排除，另一台电脑按 CMake 正常编译即可得到同版内核。

### 口径：给定阶段下料长度 + 阶段荷载，从扰动几何独立求解

FEM 的 `findErectionState` 是「主缆**总**无应力长度守恒、塔顶鞍座可滑」：
二分 H 使三跨总 `L0(H)` 等于成桥下料总长，各段 L0 随 H 重新分配。
ANSYS LINK180 则是「各段下料长度各自固定」。两者在**部分荷载阶段不是
同一个力学问题**。实测从成桥线形直接卸载：

| 阶段 | 铺装 5 | 合拢 4 | 吊装 3 | 吊装 2 / 1 |
|---|---|---|---|---|
| 线形差 | 3.4 mm | 15 mm | 2.18 m | ANSYS 低张力不收敛 |

所以本 benchmark 采用可严格闭环的口径：对每个阶段，把 FEM 找形给出的
阶段构形 `z_stage` 与阶段轴力 `T_stage` 换算成阶段等效下料长度

```
L0_stage[k] = L_stage[k] / (1 + T_stage[k]/EA)
```

逐段注入 `INISTATE`，施加该阶段 `[ERECTION_LOAD]` 的节点荷载，并把 ANSYS
初始几何加 **0.20 m 正弦扰动**（两端塔顶不动）。ANSYS 必须靠 `NLGEOM,ON`
自己把 LINK180 索链收敛回平衡构形——不是从答案出发。比较收敛后的线形、
逐段索力与水平力 H，即验证「同一组下料长度 + 荷载，两套软件给出同一
平衡态」。

阶段索链用 `KEYOPT(3)=0`（各阶段均受拉，避免 tension-only 反复开关）；
`KEYOPT(2)=0` Almansi 应变；`NSUBST,60,300,10 / NEQIT,100`。

| 组 | 项数 | 最大绝对差 | 量级 | 相对 |
|---|---|---|---|---|
| 主缆线形 z | 65 | 8.606e-5 m | 17.0 | 5.062e-6 |
| 主缆段轴力 T | 60 | 21.45 N | 3.195e6 | 6.714e-6 |
| 水平力 H | 5 | 21.47 N | 3.0e6 | 7.156e-6 |

5 个阶段全部 0 error / 0 warning。空缆 stage 0 不参与对标：内核把 H
二分到数值下界 100 N、矢高恒 0，是退化直线态；无横向荷载的 LINK180
索链存在多解/不稳定，报告与 JSON 里如实列出、不隐藏。

### 阶段序列（FEM `[ERECTION]`）

| stage | 名称 | 荷载合计 | H (kN) | 矢高 (m) | 跨中高程 (m) |
|---|---|---|---|---|---|
| 0 | 空缆 | 0 | 0.1 | 0 | 17.0 |
| 1 | 吊装 | 150 kN | 281.5 | 15.986 | 1.014 |
| 2 | 吊装 | 750 kN | 1196.0 | 15.050 | 1.950 |
| 3 | 吊装 | 1350 kN | 1777.7 | 14.344 | 2.656 |
| 4 | 合拢 | 1650 kN | 2259.5 | 11.950 | 5.050 |
| 5 | 铺装 | 2200 kN | 3000.0 | 12.000 | 5.000 |

铺装态精确回到成桥矢高 12 m，是施工阶段找形闭环的关键自检。

## 归档文件

```
cases/
  cb3d13_balanced_erection.dat / .model.json / _cable.txt
  cb3d13_balanced_erection_nl.inp                 [A] 全桥 APDL
  cb3d13_balanced_erection_nl_ansys.out
  cb3d13_balanced_erection_nl_ansys_{disp,reac,link}.csv
  cb3d13_balanced_erection_nl_compare_report.txt  [A]
  cb3d13_balanced_erection_balance_report.txt    [B]
  cb3d13_balanced_erection_stage{1..5}.inp        [C] 阶段 APDL
  cb3d13_balanced_erection_stage{1..5}_ansys.out
  cb3d13_balanced_erection_stage{1..5}_{disp,link}.csv
  cb3d13_balanced_erection_erection_compare_report.txt  [C]
  cb3d13_balanced_erection_manifest.txt / _run_log.txt
reference/
  cb3d13_balanced_erection_ansys.json   478 项黄金值
build_reference_erection.py            从冻结 CSV 生成 / --check 校验
compare_suspension_erection_3d.py      CTest 入口（离线 / --live-ansys）
README.md
```

478 = [A] 291 + 49 + 8 = 348，[C] 65 + 60 + 5 = 130。
每项 `ref` 来自冻结 ANSYS CSV；[C] 的线形与 H 是从 CSV 派生量
（`z_init + UZ`、`mean(T·dx/L)`），派生口径写在 JSON 的 `derived` 字段，
由回归脚本从冻结 CSV 重算校验，防止黄金值漂移。

## 运行

```powershell
# 1) 重生成算例（需要 ANSYS 许可 + BridgeStudio）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_erection.py

# 2) 生成 / 校验黄金值
C:\Users\user\anaconda3\python.exe `
  D:\AI-work\FEM-Eigen\tests\benchmark\suspension_erection_3d\build_reference_erection.py --check

# 3) CTest 离线回归（只需 cable_solver.exe，不需要 ANSYS）
ctest --test-dir D:\AI-work\FEM-Eigen\build -R suspension_erection_3d_benchmark --output-on-failure

# 4) 实时 ANSYS 复核 6 个 job
C:\Users\user\anaconda3\python.exe `
  D:\AI-work\FEM-Eigen\tests\benchmark\suspension_erection_3d\compare_suspension_erection_3d.py `
  --bin D:\AI-work\FEM-Eigen\build --live-ansys

# 5) 全量回归
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1
```

## 实测验证结果

| 项目 | 结果 |
|---|---|
| 离线回归 | **478/478 checks passed**，`worst d/tol=0.351` |
| FEM 重跑复现 | 632 个值（含全部 `[ERECTION_*]` 段），`worst rel=0` |
| [B] 平衡判据 | 4.064e-4 <= 1e-3 PASS |
| `--live-ansys` | 6 个 job，1024 个 CSV 数值，max abs diff 1.589e-7；0 error |
| `build_reference_erection.py --check` | PASS |
| 既有链路 | 新旧 `cable_solver.exe` 对 17 个冻结 `.dat` 数值结构零差异 |

## 已知边界

- **口径差异不是错误**：FEM 施工阶段「总长守恒、鞍座可滑」vs ANSYS
  「各段下料长度固定」。本 benchmark 对标的是给定阶段下料长度 + 荷载后
  两套软件是否给出同一平衡态，**不是** ANSYS 独立预测 FEM 的鞍座滑移口径。
- 空缆态、带索垂度/Ernst、大转动梁柱效应、`*STAGE` 施工阶段（与本
  `*ERECTION` 是两套机制）均未验证。
- `cable_dens=0` 关闭垂度；ANSYS LINK180 无垂度模型。
- FEM 共转列式索用 Lagrangian 工程应变，ANSYS NLGEOM 用 Almansi，
  一阶一致、差异 `O(eps²)`；本算例 eps 控制在 1e-3 量级。
- ANSYS 18.0 的 `LINK180 KEYOPT(3)=1` undocumented 警告是版本行为，
  升级后需重跑 `probes/link180_nl_probe.py` 并复核 `ansys_log_expect`。
- `build/CTestTestfile.cmake` 由 CMake 生成、被 `.gitignore` 排除，不入库。
  本机旧 `build/` 缓存的生成器是 `MinGW Makefiles` 却配了 MSYS 的
  `C:/msys64/usr/bin/make.exe`（两者不兼容，试编译阶段路径被反斜杠截断），
  因此历史上只能手工同步新增的 `add_test` 行。**干净重配已实测可行**，
  换机或重建时按 `tests/benchmark/README.md` 的「从零复现」用 Ninja 或
  `MSYS Makefiles` 生成器即可自动生成同样 16 条注册，无需手工同步。
- `--live-ansys` 在临时工作目录重跑 6 个 ANSYS job，把新生成的 CSV 与
  `cases/` 下的冻结 CSV 逐值比对（判据 `1e-6 + 1e-3|y|`），**不覆写冻结
  文件**；`build_reference_erection.py --check` 随时可确认黄金值与冻结 CSV
  仍逐项一致。
