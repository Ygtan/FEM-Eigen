# ANSYS 18.0 语义探测脚本归档（A~G 阶段）

这里保存的是打通「BridgeStudio 前处理 → FEM-Eigen → ANSYS 18.0 对标」过程中
**逐条实测出来的单元 / 荷载 / 求解语义**所用到的一次性探测脚本。它们是结论的
原始证据，不是产品代码：

- **不进 CTest、不参与任何黄金参考值生成**，可独立运行，也可只当文档查结论；
- 脚本内容按本机实测时的原样归档（保真），因此**仍含本机绝对路径常量**；
- 探测结论已固化到各 benchmark 的 `README.md` 与同源导出模块里，日常回归
  完全离线，不依赖这些脚本，也不依赖 ANSYS 许可。

归档位置选在 `nonlinear_ansys_3d/probes/` 是因为 G 阶段首次把这些结论产品化；
但脚本覆盖 A~G 全部阶段，`cable_arch_3d`、`continuous_beam_3d`、
`element_loads_3d` 以及 BridgeStudio 的多个导出模块都引用了它们。

## 1. 脚本清单与实测结论

### C 阶段：3D 梁单元选型 / 局部坐标

| 脚本 | 探测内容 | 结论 |
|---|---|---|
| `beam188_probe.py` | `SECTYPE/SECDATA` 参数顺序、方向节点 K、SMISC 1..30 映射、`MP,GXY` 与扭转刚度、`SECCONTROL` 剪切开关 | BEAM188 是 Timoshenko 一阶剪切梁，SECCONTROL 无法干净关掉剪切 |
| `beam188_probe2.py` | 网格细化收敛 + 剪切开关 + `PRESOL` 内力标签 | 单单元挠度 0.859×EB，细分后收敛到 **1.0678×EB = EB + PL/(GA)** |
| `beam188_probe3.py` | `KEYOPT(3)` 形函数阶次、剪切面积量级、SMISC 1..19 标定 | 找不到与 EB 严格一致的稳定配置 |
| `beam188_probe4.py` | 放大 `MP,GXY` 让 `G·As→inf` 看能否退化为 EB | 不能干净退化 |
| `beam188_probe5.py` | 正式算例配置筛选 + 加扭矩验证 `SECCONTROL` 不影响 `G·J` | 决策：**BEAM4 作严格黄金基准，BEAM188 只作主交付与差异量化、不判 FAIL** |
| `beam4_probe.py` | `BEAM4`（`SHEARZ/SHEARY=0`）与 FEM-Eigen 3D 12×12 EB 逐项对拍 | 位移 `UZ=-3.9009715e-05`、`ROTY=5.8514572e-06`、i 端内力全部一致 |
| `arch_probe.py` | X-Z 平面斜梁 + 竖向 TRUSS 吊杆的 ANSYS 默认局部轴 | 斜梁局部轴与 FEM 一致；**纯竖直梁差 90°** 的问题首次暴露 |

### D 阶段：混合单元 / 竖向荷载 / 弹簧

| 脚本 | 探测内容 | 结论 |
|---|---|---|
| `link_probe.py` | `LINK8 / LINK180` vs FEM-Eigen 3D `TRUSS` 最小混合模型 | 轴力按 `EA/L·c⊗c` 装配，逐位一致（初次全零是漏写 `F,2,FZ` 所致） |
| `spring_probe.py` | `COMBIN14` 的 `KEYOPT(2)=0/1/2/3` 哪个给出 2D 竖向接地弹簧 | `KEYOPT(2)=2`，端点 `v=-0.0067567568 m` 与闭式解一致 |
| `sfbeam_probe.py` | `SFBEAM LKEY=1` 符号 + `ACEL` 自重对应全局 −Z 荷载 | `PRES,+q,+q` 产生向下均布荷载；`MP,DENS + ACEL,0,0,9.81` ↔ `*GRAVITY 0 0 -9.81` |
| `sfbeam_probe2.py` | `SFBEAM LKEY=2`（局部 y / 全局 Y）横向荷载符号 | ANSYS 横向 SFBEAM 语义不可靠，GY / LZ 改由内核闭式解验证 |

### F 阶段：真实竖直塔柱

| 脚本 | 探测内容 | 结论 |
|---|---|---|
| `vertical_beam_probe.py` | IY≠IZ 竖直悬臂梁同时受 FX/FY，比较 ANSYS 默认 / `K(+Y)` / `K(+X)` / `K(+Y)` 且约束 K | `default` 与 `kX` 的 UX/UY 与 FEM 对调（IY/IZ 互换）；`kY` / `kY_fixed` = `MATCH_FEM`（`UX=4.04530744e-03`、`UY=1.61812298e-02`）；K 节点无需约束，0 warning / 0 error |

### G 阶段：非线性成桥态 vs ANSYS NLGEOM

| 脚本 | 探测内容 | 结论 |
|---|---|---|
| `link180_nl_probe.py` | `R,id,A,ISTRN`、`KEYOPT(3)`、NLGEOM 下的应变度量 | `ISTRN` 对 LINK180 **静默无效**；`KEYOPT(3)=1` = tension-only（压缩时轴力归零）、`2` = compression-only；`KEYOPT(2)=0` Almansi、`1` 对数应变，与 FEM 共转列式 Lagrangian 工程应变一阶一致、差 `O(eps²)` |
| `link180_inistrain_probe.py` | 零位移下产生 `N0=EA·0.01=2.0e6 N` 的可行机制（INISTATE / BFE TEMP / BF TEMP，NLGEOM OFF/ON） | 必须 `INISTATE,SET,DTYP,EPEL` + `INISTATE,DEFINE,<e>,,EPEL,X,eps`；四种组合均给出 `N=2.0e6` / `RF_FX=2.0e6` |
| `trial.py` | 小规模 `MODE=COMPLETION METHOD=DECK` 成桥态可行性 | 50 节点 / 63 单元 / 16 索，4 轮收敛，`closure_err=2.006e-10`、`tension_err=3.530e-09`，输出 `[ERECTED]/[CABLE]/[U]/[REACTION]` 可直接使用 |
| `diag_reaction.py` | CB3D10 中 ANSYS 支座反力与外荷载不闭合的来源 | 逐项剥离 `INISTATE`、`KEYOPT(3)`、方向节点 K，确认荷载丢失源于 `*CLOAD` 累加 vs ANSYS `F` 覆盖，而非初应变或方向节点 |

## 2. 跨设备运行前必须调整的本机常量

脚本按实测当时的本机环境原样保留，因此含以下硬编码路径。换机器时按注释修改
（或直接设同名环境变量后自行替换常量），**不需要改动任何探测逻辑**：

| 常量 | 出现位置 | 本机值 |
|---|---|---|
| ANSYS 可执行文件 | 多数脚本的 `ANSYS` / `ANSYS_EXE` | `C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe` |
| BridgeStudio 源码目录 | `arch_probe.py`、`link_probe.py`、`vertical_beam_probe.py`、`trial.py` 的 `sys.path.insert` | `D:\AI-work\BridgeStudio` |
| FEM-Eigen 求解器 | `arch_probe.py`、`link_probe.py`、`vertical_beam_probe.py` 的 `FEM_EXE` | `D:\AI-work\FEM-Eigen\build\example_3d_report.exe` |
| CB3D10 算例目录 | `diag_reaction.py` 的 `CASE` | `D:\AI-work\BridgeStudio\cases\CB3D10_completion_vs_ansys_nl` |
| MinGW 运行库 | 含 `MSYS_BIN` 的脚本（FEM-Eigen exe 需要 `libwinpthread-1.dll`） | `C:\msys64\mingw64\bin` |

运行方式（需要本机 ANSYS 许可；`trial.py` 只需 BridgeStudio + `cable_solver.exe`；
Python 需带 `pydantic`，本机为 `C:\Users\user\anaconda3\python.exe`）：

```powershell
cd FEM-Eigen\tests\benchmark\nonlinear_ansys_3d\probes
python vertical_beam_probe.py
python link180_nl_probe.py
python link180_inistrain_probe.py
python trial.py
```

脚本会把 `.inp / .out / .csv / .dat / .report.txt` 以及 ANSYS 二进制产物写到
**脚本所在目录**（`HERE`）或 `diag_*` 子目录。仓库 `.gitignore` 已把本目录下
除 `*.py` / `*.md` 外的所有文件排除，运行不会污染 `git status`；若想彻底隔离，
可先把脚本复制到临时目录再运行。

## 3. 与 benchmark 的关系

| 结论被用在 | 位置 |
|---|---|
| BEAM4 作黄金基准 | `continuous_beam_3d/README.md`、`BridgeStudio/ansys_bridge_3d.py` |
| SFBEAM / ACEL / COMBIN14 语义 | `continuous_beam/README.md`、`element_loads_3d/README.md`、`BridgeStudio/ansys_bridge_3d_loads.py` |
| LINK180 线性轴向杆 | `cable_arch_3d/README.md`、`BridgeStudio/ansys_bridge_cable.py` |
| 方向节点 K | `cable_arch_3d/README.md`、`BridgeStudio/ansys_bridge_cable.py`、`ansys_bridge_nl.py` |
| INISTATE / tension-only / `*CLOAD` 合并 | 上级 `nonlinear_ansys_3d/README.md`、`BridgeStudio/ansys_bridge_nl.py` |

ANSYS 18.0 对 `LINK180 KEYOPT(3)=1` 会给出一条 `is an undocumented option` 警告，
这是本机实测到的版本相关行为；升级 ANSYS 后应重跑 `link180_nl_probe.py` 复核，
再决定 `nonlinear_ansys_3d` 的 `ansys_log_expect` 是否调整。
