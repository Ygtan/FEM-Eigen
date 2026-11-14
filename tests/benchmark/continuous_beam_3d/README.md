# 3D 连续箱梁 benchmark（FEM-Eigen 3D / ANSYS BEAM4 / ANSYS BEAM188）

把 `BridgeStudio` 生成的**两个 3D 连续箱梁算例**接入 FEM-Eigen 回归体系：
每次构建后用当前 `example_3d_report.exe` 重算归档的前处理输出 `.dat`，
与**冻结的 ANSYS 18.0 批处理实测值**比对。默认离线运行，不需要 ANSYS 许可。

本类目与 `tests/benchmark/continuous_beam/`（2D BEAM3）**正交、互不冒充**：
2D 链用 `BridgeSpec` + `example_mixed.exe` + `.result`；
3D 链用 `ModelSpec3D` + `example_3d_report.exe` + `.report.txt`。
不改动任何 2D 主干，也不触碰 2D「默认 `.dat` 逐位不变」的黄金哈希，
以及 `verify_2d_export_import.py` 中「ANSYS 不含 BEAM4/BEAM188」的负例守卫。

---

## 1. 为什么要同时跑 BEAM4 和 BEAM188

这是本轮最重要的实测结论（探测脚本已归档在
`../nonlinear_ansys_3d/probes/beam188_probe*.py`、
`../nonlinear_ansys_3d/probes/beam4_probe.py`，
说明见 `../nonlinear_ansys_3d/probes/README.md`）。

| 内核 / 单元 | 梁理论 | 与自研 3D 内核的关系 |
|---|---|---|
| FEM-Eigen 3D `fem_3d.hpp` `BEAM` | 12×12 **Euler-Bernoulli** 空间梁，含 GJ 扭转，三次 Hermitian 形函数 | 基准 |
| ANSYS **BEAM4** | 6DOF **Euler-Bernoulli**，三次形函数；`SHEARZ/SHEARY=0` 关闭剪切 | **逐位一致** |
| ANSYS **BEAM188** | 一阶剪切变形 **Timoshenko** + **线性**形函数 | 系统性差异，可解释 |

BEAM188 实测（L=10 m 悬臂梁、`E=2.06e11`、`IYY=0.4148`、`A=0.478`）：

- 单单元端部挠度 = **0.859 × EB**（剪切锁定）；
- 细分到 40 单元收敛到 **1.0678 × EB** = `EB + PL/(GA)`（标准 Timoshenko）；
- `KEYOPT(3)=3`（三次形函数）单单元即达 1.0679，说明离散误差可消除，
  但剪切项仍在；
- `SECCONTROL,<大值>,<大值>` 或放大 `MP,GXY` 能抑制剪切，但**只作用于一个
  横向方向**，且端内力与节点平衡值仍差约 1.7%，强行数值退化并不干净；
- 官方明确「BEAM188 每根构件要用多个单元」，与 BEAM4/BEAM44 的三次插值不同。

因此本目录的做法是：

- **BEAM4 作严格黄金基准**（与 FEM-Eigen 3D 同为 EB），`rtol=1e-3`；
- **BEAM188 作主交付与商业软件对照**，结果全部归档，差异在
  `*_compare_report.txt` 的 `[4]` 节量化，**不判 FAIL、不进参考值**。

CB3D01 实测：竖向挠度峰值 BEAM4 `-1.2808525e-3 m`、
BEAM188 `-1.3194055e-3 m`，差 **3.01%**；位移最大绝对差 / 量级 = **3.23%**。
差异来自剪切柔化 + 连续梁超静定反力重分布，属预期物理差异，不是建模错误。

---

## 2. 目录结构

```
continuous_beam_3d/
├── cases/                     冻结的算例（每个 15 个文件）
│   ├── <stem>.dat              FEM-Eigen 3D 前处理输出（ModelSpec3D.to_dat）
│   ├── <stem>.model.json       模型快照
│   ├── <stem>.inp              ANSYS BEAM188 APDL（SECTYPE/SECDATA）
│   ├── <stem>_b4.inp           ANSYS BEAM4 APDL（R/RMORE，EB 严格基准）
│   ├── <stem>.report.txt       FEM-Eigen 3D 结果（6DOF 位移 + i 端 6 分量内力）
│   ├── <stem>_b4_{disp,reac,elem}.csv   BEAM4 机器可读结果（黄金参考来源）
│   ├── <stem>_{disp,reac,elem}.csv     BEAM188 机器可读结果
│   ├── <stem>_b4_ansys.out / <stem>_ansys.out   ANSYS 作业日志
│   └── <stem>_compare_report.txt / <stem>_manifest.txt
├── reference/
│   ├── cb3d01_3span_ansys.json    36 个检查项
│   └── cb3d02_space_ansys.json   77 个检查项
├── build_reference_3d.py          从冻结 CSV 抽取 ref（--check 校验脱节）
├── compare_continuous_beam_3d.py  CTest 入口
└── README.md
```

**参考值不手工填数**：`build_reference_3d.py` 从冻结的 BEAM4 CSV 抽值，
每项记录 `source=文件名:行键:列号`；`--check` 可随时验证 JSON 与 CSV 未脱节。

检查项三类：

| kind | 数据来源 | 用途 |
|---|---|---|
| `disp` | FEM-Eigen 3D `.report.txt` 的 `node/dof` | 6 自由度位移（含支座零位移） |
| `force` | `.report.txt` 单元 i 端 `N/Vy/Vz/T/My/Mz` | 与 BEAM4 SMISC 比**幅值**（符号约定不同） |
| `reaction` | 冻结的 `_b4_reac.csv` | ANSYS 支座反力自校验；FEM-Eigen 3D 无反力输出段 |

容差：`rtol=1e-3`；`force` 项 `atol=1e-6 N·m`，用于吸收理论零分量的
浮点舍入（ANSYS 会给出 ~1e-9 的舍入值）。

---

## 3. 算例

| 编号 | 模型 | 工况 | 关键结果 |
|---|---|---|---|
| CB3D01 | 3×30 m 三跨连续箱梁，单箱单室 `B=8.0 m / h=2.0 m / tw=0.02 / tf=0.025`，31 节点 / 30 单元 | 均布荷载 20 kN/m（节点化） | 反力 240.6 / 659.4 / 659.4 / 240.6 kN，中支座弯矩 −1.782 MN·m；FEM-Eigen 3D 与 BEAM4 位移相对差 3.7e-8、内力 1.0e-15 |
| CB3D02 | 同几何 | 均布 10 kN/m + 中跨跨中 `FY=−100 kN` + `MX=+50 kN·m` | 激活 `UY/ROTX/Vy/Mz/T`；与 BEAM4 位移相对差 2.0e-8、内力 3.3e-15 |

两算例竖向力平衡残差 **0**（`ΣFZ_反力 = −ΣFZ_荷载`），ANSYS 四次运行
**0 warning / 0 error**。

### 截面属性（同一几何，`ansys_bridge_3d.box_section_3d`）

```
A  = B·h − bw·hw                       = 0.478 m²
IY = (B·h³ − bw·hw³)/12                = 0.4148 m⁴  绕 y 轴，竖向(z)弯曲
IZ = (h·B³ − hw·bw³)/12                = 3.3751 m⁴  绕 z 轴，横向(y)弯曲
J  = 4·A0² / ∮(ds/t)  (Bredt)          = 1.1886 m⁴
G  = E/(2(1+ν))                        = 8.5833e10 Pa
```

坐标约定：**X 纵向 / Y 横向 / Z 竖向**；竖向荷载 = DOF 3 (`UZ`) = ANSYS `FZ`；
`IY ↔ ANSYS IYY`、`IZ ↔ IZZ`、`J ↔ Torsion Constant`。

---

## 4. 3D 内核的既有缺口（未改 C++ 源码）

`fem_3d.hpp::readModel3D` / `assembleModel3D` 的实际情况：

- 只认 APDL 风格 `*ET/*MP/*R/*TYPE/*MAT/*REAL/*EN/*BOUNDARY/*CLOAD`；
- **不认识** `*ELEMENT_LOAD` / `*GRAVITY` / `*SELF_WEIGHT`，且**静默忽略**
  （不报错、不警告），写进去会让求解"成功"但荷载全丢；
- `DENS` 只进一致质量矩阵，不产生静力自重荷载；
- 不输出 `.result` 与支座反力段，只有 `.report.txt`（6DOF 位移 + **i 端**端内力）。

因此本批算例把均布荷载**节点化**为 `*CLOAD`（梯形一致节点力：端节点 `q·Le/2`、
内节点 `q·Le`，DOF=3），ANSYS 侧用**同一组** `F,node,FZ` 施加，严格同源。
若后续需要真实分布荷载 / 自重，须在 `readModel3D` 增加关键字解析、在
`assembleModel3D` 增加固端力项，并在 `beam3DLocalEndForces` 中扣除固端力——
这是独立的 C++ 内核改造，不在本轮范围。

另外 3D 需要额外约束横向/扭转自由度消除 6 个刚体模式：
`node1` 锁 DOF 1-4，各中支座锁 DOF 3，末端锁 DOF 2-3。

---

## 5. 运行

```powershell
# CTest（离线，无需 ANSYS 许可）
ctest --test-dir D:\AI-work\FEM-Eigen\build -R continuous_beam_3d_benchmark --output-on-failure
ctest --test-dir D:\AI-work\FEM-Eigen\build -L benchmark --output-on-failure

# 全量回归
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1

# 参考值一致性 / 重生成
python tests\benchmark\continuous_beam_3d\build_reference_3d.py --check
python tests\benchmark\continuous_beam_3d\build_reference_3d.py

# 可选的实时 ANSYS 端到端复核（需许可，约 2 分钟；重跑两份 APDL 与冻结 CSV 对拍）
python tests\benchmark\continuous_beam_3d\compare_continuous_beam_3d.py ^
  --bin D:\AI-work\FEM-Eigen\build ^
  --cases tests\benchmark\continuous_beam_3d\cases ^
  --reference tests\benchmark\continuous_beam_3d\reference ^
  --live-ansys

# 重新建模 + 重新求解（会刷新归档）
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases_3d.py --case all
```

必须用 `C:\Users\user\anaconda3\python.exe`（默认 python 缺 pydantic）；
FEM-Eigen 3D exe 需要 `C:\msys64\mingw64\bin` 在 PATH。
脚本运行输出统一 ASCII，CTest 日志 / PowerShell / 重定向都不会中文乱码。

---

## 6. 验证结果

| 验证 | 结果 |
|---|---|
| `ctest -L benchmark` | 4/4 Passed，新增用例 0.41 s |
| `run_tests.ps1` 全量 | 33 个 exe + 5 组基准/回归，`ALL TESTS PASSED` |
| 离线比对 | **113/113** 通过（CB3D01 36 + CB3D02 77） |
| `build_reference_3d.py --check` | 两份 JSON 与冻结 CSV 一致 |
| `--live-ansys` 实时复核 | 12 组 CSV 与冻结值 **maxdiff = 0**，113/113 通过 |
| 同源校核（`compare_report.txt`） | FEM-Eigen 3D vs BEAM4：位移 3.7e-8 / 2.0e-8，内力 1.0e-15 / 3.3e-15，力平衡残差 0 |
| ANSYS 日志 | 4 次运行 0 warning / 0 error |

`CMakeLists.txt` 已正式登记 `continuous_beam_3d_benchmark`（LABELS benchmark），
`run_tests.ps1` 增加对应段。本机 `build/CMakeCache.txt` 的
`CMAKE_MAKE_PROGRAM` 仍是损坏路径、无法直接 `cmake -S . -B build`，
故用临时目录 + `MSYS Makefiles` 验证过语法后同步到
`build/CTestTestfile.cmake`；环境修好后重新 configure 会自动生成同样条目。

---

## 7. 后续

- **D**：斜拉桥 / 拱桥 3D 链路（可复用 `ModelSpec3D` + `cable_solver`）。
- 补齐 3D 内核的 `*ELEMENT_LOAD` / `*GRAVITY`，让 3D 算例也能用真实分布荷载与自重。
- 若需要与 BEAM188 严格数值对齐，可在 ANSYS 侧用 `KEYOPT(3)=3` + 截面级
  剪切控制做受控对照，并单独建一套 Timoshenko 参考值。
