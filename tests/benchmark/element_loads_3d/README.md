# element_loads_3d — 3D 分布荷载 / 自重 ANSYS 对标基准

E 阶段产物。目标：让 FEM-Eigen 的 **3D 线性内核**真正识别并计算
`*ELEMENT_LOAD` / `*GRAVITY` / `*SELF_WEIGHT`，并用 ANSYS 18.0 的
**BEAM4**（同为 Euler-Bernoulli 空间梁）做严格数值对标。

## 1. 背景：原内核的静默缺陷

改造前 `fem_3d.hpp::readModel3D` 只认
`*ET / *MP / *R / *TYPE / *MAT / *REAL / *EN / *BOUNDARY / *CLOAD`，
对 `*ELEMENT_LOAD`、`*GRAVITY`、`*SELF_WEIGHT` **既不解析也不报错**——
写进 `.dat` 后求解"成功"，但荷载被全部丢弃（最小简支梁算例实测：
旧 exe 输出全零，新 exe 与 Euler-Bernoulli 解析解一致）。

因此 C 阶段只能用节点化 `*CLOAD` 近似均布荷载。本阶段补齐真实分布荷载后，
C/D 阶段的既有算例与基线**均未改动**，新算例走独立目录。

## 2. 数据链路（单一数据源）

```
BridgeStudio: build_box_beam_3d_loads()  →  ModelSpec3D
        ├── model.to_dat()          →  cases/<stem>.dat        (FEM-Eigen 3D)
        └── spec3d_loads_to_apdl()  →  cases/<stem>.inp        (ANSYS BEAM4 APDL)
                                          │
   example_3d_report.exe  →  <stem>.report.txt    ANSYS180.exe -b  →  *_ansys_{disp,reac,elem}.csv
                     └──────────── compare_3d_loads() ────────────┘
                                          ↓
                       build_reference_loads.py  →  reference/<stem>_ansys.json
                                          ↓
                       compare_element_loads_3d.py  (CTest, 默认离线)
```

`.dat` 与 `.inp` 由**同一个 `ModelSpec3D`** 生成，杜绝两边手工建模偏差。
`spec3d_loads_to_apdl()` 复用 C 阶段的 `spec3d_to_apdl(element="BEAM4")`，
再注入 `MP,DENS`、`RMORE ... ADDMAS`、合并后的 `SFBEAM` 与 `ACEL`，
没有复制一份 APDL 生成逻辑。

## 3. 荷载语义对照（本机 ANSYS v180 实测）

| FEM-Eigen `.dat` | ANSYS APDL | 说明 |
|---|---|---|
| `*ELEMENT_LOAD 1-30 GZ -q1 -q2` | `SFBEAM,ELEM,1,PRES,+q1,+q2` | 全局 −Z 均布 / 梯形；`LKEY=1` |
| `*GRAVITY 0 0 -9.81` + `DENS` | `MP,DENS` + `ACEL,0,0,9.81` | ACEL 是**惯性加速度**，与重力反号 |
| `*SELF_WEIGHT ALL MASS 100` | `RMORE,...,ADDMAS` | 每单元附加质量（kg/m） |
| `*CLOAD 16 3 -100000` | `F,16,FZ,-100000` | 节点集中力，DOF3 = UZ |

方向关键字支持 `GZ`（全局 Z）、`GY`（全局 Y）、`LZ`（局部 z）。
`SFBEAM` 的横向面号（LKEY=2/3）在本机 BEAM4 竖向设置下实测不可靠
（只产生固端力、不产生位移），所以 **APDL 对标只覆盖 GZ 与纯 −Z 重力**；
GY / LZ 两个弯曲平面交给 `verify_kernel_loads.py` 的闭式解验证。

注意：`SFBEAM` 与 `ACEL` 同时使用时必须避免自重双计；`ACEL` 零分量已做
正规化（输出 `ACEL,0,0,9.81`，不出现 `-0`）。

## 4. 算例

三例均为 3×30 m 箱梁，31 节点 / 30 单元，支座节点 1 / 11 / 21 / 31。

| 算例 | 工况 | 关键结果 |
|---|---|---|
| CB3D06 `cb3d06_udl` | 真实 `*ELEMENT_LOAD` 均布荷载 q = 20 kN/m | 位移相对差 3.413e-08；端内力相对 6.209e-15；竖向力平衡残差 0 |
| CB3D07 `cb3d07_selfweight` | 叠加：UDL 20 kN/m + `DENS=2500` 与 `*GRAVITY 0 0 -9.81` + `*SELF_WEIGHT ALL MASS 100` + node16 `FZ=-100 kN` | 总外荷载 −3 043 355.5 N，ANSYS 支座反力 3 043 355.5 N；残差 9.313e-10 N（相对 3.060e-16） |
| CB3D08 `cb3d08_trapezoid` | 每单元梯形 q1=20 / q2=40 kN/m（`SFBEAM,1,1,PRES,20000,40000`） | 验证固端力 `L/20(7q1+3q2)` 等 q1/q2 系数；力平衡残差 0 |

三例 ANSYS 日志均为 `0 warning / 0 error`。

## 5. 冻结参考值

`reference/*_ansys.json` 全部由 `build_reference_loads.py` 从归档的真实 ANSYS
CSV 抽取，**没有任何手工录入的数值**。每个检查项带
`source = <csv 文件>:<行键>:<列号>`，`--check` 可检测 JSON 与 CSV 是否脱节。

每例 366 项检查：

- 位移 `31 节点 × 6 DOF = 186` 项（`UX/UY/UZ/ROTX/ROTY/ROTZ`，严格带符号）
- 梁端内力 `30 单元 × 6 分量 = 180` 项（`N/Vy/Vz/T/My/Mz`，i 端，比幅值以规避
  FEM 与 ANSYS 的剪力 / 弯矩符号约定差异，近零值 `atol=1e-6`）

三例合计 **1098 项**。`support_fz_sum` 另存支座竖向反力合计，供力平衡说明使用。

> 实现备忘：检查项的 `kind` 与 CSV 表是两个命名空间。`CSV_PATTERN` 只描述
> 三张冻结表（同时驱动 `doc["csv"]`、归档行数审计和 `--live-ansys` 文件名），
> 而 `kind="force"` 的检查来自 FEM report、其数值源头是 elem CSV。
> 因此单独设 `CHECK_SOURCE_CSV = dict(CSV_PATTERN, force=CSV_PATTERN["elem"])`
> 仅用于生成 `source` 字段，不污染三键语义。

## 6. 目录结构

```
element_loads_3d/
├── README.md                          本文件
├── build_reference_loads.py           冻结 CSV → reference JSON（--check 校验）
├── compare_element_loads_3d.py        CTest 入口（离线；--live-ansys 可选）
├── verify_kernel_loads.py             内核级闭式解 / 叠加恒等式验证
├── cases/                             30 个归档文件（3 例 × 10）
│   ├── cb3d0{6,7,8}_*.dat/.inp/.model.json/.report.txt
│   ├── cb3d0{6,7,8}_*_ansys_{disp,reac,elem}.csv
│   └── cb3d0{6,7,8}_*_ansys.out / _compare_report.txt / _manifest.txt
└── reference/
    ├── cb3d06_udl_ansys.json          366 checks
    ├── cb3d07_selfweight_ansys.json   366 checks
    └── cb3d08_trapezoid_ansys.json    366 checks
```

每个算例目录同时保留了「前处理输出 `.dat`」与「对应 ANSYS APDL `.inp`」，
符合最初打通链路时的归档约定。

## 7. 运行方法

Python 用 `C:\Users\user\anaconda3\python.exe`（默认 `python` 亦可，脚本只用标准库）。
所有控制台输出为 ASCII，避免 PowerShell / CTest / 重定向的编码差异。

```bash
cd D:/AI-work/FEM-Eigen/tests/benchmark/element_loads_3d

# 重新生成 / 校验冻结参考值
python build_reference_loads.py
python build_reference_loads.py --check

# 离线对标（默认 CTest 路径，不需要 ANSYS 许可）
python compare_element_loads_3d.py \
  --bin D:/AI-work/FEM-Eigen/build --cases cases --reference reference

# 可选：现场重跑归档 APDL，逐格比对冻结 CSV（需要 ANSYS 许可）
python compare_element_loads_3d.py --bin ... --cases ... --reference ... --live-ansys

# 内核级闭式解 / 叠加验证
python verify_kernel_loads.py --bin D:/AI-work/FEM-Eigen/build

# CTest
ctest --test-dir D:/AI-work/FEM-Eigen/build -L benchmark --output-on-failure
powershell -ExecutionPolicy Bypass -File D:/AI-work/FEM-Eigen/run_tests.ps1
```

`--live-ansys` 会重新运行归档 `.inp` 并检查 `warning/error` 计数，再把
`disp / reac / elem` 三张新 CSV 与冻结值逐格比对（判据 `1e-3*scale + 1e-6`）。
实测三例三类 CSV 全部 `maxdiff=0`。

## 8. 内核实现要点

`D:\AI-work\FEM-Eigen\fem_3d.hpp`（备份 `fem_3d.hpp.bak_e`）：

- `Model3D` 新增 `elemLoads`（`ElemLoad3D`: e1/e2/dir/q1/q2）、
  `selfWeights`（`SelfWeight3D`: all/elem/rule/value）、`g`（`Vec3`，默认 0）。
- 新增辅助函数 `upperAscii3`、`parseElemRange3`、`resolveLoadDir3`、
  `distributedLoadEquivalent3D`、`elementEquivalentLoad3D`。
- `assembleModel3D()` 在 `a.F0 = F` 之前逐单元累加一致节点荷载：
  ```cpp
  for (const auto& e : m.elements) {
      const Eigen::VectorXd fe = elementEquivalentLoad3D(m, e);
      if (fe.cwiseAbs().maxCoeff() == 0.0) continue;
      const std::vector<int> g = elementGdofs3D(a.dm, e);
      for (int r = 0; r < (int)g.size() && r < fe.size(); ++r)
          if (g[r] >= 0) F(g[r]) += fe(r);
  }
  ```
- `beam3DLocalEndForces()` 扣除固端力：`return fl - T * fe;`，否则端内力
  会漏掉分布荷载的贡献（零荷载时直接返回 `fl`，保证旧路径逐字节不变）。
- **两个弯曲平面的 FEF 弯矩符号相反**：局部 z 弯曲 `w-θy` 与局部 y 弯曲
  `v-θz` 不能共用一套符号，这是本阶段最容易写错的地方。
- `readModel3D()` 末尾显式拒绝「`*STAGE` + 分布荷载 / 自重」组合并抛
  `runtime_error`，杜绝再次出现静默忽略；阶段模型请在 `*STAGE_LOAD` 中按增量表达。

`D:\AI-work\FEM-Eigen\example_3d_report.cpp`（备份 `.bak_e`）：TRUSS 轴力输出
补上 i 端等效荷载分量 `fe_i = fe(0)*c(0)+fe(1)*c(1)+fe(2)*c(2)`，
使自重作用下的杆轴力不再恒为 0。

`build/example_3d_report.exe` 按原 CMake 参数重编译
（`-O3 -DNDEBUG -std=gnu++17 -static-libgcc -static-libstdc++`，含
amgcl / spectra / eigen3 的 `-I`），备份为 `.bak_e`。

## 9. 验证结果

| 项目 | 结果 |
|---|---|
| `build_reference_loads.py --check` | 3 份 JSON 各 366 checks 与冻结 CSV 一致 |
| 离线对标 | `366/366` × 3 = **1098/1098**，0 failure |
| `--live-ansys` | 3 例 × 3 张 CSV，`maxdiff=0`，0 warning / 0 error |
| 内核闭式解 / 叠加 | **67/67 passed**（udl_gz 5、udl_gy 4、udl_lz 5、selfweight_dens/addmass/mass/rho 各 5、trapezoid_superposition 31、truss_selfweight 2） |
| `ctest -L benchmark` | **8/8** 通过（含新增 2 条） |
| `run_tests.ps1` 全量 | `ALL TESTS PASSED` |
| 全量 CMake 构建 | `BUILD_EXIT=0`，生成 36 个 exe |

叠加恒等式的口径是严格线性叠加 `fe(q1,q2) = fe(q1,0) + fe(0,q2)`
（每个单元各自独立梯形，不是整梁连续梯形），容差对应 report.txt 的 8 位有效数字。
TRUSS 自重算例额外约束了自由端水平自由度，避免零刚度自由度。

### 既有基线无回归

2D 连续梁 84/84、3D 纯梁 113/113、索/拱混合单元 461/461、非线性 smoke 15/15、
`.dat` 逐字节回归 37/37 全部照常通过。新旧 `example_3d_report.exe` 对
`verify_3d_beam7.dat` 与 `verify_3d_variable_section.dat` 的输出
经 `diff --strip-trailing-cr` 确认 `BYTE-IDENTICAL`（零荷载路径未变）。

## 10. 已知边界

- 只支持**线性、无阶段**模型；阶段 + 分布荷载组合会显式报错，不是"不支持但不吭声"。
- 黄金基准是 **BEAM4**。BEAM188 是 Timoshenko 一阶剪切梁（单单元剪切锁定、
  细分后收敛到 `EB + PL/(GA)`），不能干净退化为 EB，故只作主交付与差异量化，不判 FAIL。
- ANSYS 侧横向 `SFBEAM` 语义不可靠；GY / LZ 方向由内核闭式解覆盖，不声称 ANSYS 对标。
- **纯竖直 BEAM** 的局部 y 轴默认规则在 FEM-Eigen 与 ANSYS 之间差 90°
  （IY/IZ 互换，D 阶段已实测）；本批算例均为水平梁，不覆盖该情形。
  F 阶段已在 `ansys_bridge_cable.py` 用方向节点 `K(+Y)` 解决，并由 CB3D09
  真实竖直塔柱算例验证，见 `cable_arch_3d/README.md`。
- `--live-ansys` 依赖本机 ANSYS 许可，因此**不进默认 CTest**。

## 11. 相关文件

BridgeStudio 侧：

- `model_spec_3d.py`：新增 `ElemLoad3D` / `SelfWeight3D` 与
  `elem_loads` / `self_weights` / `gravity` 字段；默认空时旧 `.dat` 逐字节不变。
- `ansys_bridge_3d_loads.py`：`build_box_beam_3d_loads`、
  `spec3d_loads_to_apdl`、`total_vertical_load`、`compare_3d_loads`。
- `cases/run_cases_3d_loads.py`：三例一键编排。
- `ansys_bridge_3d.py`：C 阶段的同源导出模块，**本阶段未改动**。
- `cases/CB3D06_distributed_udl/`、`CB3D07_distributed_selfweight/`、
  `CB3D08_trapezoidal_loads/`：BridgeStudio 侧原始交付目录。

FEM-Eigen 侧：

- `fem_3d.hpp`、`example_3d_report.cpp`、`build/example_3d_report.exe`（含 `.bak_e` 备份）。
- `CMakeLists.txt` 注册 `element_loads_3d_benchmark` 与
  `element_loads_3d_kernel_validation`（均 `LABELS benchmark`，`TIMEOUT 600`）。
- `run_tests.ps1` 增加两段；`build/CTestTestfile.cmake` 手动同步
  （本机 `CMakeCache.txt` 的 `CMAKE_MAKE_PROGRAM` 仍损坏，无法直接重新 configure，
  用临时目录 + `MSYS Makefiles` 验证语法后同步；该文件被 git 忽略）。
