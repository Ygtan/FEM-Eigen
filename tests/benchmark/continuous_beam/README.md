# 连续梁 ANSYS 对标 benchmark

把 `BridgeStudio` 生成的**五个**连续梁算例接入 FEM-Eigen 的回归体系：每次构建后用当前
`example_mixed.exe` 重算前处理输出 `.dat`，与**冻结的 ANSYS 18.0 批处理实测值**比对。
已注册为 CTest 用例 `continuous_beam_benchmark`（LABEL `benchmark`），并挂到
`run_tests.ps1` 全量脚本。CTest 默认离线运行，不需要 ANSYS 许可。

覆盖：两跨/三跨满布均布荷载（CB01/CB02，另有闭式解三方校核），以及自重（CB03）、
节点集中力（CB04）、弹性支承（CB05）三类扩展工况。

## 与 `../ansys` 的区别（重要）

| | `tests/benchmark/ansys` | 本目录 |
|---|---|---|
| 参考值来源 | 欧拉-伯努利**闭式解**（当时仓库没有 ANSYS 数值结果文件，用闭式解代表 BEAM4） | **真实 ANSYS 18.0 批处理导出的 CSV** |
| 模型来源 | 仓库内既有 `.dat` | `BridgeStudio` 的 `BridgeSpec`，`.dat` 与 `.inp` 由**同一个对象**生成 |
| 是否留档 APDL | 无 | `cases/*.inp` 完整留档，可随时重跑 |
| 实时复核 | 无 | `--live-ansys` 可重跑 APDL 做端到端比对 |

这样可防止"程序整体一起漂移、但自回归仍绿"的问题，且算例的前处理输出与 ANSYS 命令流
都留在仓库里，建模—计算—商业软件校核三步都可复现。

## 目录

```
continuous_beam/
├── README.md                       本文件
├── compare_continuous_beam.py      CTest / run_tests.ps1 入口（离线比对 + 可选实时 ANSYS）
├── build_reference.py              从冻结 CSV 生成 reference/*.json
├── cases/                          每个算例 <stem> 一组，共 5 组
│   ├── <stem>.dat                  BridgeStudio 前处理输出（FEM-Eigen new format）
│   ├── <stem>.inp                  同源 ANSYS MAPDL 命令流
│   ├── <stem>.model.json           BridgeSpec 模型快照
│   ├── <stem>.result               归档时 FEM-Eigen 求解结果（参考，不参与判定）
│   ├── <stem>_ansys_disp.csv       ANSYS: node UX UY ROTZ
│   ├── <stem>_ansys_reac.csv       ANSYS: node FX FY
│   ├── <stem>_ansys_elem.csv       ANSYS: elem NI VI MI NJ VJ MJ
│   ├── <stem>_ansys_spring.csv     ANSYS: node SMISC1（仅 CB05 弹性支承）
│   ├── <stem>_ansys.out            ANSYS 批处理日志（许可/求解留痕）
│   └── <stem>_compare_report.txt   首次端到端校核报告
└── reference/
    ├── cb01_2span_ansys.json            18 项
    ├── cb02_3span_ansys.json            22 项
    ├── cb03_selfweight_ansys.json       16 项
    ├── cb04_point_loads_ansys.json      15 项
    └── cb05_elastic_support_ansys.json  13 项   （共 84 项）
```

`<stem>` = `cb01_2span` / `cb02_3span` / `cb03_selfweight` /
`cb04_point_loads` / `cb05_elastic_support`。

## 算例

| 编号 | 结构 | 离散 | 荷载工况 | 支座反力 (kN) | 关键结果 |
|---|---|---|---|---|---|
| `cb01_2span` | 2×20 m 矩形 0.6×1.2 m，A=0.72 m²，I=0.0864 m⁴ | 20 单元 / 21 节点 | q=20 kN/m 全桥 | 150 / 500 / 150 | 中支座 −1.000e6 N·m；闭式解逐点核对偏差 4.6e-12 m |
| `cb02_3span` | 3×30 m 箱形 h=2.0 b=8.0 m，A=0.478 m²，I=0.4148 m⁴ | 30 单元 / 31 节点 | q=20 kN/m 全桥 | 240 / 660 / 660 / 240 | 中支座 −1.800e6 N·m |
| `cb03_selfweight` | CB01 几何 | 20 单元 / 21 节点 | q=20 kN/m **+ 自重** 55.446 kN/m（DENS=7850, g=9.81） | 565.85 / 1886.15 / 565.85 | 总荷载 3.01784 MN，中支座 −3.772e6 N·m |
| `cb04_point_loads` | CB01 几何 | 20 单元 / 21 节点 | 节点 6 / 16 各 **−200 kN** 集中力，无均布 | 62.5 / 275 / 62.5 | 中支座 −7.500e5 N·m |
| `cb05_elastic_support` | CB01 几何 | 20 单元 / 21 节点 | q=20 kN/m，中支座改 **K=2e6 N/m 弹簧** | 端支座 367.42 / 367.42，弹簧 65.15 | 中支座下沉 32.58 mm |

全部 E=2.06e11 Pa。CB01/CB02/CB04/CB05 自重关闭（`*GRAVITY 0 0` / 无 `ACEL`）；
CB03 开启。CB01 另用闭式解 `v(x)=q(L³x−3Lx³+2x⁴)/(48EI)` 逐节点核对。

## 工况语法对照（`.dat` ↔ APDL，均已实测）

| 语义 | FEM-Eigen `.dat` | ANSYS APDL | 备注 |
|---|---|---|---|
| 单元均布荷载 | `*ELEMENT_LOAD / 1-20 GY -q -q` | `SFBEAM,EID,1,PRES,q,q` | **符号相反**，正 PRES 向下 |
| 节点集中力 | `*CLOAD / 6 2 -200000` | `F,6,FY,-200000` | 真实外力，符号一致 |
| 自重（开） | `*MP 1 E=... DENS=7850` + `*GRAVITY / 0 -9.81` | `MP,DENS,1,7850` + `ACEL,0,9.81,0` | **符号相反**：FEM 存重力方向，ANSYS `ACEL` 是惯性加速度；两边一致荷载等价，实测相对差 1.3e-8 |
| 自重（关） | `*GRAVITY / 0 0` | 不写 `ACEL`（`MP,DENS` 可保留） | 注意求解器默认 `gy=-9.81`，必须显式关闭 |
| 弹性支承 | `*BOUNDARY / 11 2 0 K=2e6` | `ET,2,COMBIN14` + `KEYOPT,2,2,2` + `R,2,2e6` + 重合接地节点 `D,<gn>,ALL,0` + `E,<gn>,11` | `KEYOPT(2)=2` = 2D UY 弹簧，本机实测与 FEM-Eigen 对角加 `K` 严格等价 |
| 弹簧力 | **不输出**（弹簧自由度是自由自由度，残差为 0 被阈值过滤） | `ETABLE` / `*GET,ELEM,e,SMISC,1` | 由 `F = −K·UY` 反算；SMISC1 反号，比幅值 |

易踩的坑：

- `*GRAVITY` 是**开关**：`*MP` 给了 `DENS` 且默认 `gy=-9.81` 就会自动计自重；
  不要同时写 `*SELF_WEIGHT`（会双计）。
- ANSYS `COMBIN14` 的 `KEYOPT(2)` 不是方向编号随便填：实测 `0` 直接报错，
  `1`/`3` 弹簧不起作用（得到纯悬臂梁挠度），只有 `2` 给出正确 UY 弹簧。
- 弹簧单元/接地节点编号排在梁之后，后处理 `1..nnode` / `1..nelem` 循环才不会读到它们。
- `/POST1` 的 `ETABLE` 作用于**所有选中单元**：混用 COMBIN14 时必须先
  `ESEL,S,ELEM,,1,nelem` 只选梁单元（否则 ANSYS 会对弹簧单元请求不存在的 SMISC 2/6/12
  并报 warning），导出弹簧力前再 `ESEL,ALL`。五个算例现均 0 warning / 0 error。
- 3D 内核不支持 `*GRAVITY` / `*SELF_WEIGHT`（`DENS` 只进质量矩阵），本目录都是 2D 算例。

## 判定

- 判据 `|got − ref| <= atol + rtol*|ref|`，统一 `rtol = 1e-3`（沿用项目既有容差）。
- 特征量：跨内/跨中竖向挠度、支座零位移、全部支座竖向反力、中支座负弯矩、
  跨内最大正弯矩、端弯矩归零、端剪力幅值，以及弹性支承的**弹簧反力**。共 **84 项**。
- 支座位移／端弯矩参考值为 0，用 `atol`（1e-9 m / 1e-6 N·m）判定。
- **剪力只比幅值**：FEM-Eigen 的 `V1/V2` 与 ANSYS SMISC 2/8 符号约定相反，JSON 中
  这些项带 `"abs": true`。
- **弹簧反力**用 `kind: "spring_force"`：FEM-Eigen 的 `.result` 不列弹簧自由度反力，
  比对脚本由 `−K·(UY − 基准位移)` 反算，再与 ANSYS COMBIN14 SMISC1 的幅值比。
- 弯矩符号约定两边相同，直接比带符号值（BEAM3 SMISC 6 = i 端 M，12 = j 端 M）。
- 每次运行还做归档完整性检查：`.dat`/`.inp` 齐全，ANSYS 三个 CSV 行数与 `.dat`
  解析出的节点/单元数一致；有弹簧检查项时额外要求 `_ansys_spring.csv` 存在且非空。
- 各算例 `compare_report.txt` 另做竖向力平衡自查（反力合计 vs 外荷载合计），
  五个算例相对差均 ≤ 7.5e-9。

## 运行

CTest（默认离线，**不需要 ANSYS 许可**，约 0.5 s）：

```powershell
cd D:\AI-work\FEM-Eigen
ctest --test-dir build -R continuous_beam_benchmark --output-on-failure
ctest --test-dir build -L benchmark --output-on-failure
```

`run_tests.ps1` 全量（已包含 "Continuous beam ANSYS benchmark" 一段）：

```powershell
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1
```

手动：

```powershell
python tests\benchmark\continuous_beam\compare_continuous_beam.py `
  --bin build `
  --cases tests\benchmark\continuous_beam\cases `
  --reference tests\benchmark\continuous_beam\reference
```

**可选实时 ANSYS 复核**（需本机 v180 许可，会为每个算例执行归档的 `.inp`，5 个约 2–3 分钟；
不进 CTest，避免让日常回归依赖商业软件）：

```powershell
python tests\benchmark\continuous_beam\compare_continuous_beam.py `
  --bin build --cases tests\benchmark\continuous_beam\cases `
  --reference tests\benchmark\continuous_beam\reference --live-ansys
```

实时模式把 ANSYS 现场导出的 disp/reac/elem（以及 CB05 的 spring）CSV 与本次
FEM-Eigen 结果全量比对，五个算例全部 PASS：

```
cb01_2span LIVE UY (21 nodes): maxdiff 4.6e-12 m / scale 0.0009709 m = 4.74e-09  PASS
cb03_selfweight LIVE UY (21 nodes): maxdiff 4.9e-11 m / scale 0.003662 m = 1.34e-08  PASS
cb05_elastic_support LIVE spring force (1): max rel diff 3.25e-08  PASS
```

## 更新参考值（保持可审计）

`reference/*.json` 里每个 `ref` 都由 `build_reference.py` 从 `cases/*_ansys_*.csv`
抽取，并记录 `source`（精确到 CSV 文件 + 行键 + 列号），不手工填数。

```powershell
# 校验现有 JSON 是否仍与冻结 CSV 一致（不写文件）
python build_reference.py --check

# 重跑 ANSYS 并刷新 CSV 后，重新生成 JSON
python build_reference.py
```

完整重建链路（在 BridgeStudio 侧）：

```powershell
C:\Users\user\anaconda3\python.exe D:\AI-work\BridgeStudio\cases\run_cases.py --case all
```

它会在 `D:\AI-work\BridgeStudio\cases\CB01_*` ~ `CB05_*` 重新生成 `.dat`、`.inp`、
`.model.json`、FEM `.result`、ANSYS CSV 与 `compare_report.txt`（CB05 另有 `_spring.csv`），
把这些文件拷回 `cases/` 覆盖后再执行 `build_reference.py`。

## ANSYS 语义（本机 v180 实测）

- `ET,1,BEAM3` 与 FEM-Eigen 2D `BEAM` 同为平面 Euler-Bernoulli 梁；`R,1,A,I,H` 直接对应。
- `SFBEAM,ELEM,1,PRES,+q,+q` 的正 q 向下，与 `.dat` 的 `GY -q -q` 一致。
- SMISC 1/2/6 = i 端 N/V/M；7/8/12 = j 端 N/V/M；COMBIN14 SMISC 1 = 弹簧力。
- `*VWRITE,ARR(1,1),...` 一次写出数组全部剩余行，**外面不能再套 `*DO`**（否则重复写出行）。
- `*GET,NODE,I,RF,FY` 对非约束节点返回值无物理意义，参考值只取支座节点。
- 运行输出统一 ASCII，避免 PowerShell / CTest / 重定向各宿主解码不一致导致中文乱码。

## 本机环境备注

- `build/CMakeCache.txt` 中 `CMAKE_MAKE_PROGRAM` 指向已不存在的
  `D:/Program Files/CodeBlocks/MinGW/bin/mingw32-make.exe/mingw32-make.exe`，
  且本机只有 MSYS 版 `/usr/bin/make.exe`（与 MinGW Makefiles 生成器不兼容），
  因此当前 `cmake -S . -B build` 会在编译器 ABI 检查处失败。
- 已用临时目录 + `MSYS Makefiles` 生成器验证 `CMakeLists.txt` 语法并生成标准 CTest 条目，
  再同步进 `build/CTestTestfile.cmake`；`CMakeLists.txt` 已正式登记该测试，
  环境修复后任何一次正常 configure 都会自动生成同样的条目。
- 顺带把 CTest 里三个失效的 `C:/msys64/mingw64/bin/python3.14.exe` 引用修正为本机
  `C:/Users/user/anaconda3/python.exe`，`dat_byte_regression` / `ansys_benchmark` /
  `bridge_benchmark` 现均可运行。
- 比对脚本复用了 `../ansys/compare_ansys.py` 的 `parse_result` / `extract`，
  保证两个 benchmark 不会各自维护一套 `.result` 列映射；`spring_force` 是本目录扩展的
  唯一特殊 kind。
