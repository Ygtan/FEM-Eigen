# `tests/benchmark/` — 外部参考值回归总索引

本目录是 **BridgeStudio 前处理 → FEM-Eigen 计算 → ANSYS 18.0 校核** 链路的
回归基线集合。与 `tests/regression/`（冻结自身输出、逐字节比对）不同，这里
冻结的是**外部软件或闭式解的参考值**，用来防「程序整体一起漂移但自身回归
仍绿」的假绿。

所有算例遵守同一套约定：

1. **单一数据源**：`.dat` / `.inp` 都由 `BridgeStudio` 的 `BridgeSpec`（2D）
   或 `ModelSpec3D`（3D）同源导出，不手工抄写两份输入。
2. **冻结证据**：真实 ANSYS 18.0 批处理跑出的 CSV / `.out` 归档在
   `cases/`，由 `build_reference*.py` 提炼成 `reference/*.json` 黄金值；
   `--check` 可随时验证黄金值与冻结 CSV 逐项一致。
3. **默认离线**：CTest 与 `run_tests.ps1` 只需本机构建的 exe，不需要
   ANSYS 许可；`--live-ansys` 在临时目录重跑 APDL 并与冻结 CSV 逐值复核，
   **不覆写**任何归档文件。
4. **每次回归都重算 FEM**：`compare_*.py` 用当前 `--bin` 下的 exe 重跑
   归档 `.dat`，先验证 FEM 能逐位复现冻结结果，再与 ANSYS 黄金值比对。

## 目录与 CTest 注册

| 目录 | CTest 名称 | 内核 / 单元 | 覆盖算例 | 黄金值项数 |
|---|---|---|---|---|
| `ansys/` | `ansys_benchmark` | 2D `example_mixed.exe` / BEAM4 | 100 m 简支梁 | 外部闭式/ANSYS 特征量 |
| `bridge/` | `bridge_benchmark` | 2D 闭式解 | 解析梁基准 | 冻结解析值 |
| `continuous_beam/` | `continuous_beam_benchmark` | 2D BEAM3 | CB01~CB05 | **84** |
| `continuous_beam_3d/` | `continuous_beam_3d_benchmark` | 3D BEAM4 / BEAM188 | CB3D01~02 | **113** |
| `cable_arch_3d/` | `cable_arch_3d_benchmark` | 3D BEAM + TRUSS 混合 | CB3D03/04/09 | **695** |
| `cable_arch_3d/` | `cable_solver_smoke` | 非线性成桥态自洽 | CB3D05 | **15** |
| `element_loads_3d/` | `element_loads_3d_benchmark` | 3D 分布荷载 / 自重 | CB3D06~08 | **3 × 366** |
| `element_loads_3d/` | `element_loads_3d_kernel_validation` | 3D 荷载内核解析守卫 | 67 组 | **67** |
| `nonlinear_ansys_3d/` | `nonlinear_ansys_3d_benchmark` | NLGEOM + LINK180 成桥态 | CB3D10 | **325** |
| `suspension_ansys_3d/` | `suspension_ansys_3d_benchmark` | 悬索桥找形（地锚/自锚） | CB3D11~12 | **292 + 282 = 574** |
| `suspension_erection_3d/` | `suspension_erection_3d_benchmark` | 自动平衡边跨 + 施工阶段 | CB3D13 | **478** |
| `cable_sag_ansys_3d/` | `cable_sag_ansys_3d_benchmark` | 索垂度 / Ernst 弹模 + `W>0` 空缆态 | CB3D14~15 | **573** |
| `beam_nl_ansys_3d/` | `beam_nl_ansys_3d_benchmark` | 共转 CR 梁大转动 / 解析几何刚度 vs `BEAM188 + NLGEOM` | CB3D16 | **1776** |
| `snap_through_ansys_3d/` | `snap_through_ansys_3d_benchmark` | 球面弧长法后屈曲 / snap-through 全过程 vs `LINK180 + ARCLEN` | CB3D17 | **1686** |
| `stable_path_3d/` | `stable_path_3d_benchmark` | B3 主链路 `STABILITY=1` opt-in vs CB3D17 冻结路径 / von Mises 闭式解析锚点 | CB3D18 | **3311** |
| `frame_postbuckling_3d/` | `frame_postbuckling_3d_benchmark` | 共转 CR 梁刚接框架（Williams toggle 型）后屈曲全过程 vs `BEAM188 + ARCLEN` / 三套独立线性刚度锚点 | CB3D19 | **12513** |
| `completion_stability_3d/` | `completion_stability_3d_benchmark` | **成桥态起点**全过程稳定：`STABASIS=COMPLETION` / `R = F0 + λ·Finc − Fint(u)` / V 形双索独立闭式 `λ(v)` 对拍 `STATIC` / `TENSION` | CB3D20（纯 FEM 侧） | **2963** |
| `completion_bridge_3d/` | `completion_bridge_3d_benchmark` | **真实双塔三跨斜拉桥**成桥态两段 vs ANSYS（step1 恒载 + step2 增量荷载；荷载控制逐子步 + `ARCLEN` 终点交叉验证）/ 同源 `L0` + `INISTATE eps` 独立锚点 | CB3D21 | **2392** |
| `completion_basis_3d/` | `completion_basis_3d_benchmark` | **第十九次 ③a**：`MODE=COMPLETION` 成桥态基准同源修复（`basisModel`）/ 跨 `MODE` 独立锚点 / 默认零污染 / CB3D13 五个 `[ERECTION*]` 段共存核对 | CB3D10、CB3D13 | **110** |
| `bridge_limit_point_3d/` | `bridge_limit_point_3d_benchmark` | **第十九次 ③a**：真实斜拉桥**越过极值点 / 下降段** FEM 球面弧长 vs ANSYS 两段 `ARCLEN`；按 `ΔUX` 插值 `λ` 逐点对拍 / 16 索 `L0` + `INISTATE` + 弦长同源锚点 | CB3D22 | **140** |

另有 4 个解析/内核单测（`analytical_exact_solutions_2d/3d`、
`nonlinear_co_rotational_3d`、`shape_sensitivity_3d`）与 1 个
`dat_byte_regression`：② 阶段为 23 个 CTest，③a 新增第 24 / 25 项后
**合计 25 个 CTest**。当前实测全量 `ctest` 25/25、
`run_tests.ps1 -BuildDir _d_build` `ALL TESTS PASSED`。

> **构建工具链必须与冻结证据一致**：历史冻结输出全部由
> `C:/msys64/mingw64/bin/g++.exe`（`-O3 -DNDEBUG`）产出。若改用
> `C:/msys64/ucrt64/bin/c++.exe`，牛顿残差 / 模态向量的末位浮点会漂移，
> 逐字节守卫会红（不是力学回归）。干净 configure 请显式指定编译器：
> `cmake -S . -B _d_build -G Ninja -DCMAKE_BUILD_TYPE=Release
> -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe`。
> 第 24 项冻结证据已用该规范工具链重新生成、保持严格逐字节重放；
> 第 25 项的 FEM 证据与 BS 侧算例目录共享，改用
> `tests/benchmark/replay_tolerance.py` 的结构化跨工具链比较
> （`λ/CTRL/DS/PIVOT` 相对 `<1e-9`、`RES` 绝对 `<0.1 N`，其余行逐位一致）。

各目录的 `README.md` 记录该阶段的建模口径、ANSYS APDL 语义、容差来源与
已知边界；`nonlinear_ansys_3d/probes/` 归档 A~G 阶段全部 16 个 ANSYS 语义
探测脚本，`beam_nl_ansys_3d/probes/`、`snap_through_ansys_3d/probes/`
与 `frame_postbuckling_3d/probes/` 归档 J/K/L/N 阶段物理诊断脚本，
`completion_stability_3d/probes/` 归档第十九次 ① 的冻结 cases 生成器
与 ② 阶段的四个 ANSYS 两段求解物理诊断探针
（只入库 `.py` / `.md`，不进 CTest）。

## 从零复现（换机 / 重建）

`build/` 及其 `CTestTestfile.cmake` 由 CMake 生成、被 `.gitignore` 排除，
**不入库**。干净的 configure 会自动生成全部 25 条注册，不需要手工同步。

```bash
# MSYS2 / MinGW-w64 环境。**必须显式锁定 mingw64 编译器**（见上方工具链说明），
# 否则 CMake 可能挑到 ucrt64 的 c++.exe，逐字节守卫会红。
export PATH="/c/msys64/usr/bin:/c/msys64/mingw64/bin:$PATH"

# 方案 A：Ninja（推荐，路径处理最稳）——③a 全量回归实测用 _d_build
cmake -G Ninja \
  -DCMAKE_MAKE_PROGRAM="C:/msys64/mingw64/bin/ninja.exe" \
  -DCMAKE_CXX_COMPILER="C:/msys64/mingw64/bin/g++.exe" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython3_EXECUTABLE="C:/Users/user/anaconda3/python.exe" \
  -S D:/AI-work/FEM-Eigen -B D:/AI-work/FEM-Eigen/_d_build
cmake --build D:/AI-work/FEM-Eigen/_d_build -j2
ctest --test-dir D:/AI-work/FEM-Eigen/_d_build --output-on-failure

# 方案 B：MSYS Makefiles
cmake -G "MSYS Makefiles" \
  -DCMAKE_MAKE_PROGRAM="C:/msys64/usr/bin/make.exe" \
  -DCMAKE_CXX_COMPILER="C:/msys64/mingw64/bin/g++.exe" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython3_EXECUTABLE="C:/Users/user/anaconda3/python.exe" \
  -S D:/AI-work/FEM-Eigen -B D:/AI-work/FEM-Eigen/build_msys
```

**实测（2026-11-15，本机 MSYS2 / GCC 16.2 / mingw64）**：方案 A 干净
configure 自动生成 **25 条 `add_test`**，`cmake --build -j2` 编译产出
**全部 40 个 exe**（K 阶段 `beam_nl_probe`、L 阶段 `snap_probe`、N 阶段
`frame_probe` 之后，第十九次 ① / ② / ③a 都只改 `cable_solver.exe` /
头文件 / 冻结证据、不新增 exe），全新目录
`ctest` **25/25 Passed**。即换机后无需手工同步
`CTestTestfile.cmake`，也无需 ANSYS 许可即可全量回归；M 阶段主链路
口径与数字见 `stable_path_3d/README.md`，N 阶段 CR 梁刚接框架
后屈曲口径见 `frame_postbuckling_3d/README.md`，第十九次 ①
成桥态起点口径见 `completion_stability_3d/README.md`，第十九次 ②
真实斜拉桥两段对标口径见 `completion_bridge_3d/README.md`，③a 的
同源修复与极值点 / 下降段口径见 `completion_basis_3d/README.md`、
`bridge_limit_point_3d/README.md`。
全量重编译在本机
约需 10 分钟（头文件 `fem_3d_nl.hpp` 较大）；内存有限时建议 `-j 2`，
实测 `-j 8` 与其它任务并发会触发 `cc1plus.exe: out of memory`，与代码无关。

> **本机历史坑**：旧 `build/` 缓存的生成器是 `MinGW Makefiles`，却把
> `CMAKE_MAKE_PROGRAM` 配成 MSYS 的 `/usr/bin/make.exe`，两者不兼容
> （试编译阶段 Makefile 路径被反斜杠截断、报乱码）。不要再对该缓存直接
> 重配；用上面任一方案新建构建目录即可。`build/` 内已编译的
> `cable_solver.exe` 含 `[ERECTION_*]` 新段，日常回归仍可直接使用。

```powershell
# 全量回归（PowerShell，输出 ASCII）。-BuildDir 需指向与冻结证据同工具链的构建目录，
# 默认 build；干净 Ninja 构建用 _d_build。
powershell -ExecutionPolicy Bypass -File D:\AI-work\FEM-Eigen\run_tests.ps1 -BuildDir _d_build
```

## 新增一个 benchmark 的检查清单

1. `cases/` 放归档 `.dat` / `.inp` / `.out` / CSV / 报告；`.csv` 需被
   `.gitignore` 的 `!tests/benchmark/*/cases/*.csv` 放行。
2. `reference/*.json` 由 `build_reference*.py` 从冻结 CSV 生成，派生量
   要在 JSON 里写清 `derived` 口径并由脚本重算校验，防黄金值漂移。
3. `compare_*.py` 提供 `--bin/--cases/--reference`，离线为默认，
   `--live-ansys` 可选；退出码 0 才算 PASS。
4. `CMakeLists.txt` 里 `add_test(NAME ...)` + `set_tests_properties(...
   TIMEOUT 900 LABELS "benchmark")`，并在 `run_tests.ps1` 末尾追加一段。
5. 本文件与对应目录 `README.md` 同步更新；跑一次 `ctest` 与
   `run_tests.ps1 -BuildDir <dir>` 确认全绿；若新增冻结输出，须用与历史
   证据相同的 mingw64 工具链生成（见上方「构建工具链」说明）。
