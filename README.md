# FEM-Eigen 有限元框架使用与开发说明

这是一个基于 [Eigen](https://eigen.tuxfamily.org/) 的轻量级有限元（FEM）框架，
用 C++17 编写，目标是：**单元组装清晰、求解器可切换、便于扩展**。

当前已支持：

- 单元：一维/二维/三维桁架（杆）、二维平面应力三角形（CST）、三维四面体（T4）、二维欧拉梁、
  只受拉索（CABLE）、只受压杆（STRUT）。
- 三维梁格：三维桁架/梁、混合自由度、施工阶段、预应力、徐变收缩、影响线/面、多车道活载和包络。
- 三维输入格式：`input_format_3d.md`。
- 自动测试脚本：`run_tests.ps1`。
- 多自由度组装：每个节点可带 1~3 个自由度，自动映射到全局自由度编号。
- 可切换求解后端：Eigen 直接法/迭代法、AMGCL（代数多重网格）、CHOLMOD（SuiteSparse）。
- CHOLMOD 分阶段求解：符号分解复用、数值分解复用、单/多右端项。
- 输入文件：支持 `*.dat` 关键字格式，读取二维梁/桁架/索/压杆混合结构并求解；
  荷载包括集中荷载、单元均布/梯形荷载、自重、体系温升、梁温度梯度、初应变预应力。
- 分析类型：`*ANALYSIS LINEAR`（方案 A）和 `*ANALYSIS NONLINEAR`（方案 C，轴向单元几何非线性）。
- 混凝土梁预应力：`*PRESTRESS` 可自动生成等效荷载，或生成 `CABLE` 预应力筋并接入非线性分析。
- 影响线与车道荷载：`*INFLUENCE` 输出指定响应影响线；`*LANE_LOAD` 支持 `LINEAR` 影响线叠加和 `NONLINEAR` 直接非线性扫描包络。
- 约束方程：`*CE` 支持 APDL 风格多点约束方程，采用罚函数法兼容现有对称求解器。
- 混凝土徐变与收缩：`*CREEP`、`*SHRINKAGE` 第一阶段采用 AEMM 龄期调整有效模量法。
- 日照温度梯度：`*TEMP_PROFILE` 按规范温度曲线计算梁截面 `N_T`、`M_T`。
- 施工阶段正装分析：`*STAGE`、`*STAGE_LOAD`、`*STAGE_ELEMENT_LOAD`，支持单元激活/钝化；
  阶段预应力 `*STAGE_PRESTRESS`，支座激活/钝化 `*STAGE_BOUNDARY` / `*STAGE_UNBOUNDARY`；
  阶段徐变收缩 `*STAGE_CREEP` / `*STAGE_SHRINKAGE`；
  阶段温度梯度 `*STAGE_TEMP_PROFILE`；
  施工临时荷载卸载 `*STAGE_UNLOAD`。
- 阶段时变曲线：`*CREEP_CURVE` / `*SHRINKAGE_CURVE`，按阶段 `TIME` 插值徐变、收缩。
- 规范徐变收缩模型：`*CREEP_MODEL` / `*SHRINKAGE_MODEL`，内置简化 CEB-FIP 1990 公式。
- 倒拆分析：`example_stage --backward`；按施工逆序移除荷载、反向激活/钝化单元，并继承非线性状态进行倒拆近似。
- 正装-倒拆迭代：`*STAGE_TARGET` + `*STAGE_CONTROL`，支持多变量/多目标、变量上下界和 Levenberg-Marquardt 求解。
- 迭代目标支持节点位移 `DISP`、支座反力 `REACTION`、单元内力 `FORCE`、单元应力 `STRESS`。
- 多目标支持权重 `WEIGHT`，用于 Levenberg-Marquardt 加权最小二乘。
- 支持不等式约束 `*STAGE_CONSTRAINT`，可限制位移、反力、内力、应力上下界。
- 新增 `*CAMBER`，按最终施工位移计算节点预拱度并输出。
- 新增非线性预拱度迭代：`example_stage --camber-iterate`。
- 新增施工反馈修正：`*FEEDBACK` + `example_stage --feedback`。
- 新增施工阶段可视化：`example_stage_viz`，输出每阶段模型、位移、弯矩、剪力、轴力、上下翼缘应力 SVG。
- 新增施工阶段包络：`example_stage_envelope`，输出各单元施工阶段内力/应力最大最小值。
- 施工阶段包络图：轴力、剪力、弯矩、上翼缘应力、下翼缘应力 SVG。
- 单元应力：桁架/索/压杆轴向应力；梁上下翼缘正应力（需输入抗弯模量 `WT`/`WB`）。
- 可视化：`example_viz` 输出 SVG（模型 + 变形 + 弯矩图 + 应力图）。

---

## 目录

1. [环境与编译](#1-环境与编译)
2. [快速开始：运行示例](#2-快速开始运行示例)
3. [使用说明](#3-使用说明)
4. [示例程序清单](#4-示例程序清单)
5. [开发说明](#5-开发说明)
6. [求解器选型建议](#6-求解器选型建议)
7. [已知问题与注意事项](#7-已知问题与注意事项)

---

## 1. 环境与编译

### 1.1 工具链

本工程在 Windows + MSYS2 上开发，使用 `mingw64` 工具链：

- 编译器：`g++ 16.x`（`C:\msys64\mingw64\bin\g++.exe`）
- 构建工具：`CMake` + `Ninja`（`C:\msys64\mingw64\bin\cmake.exe`、`ninja.exe`）
- C++ 标准：`C++17`

### 1.2 依赖

| 依赖 | 用途 | 来源 |
|------|------|------|
| Eigen 3.3.9 | 矩阵/稀疏矩阵与基本求解器 | 项目内置 `external/eigen3` |
| AMGCL 1.5.0 | AMG 预条件 + BiCGStab 后端 | 项目内置 `external/amgcl` |
| SuiteSparse CHOLMOD 7.x | 高性能稀疏 Cholesky 直接法 | 需用 MSYS2 安装（见下） |

Eigen 和 AMGCL 是**纯头文件库**，已经内置到 `external/` 目录，无需额外安装。
CHOLMOD 需要链接动态库，用 MSYS2 安装一次：

```powershell
# 先设置 MSYS2 环境变量
$env:PATH = "C:\msys64\usr\bin;C:\msys64\mingw64\bin;$env:PATH"

# 同步数据库并安装 SuiteSparse（含 CHOLMOD）
pacman -Sy --noconfirm
pacman -S --noconfirm --needed mingw-w64-x86_64-suitesparse
```

> 说明：如果没装 CHOLMOD，框架仍然能编译运行，只是 `SolverType::CHOLMOD`
> 会不可用（运行时抛异常）。CMake 会自动检测 CHOLMOD 头文件是否存在。

### 1.3 编译命令

```powershell
$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"

cmake -S "D:\AI-work\FEM-Eigen" -B "D:\AI-work\FEM-Eigen\build" `
      -G Ninja -DCMAKE_BUILD_TYPE=Release

cmake --build "D:\AI-work\FEM-Eigen\build"
```

编译产物在 `build/` 目录下，包括 `fem_solver.exe` 和各个 `example_*.exe`。

> **已知小问题**：这台机器上 Windows Defender 偶尔会锁住 Ninja 的 `.d` 依赖文件，
> 导致 `cmake --build` 偶发报 `file is being used by another process`。
> 这是环境问题，与代码无关，**重新执行一次 `cmake --build` 即可**。
> 也可以换个构建目录或用 `-DCMAKE_DEPENDS_USE_COMPILER=FALSE` 重新配置绕开。

---

## 2. 快速开始：运行示例

所有示例都是控制台程序，运行后打印结果并自检 `[PASS]` / `[FAIL]`。

```powershell
cd D:\AI-work\FEM-Eigen

.\build\fem_solver.exe          # 最初的一维杆例子（3 单元）
.\build\example_1d_bar.exe      # 一维杆（框架版）
.\build\example_2d_truss.exe    # 二维平面桁架（每节点 2 自由度）
.\build\example_2d_plane.exe    # 二维平面应力 CST 三角形
.\build\example_3d_tet.exe      # 三维四面体
.\build\example_solvers.exe     # 演示各求解后端
.\build\example_amgcl.exe       # 2 万自由度基准对比（LDLT/CHOLMOD/AMGCL）
.\build\example_cholmod.exe     # CHOLMOD 分阶段 + 多右端项
```

如果中文显示为乱码，先执行 `chcp 65001`（程序本身也会尝试把控制台切到 UTF-8）。

---

## 3. 使用说明

框架的公共接口都在头文件 `fem.hpp` 里，命名空间为 `fem`。使用时：

```cpp
#include "fem.hpp"
using namespace fem;
```

### 3.1 核心类型

```cpp
using SparseMat = Eigen::SparseMatrix<double>;  // 稀疏矩阵（列压缩 CSC）
using DenseMat  = Eigen::MatrixXd;              // 稠密矩阵
using DenseVec  = Eigen::VectorXd;              // 稠密向量
```

### 3.2 单元类型

所有单元都有两个统一接口：

```cpp
int dofsPerNode() const;  // 每个节点的自由度数
int numNodes() const;     // 单元节点数
DenseMat stiffness() const; // 单元刚度矩阵
```

当前内置四个单元：

| 单元 | 字段 | 每节点自由度 | 说明 |
|------|------|------------|------|
| `TrussElement` | `E, A, dim, x0, x1` | `dim`（1/2/3） | 两节点杆/桁架，自动做方向余弦变换 |
| `TriangleCST` | `E, nu, t, x0, x1, x2` | `2` | 常应变三角形，平面应力 |
| `Tetra4` | `E, nu, x0, x1, x2, x3` | `3` | 常应变四面体，各向同性三维实体 |
| `EulerBeam2D` | `E, A, I, x0, x1` | `3` | 二维欧拉-伯努利梁（u/v/θ，含轴向+弯曲） |

`TriangleCST` 和 `Tetra4` 还额外提供：

```cpp
DenseMat B() const;   // 应变-位移矩阵
DenseMat D() const;   // 本构矩阵
double area() const;  // 或 volume()
```

这些可以用于后处理（由节点位移算应变/应力：`σ = D B u_e`）。

### 3.3 自由度编号与组装

全局自由度编号约定：**节点号 × 每节点自由度 + 本节点内自由度序号**。

```cpp
int dofIndex(int node, int localDof, int dofsPerNode);
std::vector<int> nodeDofs(const std::vector<int>& nodes, int dofsPerNode);
```

组装流程：先对每个单元算出 `stiffness()` 和它对应的全局自由度编号，
再统一组装成全局稀疏刚度矩阵：

```cpp
std::vector<DenseMat> kes;              // 各单元刚度矩阵
std::vector<std::vector<int>> gdofs;    // 各单元全局自由度编号

// ... 填充 kes 与 gdofs ...

SparseMat K = assemble(kes, gdofs, nDofs);
```

底层还提供 `scatter(Ke, gdofs, triplets)`，把单个单元矩阵散入三元组列表，
一般不需要直接调用。

### 3.4 载荷与边界条件

载荷向量是一个 `DenseVec F`（长度 = 总自由度数），按需给对应自由度赋值即可。

Dirichlet 边界条件（指定位移）用行/列消去法施加：

```cpp
void applyDirichlet(SparseMat& K, DenseVec& F,
                    const std::vector<int>& constrainedDofs,
                    const std::vector<double>& prescribedValues);
```

例：固定节点 0（`u0 = 0`，1 自由度/节点）：

```cpp
applyDirichlet(K, F, {0}, {0.0});
```

例：二维问题里固定节点 0 的 x、y 两个自由度，并给节点 1 的 y 方向一个非零位移：

```cpp
applyDirichlet(K, F,
               {2*0+0, 2*0+1, 2*1+1},
               {0.0,   0.0,   0.01});
```

该函数会**同时清零被约束自由度所在的行和列**并置对角线为 1，
因此矩阵仍保持对称，适合对称求解器。

### 3.5 求解器接口

统一求解入口：

```cpp
DenseVec solve(const SparseMat& K, const DenseVec& F,
               SolverType type = SolverType::SimplicialLDLT);
```

`SolverType` 可选项：

| 枚举值 | 类型 | 适用 |
|--------|------|------|
| `SimplicialLDLT` | 直接法 | **默认**，对称正定（SPD） |
| `SimplicialLLT` | 直接法 | SPD（Cholesky） |
| `SparseLU` | 直接法 | 一般方阵 |
| `SparseQR` | 直接法 | 最小二乘 / 欠定 / 满秩 |
| `ConjugateGradient` | 迭代法 | SPD（Jacobi 预条件） |
| `BiCGSTAB` | 迭代法 | 一般矩阵（Jacobi 预条件） |
| `AMGCL` | 迭代法 | 大 SPD（AMG + BiCGStab），需 AMGCL 头文件 |
| `CHOLMOD` | 直接法 | SPD 高性能（Cholesky），需 SuiteSparse |

例：

```cpp
DenseVec U = solve(K, F);                              // 默认 SimplicialLDLT
DenseVec U = solve(K, F, SolverType::CHOLMOD);         // 换 CHOLMOD
DenseVec U = solve(K, F, SolverType::ConjugateGradient);
```

`solveSPD(K, F)` 是旧接口的兼容别名，等价于 `solve(K, F, SimplicialLDLT)`。

### 3.6 CHOLMOD 分阶段求解（复用 + 多右端项）

对“同一结构、多次数值求解”的场景（例如非线性迭代、时间步进、多载荷工况），
推荐使用 `CholmodSolver`：

```cpp
CholmodSolver solver;

solver.analyzePattern(K0);   // 符号分解：同一非零结构只做一次

solver.factorize(K1);        // 数值分解：复用符号分解
DenseVec u1 = solver.solve(f1);   // 单右端项

solver.factorize(K2);        // K2 与 K1 结构相同，仅数值不同
DenseVec u2 = solver.solve(f2);

solver.factorize(K1);
DenseMat U = solver.solve(B); // 多右端项：B 是 n×k，一次分解解 k 个载荷
```

`CholmodSolver` 的完整接口：

| 方法 | 说明 |
|------|------|
| `analyzePattern(K)` | 符号分解，只依赖非零结构，开销较大，同一结构调用一次 |
| `factorize(K)` | 数值分解，复用已有符号分解 |
| `solve(F)` | 单右端项，返回 `DenseVec` |
| `solve(B)` | 多右端项，`B` 为 `n×k`，返回 `n×k` 的 `DenseMat` |
| `clear()` | 释放内部矩阵/因子，之后可重新 `analyzePattern` |

类遵循 RAII，析构时自动释放 CHOLMOD 内存，不可拷贝。
`solve(K, F, SolverType::CHOLMOD)` 内部就是 `analyzePattern → factorize → solve`
的一次性便捷封装。

### 3.7 一个最小完整例子

下面是一维杆（3 单元、4 节点、左端固定、节点 2 受 1000 N）的完整程序：

```cpp
#include "fem.hpp"
#include <iostream>
#include <vector>

int main() {
    using namespace fem;
    setupConsoleUTF8();

    const double E = 200e9, A = 0.01;
    const int dim = 1, nNodes = 4, dofsPerNode = dim, nDofs = nNodes;

    std::vector<Eigen::Vector3d> x = {{0,0,0},{1,0,0},{2,0,0},{3,0,0}};
    std::vector<std::pair<int,int>> conn = {{0,1},{1,2},{2,3}};

    std::vector<DenseMat> kes;
    std::vector<std::vector<int>> gdofs;
    for (auto [a,b] : conn) {
        TrussElement e{E, A, dim, x[a], x[b]};
        kes.push_back(e.stiffness());
        gdofs.push_back(nodeDofs({a,b}, dofsPerNode));
    }

    SparseMat K = assemble(kes, gdofs, nDofs);
    DenseVec F = DenseVec::Zero(nDofs);
    F(2) = 1000.0;
    applyDirichlet(K, F, {0}, {0.0});

    DenseVec U = solve(K, F, SolverType::CHOLMOD);
    for (int i = 0; i < nNodes; ++i)
        std::cout << "u" << i << " = " << U(i) << "\n";
    return 0;
}
```

输出应为 `u = [0, 5e-7, 1e-6, 1e-6]`。

---

### 3.8 输入文件（`*.dat`）与混合结构求解

二维梁/桁架/索/压杆结构可以直接用关键字格式的文本文件描述，然后用 `fem_io.hpp` 读取求解。
自由度约定：**每节点 3 个自由度 `[u, v, θ]`**，文件里 `dof` 编号为 `1=u`、`2=v`、`3=θ`。
桁架/索/压杆只贡献 u/v，梁贡献 u/v/θ，各类型可以混合。

推荐使用新格式（`*ET/*MP/*R/*EN`）编写复杂模型，完整约定见 `input_format.md`；
下面的旧式内联写法仍兼容，便于快速建立小模型。

文件格式（`#` 或 `!` 开头为注释，关键字大小写不敏感）：

```text
*TITLE
<标题>

*ANALYSIS LINEAR                 ! 可选；LINEAR=方案A，NONLINEAR=方案C

*NODE
<节点号> <x> <y>

*ELEMENT, TYPE=TRUSS, E=<值>, A=<值>[, ALPHA=<α>]
<单元号> <节点i> <节点j>

*ELEMENT, TYPE=CABLE, E=<值>, A=<值>[, ISTRN=<初应变>]
<单元号> <节点i> <节点j>

*ELEMENT, TYPE=STRUT, E=<值>, A=<值>[, ISTRN=<初应变>]
<单元号> <节点i> <节点j>

*ELEMENT, TYPE=BEAM, E=<值>, A=<值>, I=<值>[, ALPHA=<α>, H=<h>, WT=<Wt>, WB=<Wb>]
<单元号> <节点i> <节点j>

*BOUNDARY
<节点号> <dof> <值>

*CLOAD
<节点号> <dof> <值>

*GRAVITY
<gx> <gy>                        # 重力加速度，默认 0 -9.81（+Y 向上）

*ELEMENT_LOAD
<单元号> <方向> <q1> <q2>          # 单元均布/梯形荷载（N/m）

*SELF_WEIGHT
<单元号或ALL> <RHO|MASS> <参数>    # 自重

*TEMPERATURE
<ΔT>                             # 体系温度升降（所有单元）

*ELEMENT_TEMPERATURE
<单元号> UNIFORM  <ΔT>            # 指定单元整体升温
<单元号> GRADIENT <ΔT>            # 梁温度梯度（需 H）

*END
```

**荷载方向**（`*ELEMENT_LOAD` 的 `<方向>`）：

| 写法 | 含义 |
|------|------|
| `GX` / `GY` | 全局 X / Y |
| `LX` / `LY` | 局部轴向 / 横向（`LX`=沿 i→j，`LY`=垂直杆轴） |
| `G<θ>` | 与全局 +X 成 θ 角（逆时针为正，支持小数） |
| `L<θ>` | 与杆轴 i→j 成 θ 角（逆时针为正，支持小数） |

例：`G30`、`G-90`、`L25.89`、`L-90`。

**自重规则**：

- `RHO`：荷载强度 = `ρ·g·A`（`A` 取该单元截面积）。
- `MASS`：荷载强度 = `m·g`（`m` 为线质量）。
- `g = |重力加速度|`，方向取 `*GRAVITY` 给出的方向。

**温度**：

- `*TEMPERATURE <ΔT>`：所有单元统一温升，等效轴向荷载 `EA·α·ΔT`。
- `*ELEMENT_TEMPERATURE UNIFORM <ΔT>`：指定单元整体温升。
- `*ELEMENT_TEMPERATURE GRADIENT <ΔT>`：指定梁单元上下表面温差，弯曲等效弯矩 `EI·α·ΔT/h`（需 `H`）。
- 温度荷载需要 `*ELEMENT` 里的 `ALPHA`（热膨胀系数）。

**单元应力**：

- 桁架：轴向应力 `σ = N/A`，输出两端 `σ1`、`σ2`。
- 梁：上下翼缘正应力 `σ_top = N/A − M/WT`、`σ_bottom = N/A + M/WB`，输出两端共 4 个值。
- 应力符号：受拉为正；“上翼缘”=局部 +y 侧（水平杆 i→j 即全局上方），“下翼缘”=局部 −y 侧。
- 梁应力需要 `*ELEMENT` 里的 `WT`（上翼缘抗弯模量）和 `WB`（下翼缘抗弯模量）。

**CABLE / STRUT 与初应变**：

- `CABLE` 只受拉：计算轴力 `N ≤ 0` 时退出工作。
- `STRUT` 只受压：计算轴力 `N ≥ 0` 时退出工作。
- 轴向机械轴力统一为 `N = E·A/L·ΔL + E·A·ISTRN`（受拉为正）。
- `ISTRN > 0` 表示预拉应变；`ISTRN < 0` 表示预压应变，也是 cable/strut 的初始松弛/预紧状态。
- 求解默认采用**方案 A**：线性系统 + cable/strut 开关迭代。
- `*ANALYSIS TYPE=NONLINEAR` 启用**方案 C**：Newton-Raphson 几何非线性（轴向单元） + cable/strut 开关；
  `TOL` 和 `MAXITER` 可分别控制收敛容差和最大迭代次数。
- 方案 C 中 `SimplicialLDLT`、`CHOLMOD`、`AMGCL` 仍作为 Newton 迭代内部的线性求解器使用；切线非正定时可改用 `--solver LU` 或 `--solver BICGSTAB`。

**混凝土梁预应力 `*PRESTRESS`**：

```text
*PRESTRESS
1 ELEM=1-10 P=1500e3 F=0.20 E1=0.0 E2=0.0 LOSS=0.85 EPS=195e9 APS=0.001 METHOD=AUTO
```

- `P` 为张拉控制力，`LOSS` 为总损失系数，有效预应力 `Pe=P·LOSS`。
- `F` 为跨中下垂量，`E1/E2` 为两端偏心距，向下为正。
- `METHOD=EQUIV`：生成等效均布荷载和端部轴力/弯矩，适合线性分析。
- `METHOD=TENDON`：自动生成 `CABLE` 预应力筋和 `RIGID_PIN` 销接刚臂，适合非线性分析。
- `METHOD=AUTO`：`LINEAR` 使用 `EQUIV`，`NONLINEAR` 使用 `TENDON`。
- 不需要特别注明端节点/内部节点；程序按 `ELEM` 顺序自动识别首端节点、末端节点和中间内部节点。

**完整示例**（悬臂梁：自重 + 均布荷载 + 温度）：

```text
*NODE
1 0.0 0.0
2 4.0 0.0

*ELEMENT, TYPE=BEAM, E=200e9, A=0.01, I=1e-4, ALPHA=1.2e-5, H=0.5
1 1 2

*BOUNDARY
1 1 0.0
1 2 0.0
1 3 0.0

*GRAVITY
0 -9.81

*SELF_WEIGHT
ALL RHO 7850

*ELEMENT_LOAD
1 LY -20000 -20000

*TEMPERATURE
15

*CLOAD
2 2 -10000

*END
```

读取、求解与结果输出：

```cpp
#include "fem_io.hpp"

Model2D m = readModel2D("model.dat");
ModelResult r = solveModel2D(m, SolverType::SimplicialLDLT);
std::vector<ElementForce> forces = elementForces2D(m, r.U);
std::vector<ElementStress> stresses = elementStress2D(m, forces);
writeResult("model.result", m, r, forces);
// r.U 位移；r.R 反力；forces 单元内力；stresses 单元应力（梁需 WT/WB）
```

命令行示例：

```powershell
.\build\example_mixed.exe .\cantilever.dat
.\build\example_mixed.exe .\mixed_truss_beam.dat
.\build\example_mixed.exe .\arch_bridge_100m.dat
.\build\example_mixed.exe .\model.dat --solver CHOLMOD   # 运行时选择求解器
.\build\example_mixed.exe .\model.dat --compare          # 对比 LDLT/CHOLMOD/AMGCL
.\build\example_mixed.exe .\verify_cable_nonlinear.dat   # 方案 C 几何非线性
.\build\example_mixed.exe .\verify_prestress_equiv.dat   # 预应力等效荷载法
.\build\example_mixed.exe .\verify_prestress_tendon.dat  # 自动生成预应力筋
.\build\example_influence.exe .\verify_influence.dat     # 影响线与车道荷载包络
.\build\example_mixed.exe .\verify_ce.dat                # 约束方程 CE
.\build\example_mixed.exe .\verify_creep_shrinkage.dat   # 徐变收缩 AEMM
.\build\example_mixed.exe .\verify_temp_profile.dat      # 规范温度梯度
.\build\example_mixed.exe .\bridge_tied_arch_80m.dat     # 80m 综合系杆拱桥
.\build\example_influence.exe .\bridge_tied_arch_80m.dat # 综合算例影响线/包络
.\build\example_influence.exe .\verify_lane_nonlinear.dat # 非线性车道荷载直接扫描
.\build\example_stage.exe .\verify_stage.dat             # 施工阶段正装分析
.\build\example_stage.exe .\verify_stage_advanced.dat    # 阶段预应力/支座变化
.\build\example_stage.exe .\verify_stage.dat --backward  # 倒拆分析
.\build\example_stage.exe .\verify_stage_iterate.dat --iterate  # 正装-倒拆迭代
.\build\example_stage.exe .\verify_stage_multi_iterate.dat --iterate  # 多变量迭代
.\build\example_stage.exe .\verify_stage.dat --camber-iterate       # 非线性预拱度迭代
.\build\example_stage_viz.exe .\verify_stage.dat                    # 施工阶段 SVG
.\build\example_stage_envelope.exe .\verify_stage.dat              # 施工阶段包络
.\build\example_3d_truss.exe .\verify_3d_truss.dat                # 三维桁架基础
.\build\example_3d_beam.exe .\verify_3d_beam.dat                  # 三维梁基础
.\build\example_3d_stage.exe .\verify_3d_stage.dat                # 三维施工阶段
.\build\example_3d_influence_surface.exe .\verify_3d_influence_surface.dat  # 三维影响面
.\build\example_viz.exe .\bridge_tied_arch_80m.dat       # 综合算例 SVG
.\build\example_viz.exe .\model.dat       # 生成 7 张独立 SVG
```

`example_mixed` 支持 `--solver` 选择求解器（`LDLT`/`LLT`/`LU`/`QR`/`CG`/`BICGSTAB`/`AMGCL`/`CHOLMOD`），
并打印求解耗时；`--compare` 用于对比不同求解器在小/大模型上的耗时。

`example_viz` 会生成 7 张独立 SVG：`_model`（模型）、`_disp`（位移）、`_moment`（弯矩）、
`_shear`（剪力）、`_axial`（轴力）、`_stress_top`（上翼缘应力）、`_stress_bottom`（下翼缘应力）。
每张图都带节点编号、自动居中、四边留边距。模型里梁缺 `WT/WB` 时会跳过上下翼缘应力图。

解析校核算例：`verify_udl.dat`、`verify_selfweight.dat`、`verify_temp_uniform.dat`、
`verify_temp_gradient.dat`、`verify_stress.dat`；cable/strut 校验：`verify_cable_tension.dat`、
`verify_cable_slack.dat`、`verify_strut.dat`；方案 C 几何非线性校验：`verify_cable_nonlinear.dat`；
预应力校验：`verify_prestress_equiv.dat`、`verify_prestress_tendon.dat`、
`verify_prestress_ecc_tendon.dat`（端节点偏心弯矩校验）；影响线/车道荷载校验：`verify_influence.dat`。
约束方程校验：`verify_ce.dat`；徐变/收缩校验：`verify_creep_shrinkage.dat`；
温度梯度校验：`verify_temp_profile.dat`。
非线性车道荷载扫描校验：`verify_lane_nonlinear.dat`。
施工阶段校验：`verify_stage.dat`。
阶段预应力/支座变化校验：`verify_stage_advanced.dat`。
正装-倒拆迭代校验：`verify_stage_iterate.dat`。
多变量正装-倒拆迭代校验：`verify_stage_multi_iterate.dat`。

---

## 4. 示例程序清单

| 文件 | 演示内容 |
|------|---------|
| `One-dimension.cpp` | 最初的独立一维杆例子（不依赖 `fem.hpp`） |
| `example_1d_bar.cpp` | 框架版一维杆，1 自由度/节点 |
| `example_2d_truss.cpp` | 两杆三角形平面桁架，2 自由度/节点，坐标变换 |
| `example_2d_plane.cpp` | 二维平面应力 CST，单位正方形 patch test |
| `example_3d_tet.cpp` | 三维四面体 patch test，3 自由度/节点 |
| `example_3d_truss.cpp` | 三维桁架基础 |
| `example_3d_beam.cpp` | 三维梁基础 |
| `example_3d_stage.cpp` | 三维施工阶段 |
| `example_3d_influence.cpp` | 三维影响线/车道荷载 |
| `example_3d_influence_surface.cpp` | 三维影响面和 SVG |
| `example_3d_beam_live.cpp` | 三维梁活载内力/应力包络 |
| `example_3d_grillage.cpp` | 三维梁格综合算例 |
| `example_solvers.cpp` | 同一个问题跑遍各求解后端，校验结果一致 |
| `example_amgcl.cpp` | 2 万自由度二维网格，对比 LDLT / CHOLMOD / AMGCL 耗时 |
| `example_cholmod.cpp` | CHOLMOD 分阶段复用 + 多右端项 |
| `example_mixed.cpp` | 读取 `*.dat` 输入文件，求解二维梁/桁架混合结构 |
| `example_influence.cpp` | 影响线计算与车道荷载内力/应力包络 |
| `example_stage.cpp` | 施工阶段正装分析 |
| `example_stage_viz.cpp` | 施工阶段模型与位移 SVG |
| `example_viz.cpp` | 生成 SVG 可视化（模型/变形/弯矩图/应力图） |

示例数据文件（`*.dat`）：`cantilever.dat`、`mixed_truss_beam.dat`、`arch_bridge_100m.dat`，
以及五个解析校核 `verify_udl.dat` / `verify_selfweight.dat` / `verify_temp_uniform.dat` /
`verify_temp_gradient.dat` / `verify_stress.dat`。
综合算例：`bridge_tied_arch_80m.dat`。
APDL 基本模型：`bridge_tied_arch_80m_apdl.inp`。

---

## 5. 开发说明

### 5.1 整体架构

```
fem.hpp
├── 类型别名（SparseMat / DenseMat / DenseVec）
├── 自由度与组装：dofIndex / nodeDofs / scatter / assemble
├── 边界条件：applyDirichlet
├── 求解器：SolverType / solve / solveSPD
│     ├── Eigen 后端（直接法 + 迭代法）
│     ├── AMGCL（可选，__has_include 自动探测）
│     └── CHOLMOD（可选，CholmodSolver 类 + solveCHOLMOD 封装）
└── 单元：TrussElement / TriangleCST / Tetra4 / EulerBeam2D

fem_io.hpp（二维输入与求解层）
├── Model2D / readModel2D：*.dat 关键字读取
├── solveModel2D：组装 + 边界 + 求解 + 反力
├── 荷载等效节点荷载：均布/自重/温度
├── elementForces2D / elementStress2D / writeResult：内力、应力与结果输出
└── 方向解析：GX/GY/LX/LY/G<θ>/L<θ>
```

设计要点：

- **纯头文件**：除 CHOLMOD 外无需单独编译，直接 `#include "fem.hpp"`。
- **可扩展单元**：新增单元只要实现 `dofsPerNode() / numNodes() / stiffness()`。
- **可扩展求解器**：在 `SolverType` 加一项、在 `solve()` 加一个分支即可。
- **可选依赖**：AMGCL / CHOLMOD 用 `__has_include` 自动探测，缺依赖也能编译。

### 5.2 如何新增一个单元类型

新增单元的约定：

1. 定义一个结构体，保存材料参数与节点坐标。
2. 实现 `dofsPerNode()`、`numNodes()`、`stiffness()`。
3. （可选）实现 `B()`、`D()` 供后处理。
4. 用户侧用 `nodeDofs({...节点...}, dofsPerNode)` 得到全局自由度编号。

以新增“二维四节点双线性四边形（Q4，平面应力）”为例的骨架：

```cpp
struct Quad4PlaneStress {
    double E, nu, t;
    Eigen::Vector2d x0, x1, x2, x3;

    int dofsPerNode() const { return 2; }
    int numNodes() const { return 4; }

    DenseMat stiffness() const {
        // 用 2×2 高斯积分：Ke = Σ w * t * B(ξ,η)^T D B(ξ,η) * |J(ξ,η)|
        // B 由形函数对物理坐标的导数构造，形状函数：
        // N_a = 0.25*(1±ξ)(1±η)
        // ...
    }
};
```

要点：`stiffness()` 返回的矩阵行数必须等于 `numNodes() * dofsPerNode()`，
并且单元内自由度的排列顺序要和 `nodeDofs(...)` 生成的顺序一致
（默认按节点顺序、每节点内按自由度顺序展开）。

### 5.3 如何新增一个求解器后端

以接入一个新后端为例，需要改三处：

1. 在 `enum class SolverType` 中加一项：

```cpp
enum class SolverType {
    // ...
    MyNewSolver
};
```

2. 实现求解函数（或直接内联到 `solve()`）：

```cpp
inline DenseVec solveMySolver(const SparseMat& K, const DenseVec& F) {
    // 调用你的求解器，返回位移向量
}
```

3. 在 `solve()` 的 `switch` 里加分支：

```cpp
case SolverType::MyNewSolver:
    return solveMySolver(K, F);
```

如果该后端依赖外部库，建议参照 AMGCL / CHOLMOD 的做法：

- 用 `__has_include(<某头文件>)` 探测，定义 `FEM_HAS_XXX` 宏；
- 头文件不存在时在分支里抛“后端未启用”的异常；
- 在 `CMakeLists.txt` 的 `fem_add_executable` 里加上对应的头文件路径和链接库。

### 5.4 代码规范与注意事项

- 全局自由度编号统一为 `node * dofsPerNode + localDof`。
- `applyDirichlet` 会修改传入的 `K`、`F`，之后不要再复用原矩阵做别的分解。
- 单元坐标系：`TriangleCST` / `Tetra4` 的 `B()` 基于参考单元 + Jacobian 逆，
  不要求特定节点绕向，但面积/体积取绝对值。
- 求解器异常统一抛 `std::runtime_error`，上层可 `try/catch`。
- Windows 下控制台中文输出建议在 `main()` 开头调用 `setupConsoleUTF8()`。

---

## 6. 求解器选型建议

按“矩阵性质 + 规模”选择：

| 场景 | 推荐 |
|------|------|
| 小规模（<1e4 DOF）SPD | `SimplicialLDLT`（默认，简单够快） |
| 中等规模 SPD，追求稳定 | `CHOLMOD` |
| 大三维 SPD（>1e5，内存敏感） | `AMGCL`（AMG，O(n) 内存） |
| 非对称 / 一般矩阵 | `SparseLU` 或 `BiCGSTAB` |
| 欠定 / 最小二乘 | `SparseQR` |
| 同一结构多次求解、多载荷工况 | `CholmodSolver`（分阶段复用） |

说明：在本机 2 万自由度二维问题的实测中，`SimplicialLDLT` 反而比 `CHOLMOD`、`AMGCL`
更快（规模还太小，直接法开销分摊不足）；AMG 的优势要到更大三维问题才体现。
选择时应以实际规模和矩阵特性为准。

---

## 7. 已知问题与注意事项

1. **Ninja 依赖文件偶发被锁**（Windows Defender）：重跑 `cmake --build` 即可。
2. **MKL PARDISO 未接入**：conda 的 MKL 为 MSVC 编译，在 mingw64 进程里
   PARDISO 对非平凡矩阵会卡死；如需 PARDISO，建议改用 MSVC（`cl.exe`）编译。
3. **AMGCL 需要 `AMGCL_NO_BOOST`**：框架已自动定义该宏，避免 Boost 依赖。
4. **CHOLMOD 是外部动态库**：运行 `example_amgcl.exe` / `example_cholmod.exe` 前，
   需保证 `C:\msys64\mingw64\bin` 在 `PATH` 里（能找到 `cholmod.dll` 等）。
