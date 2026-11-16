# CB3D28 — 分岔点检测与分支切换（第二十一次力学④）

## 目的

补齐 B3 五项力学中的第④项：**利用路径跟踪已逐点输出但此前未使用的
`pivotRatio` 与一致切线特征值检测分岔临界点，在临界点施加屈曲模态
扰动实现分支切换，跟踪分岔后分支**。全部 opt-in，默认关闭、逐位不变。

## 物理模型

高拱 von Mises 两杆桁架（`B=1, H=3, L0=√10`），两端固定、顶点 Z 约束，
顶点受竖向参考荷载。与力学①/②的低拱（`B=2, H=0.5`）不同，此几何下
**反对称分岔先于对称极值点**发生：

| 临界类型 | 闭式条件 | 位置 |
|---|---|---|
| 反对称分岔 | `L³ = L0·h²`（`(B²+h²)^1.5 = L0·h²`） | `v_bif ≈ 0.44735`，`λ_bif ≈ 26.0152` |
| 对称极值点 | `L³ = B²·L0` | `v_sym ≈ 1.92555`，`λ_sym ≈ 82.3708` |

反对称方程有两个正根 `h ≈ 0.838 / 2.553`；路径上先遇到的是较大
`h`（较小 `v`）的根。闭式由 `build_reference.py` 独立求根生成，不读
FEM CSV 的物理数值。

## 实现挂点

- `stability_path_3d.hpp::ArcLengthOptions3D`：新增 `bifurcate / bifPivot /
  bifMode / bifAmp / bifEigTol`，默认 `0 / 0.05 / 1 / 0 / 1e-3`。
- 每个收敛点先用 `pivotRatio` 做轻量入口（实测分岔点 `≈8.9e-9`），再用
  **当前一致切线自由子块的最小特征值** `|μ_min|/max|diag(Kff)| < 1e-3`
  确认临界点（实测 `≈2.9e-4`）。
- 确认后沿 `computeBucklingModeShape3D(m, bifMode)` 的模态施加 `bifAmp`
  位移扰动（只改位移、不改 `L0`/荷载/求解器），后续弧长步从扰动构形继续。
- `ArcLengthPoint3D` / `ArcLengthResult3D` 新增 `bif*` 诊断字段，默认全 0；
  `snap_probe.cpp` 仅在 `BIFURC=1` 时输出 `S,...bif_*` 与 `M` 行。

## 关键实测

| 量 | 主路径 `bif_off` | 反对称切换 `bif_anti` |
|---|---|---|
| `λ_max` | 82.3675（≈闭式 82.3708） | 26.2353（分支承载力降 68%） |
| 切换点 `v` | — | 0.44630（闭式 0.44735） |
| 切换点 `λ` | — | 25.9567（闭式 26.0152） |
| 模态投影幅值 | 无 `M` 行 | 0.443 → 1.501（显著破缺） |
| `v_final` | +2.926（继续下降段） | −1.501（分支反转） |

负例：`bif_sym`（对称模态）模态幅值 `≈5e-17`、仍满足对称解析；
`bif_zero`（零幅值）与 `bif_badmode`（模态越界）检测到但 `bif_switched=0`，
路径 P 行与主路径去 id 后逐位一致。

## 复现

```bash
# 从源码重建并冻结（仅证据维护时需要）
_d_build/snap_probe.exe tests/benchmark/bifurcation_switch_3d/cases/cb3d28_bifurcation_switch.snap \
                        tests/benchmark/bifurcation_switch_3d/cases/cb3d28_bifurcation_switch_fem.csv
python tests/benchmark/bifurcation_switch_3d/build_reference.py

# 纯离线守卫（ctest 第 31 项同款）
python tests/benchmark/bifurcation_switch_3d/compare_bifurcation_3d.py \
  --bin _d_build \
  --cases tests/benchmark/bifurcation_switch_3d/cases \
  --reference tests/benchmark/bifurcation_switch_3d/reference
```

结果：`42/42 checks passed`。

## 边界（如实登记）

- 这是 **FEM 自洽深度**：分岔检测用一致切线特征值、切换用线性屈曲模态
  扰动，均在本项目 CR 桁架内核上闭式锚定；**不跑 ANSYS**。
- 本项目 RIKS 是 Wempner–Riks 正交型，与 ANSYS `ARCLEN` 球面型不同源，
  本项不宣称与 ANSYS 分岔后处理同源对拍。
- 多阶分岔按首次候选锁定只切换一次；`bifAmp=0` 只检测不切换。
- 全部改动 opt-in，B3 仍 `P` 不 `Y`（待五项力学全部收口后统一评估）。
