# 第 37 项 CTest：成桥态无应力曲率 benchmark（B5 第二十六次，opt-in）

## 这是什么

`B5 成桥态「无应力曲率」` 的 FEM 自洽验证。**κ0 是制造 / 设计给定的无应力曲率
（1/m），程序不反算**；本轮只验证「给一个 κ0，FEM 是否自洽」，**不接真实制造
参数、不做 ANSYS 对拍**。能力边界见 `BridgeStudio/能力边界矩阵.md` 与
`BridgeStudio/B类后续项接入规格.md` §4.1 / §6 第 28 条。

## 力学口径

常曲率 κ0 的等效固端力（局部 12 分量）：

```text
f_eq(θz1) = −EIz·κ0z     f_eq(θz2) = +EIz·κ0z
f_eq(θy1) = −EIy·κ0y     f_eq(θy2) = +EIy·κ0y      （θy = −dw/dx，同式）
```

经 `Tᵀ` 转全局后并入 `fem_3d.hpp::elementEquivalentLoad3D()`。该函数同时供
「组装 `a.F0`」与「端内力 `fl − T·fe`」使用，因此**一处接入两端自洽**：

```text
K·u = f_ext + f_eq        端内力 = K·u − f_eq = EI(κ − κ0)
```

曲率模态用 `Eeff = E/(1+χφ)` 与刚度同源，保证零外荷载时 `K·u0 = f_eq` 精确成立。

## 算例

| 算例 | 内容 | 锚点 |
|---|---|---|
| `uc01_cantilever_z` | 悬臂 + K0Z，零荷载 | `v2 = κ0L²/2`、`θz2 = κ0L`、端内力恒零 |
| `uc02_cantilever_y` | 悬臂 + K0Y | `w2 = −κ0L²/2`、`θy2 = +κ0L`、端内力恒零 |
| `uc03_fixed_z` | 两端固结 + K0Z | `U = 0`、`Mz = −EIz·κ0` |
| `uc04_fixed_y` | 两端固结 + K0Y | `U = 0`、`My = −EIy·κ0`、`T = 0` |
| `uc05_span_k0` | 简支两单元 + K0Z + 跨中 P | `v_mid = −PL³/48EI − κ0L²/8` |
| `uc06_span_base` | 同 uc05 但**无 κ0** | 默认路径零污染基线 |
| `uc09_torsion_j` | J/G 显著不同 + K0Z | `T = 0`（曲率不进扭转，J 不变） |
| `uc07_bad_truss` | TRUSS 给 K0Z | **硬负例**：显式拒绝 |
| `uc08_bad_nan` | K0Z=nan | **硬负例**：显式拒绝 |

材料：`E = 3.45e10 Pa`、`A = Iy = Iz = J = 1.0`、`G = 1e10`、`L = 10 m`、
`κ0 = 1e-3 1/m`、`P = 1000 N`。

## 关键断言（111 checks）

* 13 份冻结证据 sha256 未漂移，当前二进制逐字节复现；
* 逐节点位移 / 逐单元端内力对齐独立 Python 闭式 oracle（偏差 < 1e-9）；
* **自平衡**：悬臂零外荷载端内力恒零；
* **刚度不变**：`uc05` 与 `uc06` 端内力逐键相同（κ0 只改无应力构形）；
* **位移叠加**：跨中差恰为 `κ0L²/8`、支座转角差恰为 `κ0L/2`；
* **J 不变**：`uc09` 的 `T` 恒为 0；
* **默认零污染**：不带 `--curvature-detail` 不生成 `.curvature.txt`，无 κ0 的算例
  报告逐字节等于冻结基线，其它 benchmark 的 `reference/` 无多出的 curvature 文件；
* 两个硬负例显式拒绝且错误信息可核对；
* FE / BridgeStudio 源码契约守卫。

## 如实边界（不因本项涂改）

* **未覆盖共转列式** `CRElemData` 的初始曲率（需在其局部变形里减 `d0`，另立条目）；
* **未接真实制造参数**：κ0 的取值依据由设计 / 制造给定，程序不反算；
* **未做成桥态 `My` 全桥 ANSYS 对拍**；
* **未研究**徐变折减（`χφ`）与 κ0 的耦合口径；
* 不得把「节点预抬高 / 制造线形 `Xm`」写成「无应力曲率已实现」。

## 复现

```bash
python build_reference.py --bin <build_dir>          # 重新冻结（会校验 oracle）
python compare_unstressed_curvature_3d.py --bin <build_dir> \
    --cases ./cases --reference ./reference
```
