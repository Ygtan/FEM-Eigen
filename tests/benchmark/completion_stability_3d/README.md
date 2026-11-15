# 第 22 项 CTest：`completion_stability_3d`

**第十九次（2026-11-15，a 阶段 ①）** —— 把球面弧长法平衡路径的起点从
`readModel3D()` 参考构形换成**成桥态收敛位移 `u0` + 恒定荷载 `F0`**，
并在 `cable_solver.exe` 主链路结果文件里输出
`[STABILITY] BASIS=completion_state`。本项**不跑 ANSYS**：真实桥梁算例的
ANSYS 两段对标属阶段 ②，能力矩阵 B3 该行**仍 `P` 不标 `Y`**。

## 与第 19 / 20 / 21 项的关系

| 项 | 基准 | 起点 | λ 放大的量 | 求解入口 |
|---|---|---|---|---|
| 19（CB3D17） | `reference_configuration` | `u=0` | 总参考荷载 | `snap_probe.exe`（独立探针） |
| 20（CB3D18） | `reference_configuration` | `u=0` | 总参考荷载 | `cable_solver.exe STABILITY=1` |
| 21（CB3D19） | `reference_configuration` | `u=0` | 总参考荷载 | `frame_probe.exe`（CR 梁独立探针） |
| **22（本项）** | **`completion_state`** | **成桥态 `u0`** | **只放大 `*STABILITY_LOAD` 增量荷载** | **`cable_solver.exe STABASIS=COMPLETION`** |

两种 `BASIS` 的 λ **不是同一个量**：参考构形基准 `R = λ·Fref − Fint(u)`；
成桥态基准 `R = F0 + λ·Finc − Fint(u)`，`F0` 不随 λ 缩放。报告与计算书
必须分列，不得互相冒充。

## 力学模型与独立闭式锚点

V 形双索单自由度体系：两根 `CABLE` 从 `(±B, H)` 连到低点 `(0,0,0)`，
低点只留 dof 2（y）。参考构形带预拉（`L0 = 0.995·hypot(B,H) < L_ref`），
避免 `CABLE` tension-only 在起点松弛导致零切线。恒载 `F0` 向下，
成桥态低点下移 `V0`。

**不读求解器的闭式解**（`build_reference_completion.py` 独立实现）：

```
L(v)  = hypot(B, H+v)
N(v)  = EA/L0 · (L(v) − L0)            工程应变，与 fem_3d_nl.hpp 同口径
P(v)  = 2 N(v) (H+v) / L(v)            两根对称索的竖向合力
F0    = P(V0)                          成桥态平衡
T*    = N(V0)                          MODE=TENSION 的目标索力
λ(v)  = [P(v) − F0] / PINC             成桥态基准的增量荷载因子

K_yy(v) = 2[(EA/L0)c_y² + (N/L)(1−c_y²)],  c_y = (H+v)/L
dλ/dv|V0 = K_yy(V0)/PINC
```

`K_yy` 与中心差分 `dP/dv` 是**两套独立推导**，实测相对差
`3.964e-10`；路径第一点割线斜率 `51.156` 对解析 `51.150`（差 `1.3e-4`）。

参数：`B=2.0 m`、`H=0.6 m`、`E=2.1e11 Pa`、`A=1e-3 m²`、
`L0=2.07762099527 m`、`V0=0.10 m`、`F0=2.76082556144e6 N`、
`T*=4.17863177218e6 N`、`PINC=5.0e5 N`、
`CTRLNODE=2 CTRLDOF=2 CTRLTARGET=0.05`（**相对 `u0` 的增量**）、
`DS0=0.002 DSMAX=0.006 PSI=0.5 STABTOL=1e-9 STRAIN=0`。

## 冻结 cases（6 个）

| stem | MODE | 输出 | 预期 |
|---|---|---|---|
| `vshape_static` | `STATIC`（`.dat` 显式给 `L0`） | `_cable.txt` | `BASIS=completion_state`，241 个平衡点 |
| `vshape_tension` | `TENSION`（索长迭代命中 `T*`） | `_cable.txt` | 同上，另含 `[LENGTH] KIND=cable_length STATUS=ok` |
| `vshape_no_load` | `STATIC` | `_cable.txt` | `ERROR=zero_incremental_load`，无平衡点 |
| `vshape_bad_basis` | `STABASIS=FOO` | `_stability.txt` | `ERROR=unsupported_stability_basis`，`BASIS=unknown` 不伪装成参考构形 |
| `vshape_unconverged` | `COMPLETION`（纯索无梁，几何闭合不收敛） | `_cable.txt` | `ERROR=completion_not_converged` |
| `vshape_default` | `STATIC`，有 `*STABILITY_LOAD` 但无 `STABILITY` 键 | `_cable.txt` | 无 `[STABILITY]` 段，`[U](2,2) = −V0` |

`vshape_unconverged` 是**负例**：纯索体系没有梁，`MODE=COMPLETION` 的
制造线形几何闭合本来就不收敛；它只用来验证「未收敛位移不得当路径起点」，
不是能力声明。真实斜拉 / 悬索桥的 `COMPLETION` 闭环另有 CB3D10~CB3D13。

## 黄金值与实测

```
2963 checks + 18 strings + 18 criteria + 6 cases
FEM CLI 复现 6/6（_cable.txt 除 elapsed_ms 逐行一致；_stability.txt 逐字节）
独立闭式 worst|λ−λ_closed| = 1.505e-13 (STATIC) / 3.515e-08 (TENSION)
          worst|λ·PINC−(P(v)−F0)| = 7.5e-08 N / 1.76e-02 N
起点自洽 CTRL0 ≡ [U](2,2)；F0_NORM ≡ |F0|；FINC_NORM ≡ PINC
          RES0/|F0| = 1.2e-14 (STATIC) / 5.5e-09 (TENSION)
默认零污染 既有 CB3D10 1257 行逐位一致（除 elapsed_ms）
          M 阶段 von Mises REFERENCE 基准 _stability.txt 逐字节一致
          *STABILITY_LOAD 不进 m.loads：默认态 U(2,2) 仍是 −0.1
```

`_cable.txt` 头部有 `elapsed_ms` 运行耗时，按项目既有约定归一化后比对；
`[STABILITY]` 全段仍逐字节一致。`_stability.txt` 无计时字段，直接逐字节。

## 复现

```bash
# 离线 CTest（不需要 ANSYS 许可）
python compare_completion_stability_3d.py \
    --bin <build_dir> --cases cases --reference reference

# 重新生成黄金值（会覆写 reference/*.json，正常回归不要跑）
python build_reference_completion.py
python build_reference_completion.py --check
```

`cases/*.dat` 与 `cases/*_cable.txt` 由同一份 `cable_solver.exe` 生成，
没有第二份手工输入；本项不涉及 BridgeStudio 导出与 ANSYS `.inp`。

## 如实边界

- ① 阶段只交付 **FEM 侧自洽 + 独立闭式锚点**；没有 ANSYS
  「成桥态等效初应变 + 两段求解」对标，没有施工阶段链路接入，
  没有 CR 梁后屈曲的主链路版本。
- `BASIS=reference_configuration`（第十七次）与 `STABILITY=0`
  的输出**逐位不变**；`fem_3d_nl.hpp` 仍是荷载控制 Newton、
  `NonlinearOptions3D` 默认值零改动，
  `[BUCKLING] ERROR=nonlinear_state_not_supported` 原样保留。
- `[STABILITY] STATE` 行的 `RES0` 是起点闭合残差 `|F0−Fint(u0)|`，
  **不是补 0 的占位**；`CTRL` 是相对成桥态的位移增量，不是绝对位移。
- 验收三阶段表与四条硬口径见
  `BridgeStudio/B类后续项接入规格.md` §2.6；纯离线守卫
  `BridgeStudio/verify_completion_stability_3d.py`（32/32）。
