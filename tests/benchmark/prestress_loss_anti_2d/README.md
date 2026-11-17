# 第 33 项 CTest：σl2 附录 G 反摩擦（B4 第二十三次）

> 建立：2026-11-16。用户「同意建议的顺序」后，B4 第一项 σl2 的冻结证据与离线守卫。

## 能力口径

`prestress_loss.hpp::PrestressLossInput` 新增两个 opt-in 键，默认 `L2ANTI=0` / `L2REPX=0`：

- `L2ANTI=1` 时 σl2 改按 JTG 3362-2018 **附录 G** 曲线筋反向摩擦计算；
- 代表截面只支持 **张拉端 `x=0` 峰值**（用户 2026-11-16 拍板：规范无统一代表值规定，
  张拉端保守且与 G.0.2-4 闭式一致），`L2REPX≠0` 一律显式拒绝；
- `L2ANTI=0` 时仍为原式 `σl2 = ΣΔl/l·Ep`，**既有数值行逐位不变**。

内核分流（`computePrestressLosses`）：

```
σl1(l) = σcon·[1 − exp(−(μθ + k·l))]        # 锚固端全长 l，不是计算截面 x
Δσd    = σl1(l)/l                            # Pa/m
lf     = sqrt(ΣΔl·Ep / Δσd)                  # SI 同量纲：sqrt(m·Pa/(Pa/m)) = m

lf ≤ l :  Δσ = 2·Δσd·lf                      # G.0.2-4，BRANCH=1
lf > l :  Δσ′= Δσd·(lf²/l + l)               # G.0.2-5 等腰梯形面积协调，BRANCH=2
```

`lf=l` 两支连续（独立验算 `175581777.17704442 Pa`）。`Δσd=0`（无摩擦）时反向摩擦
无定义，退化为直线式并把 `BRANCH` 记为 0 —— 这是物理边界，不是缺参兜底。
`ΣΔl≤0` / `Ep≤0` / `L2REPX≠0` 一律 throw，不静默补值。

`PrestressLossResult` 新增 `lf` / `deltaSigmaD` / `deltaSigmaPeak` / `l2Branch`；
2D 明细段 `prestress_loss_output_2d.hpp` 仅在 `l2Branch≥1` 时向数据行追加
`LF=` / `DSL2D=` / `DSL2PEAK=` / `BRANCH=`，12 位有效数字不变。
3D 侧由新增只读探针 `loss_probe_3d.exe` 暴露同一组结果（只调
`readModel3D` + `resolvePrestressLosses3D`，不自带力学）。

## 冻结算例

| case | 配置 | 期望 |
|---|---|---|
| `anti_default` | 不带 `L2ANTI`（22 键原样） | 直线式，无反摩擦字段；与仓库既有 `verify_prestress_loss.dat` 明细段逐字节一致 |
| `anti_le` | `L2ANTI=1`，`ΣΔl=0.00113` | `LF≈5.01 m < l`，`BRANCH=1`，`Δσ=2Δσd·lf` |
| `anti_gt` | `L2ANTI=1`，`ΣΔl=0.018` | `LF≈20.0 m > l`，`BRANCH=2`，`Δσ′=Δσd(lf²/l+l)` |
| `anti_nodl` | `L2ANTI=1` 但 `DELTA_L=0` | 退出码非 0，报「附录 G 反摩擦需要 ΣΔl>0」，不产出明细 |
| `anti_badrepr` | `L2ANTI=1` + `L2REPX=1` | 退出码非 0，报「仅支持代表截面 x=0」 |
| `anti_nofiric` | `L2ANTI=1` 且 `MU=K=0` | 无摩擦退化，`BRANCH=0`，σl2 仍为直线式，不追加字段 |
| `anti_3d` | 3D 模型，与 `anti_le` 同一组 `*PRESTRESS_LOSS` 参数 | 2D/3D 共 13 键逐值同值 |

## 判据（117 项）

- 7 份 `.dat`、5 份冻结证据的 sha256 与 reference 对齐；
- 当前二进制重跑逐字节复现冻结证据；
- 5 条正常路径逐键对齐**独立 Python 闭式 oracle**（从 `.dat` 读参数，不读求解器输出）：
  σl1~σl6 / σpe / RATIO / LF / DSL2D / DSL2PEAK / BRANCH；
- 单位量纲锚点：σcon、σl2、Δσd 为 Pa / Pa/m，lf 在 0.1~100 m，NOTE 写明单位；
- 支闭式交叉、`lf=l` 连续性（不读求解器）；
- 两个负例显式拒绝且不产出明细文件；
- 默认路径零污染（逐字节）与未启用时数据行无反摩擦字段；
- 2D 明细段与 3D 探针 13 键同值；3D 探针源码不自带力学。

## 复现

```bash
cd D:\AI-work\FEM-Eigen
ninja.exe -C _d_build -j 2 example_mixed.exe loss_probe_3d.exe
python tests/benchmark/prestress_loss_anti_2d/compare_prestress_loss_anti_2d.py \
  --bin _d_build \
  --cases tests/benchmark/prestress_loss_anti_2d/cases \
  --reference tests/benchmark/prestress_loss_anti_2d/reference
# 重建黄金值（会覆盖 cases/*.loss.txt 与 reference JSON）：
python tests/benchmark/prestress_loss_anti_2d/build_reference.py --bin _d_build
```

## 边界

- opt-in 研究能力，默认 `L2ANTI=0`；既有 `.dat`、`.loss.txt` 数值行与冻结证据逐位不变。
- 只实现附录 G 的**代表截面峰值**，不输出 `σl2(x)` 沿程分布；非 x=0 代表截面显式拒绝。
- `σl6` 仍为用户显式录入规范表值，本项不接龄期自动发展、不调用拟合函数兜底。
- `anti_gt` 分支中 `σl5 = -7.37e6 Pa` 是 JTG 6.2.6-1 在低 σpe 下的数学行为
  （`0.52σpe/fpk − 0.26 < 0`），oracle 与内核同式一致；不代表容许应力或安全系数判定。
- 本项**不跑 ANSYS**：附录 G 反摩擦目前只有 FEM 闭式自洽证据，未与 ANSYS 实体单元/
  接触摩擦对拍；不把 B4 能力矩阵涂成 `Y`。
- 3D 探针只验证同一份参数在 2D/3D 共用内核算出同一组分项值，不代表 3D 用户报告层
  已新增独立的 `[PRESTRESSLOSS3D]` 段。
