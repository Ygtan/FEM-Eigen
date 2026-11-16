# 第 28 项 CTest：位移控制 / 线搜索（B3 第二十一次力学①）

> 建立：2026-11-16。第二十一次五项力学分阶段推进的第①项。

## 能力口径

`stability_path_3d.hpp::ArcLengthOptions3D` 新增 opt-in 字段
`controlMode` / `ctrlDispStep` / `lineSearch`（默认 0/0/0，既有输出逐位不变）：

- **位移控制**（`controlMode=1`）：每步把控制 DOF 相对步起点的位移增量
  锁定为 `ctrlDispStep`（符号在第一步按荷载响应方向锁定，越过极限点后
  不变），λ 由线性约束 `du(ctrl)=Δ` 唯一确定：
  `dλ = (Δ − base(ctrl)) / duF(ctrl)`。这是越过极值点 / snap-through 的
  经典位移控制，天然不依赖 λ 单调。
- **线搜索**（`lineSearch=1`）：Armijo 回溯，对整步校正缩放 `alpha·duFull`
  与 `alpha·ddFull`；因 `duAcc(ctrl)` 恒等于目标增量，约束对任意 `alpha`
  成立，线搜索只改变收敛稳健性、**不改变平衡路径**。正常步 `alpha=1`
  立即接受，0 次回溯，与不开线搜索逐位一致。

入口 `snap_probe.exe` 的 `*.snap` 在 `C` 行新增 opt-in 键
`CTRLMODE` / `LINESEARCH` / `CTRLSTEP` / `CTRLNODE` / `CTRLDOF`
（默认 0/0/0/2/2，既有 `CB3D17` `.snap` 逐位不变）。

## 物理模型与独立锚点

von Mises 两杆浅桁架（与 CB3D17 同源）：节点 1(-B,0)–2(0,H)–3(B,0)，
两端固定、顶点 Z 约束，顶点 -Y 向参考荷载 `PREF`。闭式解析（oracle 不读
求解器）：

```
L(v)  = sqrt(B² + (H−v)²)
N(v)  = EA·(L−L0)/L0,   L0 = sqrt(B²+H²)
P(v)  = −2·N·(H−v)/L
λ(v)  = P(v)/PREF
```

极限点由黄金分割搜索 `P(v)` 最负处独立求得（`λ_cr = 0.5942098220`）。

## 冻结算例

`cases/cb3d25_displacement_control.snap` 四个 case：

| case | 配置 | 期望 |
|---|---|---|
| `dc_plain` | 位移控制，Δ=0.01，无线搜索 | 收敛至 v=1，100 步，0 回溯 |
| `dc_ls` | 同 `dc_plain` + 线搜索 | 同路径（逐点 \|Δλ\|≈0），有回溯被触发 |
| `dc_ls_fine` | Δ=0.005 + 线搜索 | 200 步，λ_max 更接近 λ_cr |
| `dc_badctrl` | 控制 DOF=节点1（固定端） | `singular_control_dof`，不收敛、无路径推进 |

## 判据（25 项）

- 当前二进制重跑 `_fem.csv` 逐字节一致；冻结 CSV sha256 = reference；
- 三条正常路径：收敛、走到 v=1、闭式 `P(v)` 逐点误差 ≤ 1e-6 N、
  `worst_dP_rel ≤ 1e-12`；
- 线搜索与无线搜索逐点 `|Δλ| ≤ 1e-10`；无线搜索 0 回溯 / alpha=1，
  线搜索确有回溯；
- 粗步 / 细步 `λ_max` 与闭式 `λ_cr` 分别在 5e-4 / 1e-4 内，细步更接近；
- 负例报错且不收敛、无路径推进。

## 复现

```bash
cd D:\AI-work\FEM-Eigen
cmake --build _d_build --target snap_probe
python tests/benchmark/displacement_control_3d/compare_displacement_control_3d.py \
  --bin _d_build \
  --cases tests/benchmark/displacement_control_3d/cases \
  --reference tests/benchmark/displacement_control_3d/reference
# 重建黄金值：
python tests/benchmark/displacement_control_3d/build_reference.py
```

## 边界

- opt-in 研究能力，默认全关；`STABILITY=0` / 不开新键的既有输出逐位不变。
- RIKS / 位移控制 / 线搜索与 ANSYS `ARCLEN`（球面弧长）参数化不同源，
  本项不跑 ANSYS、不代替球面型 ANSYS 证据。
- 尚未接入 `cable_solver.exe` 的 `readModel3D` 成桥态主链路，只是
  `snap_probe` 独立探针层；这一项完成**不把 B3 能力矩阵提前涂成 `Y`**。
- 其余四项（初始缺陷 / 材料非线性 / 分支屈曲 / 动力 snap-through）仍未实现。
