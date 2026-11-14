# 第 20 项 CTest：B3 `[STABILITY]` 主链路 opt-in 对标

## 这一项验证什么

CB3D17（[`../snap_through_ansys_3d/`](../snap_through_ansys_3d/)）用
ANSYS 18.0 `ARCLEN` 验证了**独立探针** `snap_probe.exe` 的球面弧长法。
本项验证第十七次把**同一内核**以 `STABILITY=1` 接入
`cable_solver.exe` 主程序之后：

1. 平衡路径与 CB3D17 冻结结果**逐位同源**；
2. `STABILITY=0`（默认）时既有成桥态闭环**零污染**；
3. 独立闭式解析锚点成立；
4. `ERROR=` 双函数在缺参数 / 无荷载 / 不支持方法时不静默。

**本项不再跑 ANSYS**：ANSYS 对标证据已在 CB3D17 归档，本项是
同源内核的主链路接入回归。

## 算例与口径

von Mises 两杆浅桁架（与 CB3D17 完全相同）：

```text
节点 1(-B,0,0) — 2(0,H,0) — 3(B,0,0)，两端全约束、顶点 Z 约束
两根 TRUSS；参考荷载 Fref 施加在顶点 -Y，幅值 PREF
L0 = hypot(B,H) = 2.061552812809 m（参考构形即无应力构形）
E=2.1e11 Pa, A=1.0e-3 m², B=2.0 m, H=0.5 m, PREF=2.0e6 N
```

`PREF` 必须**大于解析极限荷载** `P_cr≈1.20e6 N`，否则 λ=1 就到头、
看不到下降段。

独立解析锚点（闭式，不读任何求解器输出）：

```text
L(v) = sqrt(B² + (H−v)²)
N = EA·ln(L/L0)         STRAIN=1，对数应变（ANSYS LINK180 KEYOPT(2)=1）
N = EA·(L−L0)/L0        STRAIN=0，工程应变（fem_3d_nl.hpp 口径）
P(v) = −2·N·(H−v)/L
P(H) = 0,  P(2H) = 0
极值点 dP/dv = 0 由黄金分割（200 次）求得
```

`λ·PREF` 应与 `P(v)` 相等；`v = -CTRL`（段内 CTRL 是带符号位移）。

| 文件 | 说明 |
| --- | --- |
| `cases/von_mises_log.dat` | `STRAIN=1`，与 CB3D17 `log_fine` 同参 |
| `cases/von_mises_eng.dat` | `STRAIN=0`，与 CB3D17 `eng_fine` 同参 |
| `cases/von_mises_bad_ctrl.dat` | 未给 `CTRLNODE` → `ERROR=invalid_control_dof` |
| `cases/von_mises_bad_method.dat` | `PATHMETHOD=RIKS` → `ERROR=unsupported_path_method` |
| `cases/von_mises_no_load.dat` | 无 `*CLOAD` → `ERROR=zero_reference_load` |
| `cases/*_stability.txt` | 冻结的 `<stem>_stability.txt` 输出 |
| `reference/cb3d18_stable_path.json` | **3311 checks + 5 strings + 14 criteria** |

13 个 opt-in 键与段格式见 `../../input_format_3d.md` 的
「非线性稳定全过程（`[STABILITY]`）」章。`STABTOL` / `STABITER`
带前缀是因为 `*CABLE_BRIDGE` 既有的 `TOL` / `MAXITER` 属于成桥态
几何闭环，同名复用会悄悄改掉既有算例默认参数。

## 实测结果

| 项 | 数值 |
| --- | --- |
| 冻结黄金值 | 3311 checks、5 strings、14 criteria、5 cases |
| 默认零污染 | 新 exe 跑 CB3D10 既有算例，`_cable.txt` 除 `elapsed_ms` 外**逐字节一致**，无 `[STABILITY]` 段、`[BUCKLING] ERROR` 仍在 |
| FEM CLI 复现 | 5 个算例重跑，与冻结输出**逐字节一致** |
| `von_mises_log` | 274 点，`worst\|dP\|=7.407e-08 N`，`λ_cr=0.600244123`，离散 `λ_max=0.600233139`，二次插值 `λ_max_fit=0.600245230`（误差 1.11e-6），`λ_min=-0.600004`，`v_final=1.001047`，`maxRES=1.940e-03 N` |
| `von_mises_eng` | 272 点，`worst\|dP\|=5.210e-08 N`，`λ_cr=0.594209822`，`λ_max_fit=0.594205674`（误差 4.15e-6），`λ_min=-0.594189`，`v_final=1.001380`，`maxRES=2.000e-03 N` |
| 与 CB3D17 对拍 | 2 条曲线（`log_fine` / `eng_fine`）**worst\|d\| = 0.000e+00** |
| ERROR 负例 | 3 类码正确、均不产出平衡点、`[BUCKLING]` 段保留 |

容差：`|λ·PREF − P(v)| ≤ 1e-4 N`、离散 `λ_max` 误差 ≤ `1e-3`、
二次插值极值误差 ≤ `5e-4`、`RES ≤ STABTOL·|Fref| = 2.002e-3 N`、
`|v_final − 1.0| ≤ 5e-3`、与 CB3D17 路径 `worst|d| = 0`。

`maxRES` 逼近 `STABTOL·PREF = 2e-3 N` 是收敛判据
`|R| ≤ tolForce·|Fref|` 的直接结果，不是发散；相对荷载量级
`2e6 N` 为 `1e-9`。

## PIVOT 的口径（不要误读）

段内 `PIVOT = min|diag(U)| / max|diag(U)|` 是**全矩阵** LU 的对角
量级比：约束自由度被置 1、自由自由度对角量级约 `EA/L ≈ 1e8`，
因此绝对值本来就很小（本基准题约 `5.2e-9` 且沿路径近似恒定）。
**它是量级比，不是切线奇异判据**；极限点由 `LAMBDA` 的转折判定。
这一口径与 CB3D17 `snap_probe` 完全同源，未作任何修改
（改内核会破坏 L 阶段 1686 项冻结黄金值）。

ANSYS 在极限点输出的 `One indefinite tangent matrix is detected`
是切线奇异的物理标志，不是求解失败。

## 如实边界

- `STABILITY=1` 是 `cable_solver.exe` 的**独立 opt-in 模式**，在
  `readModel3D` 参考构形上跟踪路径；**不等于**成桥态闭环 / 索长迭代 /
  施工阶段链路已具备全过程稳定分析。既有主流程一字未动。
- 只有球面弧长法。位移控制 / Riks / 分支屈曲与分支切换 / 初始缺陷
  与缺陷敏感性 / 材料非线性 / CR 梁后屈曲（Williams toggle、
  Lee frame）/ 真实桥梁全过程对标**均未覆盖**。
- λ 是参考荷载放大倍数，**不是规范稳定安全系数**；JTG 3362 无
  λ 合格限值，本项不做「通过 / 不通过」判定。
- 能力矩阵该行标 `P` 不标 `Y`，边界详见
  `BridgeStudio/B类后续项接入规格.md` §2.5。

## 复现

```bash
# 生成 / 校验黄金值
python build_reference_stable.py
python build_reference_stable.py --check

# 第 20 项 CTest 入口（--bin 指向含 cable_solver.exe 的构建目录）
python compare_stable_path_3d.py \
  --bin <build_dir> --cases cases --reference reference

# 总回归
ctest --test-dir <build_dir> -R stable_path_3d_benchmark --output-on-failure
```

`build_reference_stable.py` 的 `checks` 全部来自冻结输出，
`criteria` 全部来自上面的闭式解析；Python 不组装 K、不重算 λ。
