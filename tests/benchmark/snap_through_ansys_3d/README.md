# CB3D17 — 后屈曲 / snap-through 全过程：球面弧长法 vs ANSYS ARCLEN

`snap_through_ansys_3d_benchmark`（第 19 项 CTest）验证 **FEM-Eigen
`stability_path_3d.hpp` 球面弧长法** 与 **ANSYS 18.0 LINK180 + ARCLEN**
在 von Mises 两杆浅桁架的完整 S 形荷载–位移路径上逐点一致。

离线即可全量回归，**不需要 ANSYS 许可**；`--live-ansys` 才重跑 3 个作业。

---

## 1. 为什么需要这一步

`fem_3d_nl.hpp::solveModel3DNonlinear` 是**荷载控制 Newton**：
`Ftarget = lam * Fext`，`lam = istep/nSteps` 是预设常量。在荷载–位移
曲线的极限点处切线刚度奇异，荷载控制无法越过，下降段与 snap-through
全过程都取不到。仓库 `B类后续项接入规格.md` §2 把「路径跟踪控制、λ
作为未知量的增广系统、越过极值点的下降段、独立解析锚点」列为明确缺口。

L 阶段补的是**独立 opt-in 探针**，不是把主链路悄悄改成已支持 B3：

| | 状态 |
|---|---|
| `stability_path_3d.hpp` | 新文件：Crisfield 球面弧长法，λ 为未知量 |
| `fem_3d_nl.hpp` / `NonlinearOptions3D` | **未改动**，主求解器仍是荷载控制 Newton |
| `cable_solver.cpp` / `[STABILITY]` | **未接入**，不新增 `STABILITY=1` 或 `[STABILITY]` |
| `input_format_3d.md` 主格式 | **未改**，只在附录登记探针参数文件 |
| `verify_b_class_spec.py` 守卫断言 | **仍全部通过**（主链路确实未接入 B3） |

## 2. 内核算法

复用 `fem_3d_nl.hpp` 的共转单元内力与**解析一致切线**
（`crMakeElemData` / `crElementInternal`），新增的只有增广控制：

```
R(u,lam) = lam * Fref - Fint(u) = 0
du^T du + psi^2 * dlam^2 = ds^2          (球面约束)

预测: K*duF = Fref,  dlam = sign*ds/sqrt(duF^T duF + psi^2)
校正: K*duR = R,     K*duF = Fref,  du = duR + dlam*duF
      -> a*dlam^2 + b*dlam + c = 0
      取使增量与上一收敛步增量夹角余弦最大的根（防极限点回退）
```

`ArcLengthOptions3D::strainMeasure` 控制轴向杆本构：

- `0` 工程应变 `N = EA(L-L0)/L0` —— `fem_3d_nl.hpp` 内核口径，切线 `k1=EA/L0`
- `1` 对数应变 `N = EA·ln(L/L0)` —— ANSYS LINK180 `KEYOPT(2)=1` 口径，切线 `k1=EA/L`

两者都是共转、都用解析切线，只差本构；这样把「应变度量差异」与
「弧长算法差异」严格分离。

## 3. 独立解析锚点

von Mises 两杆浅桁架：`b=2 m`、`h=0.5 m`、`A=1e-3 m²`、`E=2.1e11 Pa`，
`L0=sqrt(b²+h²)=2.061552813 m`，两端固定、顶点施加 `-Y` 参考荷载
`P_ref=2 MN`。

```
L(v) = sqrt(b^2 + (h-v)^2)
P(v) = -2 * N * (h-v) / L
```

| 应变度量 | v_cr (m) | P_cr (N) | λ_cr = P_cr/P_ref |
|---|---|---|---|
| 工程应变 `N=EA(L-L0)/L0` | 0.214246422 | 1188419.6440 | 0.594209822 |
| 对数应变 `N=EA·ln(L/L0)` | 0.215213664 | 1200488.2463 | 0.600244123 |

另有 `P(h)=0`（桁架压平）、`P(2h)=0`（完全翻转）。完整路径：
λ 从 0 → +0.6002 → 0 → −0.6002 → 0，是一条 S 形曲线。

## 4. ANSYS 口径（三处关键诊断）

**① `P_ref` 必须大于极限荷载。** 若取 `P_ref=1 MN`，ANSYS 在
`TIME=1`（即 λ=1、P=1 MN）就结束荷载步，此时 `v=0.120 m < v_cr`，
**根本没到极限点**。改成 2 MN 后 λ 的峰值变成 0.6002，之后 λ 下降、
穿过 0 变负，完整路径才被 ARCLEN 跟踪到。

**② `KEYOPT(2)` 是应变度量开关。** 实测：

| ANSYS 配置 | λ_max | 与解析对比 |
|---|---|---|
| `KEYOPT(2)=1` | 0.6002024 | 对数应变解析 0.6002441，差 4.2e-5（离散采样） |
| `KEYOPT(2)=0` | 0.6125744 | 高约 2.06%，属另一种应变度量 |

`KEYOPT(2)=1` 曲线 231 个子步与 `P_log(v)` 的**逐点最大差仅
4.594e-6 N**（相对 3.8e-10）。所以正式对标用 `=1`，`=0` 只作
应变度量对照（`sn_eng`，在报告里标 `REF` 而非 `PASS/FAIL`）。

**③ 极限点的 warning 是物理标志。** ANSYS 会提示
`One indefinite tangent matrix is detected for TIME = 0.6126` —— 切线
刚度在极限点穿过零，正是 snap-through 的证据；日志末尾
`MAX. DISP. WAS REACHED (ARCTRM COM.)` 说明按位移终止。3 个 job 均
**0 error / 1 warning**。

APDL 骨架：`LINK180` + `NLGEOM,ON` + `AUTOTS,OFF` + `ARCLEN,ON,maxarc,1e-7`
+ `NSUBST/NEQIT` + `TIME,1.0` + `OUTRES,ALL,ALL` + `ARCTRM,U,1.0,2,UY`。
逐子步用 `/POST1` + `SET,,,,,,,I` 取顶点 `U,Y`、两支座 `RF,FY` 与
`ACTIVE/SET/TIME`，一次性 `*VWRITE` 出 5 列 CSV。反力反算
λ=(R1+R3)/P_ref，与 ANSYS 内部 `TIME` 互相校验（worst 1.86e-6）。

## 5. 实测结果

| 判据 | 实测 | 限值 |
|---|---|---|
| 黄金值 | **1686 项**（3 job × 每行 3 列） | — |
| FEM CLI 重跑复现冻结 CSV | **7707 值，worst \|d\|=0** | 逐位一致 |
| ANSYS 331 子步 vs 解析 `P_log(v)` | worst **3.36e-5 N** | 1e-3 N |
| FEM 4 条路径 vs 各自解析锚点 | worst **2.66e-8 N** | 1e-3 N |
| FEM 极值点二次插值误差 | worst **4.15e-6** | 1e-3 |
| ANSYS 极值点二次插值误差 | **4.39e-7** | 1e-3 |
| FEM 弧长曲线 vs ANSYS ARCLEN 曲线（二次插值） | worst **3.35e-6** | 5e-5 |
| λ 反力反算 vs ANSYS TIME | worst **1.86e-6** | 1e-5 |
| 下降段 | λ_min_fit=−0.6002445，v_max=1.0026 m | λ<−0.5、v≥0.99 |
| 应变度量对照 | `KEYOPT(2)=0` 高 **+2.061%** | 1.5%~3.0% |

FEM 离散平衡点未必恰好落在极值点，故用相邻三点二次插值求
λ*；`DSMAX=0.01` 时插值误差降到 1e-6 量级。FEM↔ANSYS 曲线互比同样
用二次插值：线性插值会引入 2.47e-4 的纯插值误差，二次插值把它压到
3.3e-6，判据才真正反映两套求解器的一致性。

## 6. 文件与运行

`cases/` 下 15 个文件：`.snap` / `.model.json` / 3×(`.inp`,
`_curve.csv`, `_ansys.out`) / FEM CSV / 对比报告 / manifest / run_log。
ANSYS 二进制与数据库产物由仓库 `.gitignore` 排除。`probes/` 归档
物理诊断脚本（弧长法原型、ARCLEN 参数/应变度量扫描）。

```bash
# 离线 CTest（无需 ANSYS 许可）
python tests/benchmark/snap_through_ansys_3d/compare_snap_through_ansys_3d.py \
  --bin build \
  --cases tests/benchmark/snap_through_ansys_3d/cases \
  --reference tests/benchmark/snap_through_ansys_3d/reference

# 黄金值自洽校验 / 重生成
python tests/benchmark/snap_through_ansys_3d/build_reference_snap.py --check
python tests/benchmark/snap_through_ansys_3d/build_reference_snap.py

# 实时 ANSYS 复核（临时目录重跑 3 个 job，不覆写冻结 CSV）
python .../compare_snap_through_ansys_3d.py ... --live-ansys

# 重新生成输入并重跑全部作业（BridgeStudio 侧）
python D:/AI-work/BridgeStudio/cases/run_cases_snap.py \
  --sync-to D:/AI-work/FEM-Eigen/tests/benchmark/snap_through_ansys_3d/cases
```

`--live-ansys` 在临时目录重跑并逐 CSV 单元 diff，**不覆写** `cases/`。
ANSYS 输出为 `E20.12`，复核容差 `1e-6 + 1e-10*|ref|`。

## 7. 边界

- 本 benchmark 是**独立单元级探针**，不表示 `cable_solver.cpp` 主链路
  已支持 B3 非线性稳定全过程；`[BUCKLING] ERROR=nonlinear_state_not_supported`
  保持不变。
- 弧长求解器目前只验证共转**轴向杆**；CR 梁的 snap-through
  （Williams toggle / Lee frame）、分支/分岔处理、初始缺陷敏感性、
  材料非线性均未覆盖。
- 切线奇异指标 `pivotRatio` 已逐点输出，但尚未作为分支判据接入。
