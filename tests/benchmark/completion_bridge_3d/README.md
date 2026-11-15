# 第 23 项 CTest：`completion_bridge_3d`

**第十九次（2026-11-15，a 阶段 ②）** —— 把第 22 项的
`STABASIS=COMPLETION` 从 V 形双索单自由度解析题放到**真实双塔三跨
斜拉桥**（50 节点 / 63 单元 / 16 索），并与 ANSYS 18.0
「step1 恒载成桥态 + step2 增量荷载」两段求解逐子步对拍。
CTest 本身**完全离线**（冻结 CSV + 重跑 `cable_solver.exe`），
不需要 ANSYS 许可；`--live-ansys` 才临时重跑 2 个 job。
② 完成后 B3 能力矩阵**仍 `P` 不标 `Y`**（③ 才推 `Y`）。

## 与第 20 / 21 / 22 项的关系

| 项 | 体系 | 起点 | 对标对象 | 求解入口 |
|---|---|---|---|---|
| 20（CB3D18） | von Mises 两杆浅桁架 | `u=0` 参考构形 | CB3D17 冻结路径 | `cable_solver.exe STABILITY=1` |
| 21（CB3D19） | Williams toggle CR 梁框架 | `u=0` 参考构形 | `BEAM188 + ARCLEN` | `frame_probe.exe` |
| 22（CB3D20） | V 形双索单自由度 | 成桥态 `u0` | 独立闭式 `λ(v)`（无 ANSYS） | `cable_solver.exe STABASIS=COMPLETION` |
| **23（本项）** | **真实双塔三跨斜拉桥** | **成桥态 `u0`** | **ANSYS 两段求解（荷载控制 + ARCLEN）** | **同一 `cable_solver.exe STABASIS=COMPLETION`** |

第 22 项证明「力学自洽」，本项证明「真实桥梁尺度下与商业软件一致」，
两者**不能互相替代**：第 22 项没有 ANSYS，本项路径单调硬化、无极限点。

## 同源导出与物理口径

唯一数据源：`BridgeStudio/ansys_bridge_completion.py` 的
`COMPLETION_CFG + EPS0 + PINC + LAMMAX`，由
`BridgeStudio/cases/run_cases_completion.py` 编排：

```
ModelSpec3D(mode=STATIC, 设计线形)
  ├─ *.dat   : *EN ... L0 = L_design/(1+EPS0)
  │            + STABASIS=COMPLETION 键
  │            + *STABILITY_LOAD (跨中 1016 节点 Z 向 -PINC)
  └─ *.inp   : 同一份设计线形坐标
               + INISTATE eps = (L_design-L0)/L0 = EPS0（16 根索）
               + step1 恒载静力 TIME=1
               + step2 F_orig -> F_orig - LAMMAX*PINC，TIME=2 / 6
```

| 键 | 值 |
|---|---|
| `EPS0` | `8.0e-4`（统一初始应变） |
| `PINC` | `1.0e6 N`（FEM 增量参考荷载） |
| `LAMMAX` | `4.55`（ANSYS step2 目标；FEM `λ_final=4.5530675`） |
| `CTRL` | `(1016, 3)`，跨中节点 Z 向位移 |
| `VMAX` | `0.07 m`（ARCTRM 保护阈值，两 job 都靠 TIME 终止） |

FEM 侧 `R = F0 + λ·Finc − Fint(u)`，起点是 `MODE=STATIC` 的成桥态收敛
位移 `u0`；`[STABILITY] STATE` 行给出 `CTRL0 / F0_NORM / FINC_NORM / RES0`。

ANSYS 侧**不用反力换算 λ**：完整斜拉桥 50 个 `F` / 34 个约束节点上
`ΣRF,FZ` 在 step2 随加载反而下降，反力求和不可用。改为

```
lambda_FEM = (TIME − 1) * LAMMAX        # 仅荷载控制主 job cs_load
DeltaUZ    = UZ − UZ_step1              # 与 FEM [STABILITY] 的 CTRL 对拍
```

`cs_arc` 的 `TIME` 是 ARCLEN 的上界、不是 λ，只做**终点位移交叉验证**。

**为什么用 `MODE=STATIC` 而不是 `COMPLETION`**：`MODE=COMPLETION` 的索长
迭代会改写 `m.elements[*].L0`，与收敛位移 `u0` 不同源，实测 CB3D10
`RES0=3.15 MN`（与 `F0_NORM` 同量级）。`STATIC` + 显式 `L0` 时
`RES0=7.97e-3 N`，与 `[NEWTON] RESLAST=5.96e-3 N` 同量级，起点严格闭合。
该缺陷登记为**独立后续项**，不是本算例的能力声明。

## 冻结 cases（15 个文件）

| 类别 | 文件 | 内容 |
|---|---|---|
| 唯一数据源 | `*.cfg.json` / `*.model.json` / `*.dat` | 配置常量、模型快照、FEM 输入 |
| FEM 主链路 | `*_cable.txt` | `[STABILITY] BASIS=completion_state`，246 个平衡点 |
| 默认负例 | `cs_default.dat` / `cs_default_cable.txt` | 无 opt-in 键、无 `[STABILITY]` 段 |
| ANSYS 主 job | `cs_load.inp` / `_curve.csv` / `_ansys.out` | 荷载控制 step2，40 个子步（step1 12 + step2 28） |
| ANSYS 交叉验证 | `cs_arc.inp` / `_curve.csv` / `_ansys.out` | `ARCLEN,ON`，45 个子步 |
| 结论 | `*_compare_report.txt` / `_manifest.txt` / `_run_log.txt` | 15 项判据全 PASS |

`*.inp` 里的 `INISTATE,DEFINE,...,EPEL,X,eps` 共 16 条，逐条与 `EPS0`
同源；`F,1016,FZ` 恰好 2 条（step1 原始恒载、step2 目标荷载）；
`SOLVE` 恰好 2 条 —— 这三项是「同一参考构形、两段求解」的硬证据，
CTest 逐条核查，防止拿参考构形 `ARCLEN` 冒充成桥态两段。

## 黄金值与实测

```
2392 checks + 12 strings + 35 criteria + 4 cases
FEM CLI 复现 2/2（_cable.txt 除 elapsed_ms 逐行一致）
同源锚点  16 根索 L0*(1+EPS0) - 设计弦长      : dev = 0
          16 条 INISTATE eps - EPS0            : dev = 1.830e-16
          APDL F_end - (F_orig-LAMMAX*PINC)   : dev = 3.306e-07 N
起点自洽  RES0=7.973e-3 N ; RES0/|F0|=2.696e-09 ; RES0/RESLAST=1.338
          STATE CTRL0 ≡ [U](1016,3)，差 3.06e-14
          FEM CTRL0 vs ANSYS step1 UZ         : |d|=1.736e-06 m
FEM 路径  246 点，λ 0 -> 4.5530675，单调硬化、无极限点
逐子步    28 子步 worst|DeltaUZ - FEM CTRL|   : 1.852e-05 m (tol 1e-4)
          初始斜率 k_ansys=9.0972e7 / k_fem=9.0985e7 N/m，rel=1.374e-4
终点交叉  ARCLEN vs 荷载控制 UZ               : |d|=6.901e-08 m
ANSYS 日志 cs_load 0 error / 1 warning；cs_arc 0 error / 2 warning
默认零污染 cs_default 无 [STABILITY]；既有 CB3D10、M 阶段 von Mises
          REFERENCE 输出除 elapsed_ms 外逐位不变
报告结论  RESULT: ALL CHECKS PASSED
```

`_cable.txt` 头部有 `elapsed_ms` 运行耗时，按项目既有约定归一化后比对；
`[STABILITY]` 全段仍逐字节一致。

## 复现

```bash
# 离线 CTest（不需要 ANSYS 许可）
python compare_completion_bridge_3d.py \
    --bin <build_dir> --cases cases --reference reference

# 临时目录重跑 2 个 ANSYS job，逐 CSV 单元 diff，绝不写回冻结 cases/
python compare_completion_bridge_3d.py \
    --bin <build_dir> --cases cases --reference reference --live-ansys

# 重新生成黄金值（会覆写 reference/*.json，正常回归不要跑）
python build_reference_bridge.py
python build_reference_bridge.py --check
```

冻结证据由 `BridgeStudio/cases/run_cases_completion.py` 生成，用
`--sync-to` 同步到本目录 `cases/`；`--sync-only` 只做拷贝、不重跑求解器。

## 如实边界

- ② 只交付 **真实斜拉桥 + ANSYS 两段对标（荷载控制 + ARCLEN 交叉验证）**；
  B3 推 `Y` 属阶段 ③，本项完成后能力矩阵**仍 `P` 不标 `Y`**。
- 本算例路径**单调硬化、无极限点**，ANSYS 荷载控制即可覆盖全程；
  「越过极限点 / 下降段 / 反向荷载」能力由 L（CB3D17）、
  N（CB3D19）阶段独立基准题与 FEM 球面弧长内核证明，不由本项证明。
- `BASIS=completion_state` 的 λ 只缩放 `*STABILITY_LOAD` 增量荷载，
  **不是规范稳定安全系数**；`reference_configuration` 的 λ 是另一个量，
  两者不得互相冒充（第 22 项 §2.6 已立口径）。
- `MODE=COMPLETION` 索长迭代改写 `L0` 导致 `u0` 与 `m` 不同源是
  **已知缺陷**，本算例用 `MODE=STATIC` 绕开；施工阶段链路接入、
  分支屈曲 / 初始缺陷 / 材料非线性仍未覆盖。
- 验收三阶段表见 `BridgeStudio/B类后续项接入规格.md` §2.6；
  纯离线守卫 `BridgeStudio/verify_completion_stability_3d.py`（47/47，
  ② 阶段扩展 `stage2_contract()` 后由 32/32 增至 47/47）。
