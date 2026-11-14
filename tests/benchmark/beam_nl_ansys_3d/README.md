# CB3D16 · 大转动梁柱效应：共转列式 CR 梁 vs ANSYS BEAM188 + NLGEOM

> 状态：**8/8 ANSYS job · 0 error / 0 warning · 1776 项冻结黄金值 · 离线回归 PASS**
> 注册名 `beam_nl_ansys_3d_benchmark`（第 18 项 CTest，`TIMEOUT 900`）

## 1. 这个基准解决什么

J 阶段收尾时留了一条 TODO：「CR 内核（共转 + 解析几何刚度）已具备，
但需要 BEAM188 + NLGEOM 的独立算例」。原因是 A~J 阶段的 ANSYS 对标里：

- `nonlinear_ansys_3d`（CB3D10）用 **BEAM4 + LINK180 + NLGEOM**，主梁成桥
  位移是 mm 量级，BEAM4 不含大转动，几何非线性主要体现在索上；
- `continuous_beam_3d` 里的 BEAM188 是**纯线性**（`ANTYPE,STATIC`，无
  `NLGEOM`），只用来量化 Timoshenko 与 EB 的差异；
- `test_3d_nl.exe` 虽然已有卷圆 / 梁柱 P-Δ 的解析解断言，但模型是
  C++ 硬编码、**不读文件、不导出结果**，无法与 ANSYS 端到端比对。

本基准把 `fem_3d_nl.hpp` 的求解能力通过 `beam_nl_probe.exe` 暴露成 CSV，
与 ANSYS 18.0 的 BEAM188 + NLGEOM 做同源对标，覆盖两条独立物理：

| 组 | 物理 | 解析标尺 | 工况 |
|---|---|---|---|
| **CIRCLE** | 有限转动、随动三元组、内力平衡 | 纯弯曲 elastica 圆弧 | `λ = M/(2πEI/L)` = 0.25 / 0.50 / 0.75 / 1.00 |
| **PDELTA** | 解析几何刚度 `Kσ`、梁柱效应 | Timoshenko 梁柱精确解 | `P/Pcr` = 0 / 0.25 / 0.50 / 0.75 |

### 两个闭式

```
CIRCLE   theta = M*L/(E*Iz),  R = L/theta
         x_tip = R*sin(theta),  y_tip = R*(1-cos(theta)),  rot_tip = theta
         lambda=1 -> theta=2*pi -> 尖端回到根部（整圆）

PDELTA   Pcr = pi^2*E*Iz/L^2,  k = sqrt(P/(E*Iz))
         y_ex = Q/(2*P*k) * (tan(k*L/2) - k*L/2)
         y_lin = Q*L^3/(48*E*Iz),  放大系数 = y_ex/y_lin
```

## 2. 同源数据链

```
BridgeStudio/ansys_bridge_beamnl.py
   BNL_GROUPS / BNL_CASES（唯一数据源）
        ├── cb3d16_beam_nl.beam ──> FEM-Eigen/beam_nl_probe.exe ──> *_fem.csv
        └── 8 份 *.inp ──────────> ANSYS180.exe -b ─────────────> *_disp.csv
                                                                      │
FEM-Eigen/tests/benchmark/beam_nl_ansys_3d/                          │
   build_reference_beamnl.py  <──────── 冻结 CSV ─────────────────────┘
   reference/cb3d16_beam_nl_ansys.json  (1776 checks + 4 criteria)
   compare_beam_nl_ansys_3d.py          (CTest 入口)
```

- `beam_nl_probe.cpp` **不重新推导任何力学**：模型构造、约束、荷载与
  `test_3d_nl.cpp::testCRCantileverCircle` / `testCRBeamColumn` 一致，
  平衡求解直接调 `solveModel3DNonlinear()`（共转三元组 +
  `crBeamTangent12` 解析几何刚度 + Newton/荷载步）。
- `.beam` 为纯 ASCII：`G <gid> key=...` 定义组，`C <cid> <gid> <KIND> <lambda>`
  定义工况；`M` 与 `P` 都只从 `λ + E,Iz,L` 派生，APDL 与 CLI 用同一公式。
- APDL 纯 ASCII；`*CFOPEN` + 一次 `*VWRITE` 写整个数组（外层套 `*DO`
  会重复写首行，历史坑）；输出 7 列 `[id, UX,UY,UZ,ROTX,ROTY,ROTZ]`。

## 3. BEAM188 的口径处理（实测诊断，不是猜测）

BEAM188 是 **Timoshenko** 单元，FEM-Eigen 的 `BEAM` 是 **Euler-Bernoulli**。
默认 `BEAM188 + ASEC` 对标简支梁时，线性基线就有 `-5.07e-4` 相对差；
把 `NLGEOM` 关掉后差异**完全不变**，证明与几何非线性无关。

### 3.1 偏差分解

| 诊断 | 结果 |
|---|---|
| `NLGEOM` OFF / ON（`P=0`） | `-5.072e-4` / `-5.112e-4`，几乎相同 → 线性单元刚度问题 |
| `NSUBST` 10 / 40 / 200 | `-5.072e-4` 不变 → 非线性收敛充分，非子步问题 |
| FEM `LOADSTEPS` 6 / 20 / 40 | 结果逐位相同 → 系统保守，无加载路径依赖 |
| 网格 n=16/32/64/128/256（相对 EB 理论） | `-2.96e-3 / -5.07e-4 / +1.07e-4 / +2.60e-4 / +2.98e-4` |
| `MP,GXY` 提高 1e4 倍（n=32） | 挠度单调变硬 → 剪切柔度 ∝ 1/G |

拆成两项：

1. **横向剪切锁定**（细长梁全积分过刚），按 `1/n²` 消失；
2. **真实 Timoshenko 剪切柔度** `δ_s = Q·L/(4·G·A)`，网格加密后收敛到
   `+2.98e-4`（ASEC 截面等效 `κ≈1`）。

### 3.2 最终配置

| 组 | `KEYOPT(1,3)` | `MP,GXY` | `NELEM` | `A` | `CNVTOL` |
|---|---|---|---|---|---|
| `cir` | `0`（默认线性插值） | `8.1e10`（真实值） | 8 | 0.02 m² | `CNVTOL,F,1.0E-2,0.005,2` |
| `pd` | `2`（三次插值，消除锁定） | `8.1e13`（×1000，消除柔度） | 64 | 20 m² | 不设 |

- `KEYOPT(3)=2` 让线性基线从 `-5.07e-4` 变成 `+3.11e-4`（只剩真实剪切
  柔度）；再叠加 `G×1000` 后变成 **`+2.33e-7`**，`KEYOPT(3)=3` 甚至到
  `-6.5e-8`。选 `2` 是 ANSYS 文档中标准的 cubic 插值。
- **CIRCLE 组不需要任何处理**：纯端弯矩无横向剪力，剪切两项都不出现，
  保持默认 `KEYOPT(3)=0` 与真实 G，实测 FEM vs ANSYS 逐节点差
  `4e-10 ~ 4e-9 m`、转角 `4e-10 rad`。
- `A=20 m²`（PDELTA）把轴向压缩效应压到可忽略，使 FEM 与 Timoshenko
  「轴线不可伸长」假定一致；`n=64` 让 FEM 离散误差降到 `~6e-4`。
- CIRCLE 端弯矩加载会让 ANSYS 报一条 `MINREF` 提示（计算出的参考力
  `3.87e-10 < 默认 1e-2`）。显式 `CNVTOL,F,1.0E-2,0.005,2` 给出**同一个**
  参考力、收敛判据不变，warning 归零，尖端结果逐位不变。

> 网格收敛（FEM vs ANSYS 相对差，`P/Pcr=0.75`）：
> n=32 `-2.39e-3` → n=64 `-5.97e-4` → n=128 `-1.47e-4`，标准 `1/n²`。
> 说明残余差是 FEM 侧离散误差，不是两套单元的力学分歧。

## 4. 实测结果（ANSYS 18.0，2026-11-14）

**8 个 job 全部 `0 error / 0 warning`**（`NSUBST` 100 / 60，`NEQIT,200`，
`OUTRES,ALL,LAST`）。

| 组 | 工况 | FEM vs ANSYS `max\|ΔU\|` | `max\|ΔR\|` | `worst d/tol` |
|---|---|---|---|---|
| CIRCLE | cir025 | `3.81e-10 m` | `2.05e-10` | `6.53e-8` |
| CIRCLE | cir050 | `5.04e-10 m` | `4.94e-10` | `7.86e-8` |
| CIRCLE | cir075 | `3.77e-9 m` | `3.85e-10` | `1.55e-7` |
| CIRCLE | cir100 | `3.61e-10 m` | `4.10e-10` | `4.06e-8` |
| PDELTA | pd000 | `4.81e-9 m` | `2.04e-10` | `0.084` |
| PDELTA | pd025 | `8.75e-7 m` | `2.67e-7` | `0.142` |
| PDELTA | pd050 | `3.56e-6 m` | `1.10e-6` | `0.231` |
| PDELTA | pd075 | `2.34e-5 m` | `7.34e-6` | `0.677` |

关键判据：

- **2π 整圆闭合**：尖端相对根部 `2.630e-12 m`；尖端转角 `6.283185307`
  vs `2π`，相对误差 `2.86e-11`；四个 λ 的转角误差 `≤1.31e-10`。
- **P-Δ 放大系数**：`1.000 → 1.329 → 1.986 → 3.957`（ANSYS），
  与 Timoshenko 精确解 `1.000 → 1.329 → 1.986 → 3.958` 一致；
  ANSYS `ymid` vs 精确解最大相对差 `1.29e-4`。
- **离散误差说明**：CIRCLE 组 FEM/ANSYS 与解析圆弧差 `1e-2 ~ 4e-2 m`
  （θ=π/2 时 `1.45e-2`，θ=π 时 `4.11e-2`），这是 **8 段直弦离散**
  的 `O(θ²/n²)` 误差，FEM 与 ANSYS **共同**呈现；两者互相之间仍只有
  `1e-9 m`。θ=2π 时离散误差恰好闭合，尖端回到原点。

## 5. 冻结内容与判据

`reference/cb3d16_beam_nl_ansys.json`：

- **1776 项 `checks`** = `Σ (NELEM+1)×6` = `4×9×6 + 4×65×6`，
  每项含 `ref`（ANSYS CSV 原值）、`scale`（该 case 该 DOF 列的最大值）、
  `atol=1e-8`、`rtol=2e-3`。
- 判据：`|FEM − ANSYS| ≤ atol + rtol·scale`，`worst d/tol = 0.677`。
- **4 条解析 criteria**：整圆闭合 `<1e-6 m`、尖端转角相对误差
  `<1e-6`、P-Δ 放大系数单调且终值 `>3.0`、`ymid` vs 精确解 `<5e-3`。
- `derived_rules` 写明每个派生量的计算口径；`--check` 会从冻结 CSV
  重算全部 1776 项与 4 条 criteria，防止黄金值漂移。

`compare_beam_nl_ansys_3d.py` 离线执行四层校验：

1. 重跑 `beam_nl_probe.exe`，与冻结 `*_fem.csv` 的 **2016 个值**逐位比对
   （实测 `worst abs diff = 0`）；
2. 从冻结 CSV 重算黄金值防漂移；
3. 1776 项 FEM vs ANSYS 逐项判据；
4. 4 条解析判据 + 8 份 ANSYS 日志 `0 error / 0 warning`。

`--live-ansys` 在临时目录重跑 8 个 BEAM188 job 并逐 CSV 单元 diff，
**不覆写** `cases/` 下的冻结文件。实测 1776 个 CSV 值最大绝对差
`9.63e-9`（`pd050` 转角列，相对约 `1.6e-6`），来自 NLGEOM 的 Newton
收敛容差，因此实时复核阈值为 `1e-8 + 1e-5*|ref|`；离线路径仍以冻结
CSV 的 `2e-3` 相对判据为准。

### 命令

```bash
# 离线（无需 ANSYS 许可）
ctest --test-dir build -R beam_nl_ansys_3d_benchmark --output-on-failure
python tests/benchmark/beam_nl_ansys_3d/build_reference_beamnl.py --check
python tests/benchmark/beam_nl_ansys_3d/compare_beam_nl_ansys_3d.py \
    --bin build \
    --cases tests/benchmark/beam_nl_ansys_3d/cases \
    --reference tests/benchmark/beam_nl_ansys_3d/reference

# 实时 ANSYS 复核（8 个 job）
python tests/benchmark/beam_nl_ansys_3d/compare_beam_nl_ansys_3d.py \
    --bin build --cases ... --reference ... --live-ansys

# 重新生成算例（BridgeStudio 侧，需要 ANSYS 许可）
python cases/run_cases_beamnl.py \
    --sync-to D:/AI-work/FEM-Eigen/tests/benchmark/beam_nl_ansys_3d/cases
```

## 6. 已知边界

- **不含后屈曲 / snap-through**：FEM 内核是荷载控制 Newton + 荷载步，
  无弧长法；`P/Pcr=0.75` 仍在稳定分支内。极值点全过程、Williams
  toggle frame / Lee frame 尚未覆盖。
- **不含支座反力与单元端内力对标**：CLI 目前只导出节点位移/转角与
  解析派生量；`nonlinear_ansys_3d` 的反力/索力对标仍由 CB3D10 承担。
- **Timoshenko 差异是被显式消除而非被建模**：`G×1000` 只用于同口径
  对标，不代表真实材料；若要做 Timoshenko 差异量化，应回到
  `continuous_beam_3d` 或本目录记录的诊断数据。
- **CIRCLE 的解析圆弧误差是离散误差**，不是内核误差；若需更高绝对
  精度，增大 `NELEM` 即可（FEM 与 ANSYS 会同步收敛）。
- `BEAM188` 的 `KEYOPT(3)=2/3` 在本机 ANSYS 18.0 实测均有效且
  `0 error`，但属单元插值选项，换版本需重新复核。

## 7. 归档文件

`cases/` 下 30 个文件：`.beam` / `.model.json` / 8×(`.inp`,
`_disp.csv`, `_ansys.out`) / FEM CSV / 对比报告 / manifest / run_log。
`probes/` 归档 K 阶段物理诊断原型（剪切锁定 / 剪切柔度 / 插值阶 /
CNVTOL 四组扫描）与说明，**不进 CTest**；`reference/` 为 1776 项黄金值。
ANSYS 二进制与数据库产物由仓库 `.gitignore` 排除。
