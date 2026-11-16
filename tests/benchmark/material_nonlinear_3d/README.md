# 第 30 项 CTest：材料非线性本构接入 3D CR 梁（力学③）

第二十一次分阶段推进的力学③：`section.hpp` 的 `materialStress` /
`sectionResponse`（STEEL 理想弹塑性、CONCRETE 受拉开裂）从孤岛接入
`fem_3d_nl.hpp::crElementInternal` 的 BEAM 分支，使 CR 梁在 opt-in 时
具备纤维截面材料非线性与一致材料切线。

## 边界（默认零污染）

- `Element3D` 新增 opt-in 字段 `matKind`（默认 `ELASTIC`）、`matFy`、
  `matFt`、`matNfy`、`matNfz`；`.dat` 通过 `*R ... MATKIND=STEEL FY=...`
  等键请求。
- `crMakeElemData` 仅在 `type==BEAM` 且 `matKind∈{STEEL,CONCRETE}` 时
  构造 `FiberNlData`；其他一律 `enabled=false`，逐位走原 `Kl12*d`。
- `crElementInternal` 的 BEAM 分支在 `enabled=false` 时与改动前**逐字节
  同源**；全部既有 29 项 CTest（含第 28/29 项力学①/②证据）未改动。
- 扭转 `q[3]=GJ/L·(θx1−θx2)` 保持线弹性，纤维本构不覆盖扭转。

## 纤维离散：Gauss-Legendre（关键设计）

纤维位置/面积取 `nfy×nfz` 点 Gauss-Legendre 规则而非等分中点：

```text
h = sqrt(12·Iz/A),  b = sqrt(12·Iy/A)
y_i = (h/2)·g_i,  w_i = 纤维面积权重 (h/2)·gw_i
z_j = (b/2)·g_j,  area_ij = (h/2·gw_i)·(b/2·gw_j)
```

Gauss 对 `y²,z²` 二阶多项式精确积分 ⇒ 离散 `A/Iy/Iz` 与解析矩形截面
**逐位一致**，弹性段 `sectionResponse` 与 `Kl12·d` 完全对齐。
（中点离散会有 O(1/n²) 误差：8×4 网格下 Iz 偏 1.6%，已弃用。）

`crFiberLocalTangent` 用中心差分 `∂q/∂deff`，自动覆盖 D12/D13/D23
耦合与弹塑性加载/卸载活动集；端点剪力由梁段平衡
`q[1]=(q5+q11)/L`、`q[2]=−(q4+q10)/L` 补齐，与 Kl12 严格相等。

## 输入与输出

```bash
mat_nl_probe.exe cases/cb3d27_material_nl.mat out.csv
```

- `G <gid> ...` 定义截面/材料（MATKIND、FY、FT、NFY、NFZ）；
- `S <sid> <gid> EPS=.. KZ=.. KY=..` 截面响应；
- `F <fid> <gid> <12 个 deff>` 纤维梁局部力/切线。

CSV：`SECT,...`（N/My/Mz + 3×3 截面切线）、`FBEAM,...`（q0..q11）、
`FTANG,...`（材料敏感 7×7 子块）、`FEQ,...`（弹性 vs Kl12 最大差）。

## 独立锚点（491 项检查）

`build_reference.py` 用 `numpy.polynomial.legendre.leggauss` **独立**
生成 Gauss 点（与 C++ Newton 求根无关），逐条重算全部期望值；
`compare_material_nl_3d.py` 校验：

1. 冻结 CSV / 输入 sha256 在档；
2. 当前 exe 重跑冻结 `.mat` → 与冻结 CSV **逐字节一致**；
3. 14 条 SECT + 5 条 FBEAM + 5 条 FTANG + 2 条 FEQ 逐项比对；
4. 力学语义：单轴弹性/钢屈服平台/混凝土开裂退出、矩形塑性铰逐阶段
   屈服 Mz 单调且 `My<M3<Mp=1.5My`、塑性切线 < 弹性切线、端剪力平衡、
   塑性弯曲削峰、纤维轴力开裂归零。

```text
ALL CHECKS PASSED (491/491)
```

## 回归

```bash
ctest --output-on-failure          # 30/30
powershell -File run_tests.ps1     # ALL TESTS PASSED
```

B3 能力矩阵仍如实保持 `P` 不 `Y`（五项力学尚未全部收口），不靠单项
完成提前打勾。
