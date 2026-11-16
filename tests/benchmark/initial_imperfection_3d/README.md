# `initial_imperfection_3d` — 第二十一次力学②：初始缺陷形态 / 幅值 / 敏感性

## 口径

第二十一次五项力学推进的第②项，**opt-in、默认全关、离线、B3 仍 `P` 不 `Y`**。

- **缺陷源**：线性屈曲模态 `(K+λKg)·φ=0` 的最低模态（ von Mises 两杆
  对称竖向 / 第 2 阶反对称横向），由结构自身几何刚度决定，不是人为
  给定的扰动。内核 `computeBucklingModeShape3D` 是**独立新增函数**，
  不改 `fem_3d.hpp::computeBucklingResult3D`（它只返回 factors，被多处调用）。
- **施加方式**：模态按 `IMPAMP` 缩放加到节点坐标；TRUSS 的 `L0` 同步
  改成缺陷后弦长 → **无应力几何缺陷**（制造误差型，λ=0 即平衡态）。
  若保留 L0 不变（残余应力型），大幅值（≥0.005）会把切线刚度推到
  奇异（pivotRatio≈5e-9）、弧长法失稳，故本项验收采用无应力口径。
- **敏感性**：对称模态幅值 0→0.02 时 `λ_max` 单调下降；反对称模态
  使顶点横向偏移、对称破缺。

## 输入键（`*.snap` 的 `C` 行，默认全关）

| 键 | 默认 | 含义 |
|---|---|---|
| `IMPERF` | 0 | 1=施加屈曲模态初始缺陷 |
| `IMPMODE` | 0 | 模态序号（0=对称竖向，1=反对称横向） |
| `IMPAMP` | 0.0 | 缺陷幅值（长度） |

`IMPERF=0` 时整段施加逻辑跳过，诊断 `S,` 行不输出，既有 `.snap`
输出与第 19/28 项冻结 CSV 逐位不变。

诊断 `S,` 行（仅 `IMPERF=1`）：`imp_ok` / `imp_lambda`（线性屈曲因子，
与非线性 λ_max 不同源）/ `imp_dx2` / `imp_dy2`（缺陷引起顶点坐标变化）/
`imp_maxphi`（模态最大分量，归一化=1）。

## 独立解析锚点

缺陷后几何 `h=H−a`、`L0'=hypot(B,h)`：

```
L(v)=sqrt(B²+(h−v)²),  N=EA(L−L0')/L0'
P(v)=−2N(h−v)/L,       λ_cr=P(v_cr)/PREF   (dP/dv=0 黄金分割)
```

细步（`DS=0.005`）+ 相邻三点二次插值 `λ_max_fit` 与闭式 `λ_cr`
误差 ≤1.5e-4；闭式 `P(v)` 逐点锚点走 `tol_dP_abs=1e-5 N`。

## 测试内容（35 项检查）

1. 当前 `snap_probe.exe` 重跑冻结 `.snap` → 与冻结 CSV 逐字节 / 逐值一致；
2. 对称模态 a0–a3 的 `λ_max_fit` 对齐独立闭式 `λ_cr`；
3. 闭式 `P(v)` 逐点锚点；
4. `λ_max` 随幅值严格单调下降（闭式与 FEM 同步）；
5. 反对称模态对称破缺（顶点 X 偏移 = 幅值、Y 不偏移）；
6. `IMPMODE=9` 越界负例 → `imp_ok=0`，不静默；
7. `IMPERF=0` 零污染：不输出诊断键、与零幅值路径逐位一致。

## 复现

```bash
cd tests/benchmark/initial_imperfection_3d
python build_reference.py                       # 重建 reference JSON
python compare_imperfection_3d.py --bin <build> --cases cases --reference reference
```

CTest：`ctest -R initial_imperfection_3d_benchmark --test-dir <build> -V`。

## 边界

- 只在 `snap_probe` 独立探针层，**未接 `cable_solver` 成桥态主链路**；
- 不跑 ANSYS；本项不把能力矩阵 B3 行 `P→Y`；
- 缺陷建模为无应力几何偏差（L0 同步），残余应力型未纳入本项验收。
