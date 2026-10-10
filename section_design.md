# 纤维截面模型设计文档

## 1. 目标与定位

把"截面"从手工输入的一堆标量（`A/Iy/Iz/J/WT/WB`）升级为一个**可积分的纤维集合**，
使下列能力统一退化为"对纤维循环"：

- 截面属性（A、形心、惯性矩、截面模量）
- 温度梯度 / 温度剖面的截面合力
- 组合截面（钢 + 混凝土，不同材料不同龄期）
- 变截面（沿梁长的截面变化）
- 材料非线性（纤维级 σ-ε 本构：混凝土开裂、钢筋屈服）
- 徐变收缩分层（不同纤维龄期差）

核心思想：**截面的一切响应都是 `Σ f(纤维)`**。一个实现，处处复用。

## 2. 数据结构

```cpp
struct Fiber {
    double y, z;      // 纤维形心坐标（相对参考点）
    double area;      // 面积
    double E;         // 弹性模量（组合截面时各纤维不同）
    double hy, hz;    // 半高/半宽（用于精确的截面边缘）
    int    mat;       // 材料索引（本构/徐变/收缩用，可扩展）
};

struct SectionProperties {
    double A, cy, cz;        // 面积、形心
    double Iy, Iz, Iyz;      // 惯性矩、惯性积（相对形心）
    double yTop, yBot;       // 上下缘坐标（相对形心）
    double zLeft, zRight;    // 左右缘坐标（相对形心）
    double Wtop(), Wbot(), WyLeft(), WyRight();
};
```

## 3. 截面属性积分（已实现）

`computeSectionProperties(fibers)`：

```
A   = Σ area_i
cy  = Σ y_i·area_i / A      cz = Σ z_i·area_i / A
Iy  = Σ (z_i-cz)²·area_i    Iz = Σ (y_i-cy)²·area_i
Iyz = Σ (y_i-cy)(z_i-cz)·area_i
Wtop = Iz/|yTop|  Wbot = Iz/|yBot|  Wy = Iy/|zLeft|, Iy/|zRight|
```

边缘坐标用 `y±hy`、`z±hz` 精确到截面边界（避免纤维形心≠边界的一阶离散误差）。

## 4. 纤维剖分（已实现）

- `makeRectangularFibers(b, h, nz, ny, E)`：矩形网格。
- `makeIFibers(bf, tf, tw, hw, n, E)`：I 型（两翼缘 + 腹板）。

计划补充：箱型、T 型、任意多边形（扫描线法剖分）。

## 5. 组合截面（步骤 1）

钢混组合梁：同一截面内不同纤维不同 `E`。用**等效模量法（换算截面）**：

- 选基准模量 `E_ref`（通常取混凝土 E_c）。
- 每个纤维的换算面积 `area* = area × E/E_ref`。
- 用 `area*` 积分得到换算截面属性（A*、形心、I*）。
- 真实应力恢复：`σ = (E/E_ref) × σ*`，其中 σ* 由换算截面算出。

新增 `Section { fibers, E_ref }` 与 `computeTransformedProperties(section, E_ref)`。

## 6. 变截面（步骤 2）

沿梁长截面变化（连续箱梁梁高渐变）：

- 单元级做法（推荐先做）：每个梁单元用一个**平均/两端插值的 Section** 生成
  A/Iy/Iz/W，仍用现有 EulerBeam2D 刚度，只是属性按单元取。
- 精确做法（后续）：`stiffness = Σ w_g · B(ξ_g)^T · D(ξ_g) · B(ξ_g) · |J|`，
  在 Gauss 点上按 `Section(x)` 取 `D`（Gauss 变刚度积分）。

接口：`SectionAt { double x; Section sec; }`，`sectionAt(x)` 插值。

## 7. 接入梁单元（步骤 3）

让梁单元从 `Section` 取属性，而不是手输 `A/Iy/Iz/W`：

- `EulerBeam2D` 增加"按 Section 构造"：`EulerBeam2D(section, E_ref, L)`。
- 温度梯度 / 温度剖面改为 `sectionTemperatureResultants(section, T(z))`。
- 应力改为 `sectionStress(section, N, My, Mz)`（四角点 / 任意纤维）。
- 3D 梁同理（`Element3D` 关联 `Section`）。

保持向后兼容：无 `*SECTION` 时仍用手输 `A/Iy/Iz/W`。

## 8. 材料非线性（步骤 4）

纤维级 σ-ε 本构，逐纤维判断、积分截面切线：

```
ε_i = ε_0 + y_i·κ_z + z_i·κ_y      // 平截面假设
σ_i = σ(ε_i, 材料本构)              // 混凝土受压/开裂、钢理想弹塑
N = Σ σ_i·area_i
M = Σ σ_i·y_i·area_i（及 z 项）
切线：E_t,i = dσ/dε，组装 D_t = Σ E_t,i·[1,y,z; y,y²,yz; z,yz,z²]·area_i
```

计划本构：混凝土（Hognestad 受压 + 拉裂退出）、钢筋/钢（理想弹塑性 / 双折线）。

## 9. 路线图

1. 组合截面（等效模量换算）
2. 变截面（单元级平均截面）
3. 接入梁单元（Section → A/Iy/Iz/W + 温度/应力纤维化）
4. 材料非线性（纤维本构 + 截面切线）
