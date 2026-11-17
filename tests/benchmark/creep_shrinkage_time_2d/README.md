# 结构时变徐变/收缩 JTG 3362 规范算例书锚点 benchmark（第 35 项 CTest）

第二十四次 B4 · 结构**整体时变**徐变链，离线、opt-in 无关（这组函数本来就是
`code_loads.hpp` 的连续拟合实现，非新增开关）。验证入口：

```
creep_probe_2d.exe <file>.dat [out.csv]   →  M,<KEY>,<value> / S,<age>,<KEY>,<value>
```

只读探针把内核中间量暴露成 CSV，供独立闭式 oracle 逐键核对；探针不含任何力学
求解逻辑（只调 `readModel2D` + `creepCoefficientAt` / `shrinkageStrainAt` +
`jtg3362NominalCreep` / `jtg3362ShrinkageNominal`）。

## 与 σl6 benchmark 的口径隔离（关键）

这是**两条不同的链**，刻意不相交，禁止混用：

| 维度 | 本 benchmark（结构时变链） | σl6 benchmark（预应力损失附录 C 表值链） |
|---|---|---|
| φ0/εcs0 来源 | `code_loads.hpp` **连续拟合** | 用户**显式录入**经核对的规范表值 |
| βH 的 RH | **实际 RH** 连续值 | RH **代表值** 55 / 80 |
| RH=75% 的 βH | ≈803.25 | ≈961.76 |
| εcs0（RH75/fcm40） | ≈3.674e-4（连续 βRH） | 0.310e-3（表值） |
| 用途 | 结构整体时变 AEMM 阶段分析 | 预应力损失 σl6 标量 |
| oracle | 独立 Python 复现连续拟合 | 独立 Python 复现表值发展公式 |

compare 脚本含**硬断言**锁死这条边界：本链 βH 必须是 803.25 而非 961.76，
εcs0 必须是 3.674e-4 而非 0.310e-3；探针源码不得出现硬编码代表值 55/80。

## 连续拟合公式（`code_loads.hpp` 附录 C.2.1 / C.1）

```
φRH   = 1 + (1 − RH/100)/(0.46·(h/100)^(1/3))
βF    = 5.3 / sqrt(fcm/10)
βT0   = 1/(0.1 + t0^0.2)
φ0    = φRH·βF·βT0
βH    = min(150[1+(1.2·RH/100)^18](h/100)+250, 1500)      （实际 RH）

εS    = (160 + 50(9 − fcm/10)) × 1e-6
βRH   = 1.55(1 − (RH/100)^3)
εcs0  = εS·βRH                                              （连续 βRH）

φ(t,t0) = φ0·[(t−t0)/(βH + (t−t0))]^0.3
βs(t,ts)= sqrt((t−ts)/(350(h/100)^2 + (t−ts)))
εcs(t,t0)= εcs0[βs(t,ts) − βs(t0,ts)]
E_eff   = E/(1 + χ·φ(t))                                    （AEMM，χ=0.8）
```

条文依据：JTG 3362-2018 附录 C；实现见 `code_loads.hpp`。

## 算例

| 算例 | 说明 |
|---|---|
| `creep_anchor` | C50/RH=75%/h=320.7mm/fcm=40/t0=7d/ts=3d；龄期 28d、97d；规范算例书锚点 |
| `creep_rh55` | 同锚点但 RH=55%，验证 φ0↑、βH↓ 单调性 |
| `creep_t0late` | 同锚点但 t0=28d，验证 βT0↓、φ0↓（加载越晚徐变越小） |
| `creep_norh` | 负例：缺 RH → 显式拒绝 |
| `creep_nofcm` | 负例：缺 FCM → 显式拒绝 |
| `creep_noh` | 负例：缺 H → 显式拒绝 |
| `creep_curve` | 负例：旧 CURVE 口径未给全 RH/H/FCM → 显式拒绝 |

### 规范算例书锚点（creep_anchor，RH=75%/fcm=40/h=320.7mm/t0=7d）

- φRH≈1.2296、βF≈3.3507、βT0≈0.5591 → **φ0≈2.301490**
- **βH≈803.2530**（实际 RH=75% 连续值）
- εS≈5.500e-4、βRH≈0.6680 → **εcs0≈3.673984e-4**
- 龄期 28d：φ≈0.7653、εcs 发展、E_eff≈1.240e11
- 龄期 97d：φ≈1.1561、E_eff≈1.039e11（AEMM 折减随龄期加深）

## 验证范围（75 checks）

- 冻结 .dat/.csv 的 sha256 与 reference 对齐；
- 当前二进制逐字节重放冻结证据；
- 每条正常算例逐键对齐**独立 Python 闭式 oracle**（φ0=φRH·βF·βT0、βH 连续 RH、
  εcs0=εS·βRH，不调用求解器）；
- 逐龄期 φ(t)/εcs(t)/βs/E_eff 交叉核对；AEMM `E_eff=E/(1+0.8φ)` 单调性；
- RH / t0 敏感性单调；
- 四个负例显式拒绝（缺参 / 旧 CURVE 口径）；
- 与 σl6 口径隔离硬断言；
- 探针只读内核、无硬编码代表值；
- **h 来源守卫**：内核 `m.creepH > 0` 覆盖传入单元 h、探针 φ0/βH/βs 均取
  `m.creepH`——锁死"给定 H 时探针与 AEMM 同源"，防未来改动覆盖逻辑后分叉。

## `H` 是理论厚度，不是梁高（口径说明）

`*CREEP_MODEL ... H=320.7` 解析到 `m.creepH`（`fem_2d.cpp:1002/1627`），
`code_loads.hpp:33` 明确 h 为**理论厚度 h=2A/u（mm）**；截面几何里的梁高是
**另一个字段** `fem_2d.hpp:260` 的 `H`，与徐变计算无关。**切勿把梁高填进
`*CREEP_MODEL H=`**——梁高约为理论厚度的两倍量级，会系统性抬高 βH/φ0。

给定 `H=` 时，`creepCoefficientAt`/`shrinkageStrainAt`（`fem_2d.cpp:751/771`）
优先用 `m.creepH` 覆盖任何传入的单元 `hTheory`，因此探针（全桥 H）与阶段
AEMM（`processTimeEffects` 虽传 `e.hTheory`，但被覆盖）严格同源、结果一致。
本 benchmark 只在**显式给定理论厚度 H** 的路径上成立。

## 如实边界

- 这是**结构时变拟合链**的闭式锚点 benchmark，不是 σl6 表值链，二者不可混用。
- oracle 复现的是 `code_loads.hpp` 连续拟合公式；2D 阶段求解的 AEMM 等效荷载
  （按上一阶段轴力施加、`Δφ` 增量、模量从 `baseModulus` 重算而非逐级累乘）
  不在本 benchmark 逐键展开，仅核对 `E_eff=E/(1+χφ)` 这一环。
- **单元级理论厚度路径未覆盖**：若模型只给 RH/FCM 而**不给** `*CREEP_MODEL H`，
  内核退到逐单元 `HTH`（`fem_2d.cpp:3088/3164-3170`），而探针在 :44 直接拒绝
  退出（无全桥 H 可采样）。这条路径本 benchmark 不覆盖，属如实边界。
- `processTimeEffects:807` 无条件写 `e.hTheory`，靠被调函数二次覆盖才正确；
  若将来重构去掉覆盖逻辑，h 来源守卫会立即失败。
- 3D 侧 `creepCoefficientAt3D()` 与 2D 简化曲线同源，本 benchmark 只覆盖 2D 探针。
- 不做真实制造参数接入、不与 ANSYS CEB90/蠕变选项做数值对拍。
  ANSYS 18.0 对拍可行性评估见同目录 `ANSYS_CEB90_对拍可行性评估.md`：
  内置 `TB,CREEP` 是金属蠕变本构，无 JTG 3362 的 φ0/βH/χ 对应参数，
  **不能本构级同源对拍**；伪弹性旁证口径有损，不建议纳入冻结 benchmark。
