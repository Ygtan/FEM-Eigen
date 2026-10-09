# FEM 输入数据格式说明书（新格式）

## 1. 总体约定

- 所有关键字**大小写不敏感**（`*mp` 与 `*MP` 等价）。
- 统一采用 SI 单位：长度 `m`，力 `N`，应力/弹性模量 `Pa`，质量 `kg`，重力加速度 `m/s²`，温度 `℃`；车道均布荷载用 `N/m`，集中荷载用 `N`。
- 以 `#` 或 `!` 开头的行是注释。
- `*MP`、`*R` 数据行的行尾可以用 `!` 追加注释（不含 `=` 的字段会被忽略）。
- 解析分两步：先做**预处理**（`*FILL`、`*DO … *ENDDO`），把文件展开成普通行，再做正式解析。
- 二维结构统一按**每节点 3 个自由度**编号：
  - `1` = u(x)，`2` = v(y)，`3` = θ（转角）。
  - 桁架（TRUSS）、只受拉索（CABLE）、只受压杆（STRUT）只用 1、2；梁（BEAM）用 1、2、3。

---

## 2. 预处理宏

### 2.1 `*FILL` —— 直线等分生成节点

```text
*FILL a b
```

用**已存在**的节点 a、b 的坐标，在两点连线上线性插值，生成节点 `a+1 … b-1`。

- 要求 `a < b`，且节点 a、b 已在 `*NODE` 中定义。
- 只支持直线等分（拱形等曲线节点仍需显式列出）。

### 2.2 `*DO … *ENDDO` —— 循环

```text
*DO VAR start end [step]
  ... 循环体 ...
*ENDDO
```

- `end` **包含**；`step` 缺省为 `1`；`step` 不能为 0。
- 循环变量名**大小写不敏感**。
- 循环体内支持变量替换，表达式形式：
  - `I`
  - `I+k`、`I-k`
  - `k*I`
  - `k*I+m`、`k*I-m`

  其中 `k` 为正整数，`m` 为正整数。

  示例：`I`、`I+1`、`I-1`、`2*I`、`2*I+1`、`2*I-199`、`I+100`。

- 当前版本限制：
  - 不支持**嵌套循环**。
  - 循环体内不要使用 `*FILL`。

---

## 3. 声明表（声明一次，多处复用）

### 3.1 `*ET` —— 单元类型表

```text
*ET
1 BEAM
2 TRUSS
```

类型取值：

- `BEAM`：梁，3 自由度/节点。
- `TRUSS`：普通桁架，2 自由度/节点，可受拉也可受压。
- `CABLE`：只受拉索，2 自由度/节点，受压时退出工作。
- `STRUT`：只受压杆，2 自由度/节点，受拉时退出工作。

`CABLE` / `STRUT` 用于模拟预应力索、拉压间隙、单向支承等线性开/关问题。

### 3.2 `*MP` —— 材料表

```text
*MP
1 E=2.06e11 ALPX=1.2e-5 DENS=7850
```

| 字段 | 含义 | 单位 | 是否必填 |
|---|---|---|---|
| `E` | 弹性模量 | Pa | 是 |
| `ALPX` | 热膨胀系数 | 1/℃ | 否（温度荷载用） |
| `DENS` | 密度 | kg/m³ | 否（自动自重用） |
| `NU` | 泊松比 | 无量纲 | 否；`>0` 时梁按 Timoshenko 考虑剪切变形 |

### 3.3 `*R` —— 实常数表

```text
*R
1 A=0.08 I=0.02 WT=0.04 WB=0.04 H=1.2 ADDMASS=456.98
2 A=0.001 ADDMASS=4.98
```

| 字段 | 含义 | 单位 | 适用类型 |
|---|---|---|---|
| `A` | 截面积 | m² | 必填（梁、桁架、索、压杆） |
| `I` | 截面惯性矩 | m⁴ | 梁必填 |
| `WT` | 上翼缘抗弯模量 | m³ | 梁（应力用，缺省跳过应力图） |
| `WB` | 下翼缘抗弯模量 | m³ | 梁（应力用，缺省跳过应力图） |
| `H` | 截面高度 | m | 梁（温度梯度用） |
| `ADDMASS` | 附加线质量 | kg/m | 可选 |
| `ISTRN` | 初应变 | 无量纲 | 可选；`CABLE`/`STRUT`/`TRUSS` 用，正值为预拉应变 |
| `K` | 剪切修正系数 | 无量纲 | 可选；矩形默认 5/6，Timoshenko 剪切用 |

---

## 4. 激活命令

```text
*TYPE,1
*MAT,1
*REAL,1
```

设置当前激活的类型 / 材料 / 实常数，作用于其后的 `*EN`，直到再次修改。

---

## 5. 节点

```text
*NODE
1 0.0 0.0
21 100.0 0.0
*FILL 1 21
```

每行：`<节点号> <x> <y>`。

---

## 6. 单元

```text
*EN,eid,nodei,nodej
```

用当前激活的 `*TYPE/*MAT/*REAL` 生成一个单元（`eid` 单元号，`nodei/nodej` 两端节点号）。

```text
*TYPE,1
*MAT,1
*REAL,1
*DO I 1 20
*EN,I,I,I+1
*ENDDO
```

---

### 6.1 `*ANALYSIS` —— 分析类型

```text
*ANALYSIS TYPE=LINEAR
*ANALYSIS TYPE=NONLINEAR TOL=1e-6 MAXITER=50
```

也接受旧式无 `TYPE=` 写法：

```text
*ANALYSIS LINEAR
```

| 取值 | 含义 |
|---|---|
| `LINEAR` | 方案 A：线弹性 + cable/strut 开关迭代（默认） |
| `NONLINEAR` | 方案 C：Newton-Raphson 几何非线性（轴向单元） + cable/strut 开关 |

`TOL` 为 Newton 迭代相对收敛容差，缺省 `1e-6`；`MAXITER` 为最大迭代次数，缺省 `50`。

方案 C 中 `SimplicialLDLT`、`CHOLMOD`、`AMGCL` 等仍作为 Newton 迭代内部的线性求解器使用；若切线矩阵非正定导致 LDLT/Cholesky 类分解失败，可改用 `--solver LU` 或 `--solver BICGSTAB`。

### 6.2 CABLE / STRUT 的初应变与开关规则

轴向单元统一按以下机械轴力计算（**受拉为正**）：

```text
N = E·A/L·ΔL + E·A·ISTRN
```

其中 `ΔL` 是节点位移引起的杆长改变量（沿 i→j 伸长）。

- `CABLE`：仅当 `N > 0` 时保留；`N ≤ 0` 时本单元退出，不提供刚度、不输出内力。
- `STRUT`：仅当 `N < 0` 时保留；`N ≥ 0` 时本单元退出，不提供刚度、不输出内力。
- `ISTRN > 0`：轴向单元处于预拉状态；对 `CABLE` 表示预张紧，对 `STRUT` 表示初始松弛。
- `ISTRN < 0`：轴向单元处于预压状态；对 `STRUT` 表示预压缩，对 `CABLE` 表示初始松弛。
- `ISTRN = 0`：无初应变，普通 `CABLE` 从零轴力开始，只有被拉长才参与工作；普通 `STRUT` 只有被压缩才参与工作。

迭代方式（方案 A）：

1. 初始所有 `CABLE` / `STRUT` 均激活；
2. 组装并求解当前激活单元的线性系统；
3. 按上式检查每个 `CABLE` / `STRUT` 的轴力，违反单向约束的单元置为不激活；
4. 重复第 2~3 步，直到本步没有任何单元状态改变；最多 100 次。

若最终所有 cable/strut 都退出，系统会退化为仅由支座、普通杆件或梁保持稳定；自由节点若完全失去刚度仍会报“零刚度自由度”。

迭代方式（方案 C）：

1. 初始所有 `CABLE` / `STRUT` 均激活，位移从支座指定位移开始；
2. 组装当前构形下的内部力向量和一致切线刚度矩阵，求解 `K_t·ΔU = F_ext − F_int`；
3. 更新位移，并按当前构形的轴力重新判断每个 `CABLE` / `STRUT` 是否激活；
4. 重复第 2~3 步，直到位移增量和自由自由度残差同时满足 `TOL`，且单元激活状态不再变化。

当前方案 C 的几何非线性覆盖 `TRUSS`/`CABLE`/`STRUT`；`BEAM` 仍使用线性欧拉梁切线。轴向单元采用当前构形的工程应变和一致切线，能够反映张紧索的几何刚度效应。

### 6.3 `*PRESTRESS` —— 混凝土梁预应力自动加载

```text
*PRESTRESS
<id> ELEM=<单元列表> P=<张拉控制力> F=<跨中下垂量> E1=<起点偏心> E2=<终点偏心> LOSS=<损失系数> EPS=<筋弹模> APS=<筋面积> METHOD=AUTO
```

例：

```text
*PRESTRESS
1 ELEM=1-10 P=1500e3 F=0.20 E1=0.0 E2=0.0 LOSS=0.85 EPS=195e9 APS=0.001 METHOD=AUTO
```

字段说明：

| 字段 | 含义 | 单位 | 说明 |
|---|---|---|---|
| `ELEM` | 施加预应力的梁单元号 | — | 必须按节点顺序连续，例如 `1-10` |
| `P` | 张拉控制力 | N | 未扣除损失 |
| `F` | 抛物线相对两端连线的跨中下垂量 | m | 向下为正 |
| `E1` | 起点偏心距 | m | 向下为正 |
| `E2` | 终点偏心距 | m | 向下为正 |
| `LOSS` | 总损失系数 | — | 有效预应力 `Pe=P*LOSS` |
| `EPS` | 预应力筋弹性模量 | Pa | `TENDON` 方式使用 |
| `APS`/`AP` | 预应力筋面积 | m² | `TENDON` 方式必填 |
| `METHOD` | `AUTO`/`EQUIV`/`TENDON` | — | 缺省 `AUTO` |

`METHOD` 行为：

- `EQUIV`：把有效预应力换算成梁单元均布横向荷载、端部弯矩和轴力，直接施加到 `BEAM` 单元。
  - 要求梁组至少一端 `dof1` 可自由移动；两端都固定时无法用外荷载模拟弹性缩短。
- `TENDON`：自动生成 `CABLE` 预应力筋单元和内部 `RIGID_PIN` 销接刚臂，适合接入 `*ANALYSIS TYPE=NONLINEAR`。
- `AUTO`：`LINEAR` 分析使用 `EQUIV`；`NONLINEAR` 分析使用 `TENDON`。

等效荷载公式：

```text
Pe = P·LOSS
w  = 8·Pe·F/L²
```

`RIGID_PIN` 只约束梁节点与筋节点之间的平动刚体关系，不向筋节点传递弯矩，避免无粘结筋在转向点产生虚假附加弯矩。

端节点即使存在偏心，`RIGID_PIN` 也会通过偏心力臂把 `M=Pe·e` 作用到梁端节点；筋节点本身的转角不参与受力，因此不需要使用梁式刚臂。

**不需要在输入文件中特别注明端节点或内部节点。** 程序按 `ELEM` 给出的顺序自动识别：第一个单元的首节点为起点端节点，最后一个单元的末节点为终点端节点，其余中间经过的节点为内部节点。

`*PRESTRESS` 在整份输入文件读取完成后统一处理，因此 `*ANALYSIS` 写在它前面或后面都可以。

## 7. 约束与荷载（编号支持 `a-b` 区间）

以下各段里，节点号、单元号、自由度号都可以写成单个编号，也可以写成区间 `a-b`（包含两端）。

### 7.1 `*BOUNDARY` —— 支座约束

```text
*BOUNDARY
1 1-3 0.0        ! 节点1 的 u、v、θ 全固定
21 2-3 0.0       ! 节点21 的 v、θ 固定
```

每行：`<节点|区间> <dof|区间> <值>`。

### 7.2 `*CLOAD` —— 节点集中荷载

```text
*CLOAD
6 2 -1000        ! 节点6 的 dof2 施加 -1000
2-5 2 -100       ! 节点2..5 的 dof2 各施加 -100
```

每行：`<节点|区间> <dof|区间> <值>`。

### 7.3 `*ELEMENT_LOAD` —— 单元分布荷载

```text
*ELEMENT_LOAD
1-10 LY -100 -100
```

每行：`<单元|区间> <方向> <q1> <q2>`。

方向 `dir` 取值：

- `GX` / `GY`：全局 x / y 方向；
- `LX` / `LY`：单元局部轴向 / 横向；
- `G<角度>` / `L<角度>`：与全局 X / 与杆轴成指定角度。

`q1`、`q2` 为单元 i、j 两端沿该方向的分布荷载强度。

### 7.4 `*TEMPERATURE` —— 全模型均匀温升

命令式（值写在同一行）：

```text
*TEMPERATURE 25
```

也兼容旧式（`*TEMPERATURE` 后另起一行写值）。

### 7.5 `*ELEMENT_TEMPERATURE` —— 单元温度

```text
*ELEMENT_TEMPERATURE
1-2 UNIFORM 30        ! 单元1..2 均匀温升 30
3-5 GRADIENT 10       ! 单元3..5 温度梯度 10
```

每行：`<单元|区间> UNIFORM|GRADIENT <dT>`。

- `UNIFORM`：单元整体温升。
- `GRADIENT`：梁上下表面温差（需要实常数里有 `H`）。

### 7.6 `*GRAVITY` —— 重力加速度

```text
*GRAVITY
0 -9.81
```

每行：`<gx> <gy>`，用于自动自重。

### 7.7 `*SELF_WEIGHT` —— 显式自重（旧格式兼容）

```text
*SELF_WEIGHT
ALL RHO 7850
```

旧格式文件仍可用。新格式建议改用材料的 `DENS` 自动计算自重（见下）。

### 7.8 `*INFLUENCE` —— 指定影响线响应

```text
*INFLUENCE
1 DISP 3 2
2 REACTION 1 1
3 FORCE 5 M 2
4 STRESS 5 TOP 1
```

每行格式：

- `DISP <节点> <dof>`：节点位移影响线。
- `REACTION <节点> <dof>`：支座反力影响线。
- `FORCE <单元> <N|V|M> <端号>`：单元内力影响线。
- `STRESS <单元> <AXIAL|TOP|BOT> <端号>`：单元应力影响线。

`端号` 为 `1` 或 `2`，表示单元 i 端或 j 端。

### 7.9 `*LANE_LOAD` —— 车道荷载与包络

```text
*LANE_LOAD
1 NODES=1-21 Q=-10.5 P=-360 DIR=GY
```

字段：

- `NODES`：车道纵向节点号，必须按行车方向顺序排列。
- `Q`：均布车道荷载强度，单位 `N/m`，沿 `DIR` 方向，可为负。
- `P`：集中荷载，单位 `N`，沿 `DIR` 方向，可为负。
- `DIR`：`GX` 或 `GY`，缺省 `GY`。
- `METHOD`：`LINEAR`（缺省）或 `NONLINEAR`。

`METHOD` 行为：

- `LINEAR`：保留单位荷载影响线叠加方法；
- `NONLINEAR`：对每个集中荷载位置施加实际 `Q`、`P`，调用非线性求解器直接扫描包络，不通过单位影响线线性叠加。

程序先按方案 A 计算单位移动荷载影响线，再按方案 B 对 `Q + P` 进行最不利加载：

- 均布荷载加载在影响线同号区段；
- 集中荷载放在影响线最大/最小纵标处。

### 7.10 `*CE` —— 约束方程

采用罚函数法实现，兼容现有 LDLT / CHOLMOD / AMGCL 求解器。

```text
*CE_PARAMETER PENALTY=1e12

*CE
1 0.0 2 1 1.0 3 1 -1.0
```

每行格式：

```text
<id> <CONST> <node1> <dof1> <coef1> <node2> <dof2> <coef2> ...
```

表示：

```text
coef1·u(node1,dof1) + coef2·u(node2,dof2) + ... = CONST
```

`*CE_PARAMETER` 可设置全局罚系数 `PENALTY` 或 `ALPHA`，缺省 `1e12`。罚函数形式为：

```text
K_pen = α·c·cᵀ
F_pen = α·CONST·c
```

### 7.11 `*CREEP` / `*SHRINKAGE` —— 混凝土徐变与收缩

第一阶段采用**龄期调整有效模量法 AEMM**：

```text
*CREEP
ALL PHI=2.0 CHI=0.8

*SHRINKAGE
ALL ECS=3.0e-4
```

也可以只对指定单元施加：

```text
*CREEP
1-10 PHI=2.0 CHI=0.8

*SHRINKAGE
1-10 ECS=3.0e-4
```

说明：

- `PHI`：徐变系数 `φ`
- `CHI`：老化系数 `χ`，缺省 `0.8`
- `ECS`：自由收缩应变，正值为缩短
- 徐变按 `E_eff = E/(1 + χ·φ)` 修正单元弹性模量
- 收缩按 `-ECS` 作为等效初应变节点荷载施加
- `ALL` 会作用于所有非刚臂单元；混合钢/混凝土结构建议使用单元范围，避免误把钢束、钢桁架也施加徐变收缩。

### 7.12 `*TEMP_PROFILE` —— 规范竖向温度梯度

```text
*TEMP_PROFILE
1-10 T1=14.0 T2=5.5 A=0.10 B=0.80 H=0.80
```

字段：

- `T1`：梁顶温度
- `T2`：距梁顶 `A` 处温度
- `A`：第二温度控制点距梁顶距离，必须大于 0
- `B`：温度降为 0 的深度，缺省取梁高
- `H`：梁高，缺省取实常数 `H`

程序对矩形截面近似计算：

```text
N_T = ∫ E·α·T(z)·b dz
M_T = ∫ E·α·T(z)·(H/2 − z)·b dz
```

然后转成梁单元轴向力和端弯矩等效节点荷载。当前假定矩形截面，宽度按 `b=A/H` 计算；复杂箱梁以后可扩展纤维截面。

### 7.13 `*STAGE` / `*STAGE_LOAD` / `*STAGE_ELEMENT_LOAD` —— 施工阶段

```text
*STAGE
1 ACTIVE=1-4
2 ACTIVE=5-8 TIME=28.0

*CREEP_CURVE
0.0 0.0
28.0 1.0
100.0 2.0

*SHRINKAGE_CURVE
0.0 0.0
28.0 0.0001
100.0 0.0003

*STAGE_LOAD
1 NODE=4 DOF=2 VALUE=-1000.0
2 NODE=8 DOF=2 VALUE=-2000.0

*STAGE_ELEMENT_LOAD
2 ELEM=5-8 DIR=GY Q1=-5000.0 Q2=-5000.0

*STAGE_PRESTRESS
1 ELEM=1-4 P=1500e3 F=0.0 E1=0.0 E2=0.0 LOSS=0.85

*STAGE_BOUNDARY
1 NODE=4 DOF=2 VALUE=0.0

*STAGE_UNBOUNDARY
2 NODE=4 DOF=2

*STAGE_CREEP
1 1-4 PHI=1.0 CHI=0.8

*STAGE_SHRINKAGE
1 1-4 ECS=0.0001

*STAGE_TEMP_PROFILE
1 1-4 T1=14.0 T2=5.5 A=0.10 B=1.20 H=1.20

*STAGE_UNLOAD
2 NODE=2 DOF=2

*STAGE_TARGET
1 NODE=2 DOF=1 VALUE=-0.0025

*STAGE_CONTROL
1 STAGE=1 PRESTRESS=1

*STAGE_CONSTRAINT
1 NODE=2 DOF=1 LOWER=-0.01 UPPER=0.01 TYPE=DISP WEIGHT=1.0

*CAMBER
2 2
```

可带目标位移：

```text
*CAMBER
2 2 0.005

*FEEDBACK
2 1 -0.0025
```

含义：

- `*STAGE`：定义阶段号和该阶段激活/钝化的单元。
- `TIME`：阶段龄期，单位天，用于时变曲线插值。
- `*STAGE_LOAD`：指定阶段节点荷载。
- `*STAGE_ELEMENT_LOAD`：指定阶段单元均布荷载。
- `*STAGE_PRESTRESS`：指定阶段预应力，第一版使用等效荷载法 `EQUIV`。
- `*STAGE_BOUNDARY`：指定阶段新增支座约束。
- `*STAGE_UNBOUNDARY`：指定阶段移除支座约束。
- `*STAGE_CREEP`：指定阶段徐变参数。
- `*STAGE_SHRINKAGE`：指定阶段收缩应变。
- `*STAGE_TEMP_PROFILE`：指定阶段温度梯度。
- `*CREEP_CURVE` / `*SHRINKAGE_CURVE`：定义徐变系数、收缩应变随龄期变化的插值曲线。
- `*STAGE_UNLOAD`：移除指定节点荷载或单元均布荷载。
- `*STAGE_TARGET`：正装-倒拆迭代的目标响应。
- `*STAGE_CONTROL`：指定用于迭代调整的阶段和预应力序号。
- `*STAGE_CONSTRAINT`：不等式约束，支持 `DISP/REACTION/FORCE/STRESS`。
- `*CAMBER`：指定需要计算预拱度的节点和自由度。
- `*FEEDBACK`：指定实测响应，用于反馈修正。

预拱度结果包含：

```text
原始坐标
最终位移
预拱度
安装坐标/安装标高
```

非线性预拱度控制：

```powershell
example_stage --camber-iterate
```

程序通过迭代调整节点坐标，使安装坐标 + 最终位移趋于设计位置。

`*STAGE_TARGET` 支持：

```text
TYPE=DISP     默认，节点位移目标
TYPE=REACTION 支座反力目标
TYPE=FORCE    单元内力目标
TYPE=STRESS   单元应力目标
```

`FORCE/STRESS` 目标示例：

```text
*STAGE_TARGET
1 ELEM=1 COMPONENT=M END=1 VALUE=-100.0 TYPE=FORCE
1 ELEM=1 COMPONENT=TOP END=1 VALUE=1.0e6 TYPE=STRESS
```

每个目标可加权重：

```text
*STAGE_TARGET
1 NODE=2 DOF=1 VALUE=-0.0025 WEIGHT=2.0
```

当前支持多变量、多目标，程序采用有限差分 + Levenberg-Marquardt 最小二乘求解；控制变量可设置 `LOWER`、`UPPER`。

当前施工阶段求解采用正装分析框架：

1. 从空结构开始；
2. 每阶段激活/钝化单元；
3. 施加该阶段荷载；
4. 求解当前阶段结构；
5. 输出该阶段位移、内力、应力。

阶段节点荷载和单元均布荷载按正装顺序**累积**；施工临时荷载的卸载命令尚未实现，后续可增加 `*STAGE_UNLOAD`。

施工阶段求解现在会继承上一阶段的：

- 单元弹性模量修正；
- 已施加的预应力等效荷载；
- 节点荷载和单元均布荷载；
- 支座约束变化。
- 非线性分析的位移初值。
- 非线性 cable/strut 的激活状态。

收缩时变曲线按阶段龄期计算**增量收缩应变**，即 `Δε_cs = ε_cs(t_i) − ε_cs(t_{i−1})`。

徐变时变曲线按阶段龄期计算**增量徐变系数**，并根据上一阶段轴力施加增量徐变等效荷载，当前为轴向徐变逐步积分。

梁单元还会根据上一阶段弯矩施加增量徐变弯曲等效荷载，因此轴向和弯曲徐变均已逐步积分。

规范模型：

```text
*CREEP_MODEL MODEL=CEB90 PHI0=2.0 T0=7.0 BETA_H=700.0 POWER=0.3
*SHRINKAGE_MODEL MODEL=CEB90 ECS_INF=3.0e-4 TS=7.0 ALPHA_DS=35.0 POWER=0.5
```

简化 CEB-FIP 1990 形式：

```text
φ(t,t0) = PHI0 · [(t−t0)/(BETA_H+t−t0)]^POWER
ε_cs(t,ts) = ECS_INF · [(t−ts)/(ALPHA_DS+t−ts)]^POWER
```

未激活节点暂按固定处理，避免零刚度自由度。

倒拆分析现在按施工逆序移除荷载、反向激活/钝化单元，并以非线性求解器继承上一倒拆步位移和 cable/strut 状态。线性模型下可得到倒拆状态；非线性模型下为增量反演近似。

---

### 7.14 `*RELEASE` —— 梁端释放（铰 / 剪力 / 轴力释放）

```text
*RELEASE
5  I M          ! 单元5 的 i 端释放弯矩（铰接）
8  2 M          ! 单元8 的 j 端释放弯矩
12 1 M V        ! 单元12 的 i 端释放弯矩 + 剪力（只传轴力）
```

每行：`<单元号> <端 1/I|2/J> <类型 N/V/M> [<类型>...]`。

| 类型 | 含义 | 释放的局部自由度 |
|---|---|---|
| `N` | 轴力释放 | 轴向 u |
| `V` | 剪力释放 | 横向 v |
| `M` | 弯矩释放（铰） | 转角 θ |

说明：

- 仅对 `BEAM` 单元有效；实现方式为静态凝聚（把被释放自由度的单元刚度消去，该端内力强制为 0）。
- 典型用途：连续梁中间铰、挂孔、铰接墩——铰节点两侧单元转角各自独立、弯矩为 0，剪力/轴力连续传递。
- 释放端若为自由转角（无其它单元贡献），求解时会自动约束零刚度转角自由度。
- 暂不支持 `*ANALYSIS TYPE=NONLINEAR` 与端部释放同时使用（会在解析/求解时明确报错）。

---

## 8. 自动自重

- 每个单元的**线质量** = `DENS × A + ADDMASS`（kg/m）。
- 自动自重荷载 = 线质量 × `|g|`，方向沿 `*GRAVITY` 指定的重力方向。
- `ADDMASS` 既参与静力自重，也保存在单元对象里（`Element.addmass`），供后续动力分析使用。

---

## 9. 完整示例（100m 拱桥）

```text
# Bridge finite element model (new format)
*TITLE
100m Arch Bridge

! Element types
*ET
1 BEAM
2 TRUSS

! Materials
*MP
1 E=206000000000 DENS=7850

! Real constants
*R
1 A=0.08  I=0.02  WT=0.04 WB=0.04 ADDMASS=456.98
2 A=0.06  I=0.015 WT=0.03 WB=0.03
3 A=0.001  ADDMASS=4.98

! Nodes: main girder (deck) 1..21, x = 0..100 equally spaced
*NODE
1 0.0 0.0
21 100.0 0.0
*FILL 1 21

! Nodes: arch rib 101..109 (parabola, listed explicitly)
*NODE
101 10 7.2
102 20 12.8
103 30 16.8
104 40 19.2
105 50 20
106 60 19.2
107 70 16.8
108 80 12.8
109 90 7.2

! Elements: main girder 1..20
*TYPE,1
*MAT,1
*REAL,1
*DO I 1 20
*EN,I,I,I+1
*ENDDO

! Elements: arch rib 101..110 (chain 1-101-102-...-109-21)
*TYPE,1
*MAT,1
*REAL,2
*EN,101,1,101
*DO I 102 109
*EN,I,I-1,I
*ENDDO
*EN,110,109,21

! Elements: hangers 201..209
*TYPE,2
*MAT,1
*REAL,3
*DO I 101 109
*EN,I+100,2*I-199,I
*ENDDO

*BOUNDARY
1 1-3 0.0
21 2-3 0.0

*GRAVITY
0 -9.81

*END
```

---

## 10. 校验规则

- `*FILL`：`a < b`，且端点节点必须已定义。
- `*DO`：`step ≠ 0`，方向一致，且必须有匹配的 `*ENDDO`。
- `*ANALYSIS`：`TYPE` 只能为 `LINEAR` 或 `NONLINEAR`；`TOL > 0`；`MAXITER ≥ 1`。
- `*PRESTRESS`：`P > 0`、`LOSS > 0`、`EPS > 0`；`ELEM` 必须按顺序连续且均为 `BEAM`；
  `METHOD=TENDON` 或 `AUTO` 在非线性分析下还需要 `APS/AP > 0`。
- `*INFLUENCE`：`DISP/REACTION` 需要 `node dof`；`FORCE/STRESS` 需要 `elem component end`。
- `*LANE_LOAD`：`NODES` 不能为空；`DIR` 只能为 `GX` 或 `GY`；`METHOD` 只能为 `LINEAR` 或 `NONLINEAR`。
- `*CE`：每个方程至少一个自由度项；`dof` 只能为 `1~3`；`*CE_PARAMETER PENALTY > 0`。
- `*CREEP`：`PHI ≥ 0`，`CHI > 0`；`*SHRINKAGE` 需要 `ECS`。
- `*TEMP_PROFILE`：`A > 0`；`B` 若指定则不能小于 `A`；梁高需由 `H` 或实常数 `H` 提供。
- `*STAGE`：阶段号必须存在；`ACTIVE/DEACTIVE` 引用的单元必须已定义；`*STAGE_LOAD` 需要 `NODE/DOF/VALUE`；
  `*STAGE_PRESTRESS` 需要 `ELEM/P`；`*STAGE_BOUNDARY/UNBOUNDARY` 需要 `NODE/DOF`；
  `*STAGE_CREEP` 需要 `PHI`，`*STAGE_SHRINKAGE` 需要 `ECS`；
  `*STAGE_TEMP_PROFILE` 需要 `A > 0`；
  `*CREEP_CURVE` / `*SHRINKAGE_CURVE` 每行需要 `time value`；
  `*STAGE_UNLOAD` 需要 `NODE/DOF` 或 `ELEM`；
  `*STAGE_TARGET` 需要 `NODE/DOF/VALUE`，`*STAGE_CONTROL` 需要 `STAGE/PRESTRESS`。
- `*EN`：使用前必须先激活 `*TYPE/*MAT/*REAL`；引用的类型/材料/实常数必须已声明。
- 梁单元必须有实常数 `I`；`TRUSS`/`CABLE`/`STRUT` 单元必须有实常数 `A`。
- `CABLE`/`STRUT` 的初应变 `ISTRN` 为可选项，缺省为 `0`。
- 区间 `a-b` 必须满足 `a ≤ b`。
- 节点号、单元号不能重复（重复会报错或产生歧义）。

---

## 11. 旧格式兼容

旧的内联写法仍然支持：

```text
*ELEMENT, TYPE=BEAM, E=200e9, A=0.01, I=1e-4, WT=1e-4, WB=1e-4
1 1 2
```

以及旧式的 `*SELF_WEIGHT`、`*TEMPERATURE`（值在下一行）等关键字均可继续使用。
