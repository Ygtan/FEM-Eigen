# 三维梁格输入格式说明

三维模块目前使用 `fem_3d.hpp`，支持三维桁架 `TRUSS` 和三维梁 `BEAM`。

## 单位

长度：m，力：N，弹性模量/应力：Pa。

## 节点

```text
*NODE
1 0.0 0.0 0.0
2 1.0 0.0 0.0
```

## 单元类型

```text
*ET
1 TRUSS
2 BEAM
3 BEAM7
```

- `BEAM`：6 自由度常规梁（圣维南扭转，`J`）。
- `BEAM7`：7 自由度薄壁梁（每节点多一个翘曲率 `φ=dθx/dx`），用于薄壁开口截面
  的约束扭转；需在 `*R` 中给出翘曲常数 `IW`（扇形惯性矩）。

## 材料

```text
*MP
1 E=200e9 DENS=7850
```

- `E`：弹性模量。
- `DENS`：材料密度 kg/m³（用于质量矩阵/自振频率与自动自重；缺省 0）。

## 实常数

桁架：

```text
*R
1 A=0.001 ISTRN=0.001
```

梁：

```text
*R
1 A=0.01 IY=1e-8 IZ=1e-8 J=1e-8 G=80e9 WT=1e-4 WB=1e-4 WY=1e-4 WTORS=2e-4 ISTRN=0.0001 ADDMASS=0
2 A=0.01 IY=1e-8 IZ=1e-8 J=1e-8 G=80e9 IW=1e-10
```

- `WT`/`WB`：绕 z 轴（竖向弯曲 Mz）的上/下缘截面模量。
- `WY`：绕 y 轴（水平弯曲 My）的截面模量；用于 `LEFT`/`RIGHT`/`MAX`/`MIN` 应力分量。
- `WTORS`：扭转截面模量 `J/r_max`；用于 `TORSION` 扭转剪应力 `τ=T/WTORS`。
- `ADDMASS`：附加线质量 kg/m（质量矩阵用；缺省 0）。
- `IW`：翘曲常数（扇形惯性矩 m⁶）；仅 `BEAM7` 使用，缺省 0 即退化为圣维南扭转。

索（`CABLE`）：

```text
*R
3 A=0.006 W=462.051 L0=99.923490771
```

- `W`：索单位长度重量 **N/m**（不是密度）。用于浅垂度弹性悬链线的
  垂度缩短量 `s(T)=w²Lh³/(24T²)`、Ernst 等效弹模
  `Eeq=E/(1+w²Lh³EA/(12T³L0))` 与索自重外荷载；**缺省 0 = 关闭垂度
  与 Ernst 修正**，索退化为直杆。生成器侧对应 `cable_dens × A × g`。
  注意它与 `*MP` 的 `DENS`（质量矩阵/自振用）不是同一个量。
- `L0`：无应力（弧）长度 **m**，缺省 0 时由几何弦长代替。ANSYS 对标侧
  对应 `INISTATE,SET,DTYP,EPEL` + `INISTATE,DEFINE,<e>,,EPEL,X,eps`。
  与 `*R ... ISTRN=`（初应变预张紧）含义不同，不可互相替代。
- 两类字段都支持**单元级覆盖**：`*EN,<id>,<n1>,<n2> W=<值> L0=<值>`。
  每根索的 W/L0 本来就各不相同，只靠 `*R` 共享一份会逼得同类索拆成
  多个实常数。单元级覆盖优先于实常数。
- 解析位置：`fem_3d.hpp` 的 `readModel3D`（`*R` 段与 `*EN` 后缀），
  闭式本构见 `cableTensionFromChord` / `cableUnstrainedLengthForTension`，
  Ernst 模量见 `fem_3d_nl.hpp::cableErnstModulus`。

> 端到端 ANSYS 对标见 `tests/benchmark/cable_sag_ansys_3d/`：ANSYS LINK180
> 是**无垂度**轴向杆，必须用分段索链 + 节点自重才能独立验证 `W>0` 的效应。

## 单元生成

```text
*TYPE,1
*MAT,1
*REAL,1
*EN,1,1,2
```

## 约束与荷载

```text
*BOUNDARY
1 1-6 0.0

*CLOAD
2 2 -1000.0
```

`*BOUNDARY` 行格式：`<节点> <dof 或 a-b> [值] [K=<刚度>]`。第四列 `K>` 即弹性支承：
该自由度接刚度 `K` 的接地弹簧，`<值>` 为弹簧基准位移（通常 0）。

```text
*BOUNDARY
1 1-6 0.0
10 2 0.0 K=2.0e6    ! 节点10 的 v 由竖向弹簧支承
```

## 徐变收缩

```text
*CREEP PHI=2.0 CHI=0.8
*SHRINKAGE ECS=3.0e-4
```

## 施工阶段

```text
*STAGE
1 ACTIVE=1,2
2 ACTIVE=3

*STAGE_LOAD
2 NODE=2 DOF=2 VALUE=-1000.0
```

`*STAGE` 行内除 `ACTIVE=`（本阶段激活的单元号，逗号分隔）外还有
`DEACTIVE=`（本阶段退出工作的单元号），两者可同时出现在一行：

```text
*STAGE
1 ACTIVE=1,2,3
2 ACTIVE=4,5 DEACTIVE=1
```

`DEACTIVE` 对应支架拆除、合龙段置换、挂篮前移这类"构件退出"，
不是"加反向荷载"——被置为不激活的单元不参与刚度组装。

### 阶段调索：`TENSION=` / `L0=`

同一行里还可以同时给出本阶段的索张拉规格，都是"单元号:数值"的逗号列表：

```text
*STAGE
1 ACTIVE=1,2,3
2 ACTIVE=4,5 TENSION=12:3.2e6,13:3.1e6
3 TENSION=12:3.35e6,13:3.35e6,14:3.0e6
4 L0=15:10.234,16:10.301
```

- `TENSION=<单元>:<目标索力 N>`：**本阶段末**这根索必须正好是这么多力（正为
  张拉）。内核返回 6 轮修正的内部迭代（`solveModel3DStages(m, 6, 1e-4)`，相对
  容差 1e-4），收敛不了就停在最后一轮、把残差如实写进输出，不假装成功。
- `L0=<单元>:<无应力索长 m>`：直接给下料长度。此时期末索力是**算出来的结果**，
  不是目标值，内核也不会反过头去改 `L0`。

两者对同一单元**互斥**（同时给直接报错）。其余四条硬校验同样报错不静默：
指向非 `CABLE` 单元、指向本阶段未激活的单元、`TENSION<=0`（放张请用 `L0=`）、
`L0<=0`。

`TENSION=` 的语义是**力**，不是"再加多少力"。给过一次的索若在后续阶段不再写
`TENSION=`，它的索力仍是状态量、会随结构变形变化；末阶段输出的 `devDes` 就是
拿它跟"最近一次明确给定过的目标"比出来的漂移，用来量化"先张的索被后张的带跑
多少"。工程上应对办法是补一轮二次调索——实测可把 16 根索全部拉回目标。

**为什么线性内核也要做这套迭代**：`assembleStiffness3D` 把 `CABLE` 当桁架
（只按 `EA/L` 组装），**根本不读 `e.L0`**。所以"改索长"对线性结果毫无影响；
目标索力必须以一对沿索向的初应变等效节点荷载（`istrn += ΔP/(Eeff·A)`）进入
方程。`ΔP = N目标 − N当前`，按索力状态算，不能按"参考几何上的初张力"算：结构
变形 0.5 m 时后者与真实索力能差近一倍（实测把先张的索带到 +89%）。

**阶段是增量法**：每阶段解的是"本阶段荷载**增量**在当前构型上的响应"
（`mc.loads = st.loads`，不累计），成桥态 = Σ 各阶段增量——见
`compute3DFinalStageBeamForces`。索力则相反，是**状态量、不累加**。

阶段求解输出（`example_3d_grillage.exe`）形如：

```text
STAGE2 E12 T=3.1997e+06 L0=10.2451 chord=10.3122 target=3.2000e+06 dev=1.1e-04 Tdes=3.2000e+06 devDes=1.1e-04
```

`T`/`L0`/`chord` 为阶段末状态，`target`/`dev` 只对本阶段明确调过的索非零，
`Tdes`/`devDes` 是相对设计目标的漂移。

### 阶段预应力：`PRESTRESS=`

施工阶段混凝土梁的纵向预应力，用**等效荷载法**，支持直线筋与二次抛物线曲线筋：

```text
*STAGE
2 ACTIVE=1,2,3,4 PRESTRESS=1-4:1.5e7:-0.35:0.85
3 ACTIVE=5,6 PRESTRESS=5,6:1.5e7:-0.35:0.85
4 ACTIVE=7,8,9,10,11 PRESTRESS=7-11:P=1.5e7:E=-0.35:S=-0.15:LOSS=0.85
```

两种写法（同一束只取其一，由值段是否含 `=` 自动判定）：

- **直线筋（位置模式）**：`PRESTRESS=<单元列表>:<P>:<e>[:<loss>]`
- **曲线筋（关键字模式）**：`PRESTRESS=<单元列表>:P=<P>:E=<e>:S=<矢高>[:LOSS=<loss>]`

字段：

- `<单元列表>`：连续梁单元的 id，`a-b` 范围或逗号列表。这些单元必须首尾相接
  成一条**开链**，链的首尾两端节点就是锚固端。
- `P`：张拉控制力 N（正）。
- `E` / `e`：**端部**偏心 m，**形心上方为正**（即 +z 方向）。下缘束写负值。
- `S` / `s`：**矢高** = 跨中偏心 − 端部偏心（m，向上为正）。`S=0` 即直线筋；
  `S<0` 跨中筋更低 → 产生向上均布力 → **上拱**。只在关键字模式出现。
- `LOSS` / `loss`：有效/张拉比，可省略（缺省 1.0），必须落在 `(0,1]`。
- `LOSSREF=<编号>`：引用 `*PRESTRESS_LOSS` 表按规范分项算损失，见下一节。
  **只在关键字模式可用，与 `LOSS` 必须且只能给一个**；位置模式不支持。

一束一个 `PRESTRESS=`，束号按本阶段内出现顺序自动编（1、2、…）。同一阶段可写
多束。曲线筋的抛物线定义在首尾锚固端连线上（矢高在弦中点），故**要求链共线**。

**等效荷载**：直线筋只有三样——两端一对轴压 `Pe`（沿梁轴，压短梁）+ 两端一对
偏心弯矩 `Pe·e`（绕横轴）。曲线筋在这三样之外再叠加两样，且必须同时施加才能
静力自平衡（竖向合力 = 0）：两端锚固切线竖向力 `4Pe·S/L`（两端同向，指束凹侧）
+ 跨内竖向均布 `q = −8Pe·S/L²` 节点化（首尾各半跨、中间节点全跨）。这里
`Pe = P·loss`。

**直线筋不产生上拱**——它只给轴压和端部弯矩，跨内无反弯；悬臂梁下缘直线束的
端部弯矩使自由端**下挠**（不是上拱），见 `test_3d.cpp` [25]。上拱来自**曲线筋**
（`S≠0`）的跨内反向均布力，且**依赖跨内单元细分**：本内核单元内弯矩仍是线性，
束跨内越细越接近抛物线的二次弯矩（单单元束的 `S` 效应会被两端抵消、退化为直线
筋）。见 `test_3d.cpp` [26] 的闭式解 `上拱 = −Pe·S·L²/(12EI)`。

硬校验（报错不静默）：`P<=0`、`loss` 不在 `(0,1]`、单元列表非空、指向非
`BEAM/BEAM7` 单元、指向本阶段未激活单元、单元不构成一条开链、曲线筋（`S≠0`）
链不共线。

与调索一样走**增量法**：等效荷载只在本阶段施加一次，产生的变形与内力通过累计
位移跨阶段延续；预加力是**状态量、不累加**。直线筋的轴向方向取锚固端连线的
单位向量，故要求链是直线（折线束本版按首尾连线近似轴向，文档不隐去这点）；
曲线筋链不共线直接报错。

阶段求解输出（`example_3d_grillage.exe`）形如：

```text
STAGE2 P1 Pe=1.275e+07
```

`Pe` 为有效预加力（= `P·loss`），只在张拉那一阶段出现，不累加。

### 阶段预应力分项损失：`*PRESTRESS_LOSS` + `LOSSREF=`

3D 与 2D **共用** `prestress_loss.hpp` 的公式与键名（JTG 3362-2018 §6.2）。
先声明一张参数表，再在束上引用；同一张表可被多束复用：

```text
*PRESTRESS_LOSS
1 SIGMA_CON=1395e6 EP=1.95e11 EC=3.45e10 FPK=1860e6 MU=0.25 K=0.0015 THETA=0.2 X=10 DELTA_L=0.006 L=10 T1=15 T2=35 SIGMA_PC4=8e6 SIGMA_PC6=8e6 PSI=1.0 ZETA=0.3 ECS=2e-4 PHI=2.0 RHO=0.005 RHOPS=1.0 POST=1

*STAGE
1 ACTIVE=1-4 PRESTRESS=1-4:P=15000000:E=-0.35:LOSSREF=1
```

- 表行格式 `<编号> KEY=VALUE …`，键名、单位与 2D `input_format.md` §6.4 完全一致；
  `REBAR=1` 切预应力螺纹钢筋的 σl5 口径，`POST=0` 切先张法。
- 内核**读完整份文件后**才解析引用，所以表写在 `*STAGE` 之前或之后都行。
  引用不存在的编号、`σpe/σcon` 不在 `(0,1]`、或未经 `readModel3D` 解析就直接
  求解，都明确报错；不会静默退回 `loss=1.0`。
- `σpe` 按张拉方式取 `σcon−σl1−σl2−σl4`（后张）或 `σcon−σl2`（先张）。
  **σl3 / σl5 / σl6 是后期损失，不进 σpe，也不进有效预加力。**
  有效比 `loss = σpe/σcon`，等效荷载仍用 `Pe = P·loss`。
- 只有 `LOSSREF` 路径才在结果里多写一段（12 位有效数字）：

```text
stage prestress loss
STAGE1 P1 SIGMA_CON=1395000000 SL1=87790888.5885 SL2=117000000 SL3=40000000 SL4=45217391.3043 SL5=20645991.0476 SL6=108364004.044 SPE=1144991720.11 RATIO=0.820782595059
```

  `stage prestress force` 段保持原来的 6 位有效数字、逐位不变。只给 `LOSS=`
  的束**不产生**这个段——内核没算过分项，就不拿一堆 0 冒充「各项损失为 0」。
- 边界与 2D 相同：σl2 只按直线筋 `ΣΔl/l·Ep`（曲线筋反摩擦长度未实现）；
  `εcs` / `φ` 需按附录 C 查表输入，程序不按龄期自动发展；分项损失不含
  筋-混凝土协同、次内力与截面应力验算。

## 三维影响线

```text
*INFLUENCE_3D
2 2
```

## 三维车道荷载

```text
*LANE_3D
NODES=1,2,3 Q=-10.0 P=-100.0 DOF=2 IMPACT=1.15 REDUCTION=0.9
```

`CODE=JTG_D60` 可按跨径自动生成车道荷载：`Q=10.5 kN/m`、`P=Pk(L)`
（L 取车道首末节点距离），无需手工填 `Q`/`P`。

```text
*LANE_3D
NODES=1,2,3 DOF=2 CODE=JTG_D60
```

横向布载（多车道）：单梁模型横向只有一个节点，算不出严格横向影响面，
规范上用**横向分布系数法**把多车道并排布载 + 偏载放大折进纵向影响线包络。
乘子 `m = LANE_COUNT · REDUCTION · ECCENTRICITY`：

```text
*LANE_3D
NODES=... CODE=JTG_D60 LANE_COUNT=3 REDUCTION=0.78 ECCENTRICITY=1.15
```

- `LANE_COUNT`：横向设计车道数（缺省 1）。
- `REDUCTION`：多车道横向折减系数（JTG D60-2015 表 4.3.1-5；3 道取 0.78）。
- `ECCENTRICITY`：偏载增大系数（横向偏载对纵向内力的放大，经验值，缺省 1.0）。

三个量都缺省时退化为单车道近似（旧的 `impact·reduction` 乘子）。

汽车冲击（JTG D60）：`IMPACT` 是**放大倍数**（`1+μ`），不是纯冲击系数 μ。索承
体系（`*CABLE_BRIDGE`/`*CABLE`）下，若 `IMPACT` 缺省（=1.0），`cable_solver`
自动在**成桥态切线刚度**上求竖向一阶弯曲基频 f1（识别主梁竖向位移 w 主导的
最低阶），按 JTG D60 算 μ（f<1.5→0.05，1.5~14→0.1767·ln f−0.0157，>14→0.45）
并乘 `1+μ`；显式写非 1.0 的 `IMPACT` 则尊重手动值。基频须用成桥态切线刚度
（张紧索切线刚度远大于弹性刚度），线弹性装配会把基频算低、μ 失真。meta 输出
`f1_vertical`/`impact_mu`/`impact_factor`。

## 三维影响面

```text
*LANE_SURFACE_3D
X_NODES=1,2,3 Y_NODES=1,2 Q=-10.0 P=-100.0 DOF=2 LANE_COUNT=2 IMPACT=1.15
```

`*LANE_SURFACE_3D` 同样支持 `CODE=JTG_D60`，并按 `LANE_COUNT` 自动横向折减。

## 三维单元影响线

按**单元**取影响系数（求某根杆的内力影响面时用），与按节点自由度的
`*INFLUENCE_3D` 并列：

```text
*ELEMENT_INFLUENCE_3D
37
41
```

每行一个单元号；也可一行写多个。

## 索承体系成桥态（斜拉桥 / 悬索桥）

以下五段**只由 `cable_solver.exe` 读取**。`example_3d_grillage.exe` 不认识
这些段——它只做线弹性的梁格求解。两种桥型共用同一个 exe，靠"有没有
`*SUSP_BRIDGE` 段"区分：有就是悬索桥，没有而给了 `*CABLE_BRIDGE`/`*CABLE`
就是斜拉桥。

### 斜拉桥：`*CABLE_BRIDGE` + `*CABLE`

```text
*CABLE_BRIDGE
MODE=COMPLETION NLSTEPS=3 NEWTON=30 MAXITER=20 TOL=1e-8 TOLT=1e-7 METHOD=MINENERGY MINT=0.02 VERBOSE=0

*CABLE
101 7           4000000.0
102 7           4000000.0
```

`*CABLE` 每行一根索：`<索单元号> <梁侧锚点节点号> [目标索力 N]`。
索力缺省或为 0 表示"自动定"（按 `METHOD` 定）。整段省略则自动取模型里
全部 `CABLE` 单元，锚点取 z 较低的那一端。

| 键 | 默认 | 含义 |
|---|---|---|
| `MODE` | `COMPLETION` | `COMPLETION` 成桥态闭环 / `TENSION` 定索力找形 / `STATIC` 纯静力 |
| `NLSTEPS` | `3` | 非线性荷载步数 |
| `NEWTON` | `30` | 每个荷载步的 Newton 迭代上限 |
| `MAXITER` | `20` | 外层迭代（索长 / 几何闭合）上限 |
| `MAXLENITER` | `15` | 索长迭代上限（悬索桥路径用） |
| `TOL` | `1e-8` | 几何闭合容差（m） |
| `TOLT` | `1e-7` | 索力相对容差 |
| `METHOD` | `MINENERGY` | 自动定索力方法；给 `RIGID` 则用刚性支承连续梁法 |
| `MINT` | `0.02` | 保力下限 = 该比例 × 最大索力 |
| `MINABS` | `0.0` | 保力下限的绝对值（N），给定后优先于 `MINT` |
| `DECKLOAD` | `0.0` | 桥面恒载合计（N），仅用于校核报告；0 = 未知 |
| `SENSITIVITY` | `LOWRANK` | 仅 `COMPLETION`：几何闭合敏感度 `LOWRANK`（默认，低秩 Woodbury）/ `FULL`（P3-3b 完整形状敏感度 `K_shape`，二次收敛但更贵）/ `FIXEDPOINT`（纯不动点，仅供对照）。未知值报错，不静默回退 |
| `VERBOSE` | `0` | 非 0 打印迭代细节；`FULL` 下另打印 `nt`、`|K_shape|max`、`|G_full|max`、`σ(I+G)` 与条件数 |

`SENSITIVITY` 是内核级诊断开关，BridgeStudio 的 `.dat` 导出不写它，
缺省即走既有路径，成桥态结果与之前逐位一致。三条路径只换 Newton 的
敏感度矩阵 `G`，线搜索与不动点兜底相同；口径与实测对比见
`USER_MANUAL.md §9.8`。

### 悬索桥：`*SUSP_BRIDGE` + `*MAINCABLE` + `*HANGER`

```text
*SUSP_BRIDGE
SYSTEM=EARTH TOWER1=3007 TOWER2=4007 SAG=30.0 DECKE0=1001 DECKE1=1060

*MAINCABLE
3001 3002 3003 ... 3040

*HANGER
2001 2002 2003 ... 2028
```

| 键 | 含义 |
|---|---|
| `SYSTEM` | `EARTH` 地锚式 / `SELF` 自锚式（其余取值按地锚式处理） |
| `TOWER1` `TOWER2` | 主缆链上的**塔顶节点号**，必须按桥纵向先后给出 |
| `SAG` | 主跨矢高（m），>0 时优先于 `SAGRATIO` |
| `SAGRATIO` | 矢跨比（`0.1` 表示 1/10）；两者都不给则报错 |
| `DECKE0` `DECKE1` | 加劲梁单元号区间，只影响 `[GIRDER]` 段的输出范围 |

`*MAINCABLE` 的单元必须**严格按"锚 → 塔1 → 塔2 → 锚"的链序**排列，内核
按这个顺序积分索多边形；顺序错了，找形出来的 H 与线形都是错的。
`*HANGER` 每根吊索必须一端接主缆、一端接加劲梁。

**两条硬校验**（违反直接报错退出，不静默计算）：

1. `TOWER1` / `TOWER2` 指到的节点上必须有梁单元。主缆若在塔顶旁边另建
   一个坐标重合的节点，主缆就成了"锚到锚的自由索"，塔不参与受力；这种
   错误求解照样收敛、位移倾角都正常，只有索力掉到设计值的百分之一。
2. 主缆两个端节点必须"要么被约束（锚碇）、要么接加劲梁（自锚式）"，
   否则是悬空索。

**两种体系的差别全在边界条件里**，内核不做特殊处理：

- 地锚式：主缆端锚在锚碇（节点全约束），H 传给地基；加劲梁端纵向固定；
  梁内轴力 ≈ 0。
- 自锚式：主缆端锚在加劲梁端节点上，H 压进加劲梁（梁内轴力 ≈ −H）；
  梁端纵向必须自由，纵向约束放跨中。

自锚式还有一条几何约束：锚点处主缆竖向分力 `V = H·Δz/L − qL/2`，令
`V = 0` 得平衡边跨 `L = √(2HΔz/q)`。边跨短于这个值，主缆在锚点就仍有
向上的分力，梁端会被拔起来（梁端支座不能承拉）。这不是偏好，是几何要求。

### 输出

与 `.dat` 同目录写 `<stem>_cable.txt`，分段纯文本，按段名顺序解析：

```text
mode=... nodes=... elements=... cables=... dofs=...
converged=... iterations=... elapsed_ms=...
closure_err=... tension_err=...
support_v_sum_kN=... applied_v_kN=... balance_pct=...
[CABLE]    elem anchor T L0 chord Ltarget dL dz
[MAINCABLE] elem T L0 z chordEnd
[HANGER]   elem deckNode T L0 chord
[GIRDER]   elem N My Vz
[GIRDEREND] elem end N My Vz
[GIRDER_LIVE] elem end NMin NMax MyMin MyMax VzMin VzMax
[REACTION] node dof R
[CAMBER]   node dx dy dz
[ERECTED]  node x y z
[U]        node dof value
[MODAL]    NMODES/NDOF/MASS/F/S ...（成桥态一致切线刚度上的模态，见 modal_output_3d.hpp）
[BUCKLING] ERROR=nonlinear_state_not_supported（本内核不做线弹性屈曲，见 buckling_output_3d.hpp）
[LENGTH]   KIND/STATUS/ERRKIND/.../I/E/NJ/J（外层索长迭代 / 几何闭合逐轮过程，见下）
[NEWTON]   KIND/STATUS/STAGE/NSTEP/.../S（内层共转列式 Newton 逐迭代过程，见下）
[OPT]      KIND/STATUS/OBJECTIVE/.../I/T（带约束成桥目标索力优化，见下）
[END]
```

`[GIRDER]` **只有悬索桥有**，每根加劲梁一行，`N`/`My`/`Vz` **各自**取两端中
绝对值较大者。它是给曲线和摘要看的——`max|N|` 可能来自 i 端、`max|My|` 来自
j 端，这一行描述的是梁上并不存在的截面，**不能拿它算应力**。

`[GIRDEREND]` 两种桥型都有，**每一端一行**（`end=1` 为 i 端、`2` 为 j 端），
所有 `BEAM`/`BEAM7` 单元都列。应力 `σ = N/A ± My/W` 必须用这一段：只有同一
端才是同一截面。主梁与塔柱都在里面——内核不替使用者判断"哪根是主梁"，
那个判断由前处理按 (材料, 实常数) 分组做。

内力是对着**制造线形**（`[ERECTED]`）算的：位移 `u` 是相对制造线形的量，
若拿设计构形当参考，算出来的梁内力会把预拱度当成变形。

### 索长迭代 / 制造线形几何闭合过程量（`[LENGTH]`）

段格式与口径定义在 `cable_length_output_3d.hpp`。它是**把内核已经算出来的
`IterationTrace3D` 写出来**，不改任何数值口径、不改变收敛判据：

```text
[LENGTH]
KIND=cable_length STATUS=ok ERRKIND=tension_rel
NCABLES=48 NITER=6 MAXITER=15 TOL=1e-07 CONVERGED=1
SENS JAVAIL=1 JREBUILD=1 JDIAGMIN=-3.05e7 JDIAGMAX=-9.6e6
I 1 0.6239684236 0.05139830483 0.1564928983 0.001698811455 1 0 REBUILD
I 2 0.007025579264 0.1572006999 0.0008461001838 9.169278796e-06 1 0 -
I 3 3.069754663e-08 0.1574791316 0 0 1 0 CONVERGED
E 9001 682510.3021 682510.3027 36.39508919 36.41630034 36.40054945
NJ=48
J 0 0 -30544922.38
J 0 1 1791540.618
```

- `KIND=cable_length`（`MODE=TENSION`）：未知量是各索无应力长度 L0，
  `ERRKIND=tension_rel` 的控制误差是 **max 相对索力偏差**（无量纲），
  由影响矩阵 `J = ∂T_i/∂L0_j = M·D` 的 Newton 步命中；`JAVAIL=1` 时
  给完整 `NJ` + `J <i> <j> <value>`（0 起）。`J` 对角应为负——放长 L0
  索力下降。冻结复用、残差停滞才重建，`JREBUILD` 记次数。
- `KIND=erection_geometry`（`MODE=COMPLETION`）：未知量是节点坐标偏移
  ΔX（制造线形预抬高），`ERRKIND=closure_m` 的控制误差是
  `max|ΔX + u|`（**单位 m**），走低秩 Newton / Woodbury + 回溯线搜索 +
  不动点兜底；**不产出 J**（`JAVAIL=0`），`SENS` 只给敏感度算子对角范围。
- `I` 行：`<it> <err> <geomErr> <stepMax> <stepRelMax> <alpha>
  <lineSearch> <flags>`；`geomErr = max|成桥弦长 − 设计弦长|`（m）只是
  **结果量、不作收敛判据**；`flags` 是 `REBUILD/CAPPED/FALLBACK/CONVERGED`
  的逗号串，无事件写 `-`。
- `E` 行：`<elem> <Ttarget_N> <Tfinal_N> <L0_m> <chord_m> <Ltar_m>`，
  力单位 **N**（前端换算成 kN）。
- `MODE=STATIC` 与悬索桥找形是闭式反算 L0 后一次求解，本来就没有外层
  迭代，写 `KIND=none STATUS=not_applicable` +
  `ERROR=no_iteration_static_mode` / `no_iteration_suspension_formfinding`。
  **不静默**：下游不给空表、不把它当求解失败。

**两条最容易读错的口径**：① 两种迭代的控制误差单位不同、判据不同，
不能画在同一根轴上；② 索力命中 ≠ 线形闭合——TENSION 实测索力收敛到
~3e-8 而成桥弦长偏差稳定在 0.15 m 左右，制造线形必须另跑 COMPLETION。
过程量只说明数值如何收敛，**不是规范意义上的设计合格判定**；未收敛如实
写 `STATUS=not_converged`。

### 内层共转列式 Newton 过程量（`[NEWTON]`）

段格式与口径定义在 `newton_output_3d.hpp`。它写的是
`solveModel3DNonlinear` 已经算出的 `NewtonTrace3D`，**不改任何数值口径**。
与 `[LENGTH]` 是两层不同的迭代：外层改模型参数（L0 或制造线形 ΔX），
每轮重跑一次完整非线性平衡；内层在给定模型上做荷载步 + Newton 求位移。
成桥态只写**收敛那一轮**的内层过程；`example_3d_stage_nl` 写出的
`<stem>.nlstage.txt` 则在每个 `STAGE` 块里各写一段（多一行 `STAGE=<id>`）。

```text
[NEWTON]
KIND=newton_cr STATUS=ok
NSTEP=3 NITER=15 MAXITER=30 TOLF=1e-08 TOLD=1e-09 CONVERGED=1
       REFACT=1 RESFIRST=19024291.18 RESLAST=0.04180729204
S 1 1 0.3333333333 19024291.18 0.4578740125 3059597.38 1 0
S 1 2 0.3333333333 4900282.175 0.1468978793 3059597.38 0 0
...
S 3 4 1 0.04180729204 6.325719466e-09 9178792.139 0 1
```

- `S <loadStep> <iter> <lambda> <|R|> <|du|> <|Ftarget|> <refactor> <converged>`：
  `iter` 是**该荷载步内**的迭代序号；`lambda = loadStep/NSTEP`；
  `|R|`/`|du|`/`|Ftarget|` 只统计自由自由度；末两列是 0|1，
  `converged=1` 表示该步之后 break。
- 收敛判据是 **或**：`|R| ≤ TOLF·max(|Ftarget|,1e-12)` 或
  `|du| ≤ TOLD·max(1,|u|)`，**每个荷载步都要满足**。换荷载步时
  `|Ftarget|` 变大，残差会跳升——整条 `|R|` 曲线不是单调的，不能按
  "跨步单调下降"去读，更不能据此判定发散。
- `REFACT` 只统计符号分解重建（首次 / 稀疏模式单调增长时新增位置），
  数值 `refactorize` 不计；是性能诊断量，不是力学量。
- 无内层过程（旧内核 / 未走几何非线性平衡）写
  `KIND=none STATUS=not_applicable` + `ERROR=no_newton_solution`，不摆空段。

过程量只说明"数值怎么收敛"，规范没有 Newton 迭代轮数的合格限值。

### 带约束成桥目标索力优化过程量（`[OPT]`）

段格式与口径定义在 `cable_opt_3d.hpp`。`OPTIMIZE=1` 时内核先解一个带箱
约束的二次型，得到目标索力 T\* 并写回 `*CABLE`，随后 `TENSION` /
`COMPLETION` 仍走各自既有的「命中 T\*」路径。**优化只负责定目标，不另造
闭环**。未请求优化时写 ERROR 段，既有索力口径逐位不变。

下面这份是 48 索算例 `OPTIMIZE=1 / OPTOBJ=MINENERGY / OPTBALANCE=0` 的
**守卫拒绝**实测：优化数学上收敛，但体系判据未满足，内核没有写回 T\*。

```text
[OPT]
KIND=constrained_cable_force STATUS=system_check_rejected
OBJECTIVE=MINENERGY START=GIVEN OBJECTIVE_KIND=NORMALIZED
NCABLES=48 NMOMENTS=128 NBALANCE=0 IRLS=0
NITER=196 CONVERGED=1 CONVBY=projected_gradient APPLIED=0
TMIN=28056.32906 TMAX=2104224.68 TMAXSRC=DEFAULT_1.5_RIGID
MPEAK_KIND=LINEAR_BARE
MPEAK0=35504728.63 MPEAK=29328621.44 BAL0=0 BAL=0
MPEAK_FINAL_KIND=NOT_APPLICABLE
MPEAK_FINAL=0
SYSCHECK=GUARDED AVAILABLE=1
SYSCHECK_LIMITS DECKSHARE_MIN=90 DECKSHARE_MAX=100 ADJRATIO_MAX=1.5
SYSCHECK_PRED DECKSHARE=89.04 SUPPORT_KN=3944.29 ADJRATIO=31.98 AT_X=261.25,267.5
SYSCHECK_FINAL_KIND=ORIGINAL_TENSIONS
SYSCHECK_FINAL DECKSHARE=99.31 SUPPORT_KN=250 ADJRATIO=1.08 AT_X=248.75,255
SYSCHECK_PASS=0 REJECT=system_deck_share_low,system_adjacent_jump
I 0 1 ...          ; <irls> <iter> <归一化目标> <max|M|(N·m)> <BAL(N)> <相对步长>
T 1 ...            ; <序号> <T0(N)> <T*(N)> <lo(N)> <hi(N)>
```

三条全过时则写 `STATUS=ok`（或 `not_converged`）、`APPLIED=1`、
`SYSCHECK_PASS=1 REJECT=none`、`MPEAK_FINAL_KIND=GIRDEREND_FINAL`、
`SYSCHECK_FINAL_KIND=OPTIMIZED_TENSIONS`。

- `OBJECTIVE=MINMAXM`（IRLS 逼近 min max|M主梁|）| `MINENERGY`（最小裸
  结构弯曲应变能）；`START=GIVEN`（`*CABLE` 段全给了目标索力）|
  `RIGID`（刚性支承连续梁法解析起点）。
- **三组量分属不同口径，段里各写一行 KIND**：`MPEAK0`/`MPEAK` 是去掉
  索后的裸结构线性化影响矩阵量，只用于优化目标；`MPEAK_FINAL` 由
  `cable_solver` 从本次 `[GIRDEREND]` 回填，是成桥态实际 `max|My|`；
  `I` 行的目标值是 `x = T/Tscale` 的**无量纲归一化值**。实测可能反向
  变化，**不得互相相减当作"弯矩改善量"**，也不得拿线性化目标冒充成桥态
  效应。未优化基线需另跑 `OPTIMIZE=0` 取 `[GIRDEREND] max|My|` 对拍。
- `APPLIED=0` 表示 T\* **未写回** `*CABLE`（求解失败、线性化目标不降
  反升，或体系守卫拒绝）。此时 `T` 行只是过程量，不能当设计索力，
  `MPEAK_FINAL_KIND=NOT_APPLICABLE`。`STATUS=objective_not_improved`
  / `system_check_rejected` / `not_converged` 都如实标注，不摆空表、
  不冒充通过。数学守卫 `objective_not_improved` 优先于体系守卫，但
  `REJECT` 行仍如实列出体系判据的超限原因码。
- **体系守卫 `SYSCHECK`（默认 `SYSGUARD=1`）**：优化目标只有弯矩、约束
  只有索力箱，只盯目标会得到「数学收敛、体系不可用」的 T\*。内核在
  **写回 T\* 之前**用与计算书工程校核表同一口径预估三条体系判据，
  不合格就不写回：
  - `DECKSHARE > DECKSHARE_MIN`（默认 90 %）：索承担桥面恒载比例，
    `100·Σ(T·ez)/DECKLOAD`，严格大于才算通过；
  - `DECKSHARE ≤ DECKSHARE_MAX`（默认 100 %）：超过即支座上拔；
  - `ADJRATIO ≤ ADJRATIO_MAX`（默认 1.5）：按支座断面分跨、跨内按
    锚点 x 排序取相邻 `max(Ti,Tj)/min(Ti,Tj)` 的最大值。
  `SYSCHECK_PRED` 用 T\* 预估（守卫依据），`SYSCHECK_FINAL` 用**成桥态
  实测索力**回填；`APPLIED=0` 时后者是原方案实测，故
  `SYSCHECK_FINAL_KIND=ORIGINAL_TENSIONS`，不是 T 行过程量。
  `REJECT` 原因码：`system_deck_share_low` / `system_support_uplift` /
  `system_adjacent_jump`。`SYSGUARD=0` 仅供诊断，写 `SYSCHECK=OFF`、
  `SYSCHECK_PASS=NA`，只记录不拦截；`DECKLOAD` 缺失或支座断面不足 2 个时
  写 `SYSCHECK=NOT_AVAILABLE AVAILABLE=0`，不假判合格也不假判失败。
  可用 `SYSDECKMIN` / `SYSDECKMAX` / `SYSADJMAX` 覆盖限值。
  这三条是**体系合理性判据，不是规范条文**，也不构成设计合格判定。
- `TMAX=0`（缺省）并不是"不限上界"：不给上界时优化器会把索力集中到
  少数索上（48 索实测 0.68~1.40 MN 被推到 0.14~6.67 MN），工程上不可用。
  程序默认取 `1.5 × RIGID 起点最大索力`，段里写
  `TMAXSRC=DEFAULT_1.5_RIGID`；这是**程序口径不是规范限值**，可用
  `TMAX=` 显式覆盖（`TMAXSRC=USER`）。
- 塔顶水平力平衡 `NBALANCE` 条是**软约束罚项**，与弯矩目标存在权衡。
  实测 48 索算例开启 `OPTBALANCE=1` 且用 `MINMAXM` 时，两者冲突会让
  线性化弯矩峰值反而升高，此时守卫拒绝采纳；`MINENERGY` 与
  `OPTBALANCE=0` 的组合已验证稳定收敛。
- `ERROR` 段：`no_optimization_requested`（默认 `OPTIMIZE=0`）、
  `not_applicable_suspension_formfinding`（悬索桥找形不适用）、
  `no_cables` / `bare_structure_is_mechanism` / `no_girder_moments` 等。
  旧内核无段由前端写 `no_opt_section`。

优化轮数、收敛率、目标下降幅度都**不构成设计规范意义上的合格判定**。

### 活载梁内力包络（`[GIRDER_LIVE]`）

`cable_solver` 也读 `*LANE_3D`（见「三维车道荷载」一节），在成桥态基础上
算**活载梁内力包络**：对每条车道、每根主梁两端，求 `N`/`My`/`Vz` 的
`min/max` 包络，输出为 `[GIRDER_LIVE]` 段。`[GIRDEREND]` 是恒载标准值，
`[GIRDER_LIVE]` 是活载包络，两者配对做荷载组合。

两条必须写死的能力边界：

1. **影响线在成桥态一致切线刚度上线性化**，不是线弹性装配——张紧索的切线
   刚度远大于其弹性刚度（Ernst 模量），用 `EA/L` 会把索算软、活载内力严重
   失真。这是几何非线性桥梁活载分析的标准做法。
2. **梁几何刚度（P-Δ）按一阶略去**：恒载已由非线性求解精确算入，活载增量
   用切线刚度求位移、再用线性弹性刚度提取梁端内力。梁的 P-Δ 是二阶小量。
3. **横向影响面用横向分布系数法近似**：单梁模型横向只有一个节点，严格横向
   影响面（横向单位力移动的响应面）需要梁格/板壳。这里按规范把横向布载折
   成 `LANE_COUNT·REDUCTION·ECCENTRICITY` 一个系数乘到纵向包络上，偏载增大
   系数是经验值，严格做法需梁格模型。

缺省车道：模型没写 `*LANE_3D` 时，内核自动识别主梁纵向链（判据
`|Δx| >= max(|Δy|, |Δz|)` 的水平 BEAM，排除塔柱与横梁）作单车道，并按
`CODE=JTG_D60` 生成车道荷载。识别不出则 `[GIRDER_LIVE]` 留空，调用方按
「未定义车道」提示，不静默当作无活载。

悬索桥额外多几行：`susp_system=` / `H_kN=` / `sag_m=` / `sag_target_m=` /
`span_m=` / `sag_ratio=` / `hanger_count` / `girder_N_max_kN` 等。
悬索桥的 `closure_err` 是**预变形量** `max|u|`（把设计构形放成成桥构形
需要预加的偏移），不是斜拉桥那个"几何闭合误差"。

---

## 附录：单元级探针参数文件（不是 `.dat`）

有两个 CLI 走**独立的纯 ASCII 参数文件**，不经 `readModel3D()` / `*ET`：

| CLI | 输入 | 用途 | 同源导出方 |
|---|---|---|---|
| `cable_sag_probe.exe` | `*.sag`（`E= / A= / LH= / NSEG=` + `C <id> <W> <L0> <amp>`） | 索垂度 / Ernst 闭式，对标分段 LINK180 索链 | `BridgeStudio/ansys_bridge_sag.py` |
| `beam_nl_probe.exe` | `*.beam`（`G <gid> key=...` + `C <cid> <gid> <CIRCLE\|PDELTA> <lambda>`） | 共转 CR 梁大转动 / 梁柱 P-Δ，对标 BEAM188 + NLGEOM | `BridgeStudio/ansys_bridge_beamnl.py` |

`*.beam` 的组参数键为
`E G L A IY IZ J NELEM LOADSTEPS MAXITER TOL Q`；`CIRCLE` 的端弯矩
`M = lambda·2πEI/L`、`PDELTA` 的轴压 `P = lambda·π²EI/L²` 都只由
`lambda` 与组参数派生，APDL 与 CLI 使用同一公式，不存在第二份手工输入。
两份参数文件都只做**单元级 / 单构件级**对标，不生成 `ModelSpec3D` 模型；
完整口径与冻结黄金值见
`tests/benchmark/cable_sag_ansys_3d/README.md` 与
`tests/benchmark/beam_nl_ansys_3d/README.md`。
