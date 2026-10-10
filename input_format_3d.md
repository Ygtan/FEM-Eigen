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

## 三维影响面

```text
*LANE_SURFACE_3D
X_NODES=1,2,3 Y_NODES=1,2 Q=-10.0 P=-100.0 DOF=2 LANE_COUNT=2 IMPACT=1.15
```
