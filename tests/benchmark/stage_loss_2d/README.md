# 第 36 项 CTest：2D 分项预应力损失阶段维度 + 阶段有效预加力 benchmark

B4 第二十五次接入（用户：「既然3D 侧早已是 Y，补 2D 阶段维度：让
`[PRESTRESSLOSS2D]` 支持 STAGE、把 σl2/σl6 的选项映射到 2D」），第二十七次
在同一 benchmark 内补「阶段有效预加力」段（`[PRESTRESSFORCE2D]`）后，
**该格由 `2D=P` 推 `2D=Y`，`Y` 口径 = 平面内完整**（不冒充三向空间量）。

## 验证了什么

- 2D 阶段束 `*STAGE_PRESTRESS ... LOSSREF=<表号>` 走与 3D 共用的
  `computePrestressLosses()`，σl1~σl6 / σpe / RATIO 逐键对齐独立 Python
  闭式 oracle（JTG 3362-2018 §6.2 + 附录 G 反摩擦 + 附录 C 显式表值龄期）。
- 阶段段头行 `BASIS=2d_stage`、行首 `STAGE<n> P<id>`，字段顺序、单位（Pa /
  无量纲）、12 位有效数字与 3D `stage prestress loss` 段一致。
- 只给 `LOSS=` 的阶段束不进明细：不摆零值行、不造阶段号。
- σl2 `L2ANTI/L2REPX`、σl6 `C2D/C_SKIP/EPS0/PHI0/C_RH/C_H/C_T0/C_TS/C_T`
  经**同一张** `*PRESTRESS_LOSS` 表自动映射，不另立 2D 公式。
- 软负例（无 LOSSREF）写 `ERROR=no_lossref_stage_prestress`；硬负例
  （LOSS+LOSSREF 冲突、LOSSREF 引用不存在表号）显式拒绝、不产出明细文件。
- 默认零污染：不带开关不写阶段文件；带/不带开关 `_stages.txt` 逐字节一致；
  既有第 34 项全局 `.loss.txt` 与 3D `.csv` 逐字节不变。
- 2D/3D 同参同值：`stage_same3d.dat` 与第 34 项 `creep_3d.dat` 的
  `*PRESTRESS_LOSS` 行逐字相同，9+8 个键逐键对拍 3D 探针。
- FE / BridgeStudio 源码契约守卫：旧全局段仍 `BASIS=2d_global`、旧
  `parse_loss_txt()` 仍只认 `2d_global + P<id>`。
- **阶段有效预加力（第二十七次）**：每阶段每束 `PE = P·LOSS` 逐行对齐独立
  oracle（`LOSSREF` 束的 `LOSS` 取 oracle 的 `σpe/σcon`）；多束多阶段编号
  `{(1,1),(1,2),(2,1)}` 不串阶段；只给 `LOSS=` 的束 `SIGMA_CON=0`（不反造）；
  软负例 `ERROR=no_stage_prestress`；两个 opt-in 开关互不依赖（同时开与单独开
  的两份文件逐字节一致，只开一个不会顺手写另一个）。

## 文件

| 路径 | 作用 |
|---|---|
| `cases/stage_default.dat` | 锚点：阶段 1 LOSSREF、阶段 2 仅 LOSS，验证只收录 LOSSREF 束 |
| `cases/stage_l2anti.dat` | σl2 附录 G 反摩擦（lf>l 分支，BRANCH=2） |
| `cases/stage_c2d.dat` | σl6 附录 C 表值龄期（βH≈962、σl6≈45.4MPa） |
| `cases/stage_cskip.dat` | C_SKIP=1 显式不参与（SL6=0、C_SOURCE=2） |
| `cases/stage_same3d.dat` | 2D/3D 同参算例（损失表与第 34 项 creep_3d 相同） |
| `cases/stage_noref.dat` | 软负例：无 LOSSREF，阶段段写 ERROR= |
| `cases/stage_badref.dat` | 硬负例：LOSSREF=99 引用不存在表号 |
| `cases/stage_conflict.dat` | 硬负例：LOSS 与 LOSSREF 同时出现 |
| `cases/force_multi.dat` | 阶段有效预加力锚点：阶段 1 两束（LOSS + LOSSREF）、阶段 2 一束 LOSS |
| `cases/force_none.dat` | 软负例：全程无 `*STAGE_PRESTRESS`，段落写 `ERROR=no_stage_prestress` |
| `cases/*.stage.loss.txt` | 冻结的阶段损失证据（由当前 mingw64 二进制生成） |
| `cases/*.stage.force.txt` | 冻结的阶段有效预加力证据（同上） |
| `build_reference.py` | 独立闭式 oracle + 冻结器 |
| `compare_stage_loss_2d.py` | CTest 比对脚本（269/269 checks） |
| `reference/stage_loss_2d.json` | 冻结哈希、cfg、oracle 与阈值 |

## 运行

```bash
cd D:/AI-work/FEM-Eigen
python tests/benchmark/stage_loss_2d/build_reference.py \
    --bin <build_dir> \
    --cases tests/benchmark/stage_loss_2d/cases \
    --reference tests/benchmark/stage_loss_2d/reference

python tests/benchmark/stage_loss_2d/compare_stage_loss_2d.py \
    --bin <build_dir> \
    --cases tests/benchmark/stage_loss_2d/cases \
    --reference tests/benchmark/stage_loss_2d/reference
```

`<build_dir>` 必须含 `example_stage.exe` / `example_mixed.exe` /
`loss_probe_3d.exe`。冻结证据由 `C:/msys64/mingw64/bin/g++.exe` 产出；
`ucrt64` 工具链会因 libm 末位漂移导致逐字节比对失败（见
`DEVELOPMENT_GUIDE.md §7`）。

## 如实边界

- 阶段段是独立 opt-in 文件，不是对旧全局段的改写；两个段互不混行，
  Python 侧 `parse_loss_txt()` / `parse_stage_loss_txt()` 严格隔离。
- 只给代表截面标量明细，不输出 σl2(x) / εcs(x) / φ(x) 沿程分布。
- 不跑 ANSYS、不做真实制造参数接入、不做规范合格判定。
- 第二十七次在同一 benchmark 内加「阶段有效预加力」段
  （`[PRESTRESSFORCE2D]` / `STAGE<n> P<id> SIGMA_CON= P= LOSS= PE=`，
  独立 `<stem>.stage.force.txt`，`PE` 与 `processPrestress` 里那一次
  `Pe = P·loss` 同源），checks 由 167 增至 269。
- B4 侧能力矩阵该格已由 `2D=P` 推 `2D=Y`，`Y` 口径 = **平面内完整**，
  不冒充三向空间量；真实制造参数与分项损失 ANSYS 对拍 2D/3D 同样没有。
