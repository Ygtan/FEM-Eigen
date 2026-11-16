# 第二十一次力学⑤：动力 snap-through benchmark（cb3d29）

## 口径

独立 opt-in 模块 `dynamic_newmark_3d.hpp` + 文件驱动 CLI `dyn_probe.exe`，
**不改** `snap_probe.cpp` / `stability_path_3d.hpp` / `fem_3d_nl.hpp` /
`fem_3d.hpp` 的任何既有路径。`DYN=0` 时只有静态闭式锚点 `S,` 行，
无 `T,` 行；默认不调用本 probe，静力链逐位不变。

物理模型：浅拱 von Mises 两杆桁架（`B=2,H=0.5`），两端固定、顶点 Z 约束，
顶点 Y 向阶跃荷载 `−PREF`。闭式：

```
L(v)=hypot(B,H−v),  P(v)=2EA(L0−L)(H−v)/(L0·L)
v_cr=0.2142464215,  P_cr=1.188419643963e6
一致质量顶点等效质量 m_c=2·rhoA·L0/3=10.788793
集中质量（行和）顶点等效质量 m_l=rhoA·L0=16.183190
```

`PREF=0.85·P_cr≈1.010157e6`，**低于静力极限**；无阻尼阶跃下实测
`v_max≈1.3615`（`t≈0.00885 s`），`t_snap≈0.0039 s`，即惯性驱动的动力跃越。

## 文件

- `cases/cb3d29_dynamic_snapthrough.dyn`：11 个 case（动力 7 + 负例 4）。
- `cases/cb3d29_dynamic_snapthrough_fem.csv`：冻结输出（mingw64 `-O3`，14 位精度）。
- `reference/cb3d29_dynamic_snapthrough.json`：输入/输出 sha256 + 闭式锚点 +
  独立 SDOF RK4（子步长 1/5，每 20 个 FEM 步冻结一个样本）。
- `build_reference.py`：从冻结 CSV 与独立积分生成参考 JSON。
- `compare_dynamic_3d.py`：离线对比脚本，`--bin/--cases/--reference`。

## 校验层次（88 项）

1. 冻结 sha256 + 当前 `dyn_probe.exe` 重跑逐字节一致；
2. 闭式锚点独立重算（`v_cr`、`P_cr`、两种等效质量、荷载比）；
3. 独立 RK4 时程：三种工况（一致质量 / 集中质量 / 低荷载）逐采样点
   位移、速度核验（位移 <2e-3 m，速度 <5 m/s）；
4. 力学行为：低于 `P_cr` 仍跃越、荷载时程（阶跃 vs 爬升 / 半正弦）
   差异、阻尼耗散、质量敏感性、dt 加密收敛、动能/应变能非负、
   外功 `≈P_ref·v_final`、初加速度 `-P_ref/m`；
5. 负例 / 零污染：`DYN=0`、`DT=0`、`DENS=0`、`BETA=0` 均无 `T,` 行，
   `DYN=0` 不输出任何运行时键。

## 复现

```powershell
# 重新生成冻结 CSV（需与历史证据同工具链）
D:\AI-work\FEM-Eigen\_d_build\dyn_probe.exe `
  D:\AI-work\FEM-Eigen\tests\benchmark\dynamic_snapthrough_3d\cases\cb3d29_dynamic_snapthrough.dyn `
  D:\AI-work\FEM-Eigen\tests\benchmark\dynamic_snapthrough_3d\cases\cb3d29_dynamic_snapthrough_fem.csv
python build_reference.py --cases ../cases --reference ../reference
python compare_dynamic_3d.py --bin D:\AI-work\FEM-Eigen\_d_build --cases ../cases --reference ../reference
```

`ctest -R dynamic_snapthrough_3d_benchmark --output-on-failure`。
本 benchmark 不跑 ANSYS；ANSYS 瞬态动力学（如 `ANSECT` / `HHT`）与本
Newmark-β 隐式积分不同源，不作为对拍对象，深度如实登记为 FEM 自洽。
