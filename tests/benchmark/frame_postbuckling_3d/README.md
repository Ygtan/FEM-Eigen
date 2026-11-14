# CB3D19：CR 梁刚接框架后屈曲 vs ANSYS BEAM188 + ARCLEN（第 21 项 CTest）

> 2026-11-14，N 阶段（b）。离线 CTest **不需要 ANSYS 许可**：所有 ANSYS
> 结果已冻结在 `cases/`，`--live-ansys` 只写临时目录，不覆写任何冻结文件。

## 1. 为什么做这个算例

CB3D17（L 阶段）用 von Mises 两杆**桁架**验证了球面弧长法能越过荷载极限点，
CB3D18（M 阶段）把同一内核以 `STABILITY=1` 接进 `cable_solver.exe` 主程序。
但两者都只用到 `TRUSS` 轴向单元，且有闭式解析解 `P(v)` 兜底；一旦换成
**必须用梁单元**的后屈曲问题，独立解析锚点就没有了，CR 梁的共转列式、
解析一致切线和弧长法选根策略也才真正受到考验。

本算例取 Williams toggle 型**刚接框架**：两根斜杆在顶点刚接、两端铰支，
每杆 24 个 `BEAM` 单元。顶点受竖向参考荷载，平衡路径依次出现：

```text
初始线性段 -> 第一极限点 -> 下降段（lambda 变负，反向荷载）
           -> 谷底 -> 再硬化段（弯曲刚度主导，lambda 快速回升）
```

与桁架不同，**全局 `lambda` 最大落在终点而不是极限点**，因此极限点判据
必须取「第一个 lambda 局部极大」，不能用 `max(lambda)`。

## 2. 内核与算法口径

* `frame_probe.exe`（新 CLI，纯新增文件）读 `*.frame`，构造 `Model3D` 后调用
  `stability_path_3d.hpp::solveArcLengthPath3D()`。
* 力学完全复用 `fem_3d_nl.hpp` 的共转 CR 梁：`expMapRot` 有限转动 +
  平分向量法随动三元组 + `crBeamTangent12` 解析一致切线。**没有重写单元公式，
  没有修改 `fem_3d_nl.hpp` / `stability_path_3d.hpp` / `NonlinearOptions3D`**。
* 增广系统：`R(u,lam) = lam*Fref - Fint(u) = 0`，
  `du^T du + psi^2*dlam^2 = ds^2`（Crisfield 球面约束），每个校正迭代同一次
  LU 分解解 `duR`、`duF` 两个右端项，由二次方程选根（与上一步增量夹角余弦最大）。
* 输出 `P` 行逐平衡点 `(step, lambda, v, residual, iterations, ds, pivot)`，
  `S` 行同源参数与派生量，`E` 行错误码。
* `k0_fem` 是内核在参考构形上组装 `K0`、解 `K0*du = Fref` 后由顶点柔度反算的
  **精确线性切线刚度**（不是路径首两点割线）。

几何 / 边界（FEM 与 ANSYS 同源）：

```text
左支座 (-B,0) -- NELEM 个 BEAM -- 顶点 (0,H) -- NELEM 个 BEAM -- 右支座 (B,0)
全部节点：UZ / ROTX / ROTY 约束（保持 XY 平面）
两端支座：UX / UY 约束，释放 ROTZ（铰支）
顶点：UX / UY / ROTZ 自由；参考荷载 Fref = -PREF * eY
控制自由度：顶点 UY，|v| >= VMAX 终止
```

## 3. 独立锚点：三套互不相干的实现

梁框架没有 `P(v)` 那样的全程闭式解，因此改用**初始切线刚度 `k0`** 作独立锚点，
并刻意让三套实现互相不知道对方：

| 实现 | 位置 | 作用 |
|---|---|---|
| 纯 Python 2D 线性刚架刚度矩阵 | `build_reference_frame.py::linear_frame_stiffness`（本 benchmark 自带） | 黄金值 / CTest 独立重算 |
| 第二套纯 Python 实现 | `BridgeStudio/ansys_bridge_frame.py::linear_frame_stiffness` | 同源导出侧独立复核 |
| FEM 内核参考构形切线 | `frame_probe.cpp::linearSlope` -> `S,k0_fem` | 待验证量 |

三者组装标准 2D Euler-Bernoulli 刚架单元（每节点 `UX/UY/ROTZ`），约束两端
平动后对顶点施加单位 `-Y` 荷载，`k0 = 1 / v_unit`。实测：

```text
k0_analytic = 1.211946322189e+07 N/m     (benchmark 侧)
k0_fem      = 1.211946322188e+07 N/m     (FEM 内核)
相对误差     = 1.907e-13                (容差 1e-8)
k0/PREF     = 3.029865805471
```

另外用 ANSYS 极小位移子步的 `lambda/v` 作第三方复核（相对误差 2.1e-3 ~ 3.1e-3，
容差 1e-2；ANSYS 第一子步已含微小几何非线性，故容差宽于内核对拍）。

## 4. ANSYS 口径与诊断结论

| 项目 | 处理 |
|---|---|
| 单元 | `BEAM188` + `SECTYPE,ASEC` `SECDATA`（直接给积分截面特性，不让 ANSYS 重算 I/J） |
| EB 退化 | `KEYOPT(3)=2`（三次形函数）+ `MP,GXY` 取 `G×1000 = 8.1e13`，沿用 K 阶段 CB3D16 口径：线性基线偏差由 -5.07e-4 降到 +2.3e-7 |
| 局部轴 | `IY = IZ` 对称截面。`beamRefTriad` 以 global Z 建局部轴，ANSYS 无 K 节点时缺省局部轴未必同向；`IY=IZ` 时绕杆轴的任意旋转都不改变刚度，从根上消除口径差 |
| 扭转 | FEM 侧 `e.G` 只进入 `GJ/L`，而 `ROTX/ROTY` 已约束，平面内不激活；两侧仍取同一 `G` 保持参数同源 |
| 应变度量 | BEAM188 **没有** LINK180 的 `KEYOPT(2)` 开关。实测完整路径主体一致到 1e-5 ~ 3e-4，无系统性的 2~3% 本构偏移 |
| 弧长控制 | `ARCLEN,ON,maxarc,1e-7` + `AUTOTS,OFF` + `NSUBST` + `NEQIT,60` + `TIME,1.0` + `OUTRES,ALL,ALL` + `ARCTRM,U,VMAX,apex,UY` |
| 反力反算 | `/POST1` 按数据集号 `SET` 逐子步 `*GET` 顶点 `U,Y`、两支座 `RF,FY`、`SET,TIME`；`lambda_r=(R1+R3)/PREF`，`lambda_t=TIME/TIME_END` |

诊断过的三件事：

1. **CR 梁离散收敛**：FEM 每杆 8 -> 12 -> 16 -> 24 单元时，与 ANSYS 的
   `worst|dlambda|` 从 1.43e-2 -> 4.79e-3 -> 1.57e-3 -> 1.34e-3。ANSYS 用
   8 单元与 24 单元结果几乎相同（三次形函数已收敛）。最终两边都取 24 单元。
2. **末端强硬化区**：`v -> 1` 时 `dlambda/dv` 很大，三点二次插值与两套求解器
   不同采样密度会把微小位移差放大成 ~1.5e-3 的 `lambda` 差。故 `VMAX` 取
   **0.95**，判据区间限定在 `0 <= v <= VMAX`，避开该段；曲线主体（含第一
   极限点、下降段、谷底、再硬化起始）全部保留。
3. **反力 / TIME 偏差**：`worst|lambda_r - lambda_t| = 5.66e-05`，仅出现在
   第一极限点附近（ANSYS ARCLEN 在该处子步的反力与伪时间瞬时不同步），
   其余子步 < 1e-6。判据容差 1e-4，并在报告里单列。

`ANSYS 日志 0 error / 1 warning`；该 warning 是 ARCLEN 极限点附近的常规提示，
与 L 阶段 `One indefinite tangent matrix is detected` 同类物理标志。

## 5. 实测结果

```text
参数    E=2.1e11 Pa  G=8.1e13 Pa  A=1e-3 m2  IY=IZ=1e-6 m4  J=2e-6 m4
        B=2 m  H=0.5 m  NELEM=24/slope  PREF=4e6 N  VMAX=0.95 m
        TOL=1e-9  MAXSTEPS=3000  MAXITER=40  DSMAX=0.01  DSMIN=1e-8

FEM     fine    659 点  lambda_cr_fit=0.089355796  lambda_min=-0.018654
                        v_final=0.950635  maxRES=3.534e-03 N
        coarse  655 点  lambda_cr_fit=0.089355767  lambda_min=-0.018653
        psi_half 656 点 lambda_cr_fit=0.089355789  lambda_min=-0.018653
        bad_noload      ERROR=zero_reference_load，不产出平衡点

ANSYS   fr_fine   120 子步  lambda_cr=0.089574  lambda_min=-0.018622
        fr_coarse  55 子步  lambda_cr=0.089752  lambda_min=-0.018607

对拍    FEM CLI 复现 13958 个值 worst|d| = 0.000e+00
        FEM vs ANSYS worst|dlambda| = 3.405e-04 (fr_fine) / 3.277e-04 (fr_coarse)
        第一极限点差 < 4e-4 (容差 2e-3)
        独立锚点 k0 相对误差 1.907e-13 (容差 1e-8)
        --live-ansys 临时目录复跑 2 个 job：
        fr_fine 480 值 max|d|=2.802e-10、fr_coarse 220 值 max|d|=2.142e-10,
        冻结 CSV 未覆写、无漂移

黄金值  12513 checks + 4 strings + 37 criteria + 6 cases
```

## 6. 文件与运行

```text
README.md
build_reference_frame.py          黄金值生成 + --check 自洽（自带独立锚点）
compare_frame_postbuckling_3d.py  第 21 项 CTest 离线入口
reference/cb3d19_frame_postbuckling.json
cases/
  cb3d19_frame_postbuckling.frame / .model.json / _fem.csv
  cb3d19_frame_postbuckling_compare_report.txt / _manifest.txt / _run_log.txt
  fr_fine.inp / fr_fine_ansys.out / fr_fine_curve.csv
  fr_coarse.inp / fr_coarse_ansys.out / fr_coarse_curve.csv
probes/README.md                 诊断探针索引
```

```bash
# 重新生成同源输入（只写文件，不跑 ANSYS）
python D:/AI-work/BridgeStudio/cases/run_cases_frame.py --generate-only

# 完整重跑（需要 ANSYS 18.0 许可，2 个 job）
python D:/AI-work/BridgeStudio/cases/run_cases_frame.py \
    --sync-to D:/AI-work/FEM-Eigen/tests/benchmark/frame_postbuckling_3d/cases

# 重建 / 校验黄金值
python build_reference_frame.py
python build_reference_frame.py --check

# 离线 CTest（无需 ANSYS）
python compare_frame_postbuckling_3d.py \
    --bin D:/AI-work/FEM-Eigen/build \
    --cases D:/AI-work/FEM-Eigen/tests/benchmark/frame_postbuckling_3d/cases \
    --reference D:/AI-work/FEM-Eigen/tests/benchmark/frame_postbuckling_3d/reference

# 可选：临时目录重跑 2 个 ANSYS job 并逐 CSV 单元 diff（不覆写冻结文件）
python compare_frame_postbuckling_3d.py ... --live-ansys
```

`ANSYS180.exe` 默认路径：`C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe`，
可用 `--ansys-exe` 覆盖。

## 7. 如实边界

* 本算例是**独立探针 `frame_probe.exe`**，不是 `cable_solver.exe` 成桥态闭环；
  `fem_3d_nl.hpp::solveModel3DNonlinear` 仍是荷载控制 Newton，不变。
* 能力矩阵里「非线性稳定全过程」仍标 `P` 不标 `Y`：B3 由 `P` 推 `Y`
  需要把路径跟踪接进成桥态 / 施工阶段闭环并完成真实桥梁对标，属后续 a 阶段。
* 只覆盖 Williams toggle 型平面刚接框架、单一本构、无初始缺陷，且是
  对称响应路径；**未覆盖** Lee frame 的 snap-back、分支屈曲、缺陷敏感性、
  材料非线性、空间框架、真实桥梁全过程稳定对标。
* `IY=IZ` 与 `G×1000` 是为对齐 ANSYS BEAM188 而选的 EB 口径，不代表内核
  支持任意非对称截面 / Timoshenko 剪切的同等精度。
* `PIVOT` 是全场 LU 对角量级比（约束自由度置 1 后近似恒定），**不是**切线
  奇异判据；极限点由 `lambda` 转折与二次插值判定。
* `--live-ansys` 写临时目录，不覆写 `cases/`；冻结 CSV 是唯一黄金来源。
