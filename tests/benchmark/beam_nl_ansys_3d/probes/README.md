# CB3D16 物理诊断探针（只读归档，不进 CTest）

`beamnl_diagnostic_probe.py` 是 K 阶段做物理可行性验证时用的**早期原型**：
在正式写 `BridgeStudio/ansys_bridge_beamnl.py` 之前，用它确认 ANSYS
BEAM188 + NLGEOM 能否在 2π 卷圆下收敛、并定位 P-Δ 组的系统偏差来源。

> 正式算例生成请用 `BridgeStudio/cases/run_cases_beamnl.py`，
> 同源参数以 `BridgeStudio/ansys_bridge_beamnl.py::BNL_GROUPS/BNL_CASES`
> 为准。本脚本不带 `KEYOPT(3)` / `CNVTOL` / 大 `G` 等最终口径设置。

## 脚本能力

`make_apdl(kind, tag, n, L, A, Iy, Iz, J, M, P, Q, nsub)` 生成
BEAM188 + NLGEOM 的悬臂端弯矩（CIRCLE）或简支梁柱（PDELTA）APDL，
`run()` 调 `ANSYS180.exe -b` 并读 `*_disp.csv`。改模块级 `probe.G`
即可改 `MP,GXY`；在返回的 APDL 字符串上替换
`KEYOPT,1,1,0` 可插入 `KEYOPT,1,3,<n>`，替换 `NLGEOM,ON` 可关掉
几何非线性或插入 `CNVTOL`。

## 已跑过的四组诊断（结果见上级 README 第 3 节）

1. **可行性**：`python beamnl_diagnostic_probe.py`
   → 4 个卷圆 + 4 个 P-Δ job，8 个 `0 error`；2π 尖端闭合
   `2.63e-12 m`、转角误差 `2.9e-11`。
2. **非线性开关 / 子步数**：`NLGEOM` OFF/ON、`NSUBST` 10/40/200
   → P=0 线性基线差 `-5.07e-4` 不变，证明与几何非线性、收敛子步无关。
3. **G 与网格扫描**：`G = 8.1e10 / 8.1e12 / 8.1e13 / 8.1e14`，
   `n = 16/32/64/128/256` → 分离出**剪切锁定**（`1/n²`）与**真实剪切
   柔度**（收敛到 `+2.98e-4`）两项。
4. **插值阶 + 收敛参考值**：`KEYOPT(3)=0/1/2/3` × 两种 `G`，
   以及 `CNVTOL,F,1.0E-2,0.005,2` → 找到 `KEYOPT(3)=2 + G×1000`
   使线性基线收敛到 `+2.3e-7`；CNVTOL 消除 `MINREF` warning 且结果
   逐位不变。

## 结论落点

这些诊断最终固化成正式配置：

- `cir` 组：`KEYOPT(3)=0`、真实 `G`、`CNVTOL,F,1.0E-2,0.005,2`、n=8；
- `pd` 组：`KEYOPT(3)=2`、`G=8.1e13`、`A=20`、n=64。

只入库 `.py` / `.md`；诊断运行产生的 `.inp` / `.out` / CSV 与 ANSYS
二进制不入库。
