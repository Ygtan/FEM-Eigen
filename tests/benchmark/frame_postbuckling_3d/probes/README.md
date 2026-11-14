# CB3D19 诊断探针

本目录归档 N 阶段定标过程中用过的、可独立复现的诊断脚本。它们**不参与
CTest**，也不写入 `cases/`；正式判据全部在 `compare_frame_postbuckling_3d.py`。

| 脚本 | 用途 | 复现命令 |
|---|---|---|
| `frame_nelem_convergence_probe.py` | CR 梁每杆单元数收敛性：确认 FEM 与 ANSYS 的差异主要来自 FEM 的共转梁离散，并据此把 `NELEM` 定标为 24 | `python frame_nelem_convergence_probe.py --exe <frame_probe.exe>` |

定标期另有两份**临时脚本不入库**：一份 numpy 版独立线性刚度
扫描（早期 `IZ=1e-7 / NELEM=8` 参数）和一份手写 BEAM188+ARCLEN
可行性 APDL 探针（`VMAX=1.0`）。前者已由
`build_reference_frame.py::linear_frame_stiffness`（纯 Python、
零第三方依赖）与 `BridgeStudio/ansys_bridge_frame.py` 两套正式
实现取代，后者已由 `ansys_bridge_frame.py::frame_apdl` 取代；
保留旧脚本会与正式 `FRAME_CFG`（`IZ=1e-6 / NELEM=24 / VMAX=0.95`）
口径冲突，故只归档上面这份可独立复现的收敛性探针与结论。

## 已归档的结论（无需重跑即可引用）

**FEM 离散收敛**（`IY=IZ=1e-6`、`PREF=4e6`、`VMAX=1.0`，对照同一份
ANSYS BEAM188 + ARCLEN 8 单元/杆曲线）：

```text
NELEM=8   worst|dlambda| vs ANSYS = 1.429e-02
NELEM=12                          4.793e-03
NELEM=16                          1.571e-03
NELEM=24                          1.337e-03
```

**ANSYS 侧离散**：8 单元/杆与 24 单元/杆的曲线几乎相同（BEAM188 三次形函数
已收敛），因此把两边统一取 `NELEM=24`，差异降到主体 1e-5 ~ 3e-4。

**残余系统差异**：`v -> 1` 的强硬化段 `dlambda/dv` 很大，插值与采样密度会把
微小位移差放大到 ~1.5e-3；最终 `VMAX=0.95` 避开该段，判据区间限定在
`0 <= v <= VMAX`。详见上级 `README.md` §4。

**ANSYS 反力 / 伪时间**：`worst|lambda_r - lambda_t| = 5.66e-05` 只出现在
第一极限点附近，其余子步 < 1e-6；`ARCLEN` 在极限点附近子步的反力与
`TIME` 瞬时不同步，属数值现象，判据容差 1e-4。

**ANSYS 应变度量**：BEAM188 无 LINK180 的 `KEYOPT(2)` 开关；本算例没有出现
L 阶段 `KEYOPT(2)=0` 那种 2~3% 的系统性本构偏移，故正式对标直接用默认
应变度量 + `KEYOPT(3)=2` / `G×1000` 的 EB 退化口径。
