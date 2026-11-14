# L 阶段物理诊断探针（只归档，不进 CTest）

这些脚本是 CB3D17 正式模块定型前的**一次性物理诊断**，用于确认
「ANSYS ARCLEN 真的能越过极限点」以及「3% 偏差来自哪里」。只入库
`.py` / `.md`；ANSYS 的 `.out` / `.csv` / `.rst` 等产物不入库。

| 脚本 | 结论 |
|---|---|
| `fem_arclength_prototype.py` | Python 版球面弧长法原型：λ 从 0 → +1.1866 → 0 → −1.1878 → 0，与解析 `P(v)` 逐点差 1.7e-8 N，证明增广系统与选根逻辑正确 |
| `ansys_arclen_probe.py` | ANSYS `ARCLEN,ON` 首版原型；暴露 `*GET,...,RF,F,Y` 语法错（应为 `RF,FY`），以及 `P_ref=1 MN` 时 λ 到 1 就停、只到 v=0.12 |
| `arclen_limit_point_probe.py` | 把 `P_ref` 提到 2 MN（> 极限荷载 1.188 MN）后，ANSYS 走完整 S 曲线：v 到 1.006 m、λ 从 +0.613 穿 0 到 −0.613 |
| `strain_measure_probe.py` | 定位 `KEYOPT(2)`：`=1` 的 λ_max=0.600202 与对数应变解析 0.600244 吻合（差 7e-5）；`=0` 为 0.612574，高约 2%。四种解析应变度量并列比对 |

正式口径与全部实测数字见上级 `README.md`。
