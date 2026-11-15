#!/usr/bin/env python3
"""跨工具链的 `[STABILITY]` 离线重放比较（第 24 / 25 项 CTest 共用）。

仓库历史冻结证据由 `C:/msys64/mingw64/bin/g++.exe` 产出，③a 的两项新
证据由 `C:/msys64/ucrt64/bin/c++.exe` 产出。两者都是 GCC 16.2、
`-O3 -DNDEBUG`，但 CRT / libm 不同，实测 `[STABILITY]` 逐点行里：

| 字段 | 跨工具链表现 | 判据 |
|---|---|---|
| `STEP` / `ITER` / `CONVERGED` | 逐位一致 | 字符串严格相等 |
| `LAMBDA` / `CTRL` / `DS` / `PIVOT` | 末位 ~1e-13 相对漂移 | 相对差 `< 1e-9` |
| `RES`（Newton 残差） | 近零量，相对差失真；实测绝对差 ≤ 4e-3 N | 绝对差 `< 0.1 N` |
| 其余所有行（含段头 `STATE` / `[U]` / `[NEWTON]` 等） | 逐位一致 | 掩掉 `elapsed_ms` 后严格相等 |

冻结文件本身的字节**仍由各 benchmark 的 sha256 + 黄金值自洽判据锁死**，
本模块只放宽「用另一套工具链编译出的 exe 重放」这一步，不放宽任何物理
判据，也不允许路径步数、迭代数、收敛标志或非逐点输出发生任何变化。

用法：
    ok, detail = compare_replay(frozen_text, replay_text)
"""
from __future__ import annotations

import re

POINT_RE = re.compile(
    r"^I STEP=(\d+) LAMBDA=(\S+) CTRL=(\S+) RES=(\S+) "
    r"ITER=(\d+) DS=(\S+) PIVOT=(\S+) CONVERGED=(\d+)$"
)
TOL_REL = 1.0e-9       # LAMBDA / CTRL / DS / PIVOT 的相对容差
TOL_RES_ABS = 0.1       # RES 的绝对容差 (N)；物理量级判据另由各 benchmark 负责
_EXACT = (0, 4, 7)      # STEP / ITER / CONVERGED：字符串严格相等
_REL = (1, 2, 5, 6)    # LAMBDA / CTRL / DS / PIVOT：相对容差


def _mask(txt: str) -> str:
    return re.sub(r"elapsed_ms=[0-9.eE+-]+", "elapsed_ms=*", txt)


def _split(txt: str) -> tuple[str, list[tuple[str, ...]]]:
    """返回 (去掉逐点行后的文本, 逐点字段元组列表)。"""
    lines = _mask(txt).replace("\r\n", "\n").replace("\r", "\n").splitlines()
    body = "\n".join(ln for ln in lines if not ln.startswith("I STEP="))
    pts = [m.groups() for ln in lines
           if ln.startswith("I STEP=") and (m := POINT_RE.match(ln))]
    n_raw = sum(1 for ln in lines if ln.startswith("I STEP="))
    if len(pts) != n_raw:
        # 出现格式不符的逐点行：记为 -1 条，调用方会因行数不一致而失败
        pts.append(("MALFORMED",) * 8)
    return body, pts


def _rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-30)


def compare_replay(frozen: str, got: str) -> tuple[bool, str]:
    """比较冻结输出与重放输出；返回 (是否通过, 诊断串)。"""
    fbody, fpts = _split(frozen)
    gbody, gpts = _split(got)
    if fbody != gbody:
        return False, "非逐点行（段头 / 结果段）不一致"
    if len(fpts) != len(gpts):
        return False, f"逐平衡点行数不一致: frozen={len(fpts)} replay={len(gpts)}"
    if not fpts:
        # 默认零污染 case 没有 [STABILITY] 段，走纯文本严格比较
        return True, "no [STABILITY] points; non-point text identical"
    for x, y in zip(fpts, gpts):
        for k in _EXACT:
            if x[k] != y[k]:
                return False, (f"STEP={x[0]} 字段索引 {k} 不一致: "
                               f"{x[k]} vs {y[k]}")
    worst_rel = 0.0
    for x, y in zip(fpts, gpts):
        for k in _REL:
            worst_rel = max(worst_rel, _rel(float(x[k]), float(y[k])))
    worst_res = max(abs(float(x[3]) - float(y[3])) for x, y in zip(fpts, gpts))
    detail = (f"{len(fpts)} points, worst rel={worst_rel:.3e}, "
              f"worst |dRES|={worst_res:.3e} N")
    if worst_rel >= TOL_REL:
        return False, f"路径字段超出相对容差: {detail}"
    if worst_res >= TOL_RES_ABS:
        return False, f"Newton 残差超出绝对容差: {detail}"
    return True, detail
