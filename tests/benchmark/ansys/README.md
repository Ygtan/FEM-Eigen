# ANSYS 对标基线

用 ANSYS BEAM4 的结果作为外部参考，校验 FEM-Eigen 二维梁的端到端计算（`.dat` 输入 →
`.result` 输出）。与 `tests/regression/` 的“黄金基线逐字节回归”不同，这里不是冻结自身输出，
而是冻结**外部参考值**，用来防“程序整体一起漂移但自身回归仍绿”的问题。

## 目录

- `reference/*.json`：黄金参考。每个 JSON 描述一个算例 + 一组特征量及其容差。
- `compare_ansys.py`：跑当前构建的 exe 重算 `.dat`，解析 `.result`，抽取特征量并比对。

## 当前算例

- `100m_simple`：100 m 简支梁，自重 + 均布荷载。特征量：跨中挠度、两端支座反力、
  跨中弯矩、跨中上下缘应力。

## 参考值来源（重要）

仓库里只有 ANSYS 的 APDL **输入脚本**（`100m-simple-ansys.txt`），**没有 ANSYS 数值结果文件**。
因此当前 `ref` 取的是 ANSYS BEAM4 对应的欧拉-伯努利闭式解——BEAM4 默认不计剪切变形，
对该工况与闭式解一致到机器精度，用它代表 ANSYS 结果是等价且可复现的。

若以后拿到真实 ANSYS 结果，直接改 `reference/*.json` 里的 `ref`/`rtol` 即可，解析器无需改动。

## 容差

判据：`|got - ref| <= atol + rtol * |ref|`。当前统一 `rtol = 1e-3`（0.1%），
远大于数值噪声，足以抓出刚度系数写错、荷载折算错、符号翻转这类真实错误。

## 手动运行

```powershell
python tests\benchmark\ansys\compare_ansys.py --bin build --src . --reference tests\benchmark\ansys\reference
```

已接入 CTest（`ansys_benchmark`）与 `run_tests.ps1`，随全量测试一起跑。
