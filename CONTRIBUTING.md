# 贡献与提交规范

本项目提交信息沿用既有中文风格。提交前请先跑测试，改动结果如有意变化需同步刷新回归基线。

## 提交信息格式

| 场景 | 格式 | 示例 |
|---|---|---|
| 修复评审报告中的 P 级问题 | `P<级别>-<编号>：<一句话说明>` | `P3-1: 将 fem_io.hpp 拆分为声明 + 实现（fem_io.cpp）静态库` |
| 一次修复多项 P 级问题 | `P2-7/P3-5/P2-6：<一句话说明>` | `P1-5/P1-6/P1-7：徐变逐阶段 AEMM、Newton 收敛保障、CHOLMOD 复用接入` |
| 对已有 P 级项的补充 | `P0-3补充：<一句话说明>` | `P0-3补充：将 verify_*.dat 逐字节回归比对固化为 CTest` |
| 清理 / 构建 / 纯工程改动 | 一句话中文，多项用 ` + ` 连接 | `清理无用文件 + 修复 CHOLMOD 链接 + 构建目录统一` |
| 文档 | `docs: <说明>` | `docs: fem_3d.hpp 增加一行说明注释` |
| 回滚 | `Revert "<原提交说明>"` | `Revert "docs: ..."` |

说明：

- P 级编号对应代码评审报告中的 P0/P1/P2/P3 问题编号，一项改动尽量只对应一个逻辑变更。
- 说明用一句话讲清「改了什么、为什么」，不超过一行；必要时把关键位置（文件/函数）写进括号。
- 中英文冒号 `:` / `：` 均可，但同一提交内保持一致。

## 提交前检查

1. 构建通过：

   ```bash
   cmake -S . -B build
   cmake --build build -j8
   ```

2. 测试全绿（解析解断言 + .dat 逐字节回归）：

   ```bash
   ctest --test-dir build --output-on-failure
   # 或（不走 CMake 时）
   powershell -ExecutionPolicy Bypass -File run_tests.ps1
   ```

3. 若某次改动**有意**改变了 `.dat` 算例结果，先刷新回归基线再提交：

   ```bash
   python tests/regression/compare_regression.py --freeze \
         --bin build --src . --baseline tests/regression/baselines
   ```

4. 提交：

   ```bash
   git add <改动的文件>
   git commit -m "P3-1: <一句话说明>"
   ```
