# FEM-Eigen 交付说明

## 项目位置

```text
D:\AI-work\FEM-Eigen
```

## 编译

```powershell
$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
cmake -S "D:\AI-work\FEM-Eigen" -B "D:\AI-work\FEM-Eigen\build_debug" -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build "D:\AI-work\FEM-Eigen\build_debug"
```

## 自动测试

```powershell
powershell -ExecutionPolicy Bypass -File .\run_tests.ps1
```

## 主要文档

- `README.md`
- `input_format.md`
- `input_format_3d.md`
- `CHANGELOG.md`
- `MODIFICATION_LOG.md`

## 核心程序

二维：

- `example_mixed`
- `example_influence`
- `example_viz`
- `example_stage`
- `example_stage_viz`
- `example_stage_envelope`

三维：

- `example_3d_truss`
- `example_3d_beam`
- `example_3d_stage`
- `example_3d_influence`
- `example_3d_influence_surface`
- `example_3d_beam_live`
- `example_3d_grillage`

## 导出目录

```text
C:\Users\user\Documents\Codex\2026-09-15\new-chat\FEM-Eigen-export-final
```
