$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
# 统一使用构建目录 build（清理后已移除 build_cable/build_debug/build_local/build_make）
$bin = Join-Path $root "build"

$tests = @(
  # 解析解断言：校验数值结果本身，而不只是"程序没崩"
  @{ exe="test_2d.exe"; args=@() },
  @{ exe="test_3d.exe"; args=@() },
  @{ exe="example_1d_bar.exe"; args=@() },
  @{ exe="example_2d_truss.exe"; args=@() },
  @{ exe="example_2d_plane.exe"; args=@() },
  @{ exe="example_3d_tet.exe"; args=@() },
  @{ exe="example_mixed.exe"; args=@((Join-Path $root "verify_ce.dat")) },
  @{ exe="example_mixed.exe"; args=@((Join-Path $root "verify_prestress_equiv.dat")) },
  @{ exe="example_stage.exe"; args=@((Join-Path $root "verify_stage.dat")) },
  @{ exe="example_stage.exe"; args=@((Join-Path $root "verify_stage_advanced.dat")) },
  @{ exe="example_3d_truss.exe"; args=@((Join-Path $root "verify_3d_truss.dat")) },
  @{ exe="example_3d_beam.exe"; args=@((Join-Path $root "verify_3d_beam.dat")) },
  @{ exe="example_3d_stage.exe"; args=@((Join-Path $root "verify_3d_stage.dat")) },
  @{ exe="example_3d_influence.exe"; args=@((Join-Path $root "verify_3d_influence.dat")) },
  @{ exe="example_3d_influence_surface.exe"; args=@((Join-Path $root "verify_3d_influence_surface.dat")) },
  @{ exe="example_3d_grillage.exe"; args=@((Join-Path $root "bridge_3d_grillage.dat")) },
  # 二维影响线：走"一次分解 + 多右端项"路径（P1-1 / P1-7），
  # 这条路径原先硬编码 SimplicialLDLT，改动后必须回归
  @{ exe="example_influence.exe"; args=@((Join-Path $root "verify_influence.dat")) },
  # 三维梁包络：曾因阶段 U 与全模型自由度编号错配而越界 abort
  @{ exe="example_3d_beam_envelope.exe"; args=@((Join-Path $root "verify_3d_beam_stage.dat")) },
  # ---- 补充覆盖（复审 P1-3）：之前漏测的目标 ----
  @{ exe="example_modal.exe"; args=@((Join-Path $root "verify_modal.dat")) },
  @{ exe="example_3d_modal.exe"; args=@((Join-Path $root "verify_3d_modal.dat")) },
  @{ exe="example_3d_beam_live.exe"; args=@((Join-Path $root "verify_3d_beam_live.dat")) },
  @{ exe="example_3d_stage_live.exe"; args=@((Join-Path $root "verify_3d_stage_live.dat")) },
  @{ exe="example_3d_stage_envelope.exe"; args=@((Join-Path $root "verify_3d_stage.dat")) },
  @{ exe="example_3d_beam_stage.exe"; args=@((Join-Path $root "verify_3d_beam_stage.dat")) },
  @{ exe="example_stage_envelope.exe"; args=@((Join-Path $root "verify_stage.dat")) },
  @{ exe="example_viz.exe"; args=@((Join-Path $root "mixed_truss_beam.dat")) },
  @{ exe="example_stage_viz.exe"; args=@((Join-Path $root "verify_stage.dat")) },
  @{ exe="example_solvers.exe"; args=@() },
  @{ exe="example_cholmod.exe"; args=@() },
  @{ exe="example_amgcl.exe"; args=@() },
  @{ exe="fem_solver.exe"; args=@() },
  # 性能回归：模态/屈曲耗时阈值（防稠密 O(n^3) 回潮）
  @{ exe="modal_bench.exe"; args=@("1000") }
)

$failed = 0
foreach ($t in $tests) {
    $p = Join-Path $bin $t.exe
    if (-not (Test-Path $p)) {
        Write-Host "[MISS] $($t.exe)"
        $failed++
        continue
    }
    & $p $t.args | Out-Null
    if ($LASTEXITCODE -eq 0) {
        Write-Host "[PASS] $($t.exe) $($t.args -join ' ')"
    } else {
        Write-Host "[FAIL] $($t.exe) $($t.args -join ' ')"
        $failed++
    }
}

# ---- .dat byte-for-byte regression vs golden baselines in tests/regression/baselines ----
$regScript = Join-Path $root "tests\regression\compare_regression.py"
$baseline  = Join-Path $root "tests\regression\baselines"
if (Test-Path $regScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== .dat byte-for-byte regression =="
        & python $regScript --bin $bin --src $root --baseline $baseline
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] .dat byte regression"
        } else {
            Write-Host "[FAIL] .dat byte regression"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip .dat byte regression"
        $failed++
    }
} else {
    Write-Host "[MISS] regression script missing: $regScript"
    $failed++
}

# ---- ANSYS 对标基线 vs tests/benchmark/ansys/reference ----
$ansysScript = Join-Path $root "tests\benchmark\ansys\compare_ansys.py"
$ansysRef    = Join-Path $root "tests\benchmark\ansys\reference"
if (Test-Path $ansysScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== ANSYS benchmark =="
        & python $ansysScript --bin $bin --src $root --reference $ansysRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] ANSYS benchmark"
        } else {
            Write-Host "[FAIL] ANSYS benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip ANSYS benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] ANSYS benchmark script missing: $ansysScript"
    $failed++
}

# ---- 真实桥梁 ANSYS 对标基线 vs tests/benchmark/bridge/reference ----
$bridgeScript = Join-Path $root "tests\benchmark\bridge\compare_bridge.py"
$bridgeRef    = Join-Path $root "tests\benchmark\bridge\reference"
if (Test-Path $bridgeScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Bridge ANSYS benchmark =="
        & python $bridgeScript --bin $bin --src $root --reference $bridgeRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Bridge ANSYS benchmark"
        } else {
            Write-Host "[FAIL] Bridge ANSYS benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip Bridge ANSYS benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] Bridge benchmark script missing: $bridgeScript"
    $failed++
}

if ($failed -eq 0) { Write-Host "ALL TESTS PASSED" } else { Write-Host "$failed TEST(S) FAILED"; exit 1 }
