$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
# 统一使用构建目录 build（清理后已移除 build_cable/build_debug/build_local/build_make）
$bin = Join-Path $root "build"

$tests = @(
  # 解析解断言：校验数值结果本身，而不只是"程序没崩"
  @{ exe="test_analytical.exe"; args=@() },
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
  @{ exe="example_3d_beam_envelope.exe"; args=@((Join-Path $root "verify_3d_beam_stage.dat")) }
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

if ($failed -eq 0) { Write-Host "ALL TESTS PASSED" } else { Write-Host "$failed TEST(S) FAILED"; exit 1 }
