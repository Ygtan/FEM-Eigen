$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$bin = Join-Path $root "build_debug"

$tests = @(
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
  @{ exe="example_3d_grillage.exe"; args=@((Join-Path $root "bridge_3d_grillage.dat")) }
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

if ($failed -eq 0) { Write-Host "ALL TESTS PASSED" } else { Write-Host "$failed TEST(S) FAILED"; exit 1 }
