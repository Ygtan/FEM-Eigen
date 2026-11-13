$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
# 统一使用构建目录 build（清理后已移除 build_cable/build_debug/build_local/build_make）
$bin = Join-Path $root "build"

$tests = @(
  # 解析解断言：校验数值结果本身，而不只是"程序没崩"
  @{ exe="test_2d.exe"; args=@() },
  @{ exe="test_3d.exe"; args=@() },
  # 几何非线性（共转列式）：刚体不变性 / 卷圆 / 梁柱 P-Δ / 索无应力长度
  @{ exe="test_3d_nl.exe"; args=@() },
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

# ---- 连续梁算例：ANSYS 18.0 BEAM3 实测冻结值 vs tests/benchmark/continuous_beam ----
# 与上面两个 benchmark 不同，这里 ref 来自真实 ANSYS 批处理导出的 CSV。
# 默认离线运行（只重跑 FEM-Eigen），不依赖 ANSYS 许可；
# 需要端到端商业软件复核时手动加 --live-ansys。
$cbScript = Join-Path $root "tests\benchmark\continuous_beam\compare_continuous_beam.py"
$cbCases  = Join-Path $root "tests\benchmark\continuous_beam\cases"
$cbRef    = Join-Path $root "tests\benchmark\continuous_beam\reference"
if (Test-Path $cbScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Continuous beam ANSYS benchmark =="
        & python $cbScript --bin $bin --cases $cbCases --reference $cbRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Continuous beam ANSYS benchmark"
        } else {
            Write-Host "[FAIL] Continuous beam ANSYS benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip continuous beam benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] Continuous beam benchmark script missing: $cbScript"
    $failed++
}

# ---- 3D 连续箱梁算例：FEM-Eigen 3D vs 冻结的 ANSYS BEAM4 (Euler-Bernoulli) ----
# BEAM188 (Timoshenko) 结果已归档，差异在 *_compare_report.txt [4] 节量化，
# 不作为黄金参考值。默认离线运行，不依赖 ANSYS 许可。
$cb3Script = Join-Path $root "tests\benchmark\continuous_beam_3d\compare_continuous_beam_3d.py"
$cb3Cases  = Join-Path $root "tests\benchmark\continuous_beam_3d\cases"
$cb3Ref    = Join-Path $root "tests\benchmark\continuous_beam_3d\reference"
if (Test-Path $cb3Script) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Continuous beam 3D ANSYS benchmark =="
        & python $cb3Script --bin $bin --cases $cb3Cases --reference $cb3Ref
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Continuous beam 3D ANSYS benchmark"
        } else {
            Write-Host "[FAIL] Continuous beam 3D ANSYS benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip continuous beam 3D benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] Continuous beam 3D benchmark script missing: $cb3Script"
    $failed++
}

# ---- 3D 斜拉桥/拱桥混合单元：FEM-Eigen 3D vs 冻结 ANSYS BEAM4+LINK180 ----
# 同一 ModelSpec3D 里 BEAM + CABLE/TRUSS；默认离线运行，不依赖 ANSYS 许可。
$caScript = Join-Path $root "tests\benchmark\cable_arch_3d\compare_cable_arch_3d.py"
$caCases  = Join-Path $root "tests\benchmark\cable_arch_3d\cases"
$caRef    = Join-Path $root "tests\benchmark\cable_arch_3d\reference"
if (Test-Path $caScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Cable/arch 3D mixed-element ANSYS benchmark =="
        & python $caScript --bin $bin --cases $caCases --reference $caRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Cable/arch 3D mixed-element ANSYS benchmark"
        } else {
            Write-Host "[FAIL] Cable/arch 3D mixed-element ANSYS benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip cable/arch 3D benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] Cable/arch 3D benchmark script missing: $caScript"
    $failed++
}

# ---- 非线性成桥态 / 索长迭代 smoke：cable_solver.exe (fem_3d_nl.hpp) ----
# 共转列式几何非线性链，无 ANSYS 基准；检查收敛 / 几何闭合 / 索长闭合
# 以及与冻结 *_cable.txt 的一致性。不需要 ANSYS 许可。
$csScript = Join-Path $root "tests\benchmark\cable_arch_3d\compare_cable_solver_smoke.py"
$csSmoke  = Join-Path $root "tests\benchmark\cable_arch_3d\nonlinear_smoke"
if (Test-Path $csScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== cable_solver nonlinear completion smoke =="
        & python $csScript --bin $bin --smoke $csSmoke
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] cable_solver nonlinear completion smoke"
        } else {
            Write-Host "[FAIL] cable_solver nonlinear completion smoke"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip cable_solver smoke"
        $failed++
    }
} else {
    Write-Host "[MISS] cable_solver smoke script missing: $csScript"
    $failed++
}

# ---- 3D 分布荷载/自重：FEM-Eigen 3D vs 冻结 ANSYS BEAM4 ----
# *ELEMENT_LOAD / *GRAVITY / *SELF_WEIGHT 一致节点荷载与固端力扣除；
# CB3D06/07/08 同源导出，默认离线运行，不依赖 ANSYS 许可。
$elScript = Join-Path $root "tests\benchmark\element_loads_3d\compare_element_loads_3d.py"
$elCases  = Join-Path $root "tests\benchmark\element_loads_3d\cases"
$elRef    = Join-Path $root "tests\benchmark\element_loads_3d\reference"
if (Test-Path $elScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== 3D element-load / selfweight ANSYS benchmark =="
        & python $elScript --bin $bin --cases $elCases --reference $elRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] 3D element-load / selfweight ANSYS benchmark"
        } else {
            Write-Host "[FAIL] 3D element-load / selfweight ANSYS benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip element-load 3D benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] element-load 3D benchmark script missing: $elScript"
    $failed++
}

# ---- 3D 分布荷载内核级闭式解 / 叠加恒等式验证 ----
$elkScript = Join-Path $root "tests\benchmark\element_loads_3d\verify_kernel_loads.py"
if (Test-Path $elkScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== 3D distributed-load kernel validation =="
        & python $elkScript --bin $bin
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] 3D distributed-load kernel validation"
        } else {
            Write-Host "[FAIL] 3D distributed-load kernel validation"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip 3D distributed-load kernel validation"
        $failed++
    }
} else {
    Write-Host "[MISS] 3D distributed-load kernel validation script missing: $elkScript"
    $failed++
}

if ($failed -eq 0) { Write-Host "ALL TESTS PASSED" } else { Write-Host "$failed TEST(S) FAILED"; exit 1 }
