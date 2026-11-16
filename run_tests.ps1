param(
    [string]$BuildDir = "build"
)

$env:PATH = "C:\msys64\mingw64\bin;C:\msys64\usr\bin;$env:PATH"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
# 默认仍是 build；干净 Ninja 构建在别的目录时用 -BuildDir 指定。
# 离线 benchmark 的「同二进制重放逐行一致」判据要求 exe 与冻结证据同源，
# 因此不能让本脚本悄悄指向修复前的旧 build/。
if ([System.IO.Path]::IsPathRooted($BuildDir)) {
    $bin = $BuildDir
} else {
    $bin = Join-Path $root $BuildDir
}

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

# ---- 非线性成桥态 vs ANSYS NLGEOM + LINK180（CB3D10）----
# cable_solver.exe 重跑 .dat，与冻结自真实 ANSYS 18.0 批处理的 325 项
# 黄金值逐项比对；默认离线，不需要 ANSYS 许可。
$nlScript = Join-Path $root "tests\benchmark\nonlinear_ansys_3d\compare_nonlinear_ansys_3d.py"
$nlCases  = Join-Path $root "tests\benchmark\nonlinear_ansys_3d\cases"
$nlRef    = Join-Path $root "tests\benchmark\nonlinear_ansys_3d\reference"
if (Test-Path $nlScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Nonlinear completion vs ANSYS NLGEOM benchmark =="
        & python $nlScript --bin $bin --cases $nlCases --reference $nlRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Nonlinear completion vs ANSYS NLGEOM benchmark"
        } else {
            Write-Host "[FAIL] Nonlinear completion vs ANSYS NLGEOM benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip nonlinear ANSYS benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] nonlinear ANSYS benchmark script missing: $nlScript"
    $failed++
}

# ---- 悬索桥找形 vs ANSYS NLGEOM + LINK180（CB3D11 地锚式 / CB3D12 自锚式）----
# 共用参数只改 SYSTEM；重跑 cable_solver.exe 与 574 项冻结黄金值逐项比对。
# 默认离线，不需要 ANSYS 许可。
$suScript = Join-Path $root "tests\benchmark\suspension_ansys_3d\compare_suspension_ansys_3d.py"
$suCases  = Join-Path $root "tests\benchmark\suspension_ansys_3d\cases"
$suRef    = Join-Path $root "tests\benchmark\suspension_ansys_3d\reference"
if (Test-Path $suScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Suspension bridge form-finding vs ANSYS NLGEOM benchmark =="
        & python $suScript --bin $bin --cases $suCases --reference $suRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Suspension bridge ANSYS NLGEOM benchmark"
        } else {
            Write-Host "[FAIL] Suspension bridge ANSYS NLGEOM benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip suspension ANSYS benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] suspension ANSYS benchmark script missing: $suScript"
    $failed++
}

$seScript = Join-Path $root "tests\benchmark\suspension_erection_3d\compare_suspension_erection_3d.py"
$seCases  = Join-Path $root "tests\benchmark\suspension_erection_3d\cases"
$seRef    = Join-Path $root "tests\benchmark\suspension_erection_3d\reference"
if (Test-Path $seScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Balanced self-anchored suspension + erection stages benchmark =="
        & python $seScript --bin $bin --cases $seCases --reference $seRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Balanced suspension + erection stages benchmark"
        } else {
            Write-Host "[FAIL] Balanced suspension + erection stages benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip balanced suspension benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] balanced suspension benchmark script missing: $seScript"
    $failed++
}

$sgScript = Join-Path $root "tests\benchmark\cable_sag_ansys_3d\compare_cable_sag_ansys_3d.py"
$sgCases  = Join-Path $root "tests\benchmark\cable_sag_ansys_3d\cases"
$sgRef    = Join-Path $root "tests\benchmark\cable_sag_ansys_3d\reference"
if (Test-Path $sgScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Cable sag / Ernst modulus + W>0 empty-cable state benchmark =="
        & python $sgScript --bin $bin --cases $sgCases --reference $sgRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Cable sag / Ernst benchmark"
        } else {
            Write-Host "[FAIL] Cable sag / Ernst benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip cable sag benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] cable sag benchmark script missing: $sgScript"
    $failed++
}

$bnScript = Join-Path $root "tests\benchmark\beam_nl_ansys_3d\compare_beam_nl_ansys_3d.py"
$bnCases  = Join-Path $root "tests\benchmark\beam_nl_ansys_3d\cases"
$bnRef    = Join-Path $root "tests\benchmark\beam_nl_ansys_3d\reference"
if (Test-Path $bnScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Large-rotation beam-column (CR kernel) vs ANSYS BEAM188 benchmark =="
        & python $bnScript --bin $bin --cases $bnCases --reference $bnRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Large-rotation beam-column benchmark"
        } else {
            Write-Host "[FAIL] Large-rotation beam-column benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip beam-column benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] beam-column benchmark script missing: $bnScript"
    $failed++
}

$stScript = Join-Path $root "tests\benchmark\snap_through_ansys_3d\compare_snap_through_ansys_3d.py"
$stCases  = Join-Path $root "tests\benchmark\snap_through_ansys_3d\cases"
$stRef    = Join-Path $root "tests\benchmark\snap_through_ansys_3d\reference"
if (Test-Path $stScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Post-buckling / snap-through (arc-length) vs ANSYS ARCLEN benchmark =="
        & python $stScript --bin $bin --cases $stCases --reference $stRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Snap-through arc-length benchmark"
        } else {
            Write-Host "[FAIL] Snap-through arc-length benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip snap-through benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] snap-through benchmark script missing: $stScript"
    $failed++
}

# ---- B3 主链路 [STABILITY] opt-in（第十七次）---------------------------
$spScript = Join-Path $root "tests\benchmark\stable_path_3d\compare_stable_path_3d.py"
$spCases  = Join-Path $root "tests\benchmark\stable_path_3d\cases"
$spRef    = Join-Path $root "tests\benchmark\stable_path_3d\reference"
if (Test-Path $spScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== B3 STABILITY=1 main-chain arc-length vs CB3D17 frozen path =="
        & python $spScript --bin $bin --cases $spCases --reference $spRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] STABILITY main-chain benchmark"
        } else {
            Write-Host "[FAIL] STABILITY main-chain benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip STABILITY benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] STABILITY benchmark script missing: $spScript"
    $failed++
}

# ---- CR 梁刚接框架后屈曲（CB3D19 / 第十八次）---------------------------
$fpScript = Join-Path $root "tests\benchmark\frame_postbuckling_3d\compare_frame_postbuckling_3d.py"
$fpCases  = Join-Path $root "tests\benchmark\frame_postbuckling_3d\cases"
$fpRef    = Join-Path $root "tests\benchmark\frame_postbuckling_3d\reference"
if (Test-Path $fpScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== CR beam frame post-buckling vs ANSYS BEAM188 ARCLEN benchmark =="
        & python $fpScript --bin $bin --cases $fpCases --reference $fpRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Frame post-buckling benchmark"
        } else {
            Write-Host "[FAIL] Frame post-buckling benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip frame post-buckling benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] frame post-buckling benchmark script missing: $fpScript"
    $failed++
}

# ---- 成桥态起点全过程稳定（CB3D20 / 第十九次 ①）------------------------
$csScript = Join-Path $root "tests\benchmark\completion_stability_3d\compare_completion_stability_3d.py"
$csCases  = Join-Path $root "tests\benchmark\completion_stability_3d\cases"
$csRef    = Join-Path $root "tests\benchmark\completion_stability_3d\reference"
if (Test-Path $csScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Completion-state arc-length stability benchmark (FEM self-consistency + closed-form oracle) =="
        & python $csScript --bin $bin --cases $csCases --reference $csRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Completion-state stability benchmark"
        } else {
            Write-Host "[FAIL] Completion-state stability benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip completion-state stability benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] completion-state stability benchmark script missing: $csScript"
    $failed++
}

# ---- 真实斜拉桥成桥态两段对标（CB3D21 / 第十九次 ②）----------------------
$cbScript = Join-Path $root "tests\benchmark\completion_bridge_3d\compare_completion_bridge_3d.py"
$cbCases  = Join-Path $root "tests\benchmark\completion_bridge_3d\cases"
$cbRef    = Join-Path $root "tests\benchmark\completion_bridge_3d\reference"
if (Test-Path $cbScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Completion-state bridge stability vs ANSYS two-step benchmark =="
        & python $cbScript --bin $bin --cases $cbCases --reference $cbRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Completion-state bridge two-step benchmark"
        } else {
            Write-Host "[FAIL] Completion-state bridge two-step benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip completion-state bridge benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] completion-state bridge benchmark script missing: $cbScript"
    $failed++
}

# ---- MODE=COMPLETION 成桥态基准同源修复（CB3D10 / CB3D13，第十九次 ③）-----
$bsScript = Join-Path $root "tests\benchmark\completion_basis_3d\compare_completion_basis_3d.py"
$bsCases  = Join-Path $root "tests\benchmark\completion_basis_3d\cases"
$bsRef    = Join-Path $root "tests\benchmark\completion_basis_3d\reference"
if (Test-Path $bsScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Completion-state basis-model fix (COMPLETION vs STATIC cross-MODE anchor) =="
        & python $bsScript --bin $bin --cases $bsCases --reference $bsRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Completion-state basis-model benchmark"
        } else {
            Write-Host "[FAIL] Completion-state basis-model benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip completion-state basis-model benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] completion-state basis-model benchmark script missing: $bsScript"
    $failed++
}

# ---- 真实桥梁越过极值点 / 下降段（CB3D22 / 第十九次 ③a）-----------------
$lpScript = Join-Path $root "tests\benchmark\bridge_limit_point_3d\compare_bridge_limit_point_3d.py"
$lpCases  = Join-Path $root "tests\benchmark\bridge_limit_point_3d\cases"
$lpRef    = Join-Path $root "tests\benchmark\bridge_limit_point_3d\reference"
if (Test-Path $lpScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Real-bridge limit point / descending branch (FEM vs ANSYS ARCLEN, offline) =="
        & python $lpScript --bin $bin --cases $lpCases --reference $lpRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Bridge limit-point benchmark"
        } else {
            Write-Host "[FAIL] Bridge limit-point benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip bridge limit-point benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] bridge limit-point benchmark script missing: $lpScript"
    $failed++
}

# ---- 施工阶段中间态起点的全过程稳定（CB3D23 / 第二十次 ③b）-----------
# STABASIS=STAGE STAGEID=2：起点是 solveModel3DStagesNonlinear() 的阶段 2
# 共转列式全量平衡态（u0 = u_stage、F0 = 该阶段累计恒载、Finc 只来自
# *STABILITY_LOAD），BASIS=stage_state 与 completion_state 分列。
# 15 份冻结证据 + 两个 ANSYS 18.0 job（0 error）纯离线比对，
# 含负例 ERROR=stage_not_found 与当前二进制重放逐字节一致。
$ssScript = Join-Path $root "tests\benchmark\stage_stability_3d\compare_stage_state_3d.py"
$ssCases  = Join-Path $root "tests\benchmark\stage_stability_3d\cases"
$ssRef    = Join-Path $root "tests\benchmark\stage_stability_3d\reference"
if (Test-Path $ssScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Construction-stage intermediate-state stability (FEM vs ANSYS, offline) =="
        & python $ssScript --bin $bin --cases $ssCases --reference $ssRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Stage-state stability benchmark"
        } else {
            Write-Host "[FAIL] Stage-state stability benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip stage-state stability benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] stage-state stability benchmark script missing: $ssScript"
    $failed++
}

# ---- Wempner–Riks 正交轨迹（PATHMETHOD=RIKS，仅 STAGE 起点 opt-in）----
# 每轮校正增量与上一收敛步正交，dλ 由线性方程唯一确定；与球面弧长
# 跟踪同一 R(u,λ)=0，作「同问题、双方法」交叉验证。闭式 von Mises
# 两杆 λ(v) 逐点核验 + CB3D23 跨基准终点 + 未知方法负例 + 二进制重放。
$rkScript = Join-Path $root "tests\benchmark\riks_stage_stability_3d\compare_riks_stage_3d.py"
$rkCases  = Join-Path $root "tests\benchmark\riks_stage_stability_3d\cases"
$rkRef    = Join-Path $root "tests\benchmark\riks_stage_stability_3d\reference"
if (Test-Path $rkScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Wempner-Riks orthogonal path at STAGE start (offline) =="
        & python $rkScript --bin $bin --cases $rkCases --reference $rkRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Riks stage-state stability benchmark"
        } else {
            Write-Host "[FAIL] Riks stage-state stability benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip Riks stage-state stability benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] Riks stage-state stability benchmark script missing: $rkScript"
    $failed++
}

# ---- 第二十一次力学①：位移控制 / 线搜索（opt-in，默认全关）----
# 控制 DOF 位移增量驱动、λ 由线性约束唯一确定；Armijo 回溯线搜索保持位移
# 约束、不改变平衡路径。von Mises 两杆闭式 P(v) 逐点锚点 + 线搜索/无线搜索
# 同路径 + 细步 λ_max→λ_cr + 奇异控制 DOF 负例 + 二进制逐字节重放。
$dcScript = Join-Path $root "tests\benchmark\displacement_control_3d\compare_displacement_control_3d.py"
$dcCases  = Join-Path $root "tests\benchmark\displacement_control_3d\cases"
$dcRef    = Join-Path $root "tests\benchmark\displacement_control_3d\reference"
if (Test-Path $dcScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Displacement control + line search (offline) =="
        & python $dcScript --bin $bin --cases $dcCases --reference $dcRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Displacement control + line search benchmark"
        } else {
            Write-Host "[FAIL] Displacement control + line search benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip displacement control benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] displacement control benchmark script missing: $dcScript"
    $failed++
}

# ---- 第二十一次力学②：初始缺陷形态 / 幅值 / 缺陷敏感性（opt-in，默认全关）----
# 线性屈曲模态生成无应力缺陷场；幅值→λ_max 单调 + 反对称模态对称破缺 +
# 独立闭式 P(v) 锚点 + IMPMODE 越界负例 + IMPERF=0 零污染 + 二进制逐字节重放。
$ipScript = Join-Path $root "tests\benchmark\initial_imperfection_3d\compare_imperfection_3d.py"
$ipCases  = Join-Path $root "tests\benchmark\initial_imperfection_3d\cases"
$ipRef    = Join-Path $root "tests\benchmark\initial_imperfection_3d\reference"
if (Test-Path $ipScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Initial imperfection + sensitivity (offline) =="
        & python $ipScript --bin $bin --cases $ipCases --reference $ipRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Initial imperfection + sensitivity benchmark"
        } else {
            Write-Host "[FAIL] Initial imperfection + sensitivity benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip initial imperfection benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] initial imperfection benchmark script missing: $ipScript"
    $failed++
}

# ---- 第二十一次力学③：材料非线性本构接入 3D CR 梁（opt-in，默认 ELASTIC）----
# 单轴应力-应变平台 + 矩形截面塑性铰逐阶段屈服 + 混凝土受拉开裂 +
# 弹性段纤维 vs Kl12 逐位对齐 + 独立 Python 纤维积分锚点 + 二进制逐字节重放。
$mnScript = Join-Path $root "tests\benchmark\material_nonlinear_3d\compare_material_nl_3d.py"
$mnCases  = Join-Path $root "tests\benchmark\material_nonlinear_3d\cases"
$mnRef    = Join-Path $root "tests\benchmark\material_nonlinear_3d\reference"
if (Test-Path $mnScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Material nonlinear CR beam (offline) =="
        & python $mnScript --bin $bin --cases $mnCases --reference $mnRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Material nonlinear CR beam benchmark"
        } else {
            Write-Host "[FAIL] Material nonlinear CR beam benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip material nonlinear benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] material nonlinear benchmark script missing: $mnScript"
    $failed++
}

# ---- 第二十一次力学④：分岔点检测与分支切换（opt-in，默认 BIFURC=0）----
# 高拱 von Mises 桁架反对称分岔闭式锚点 + pivotRatio/一致切线特征值确认
# + 模态扰动切换分支 + 三类负例 + 二进制逐字节重放。FEM 自洽，离线。
$bfScript = Join-Path $root "tests\benchmark\bifurcation_switch_3d\compare_bifurcation_3d.py"
$bfCases  = Join-Path $root "tests\benchmark\bifurcation_switch_3d\cases"
$bfRef    = Join-Path $root "tests\benchmark\bifurcation_switch_3d\reference"
if (Test-Path $bfScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Bifurcation detection & branch switching (offline) =="
        & python $bfScript --bin $bin --cases $bfCases --reference $bfRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Bifurcation switching benchmark"
        } else {
            Write-Host "[FAIL] Bifurcation switching benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip bifurcation benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] bifurcation benchmark script missing: $bfScript"
    $failed++
}

# ---- 第二十一次力学⑤：动力 snap-through（opt-in，默认 DYN=0）----
# 一致/集中质量矩阵 + Newmark-beta 隐式时程；独立 SDOF RK4 时程锚点
# + 低于静力极限的动力跃越 + 荷载时程/阻尼/质量敏感性 + 负例 + 逐字节重放。
# 独立新模块，不改静力链；FEM 自洽，离线。
$dyScript = Join-Path $root "tests\benchmark\dynamic_snapthrough_3d\compare_dynamic_3d.py"
$dyCases  = Join-Path $root "tests\benchmark\dynamic_snapthrough_3d\cases"
$dyRef    = Join-Path $root "tests\benchmark\dynamic_snapthrough_3d\reference"
if (Test-Path $dyScript) {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if ($py) {
        Write-Host ""
        Write-Host "== Dynamic snap-through / Newmark-beta (offline) =="
        & python $dyScript --bin $bin --cases $dyCases --reference $dyRef
        if ($LASTEXITCODE -eq 0) {
            Write-Host "[PASS] Dynamic snap-through benchmark"
        } else {
            Write-Host "[FAIL] Dynamic snap-through benchmark"
            $failed++
        }
    } else {
        Write-Host "[MISS] python not found, skip dynamic benchmark"
        $failed++
    }
} else {
    Write-Host "[MISS] dynamic benchmark script missing: $dyScript"
    $failed++
}

if ($failed -eq 0) { Write-Host "ALL TESTS PASSED" } else { Write-Host "$failed TEST(S) FAILED"; exit 1 }

