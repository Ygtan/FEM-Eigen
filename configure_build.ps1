<#
.SYNOPSIS
  Reproducible cross-machine setup of the FEM-Eigen build dir (Ninja + pinned mingw64 g++).

.DESCRIPTION
  Fixes the "stale cross-machine build/ cache makes cmake --build hang / configure fail" defect:
   1. Auto-detect cmake / ninja (PATH first, then MSYS2 standard locations) --
      never hardcode machine-specific paths like C:\Program Files\CMake.
   2. Pin compiler to C:/msys64/mingw64/bin/g++ (frozen-evidence byte-replay needs
      the same toolchain; never use ucrt64).
   3. If build/ cached a non-Ninja generator or a hardcoded dead cmake path,
      clean it and rebuild.

  Usage:
     powershell -ExecutionPolicy Bypass -File configure_build.ps1
     powershell -ExecutionPolicy Bypass -File configure_build.ps1 -Jobs 4

  After switching machines (or a clean checkout), run once; no need to remember
  build commands or hand-sync CTestTestfile.cmake.
#>
param(
    [int]$Jobs = 8
)

$ErrorActionPreference = "Stop"
$root  = Split-Path -Parent $MyInvocation.MyCommand.Path
$build = Join-Path $root "build"

function Resolve-Tool {
    param([string[]]$Names, [string[]]$ExtraPaths)
    foreach ($n in $Names) {
        $c = Get-Command $n -ErrorAction SilentlyContinue
        if ($c) { return $c.Source }
    }
    foreach ($p in $ExtraPaths) {
        if (Test-Path $p) { return $p }
    }
    return $null
}

# Put mingw64 first so PATH-based detection matches the pinned compiler.
if (Test-Path "C:\msys64\mingw64\bin") {
    $env:PATH = "C:\msys64\mingw64\bin;$env:PATH"
}
$cmake = Resolve-Tool -Names @("cmake") -ExtraPaths @(
    "C:\msys64\mingw64\bin\cmake.exe",
    "C:\msys64\usr\bin\cmake.exe"
)
$ninja = Resolve-Tool -Names @("ninja") -ExtraPaths @(
    "C:\msys64\mingw64\bin\ninja.exe"
)
$gpp = Resolve-Tool -Names @() -ExtraPaths @(
    "C:\msys64\mingw64\bin\g++.exe"
)

if (-not $cmake) { throw "cmake not found: install MSYS2/CMake or add it to PATH" }
if (-not $ninja) { throw "ninja not found: pacman -S mingw-w64-x86_64-ninja" }
if (-not $gpp)   { throw "mingw64 g++ not found: frozen evidence requires the mingw64 toolchain, not ucrt64" }

Write-Host "[env] cmake = $cmake"
Write-Host "[env] ninja = $ninja"
Write-Host "[env] g++   = $gpp"

$cache    = Join-Path $build "CMakeCache.txt"
$makefile = Join-Path $build "Makefile"
$bad = $false
if (Test-Path $cache) {
    $genLine = (Select-String -Path $cache -Pattern '^CMAKE_GENERATOR:INTERNAL=' | Select-Object -First 1).Line
    if ($genLine -and $genLine -notmatch 'Ninja') {
        $bad = $true
        Write-Host "[env] non-Ninja generator detected: $genLine"
    }
}
if ((Test-Path $makefile) -and (Select-String -Path $makefile -SimpleMatch 'Program Files' -Quiet)) {
    $bad = $true
    Write-Host "[env] Makefile has a hardcoded dead cmake path (Program Files)"
}
if ($bad) {
    Write-Host "[env] cleaning stale cross-machine build/ cache before rebuild"
    Remove-Item $cache -Force -ErrorAction SilentlyContinue
    Remove-Item $makefile -Force -ErrorAction SilentlyContinue
    Remove-Item (Join-Path $build "CMakeFiles") -Recurse -Force -ErrorAction SilentlyContinue
}

& $cmake -S $root -B $build -G Ninja `
    "-DCMAKE_MAKE_PROGRAM=$ninja" `
    "-DCMAKE_CXX_COMPILER=$gpp" `
    "-DCMAKE_BUILD_TYPE=Release" `
    "-DFEM_USE_CHOLMOD=ON"
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }

& $cmake --build $build -j $Jobs
if ($LASTEXITCODE -ne 0) { throw "build failed" }

Write-Host "[env] build/ ready: Ninja + mingw64 g++ (auto-detected, no hardcoded paths)"
Write-Host "[env] regression: ctest --test-dir $build --output-on-failure ; or run_tests.ps1 -BuildDir $build"
