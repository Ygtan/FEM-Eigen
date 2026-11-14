cmake_minimum_required(VERSION 3.10)
project(FEM_Solver LANGUAGES CXX)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 未指定构建类型时默认 Release（避免空优化级别导致性能骤降）。
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# MinGW 下静态链接 C++/GCC 运行库，使示例 exe 不依赖 libstdc++-6.dll 等
# 运行库 DLL，双击或复制到其它机器也能直接运行。
if(MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

# ---- Eigen3（纯头文件库）----
# 优先使用系统安装的 Eigen3 CMake 包；找不到时使用项目内自带的
# external/eigen3 头文件（Eigen 3.3.9）。
find_package(Eigen3 QUIET NO_MODULE)

# 若未找到系统 Eigen3，则构造一个指向自带头文件的 INTERFACE 目标。
if(NOT TARGET Eigen3::Eigen)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3/Eigen/Core")
        message(FATAL_ERROR "未找到 Eigen3，且 external/eigen3 中缺少自带头文件。")
    endif()
    add_library(eigen3_vendored INTERFACE)
    target_include_directories(eigen3_vendored INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3")
    message(STATUS "使用项目内自带的 Eigen3 头文件")
endif()

# ---- AMGCL（可选的高性能代数多重网格后端，头文件已自带）----
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/amgcl/amgcl/amg.hpp")
    set(FEM_AMGCL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/amgcl")
    message(STATUS "使用项目内自带的 AMGCL 头文件")
endif()

# ---- Spectra（可选的稀疏特征值库，头文件已自带）----
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/spectra/Spectra/SymEigsSolver.h")
    set(FEM_SPECTRA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/spectra")
    message(STATUS "使用项目内自带的 Spectra 头文件")
endif()

# ---- SuiteSparse CHOLMOD（可选，通过 MSYS2 pacman / vcpkg / apt 安装）----
# 优先尝试 CMake CONFIG 模式（ SuiteSparse 5.4+ / vcpkg 提供），
# 失败则回退到手动 find_path + find_library（支持 MSYS2 / MinGW 等）。
option(FEM_USE_CHOLMOD "启用 SuiteSparse CHOLMOD 后端（未安装时自动禁用）" ON)

set(FEM_CHOLMOD_READY FALSE)
if(FEM_USE_CHOLMOD)
    # 尝试 1：CONFIG 模式（最干净，无硬编码路径）。
    find_package(SuiteSparse CONFIG QUIET)
    if(SuiteSparse_FOUND AND TARGET SuiteSparse::CHOLMOD)
        set(FEM_CHOLMOD_READY TRUE)
        set(FEM_CHOLMOD_CONFIG_MODE TRUE)
        message(STATUS "找到 CHOLMOD：通过 SuiteSparse CONFIG（${SuiteSparse_DIR}）")
    else()
        # 尝试 2：手动定位头文件 + 库文件。
        # 保留 MSYS2 常用路径作为 HINTS（而非唯一的 PATHS），仍允许系统标准路径参与搜索。
        find_path(CHOLMOD_INCLUDE_DIR cholmod.h
            HINTS "C:/msys64/mingw64/include/suitesparse" "C:/msys64/ucrt64/include/suitesparse"
                  "$ENV{SUITESPARSE_ROOT}/include/suitesparse")
        if(CHOLMOD_INCLUDE_DIR)
            get_filename_component(_cholmod_inc_parent "${CHOLMOD_INCLUDE_DIR}" DIRECTORY)
            get_filename_component(_cholmod_prefix "${_cholmod_inc_parent}" DIRECTORY)
            set(_cholmod_lib_hints "${_cholmod_prefix}/lib")

            set(FEM_CHOLMOD_LIBS "")
            set(_cholmod_missing "")
            foreach(_lib cholmod suitesparseconfig amd colamd camd ccolamd openblas)
                find_library(CHOLMOD_LIB_${_lib} NAMES ${_lib} HINTS "${_cholmod_lib_hints}")
                if(CHOLMOD_LIB_${_lib})
                    list(APPEND FEM_CHOLMOD_LIBS "${CHOLMOD_LIB_${_lib}}")
                else()
                    list(APPEND _cholmod_missing "${_lib}")
                endif()
            endforeach()

            if(_cholmod_missing)
                message(STATUS "CHOLMOD: 缺少库 ${_cholmod_missing}，已禁用 CHOLMOD 后端")
            else()
                set(FEM_CHOLMOD_READY TRUE)
                message(STATUS "找到 CHOLMOD: 头文件 ${CHOLMOD_INCLUDE_DIR}，库目录 ${_cholmod_lib_hints}")
            endif()
        else()
            message(STATUS "未找到 CHOLMOD 头文件，已禁用 CHOLMOD 后端")
        endif()
    endif()
endif()

# 添加可执行目标的便捷函数。
function(fem_add_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE fem_2d fem_3d)
    if(FEM_AMGCL_INCLUDE_DIR)
        target_include_directories(${target} PRIVATE "${FEM_AMGCL_INCLUDE_DIR}")
        target_compile_definitions(${target} PRIVATE AMGCL_NO_BOOST)
    endif()
endfunction()

# 仅为真正使用 CHOLMOD 的目标链接 SuiteSparse（避免其它示例产生 DLL 依赖）。
# CONFIG 模式下直接用 imported target；手动模式下使用绝对路径库。
function(fem_link_cholmod target)
    if(NOT FEM_CHOLMOD_READY)
        return()
    endif()
    if(FEM_CHOLMOD_CONFIG_MODE)
        target_link_libraries(${target} PRIVATE SuiteSparse::CHOLMOD)
    else()
        target_include_directories(${target} PRIVATE "${CHOLMOD_INCLUDE_DIR}")
        target_link_libraries(${target} PRIVATE ${FEM_CHOLMOD_LIBS})
    endif()
endfunction()

# ---- 测试 ----
# enable_testing() 后可用 `ctest` 跑回归；解析解断言测试会校验数值结果，
# 而不只是检查程序有没有崩溃。
enable_testing()

# 解析解断言测试：拆成 2D 与 3D 两个可执行，避免一个维度崩溃吞掉另一个维度的结果。
fem_add_executable(test_2d test_2d.cpp)
fem_add_executable(test_3d test_3d.cpp)
fem_add_executable(test_3d_nl test_3d_nl.cpp)
# P3-3b：完整形状敏感度 K_shape 的独立验收（差分对拍 + 三路径逐轮误差）。
fem_add_executable(test_shape_sensitivity_3d test_shape_sensitivity_3d.cpp)
add_test(NAME analytical_exact_solutions_2d COMMAND test_2d)
add_test(NAME analytical_exact_solutions_3d COMMAND test_3d)
add_test(NAME nonlinear_co_rotational_3d COMMAND test_3d_nl)
add_test(NAME shape_sensitivity_3d COMMAND test_shape_sensitivity_3d)

# .dat 逐字节回归比对：用当前构建重算 verify_*.dat，并与其黄金基线逐字节比对。
# 依赖 Python3；未安装时仅跳过该测试，不影响编译与其它测试。
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME dat_byte_regression
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/compare_regression.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --src ${CMAKE_CURRENT_SOURCE_DIR}
                     --baseline ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/baselines)
    set_tests_properties(dat_byte_regression PROPERTIES
        TIMEOUT 900
        LABELS "regression")

    # ANSYS 对标基线：用 ANSYS BEAM4 对应的欧拉-伯努利闭式解做黄金参考，
    # 比对 100m 简支梁的跨中挠度/支座反力/跨中弯矩/跨中应力。
    add_test(NAME ansys_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/ansys/compare_ansys.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --src ${CMAKE_CURRENT_SOURCE_DIR}
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/ansys/reference)
    set_tests_properties(ansys_benchmark PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 真实桥梁 ANSYS 对标基线：80m 系杆拱桥的节点位移 + 单元内力。
    add_test(NAME bridge_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bridge/compare_bridge.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --src ${CMAKE_CURRENT_SOURCE_DIR}
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bridge/reference)
    set_tests_properties(bridge_benchmark PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 连续梁算例：BridgeStudio 同源 .dat/.inp + ANSYS 18.0 BEAM3 批处理实测值。
    # 与 ansys_benchmark 的区别是参考值取自真实 ANSYS 导出的 CSV（冻结在 cases/），
    # 不是用闭式解代替；默认离线运行，不需要 ANSYS 许可。
    add_test(NAME continuous_beam_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/continuous_beam/compare_continuous_beam.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/continuous_beam/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/continuous_beam/reference)
    set_tests_properties(continuous_beam_benchmark PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 3D 连续箱梁算例：FEM-Eigen 3D (example_3d_report.exe) 的 12x12
    # Euler-Bernoulli 空间梁 vs 冻结的 ANSYS 18.0 BEAM4 (同为 EB) 实测值。
    # BEAM188 是 Timoshenko 单元，差异单列在 *_compare_report.txt，不作黄金参考。
    # 默认离线运行，不需要 ANSYS 许可。
    add_test(NAME continuous_beam_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/continuous_beam_3d/compare_continuous_beam_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/continuous_beam_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/continuous_beam_3d/reference)
    set_tests_properties(continuous_beam_3d_benchmark PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 3D 斜拉桥 / 拱桥混合单元算例：FEM-Eigen 3D 线性内核
    # (BEAM + CABLE/TRUSS，轴向杆按 EA/L 装配) vs 冻结的 ANSYS 18.0
    # BEAM4 + LINK180 实测值。默认离线运行，不需要 ANSYS 许可。
    add_test(NAME cable_arch_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_arch_3d/compare_cable_arch_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_arch_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_arch_3d/reference)
    set_tests_properties(cable_arch_3d_benchmark PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 非线性成桥态 / 索长迭代端到端 smoke：cable_solver.exe (fem_3d_nl.hpp
    # 共转列式几何非线性)。无线性 ANSYS 基准，检查收敛、几何闭合、索长
    # 闭合以及重跑结果与冻结 *_cable.txt 的一致性。
    add_test(NAME cable_solver_smoke
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_arch_3d/compare_cable_solver_smoke.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --smoke ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_arch_3d/nonlinear_smoke)
    set_tests_properties(cable_solver_smoke PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 3D 分布荷载 / 自重：*ELEMENT_LOAD、*GRAVITY、*SELF_WEIGHT 的一致节点
    # 荷载与固端力扣除。CB3D06/07/08 由同一 ModelSpec3D 导出，参考值冻结
    # 自 ANSYS 18.0 BEAM4 实测 CSV。默认离线运行，不需要 ANSYS 许可。
    add_test(NAME element_loads_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/element_loads_3d/compare_element_loads_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/element_loads_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/element_loads_3d/reference)
    set_tests_properties(element_loads_3d_benchmark PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 同一荷载通道的内核级闭式解 / 叠加恒等式验证：GZ、GY、LZ 两个弯曲
    # 平面、DENS、ADDMASS、*SELF_WEIGHT RHO/MASS、TRUSS 自重、梯形 q1/q2。
    add_test(NAME element_loads_3d_kernel_validation
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/element_loads_3d/verify_kernel_loads.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(element_loads_3d_kernel_validation PROPERTIES
        TIMEOUT 600
        LABELS "benchmark")

    # 非线性成桥态 vs ANSYS NLGEOM：CB3D10 由同一 ModelSpec3D 导出
    # .dat / BEAM4+LINK180+NLGEOM+INISTATE .inp。重跑 cable_solver.exe，
    # 与冻结自真实 ANSYS 18.0 批处理的 325 项黄金值逐项比对；同时校验
    # 重跑结果与冻结 *_cable.txt 的一致性。默认离线，不需要 ANSYS 许可，
    # --live-ansys 可现场重跑归档 APDL 逐格复核 CSV。
    add_test(NAME nonlinear_ansys_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/nonlinear_ansys_3d/compare_nonlinear_ansys_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/nonlinear_ansys_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/nonlinear_ansys_3d/reference)
    set_tests_properties(nonlinear_ansys_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 悬索桥找形 vs ANSYS NLGEOM：CB3D11 地锚式 / CB3D12 自锚式，共用
    # 参数只改 *SUSP_BRIDGE SYSTEM。ANSYS 初始几何 = 设计几何 [ERECTED]+[U]，
    # 索下料长度取 FEM 绝对 L0。重跑 cable_solver.exe 与冻结自真实 ANSYS 18.0
    # 批处理的 574 项黄金值逐项比对，默认离线；--live-ansys 可现场复核 CSV。
    add_test(NAME suspension_ansys_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/suspension_ansys_3d/compare_suspension_ansys_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/suspension_ansys_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/suspension_ansys_3d/reference)
    set_tests_properties(suspension_ansys_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 自锚式自动平衡边跨 + 施工阶段找形：CB3D13。三层验证——
    # [A] 全桥成桥态 / [B] 梁端无抗拔判据 / [C] 5 个施工阶段各一份
    # LINK180 索链 APDL（从扰动几何独立 NLGEOM 求解）。重跑 cable_solver.exe
    # 复现含 [ERECTION_*] 新段的冻结结果，与冻结自真实 ANSYS 18.0 批处理的
    # 478 项黄金值逐项比对；默认离线，--live-ansys 可现场复核 6 个 job。
    add_test(NAME suspension_erection_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/suspension_erection_3d/compare_suspension_erection_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/suspension_erection_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/suspension_erection_3d/reference)
    set_tests_properties(suspension_erection_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 索垂度 / Ernst 等效弹模（CB3D14）+ 带索自重悬索桥空缆态（CB3D15）。
    # ANSYS LINK180 无垂度模型，因此把单根索离散成 80 段 LINK180 +
    # 水平均布节点自重 + NLGEOM，从扰动几何独立收敛出悬链线；FEM 侧由
    # cable_sag_probe.exe 调用 fem_3d.hpp / fem_3d_nl.hpp 的内核闭式。
    # 7 个 job（5 点 L0 扫描 + W=0 直杆 + W/2 低重度），573 项黄金值；
    # 另用 L0 扫描的中心差分验证 Ernst 解析切线 dT/dL0。
    add_test(NAME cable_sag_ansys_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_sag_ansys_3d/compare_cable_sag_ansys_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_sag_ansys_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/cable_sag_ansys_3d/reference)
    set_tests_properties(cable_sag_ansys_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 大转动梁柱效应（CB3D16）：FEM-Eigen fem_3d_nl.hpp 共转列式 CR 梁
    # （随动三元组 + 解析几何刚度）vs ANSYS BEAM188 + NLGEOM。8 个工况：
    # 悬臂端弯矩卷圆 0.25/0.5/0.75/1.0 圈（有限转动 2*pi）、简支梁轴压
    # P/Pcr = 0/0.25/0.50/0.75 的 P-Delta 放大。1776 项冻结黄金值。
    add_test(NAME beam_nl_ansys_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/beam_nl_ansys_3d/compare_beam_nl_ansys_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/beam_nl_ansys_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/beam_nl_ansys_3d/reference)
    set_tests_properties(beam_nl_ansys_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")
else()
    message(STATUS "未找到 Python3，跳过 .dat 逐字节回归 CTest")
endif()

# 原始的一维杆示例。
fem_add_executable(fem_solver One-dimension.cpp)

# 通用框架示例：1D/2D/3D 单元与多自由度组装。
fem_add_executable(example_1d_bar example_1d_bar.cpp)
fem_add_executable(example_2d_truss example_2d_truss.cpp)
fem_add_executable(example_2d_plane example_2d_plane.cpp)
fem_add_executable(example_3d_tet example_3d_tet.cpp)
fem_add_executable(example_3d_truss example_3d_truss.cpp)
fem_add_executable(example_3d_stage example_3d_stage.cpp)
fem_add_executable(example_3d_stage_nl example_3d_stage_nl.cpp)
fem_add_executable(example_3d_stage_envelope example_3d_stage_envelope.cpp)
fem_add_executable(example_3d_beam example_3d_beam.cpp)
fem_add_executable(example_3d_beam_stage example_3d_beam_stage.cpp)
fem_add_executable(example_3d_beam_envelope example_3d_beam_envelope.cpp)
fem_add_executable(example_3d_influence example_3d_influence.cpp)
fem_add_executable(example_3d_influence_surface example_3d_influence_surface.cpp)
fem_add_executable(example_3d_stage_live example_3d_stage_live.cpp)
fem_add_executable(example_3d_beam_live example_3d_beam_live.cpp)
fem_add_executable(example_3d_grillage example_3d_grillage.cpp)
fem_add_executable(example_solvers example_solvers.cpp)
fem_link_cholmod(example_solvers)
fem_add_executable(example_amgcl example_amgcl.cpp)
fem_link_cholmod(example_amgcl)
if(FEM_CHOLMOD_READY)
    fem_add_executable(example_cholmod example_cholmod.cpp)
    fem_link_cholmod(example_cholmod)
else()
    message(STATUS "CHOLMOD 不可用，跳过 example_cholmod")
endif()
fem_add_executable(example_mixed example_mixed.cpp)
fem_link_cholmod(example_mixed)
fem_add_executable(example_influence example_influence.cpp)
fem_link_cholmod(example_influence)
fem_add_executable(example_stage example_stage.cpp)
fem_link_cholmod(example_stage)
fem_add_executable(example_stage_viz example_stage_viz.cpp)
fem_link_cholmod(example_stage_viz)
fem_add_executable(example_stage_envelope example_stage_envelope.cpp)
fem_link_cholmod(example_stage_envelope)
fem_add_executable(example_viz example_viz.cpp)
fem_add_executable(example_modal example_modal.cpp)
fem_add_executable(example_3d_modal example_3d_modal.cpp)
fem_add_executable(example_3d_report example_3d_report.cpp)
fem_add_executable(modal_bench bench/modal_bench.cpp)
fem_add_executable(nl_bench bench/nl_bench.cpp)
# 斜拉桥成桥态 CLI 接缝：BridgeStudio 导出 .dat → 这里 → <stem>_cable.txt
fem_add_executable(cable_solver cable_solver.cpp)
# J 阶段：把索垂度 / Ernst 等效弹模的内核闭式（fem_3d.hpp + fem_3d_nl.hpp）
# 暴露成 CSV，供 tests/benchmark/cable_sag_ansys_3d 与 ANSYS 分段 LINK180
# 索链端到端对标。只读参数文件，不改动 cable_solver 的任何既有路径。
fem_add_executable(cable_sag_probe cable_sag_probe.cpp)
# K 阶段：把 fem_3d_nl.hpp 的共转列式 CR 梁（大转动 + 解析几何刚度）
# 暴露成文件驱动 CSV，供 tests/benchmark/beam_nl_ansys_3d 与 ANSYS
# BEAM188 + NLGEOM 端到端对标。力学完全复用内核，不新增单元公式。
fem_add_executable(beam_nl_probe beam_nl_probe.cpp)

# ---- 2D/3D 分库 ----
# 2D 层（fem_2d）编译为静态库；共享核心 fem.hpp 与 3D 层 fem_3d.hpp 为纯头文件。
# 3D 层用 INTERFACE 目标表达头文件与第三方依赖，不产生编译产物。
add_library(fem_2d STATIC fem_2d.cpp)
add_library(FEM_Solver::fem_2d ALIAS fem_2d)
target_include_directories(fem_2d PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/FEM_Solver>
    $<BUILD_INTERFACE:${FEM_SPECTRA_INCLUDE_DIR}>
)

add_library(fem_3d INTERFACE)
add_library(FEM_Solver::fem_3d ALIAS fem_3d)
target_include_directories(fem_3d INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/FEM_Solver>
    $<BUILD_INTERFACE:${FEM_SPECTRA_INCLUDE_DIR}>
)

if(TARGET Eigen3::Eigen)
    target_link_libraries(fem_2d PUBLIC Eigen3::Eigen)
    target_link_libraries(fem_3d INTERFACE Eigen3::Eigen)
else()
    # 安装场景下要求消费方自行提供 Eigen3；此处仅把 vendored 路径用于构建期示例。
    target_include_directories(fem_2d PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3>)
    target_include_directories(fem_3d INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3>)
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/FEM_Solver
    FILES_MATCHING PATTERN "*.hpp"
)
install(TARGETS fem_2d
    EXPORT FEM_SolverTargets
)
install(EXPORT FEM_SolverTargets
    FILE FEM_SolverTargets.cmake
    NAMESPACE FEM_Solver::
    DESTINATION lib/cmake/FEM_Solver
)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/FEM_SolverConfigVersion.cmake"
    VERSION 1.0.0
    COMPATIBILITY SameMajorVersion
)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/FEM_SolverConfig.cmake"
"include(\\\"\\\${CMAKE_CURRENT_LIST_DIR}/FEM_SolverTargets.cmake\\\")\n"
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FEM_SolverConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FEM_SolverConfigVersion.cmake"
    DESTINATION lib/cmake/FEM_Solver
)
