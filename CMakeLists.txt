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
    target_link_libraries(${target} PRIVATE fem_core)
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

# 解析解断言测试：用材料力学闭式解校验梁/桁架单元与求解器。
fem_add_executable(test_analytical test_analytical.cpp)
add_test(NAME analytical_exact_solutions COMMAND test_analytical)

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

# ---- 核心静态库（fem_io.cpp + 头文件）----
# 将 fem_io.hpp 中的实现提取到 fem_io.cpp，减少编译时间并实现二进制隔离。
add_library(fem_core STATIC fem_io.cpp)
add_library(FEM_Solver::fem_core ALIAS fem_core)
target_include_directories(fem_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/FEM_Solver>
    $<BUILD_INTERFACE:${FEM_SPECTRA_INCLUDE_DIR}>
)
if(TARGET Eigen3::Eigen)
    target_link_libraries(fem_core PUBLIC Eigen3::Eigen)
else()
    # 安装场景下要求消费方自行提供 Eigen3；此处仅把 vendored 路径用于构建期示例。
    target_include_directories(fem_core PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3>)
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/FEM_Solver
    FILES_MATCHING PATTERN "*.hpp"
)
install(TARGETS fem_core
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
