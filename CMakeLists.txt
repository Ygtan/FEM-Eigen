cmake_minimum_required(VERSION 3.10)
project(FEM_Solver LANGUAGES CXX)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# ---- SuiteSparse CHOLMOD（可选，通过 MSYS2 pacman 安装）----
# 注意：必须同时定位"头文件"与"库文件"。只 find_path 头文件会导致链接期
# `cannot find -lcholmod`，因此这里由头文件路径推导库目录并逐个 find_library，
# 任一库缺失即整体禁用 CHOLMOD，避免产生必然失败的链接目标。
option(FEM_USE_CHOLMOD "启用 SuiteSparse CHOLMOD 后端（未安装时自动禁用）" ON)

set(FEM_CHOLMOD_READY FALSE)
if(FEM_USE_CHOLMOD)
    find_path(CHOLMOD_INCLUDE_DIR cholmod.h
        PATHS "C:/msys64/mingw64/include/suitesparse" "C:/msys64/ucrt64/include/suitesparse")
    if(CHOLMOD_INCLUDE_DIR)
        # <prefix>/include/suitesparse  ->  <prefix>/lib
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

# 添加可执行目标的便捷函数。
function(fem_add_executable target source)
    add_executable(${target} ${source})
    if(TARGET Eigen3::Eigen)
        target_link_libraries(${target} PRIVATE Eigen3::Eigen)
    else()
        target_link_libraries(${target} PRIVATE eigen3_vendored)
    endif()
    if(FEM_AMGCL_INCLUDE_DIR)
        target_include_directories(${target} PRIVATE "${FEM_AMGCL_INCLUDE_DIR}")
        target_compile_definitions(${target} PRIVATE AMGCL_NO_BOOST)
    endif()
endfunction()

# 仅为真正使用 CHOLMOD 的目标链接 SuiteSparse（避免其它示例产生 DLL 依赖）。
# 使用 find_library 得到的绝对路径，而不是裸 -l 名称。
function(fem_link_cholmod target)
    if(FEM_CHOLMOD_READY)
        target_include_directories(${target} PRIVATE "${CHOLMOD_INCLUDE_DIR}")
        target_link_libraries(${target} PRIVATE ${FEM_CHOLMOD_LIBS})
    endif()
endfunction()

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
