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

# 第 39 项：3D 约束方程（索鞍顶推口径 C）内核验收。
# 独立构造三节点模型，覆盖线性/非线性两条求解路径的 CE 顶推命中、
# 缆鞍固结、无效 dof 拒绝、重复系数合并、空 CE 确定性。离线、不跑 ANSYS。
fem_add_executable(ce_3d_benchmark ce_3d_benchmark.cpp)
add_test(NAME ce_3d_constraint_equations COMMAND ce_3d_benchmark)
set_tests_properties(ce_3d_constraint_equations PROPERTIES
    TIMEOUT 300
    LABELS "benchmark")

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

    # 后屈曲 / snap-through（CB3D17）：FEM-Eigen stability_path_3d.hpp 球面
    # 弧长法（lambda 为未知量 + Crisfield 约束）vs ANSYS LINK180 + ARCLEN。
    # von Mises 两杆浅桁架完整 S 形荷载-位移路径：极限点、下降段、反向
    # 荷载段。1686 项冻结黄金值 + 解析 P(v) 独立锚点。
    add_test(NAME snap_through_ansys_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/snap_through_ansys_3d/compare_snap_through_ansys_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/snap_through_ansys_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/snap_through_ansys_3d/reference)
    set_tests_properties(snap_through_ansys_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # B3 主链路 [STABILITY] opt-in（第十七次）：cable_solver.exe 的
    # STABILITY=1 模式 vs CB3D17 冻结路径 + von Mises 闭式解析锚点。
    # 默认全关零污染、3 类 ERROR 负例，3311 项冻结黄金值。
    add_test(NAME stable_path_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stable_path_3d/compare_stable_path_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stable_path_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stable_path_3d/reference)
    set_tests_properties(stable_path_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # CR 梁刚接框架后屈曲（CB3D19 / 第十八次）：frame_probe.exe 的球面
    # 弧长法 + 共转 CR 梁 vs ANSYS BEAM188 + NLGEOM + ARCLEN。Williams
    # toggle 型框架：第一极限点 / 下降段 / 反向荷载 / 再硬化完整路径，
    # 12513 项冻结黄金值 + 独立 2D 线性刚架刚度锚点 k0。
    add_test(NAME frame_postbuckling_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/frame_postbuckling_3d/compare_frame_postbuckling_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/frame_postbuckling_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/frame_postbuckling_3d/reference)
    set_tests_properties(frame_postbuckling_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 成桥态起点全过程稳定（CB3D20 / 第十九次 ①）：STABASIS=COMPLETION
    # 把弧长路径起点换成成桥态收敛位移 u0 + 恒定荷载 F0，
    # R = F0 + lambda*Finc - Fint(u)。V 形双索单自由度体系的不读求解器
    # 闭式 lambda(v) 逐点对拍 MODE=STATIC / TENSION 两条成桥态闭环，
    # 2963 项冻结黄金值 + 默认零污染 + 4 类 ERROR 负例。不跑 ANSYS。
    add_test(NAME completion_stability_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_stability_3d/compare_completion_stability_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_stability_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_stability_3d/reference)
    set_tests_properties(completion_stability_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 真实斜拉桥成桥态两段对标（CB3D21 / 第十九次 ②）：50 节点 / 63 单元 /
    # 16 索，FEM STABASIS=COMPLETION 与 ANSYS 180「step1 恒载成桥态 +
    # step2 增量荷载」逐子步对拍。荷载控制主 job 用 lambda=(TIME-1)*LAMMAX，
    # ARCLEN job 只作终点交叉验证；2392 项冻结黄金值 + 同源 L0/INISTATE
    # 独立锚点 + 默认零污染。CTest 本身完全离线、不需要 ANSYS 许可。
    add_test(NAME completion_bridge_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_bridge_3d/compare_completion_bridge_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_bridge_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_bridge_3d/reference)
    set_tests_properties(completion_bridge_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第十九次 ③：MODE=COMPLETION 成桥态基准同源修复（CB3D10 / CB3D13）。② 登记的
    # 缺陷是求解器按值收 Model3D、内部改写 L0 / 制造线形，主链路却拿外部
    # 原模型算 F0 与 Fint(u0)（CB3D10 RES0=3.15 MN、首步 λ<0）。修复由
    # basisModel 回传同源模型；本项冻结两条独立成桥态闭环
    # （COMPLETION 制造线形基准 / STATIC 设计线形+冻结 L0 基准）的逐点
    # λ / ΔCTRL 对拍、起点闭合、默认零污染，并核对 CB3D13 悬索桥
    # 5 个 [ERECTION*] 施工阶段段与 [STABILITY] 共存且逐位一致。
    # 完全离线、不跑 ANSYS。
    add_test(NAME completion_basis_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_basis_3d/compare_completion_basis_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_basis_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/completion_basis_3d/reference)
    set_tests_properties(completion_basis_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第十九次 ③a：真实斜拉桥越过极值点 / 下降段（CB3D22）。
    # BS 侧 cases/CB3D22_limit_point_vs_ansys/ 的同源冻结证据：在 CB3D21
    # 同源桥上只缩塔刚度、增量荷载改为塔顶 2009 节点 +X，靠 P-Δ 软化
    # 产生 lambda_max 与下降段。本项完全离线，不重跑 ANSYS；判据覆盖
    # 冻结证据哈希、黄金值自洽、16 索 L0 / INISTATE / ANSYS 弦长同源锚点、
    # FEM 极值点物理判据、两段 ANSYS 求解证据、按 ΔUX 插值 λ 的
    # FEM↔ANSYS 逐点对拍，以及当前二进制重放两个 .dat。
    add_test(NAME bridge_limit_point_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bridge_limit_point_3d/compare_bridge_limit_point_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bridge_limit_point_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bridge_limit_point_3d/reference)
    set_tests_properties(bridge_limit_point_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第二十次 ③b：施工阶段中间态起点的全过程稳定（CB3D23）。
    # STABASIS=STAGE STAGEID=2 把弧长路径起点换成 solveModel3DStagesNonlinear()
    # 的阶段 2 共转列式全量平衡态 u_stage 与同源 NLStageResult3D::model，
    # F0 为该阶段累计恒定荷载、Finc 只来自 *STABILITY_LOAD；
    # [STAGES] + STAGE= / BASIS=stage_state 与另两种基准分列。
    # BS 侧 cases/CB3D23_stage_state_vs_ansys/ 的 15 份同源冻结证据：
    # 4 节点 3 杆 TRUSS 浅拱、两个 ANSYS 18.0 job（st_load 荷载控制 /
    # st_arc ARCLEN），step1 TIME=1 求同一阶段平衡态、step2 λ=TIME−1，
    # 按 |ΔUY| 插值 λ 逐点对拍。本项完全离线，不重跑 ANSYS；判据覆盖
    # 冻结证据哈希、黄金值自洽、起点闭合 / 非零初态 / 阶段递进、
    # APDL N/E/F/D 同源锚点、负例 ERROR=stage_not_found，
    # 以及当前二进制重放两个 .dat 的逐字节一致性（mingw64 同工具链）。
    add_test(NAME stage_stability_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stage_stability_3d/compare_stage_state_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stage_stability_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stage_stability_3d/reference)
    set_tests_properties(stage_stability_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 27 项：Wempner–Riks 正交轨迹（PATHMETHOD=RIKS，仅 STAGE 起点
    # opt-in）。von Mises 两杆闭式 λ(v) 逐点核验 + SPHERE/RIKS 同问题
    # 双方法交叉 + CB3D23 跨基准终点锚点 + 未知方法负例 + 当前二进制
    # 逐字节重放。完全离线、不跑 ANSYS；球面默认路径与第 20/26 项
    # 冻结输出保持逐位不变。
    add_test(NAME riks_stage_stability_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/riks_stage_stability_3d/compare_riks_stage_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/riks_stage_stability_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/riks_stage_stability_3d/reference)
    set_tests_properties(riks_stage_stability_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 28 项：第二十一次力学①位移控制 / 线搜索（opt-in）。
    # von Mises 两杆闭式 P(v) 逐点锚点 + 线搜索 / 无线搜索同路径交叉 +
    # 细步 λ_max→λ_cr 收敛 + 奇异控制 DOF 负例 + 当前二进制逐字节重放。
    # 默认全关、snap_probe 既有 .snap 与第 19 项冻结输出逐位不变；离线不跑 ANSYS。
    add_test(NAME displacement_control_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/displacement_control_3d/compare_displacement_control_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/displacement_control_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/displacement_control_3d/reference)
    set_tests_properties(displacement_control_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 29 项：第二十一次力学②初始缺陷形态 / 幅值 / 敏感性（opt-in）。
    # 线性屈曲模态生成无应力缺陷场；幅值→λ_max 单调 + 反对称模态对称破缺
    # + 独立闭式 P(v) 锚点 + IMPMODE 越界负例 + IMPERF=0 零污染；离线。
    add_test(NAME initial_imperfection_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/initial_imperfection_3d/compare_imperfection_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/initial_imperfection_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/initial_imperfection_3d/reference)
    set_tests_properties(initial_imperfection_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 30 项：第二十一次力学③材料非线性本构接入 3D CR 梁（opt-in）。
    # 单轴应力-应变平台 + 矩形截面塑性铰逐阶段屈服 + 混凝土受拉开裂 +
    # 弹性段纤维 vs Kl12 逐位对齐 + 独立 Python 纤维积分锚点 + 负例；离线。
    # 默认 MATKIND 缺席 → ELASTIC，全部既有 CTest 逐位不变。
    add_test(NAME material_nonlinear_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/material_nonlinear_3d/compare_material_nl_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/material_nonlinear_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/material_nonlinear_3d/reference)
    set_tests_properties(material_nonlinear_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 31 项：第二十一次力学④分岔点检测与分支切换（opt-in）。
    # 高拱 von Mises 桁架：反对称分岔 L^3=L0*h^2 闭式锚点 + pivotRatio/
    # 一致切线特征值双确认 + 模态扰动切换分支（承载力降 / 对称破缺）+
    # 对称模态 / 零幅值 / 越界模态负例 + BIFURC=0 零污染；FEM 自洽，离线。
    add_test(NAME bifurcation_switch_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bifurcation_switch_3d/compare_bifurcation_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bifurcation_switch_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/bifurcation_switch_3d/reference)
    set_tests_properties(bifurcation_switch_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 32 项：第二十一次力学⑤动力 snap-through（opt-in，独立新模块）。
    # 一致/集中质量矩阵 + Newmark-β 隐式积分 + 非线性 Newton 平衡迭代；
    # 独立 SDOF RK4 时程锚点 + 静力 Pcr 以下动力跃越 + 阶跃/爬升/半正弦
    # 荷载 + 阻尼 + 质量模式 + 零质量/零步长负例 + DYN=0 零污染；离线。
    add_test(NAME dynamic_snapthrough_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/dynamic_snapthrough_3d/compare_dynamic_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/dynamic_snapthrough_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/dynamic_snapthrough_3d/reference)
    set_tests_properties(dynamic_snapthrough_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 33 项：第二十三次 B4 · σl2 附录 G 反摩擦（opt-in，默认全关）。
    # 张拉端 x=0 代表值；lf≤l / lf>l 两支闭式 + 无摩擦退化 + 缺参负例 +
    # 默认逐位零污染 + 2D/3D 同参同值；离线，不跑 ANSYS。
    add_test(NAME prestress_loss_anti_2d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/prestress_loss_anti_2d/compare_prestress_loss_anti_2d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/prestress_loss_anti_2d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/prestress_loss_anti_2d/reference)
    set_tests_properties(prestress_loss_anti_2d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 34 项：第二十三次 B4 · σl6 收缩徐变附录 C 显式表值（opt-in，默认全关）。
    # 表值由用户显式录入，程序只按解析发展公式组装 εcs/φ，不查表不拟合；
    # C_SKIP 显式不参与 / C2D 完整表值 / 四负例 / 默认零污染 / 2D/3D 同值；
    # 离线，不跑 ANSYS。
    add_test(NAME prestress_loss_creep_2d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/prestress_loss_creep_2d/compare_prestress_loss_creep_2d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/prestress_loss_creep_2d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/prestress_loss_creep_2d/reference)
    set_tests_properties(prestress_loss_creep_2d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 35 项：第二十四次 B4 · 结构整体时变徐变/收缩 JTG 3362 规范算例书锚点。
    # 与 σl6 表值链刻意隔离：本链 φ0=φRH·βF·βT0、βH 用实际 RH 连续值、
    # εcs0=εS·βRH（code_loads.hpp 连续拟合）；独立 Python 闭式 oracle 逐键核对
    # φ/εcs/βs/E_eff，含 σl6 口径隔离硬断言 + 四负例 + 当前二进制逐字节重放。
    # 离线，不跑 ANSYS。
    add_test(NAME creep_shrinkage_time_2d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/creep_shrinkage_time_2d/compare_creep_shrinkage_time_2d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/creep_shrinkage_time_2d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/creep_shrinkage_time_2d/reference)
    set_tests_properties(creep_shrinkage_time_2d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 36 项：2D 分项预应力损失阶段维度（B4 第二十五次，opt-in）。
    # *STAGE_PRESTRESS 支持 LOSSREF，σl2 附录 G / σl6 附录 C 选项经同一张
    # *PRESTRESS_LOSS 表映射；阶段段独立写 <stem>.stage.loss.txt
    # （BASIS=2d_stage），与全局 BASIS=2d_global 段及既有 .loss.txt 隔离。
    # 独立 Python 闭式 oracle + 冻结证据 + 硬/软负例 + 默认零污染 +
    # 2D/3D 同参同值 + FE/BS 源码契约守卫。离线，不跑 ANSYS。
    add_test(NAME stage_loss_2d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stage_loss_2d/compare_stage_loss_2d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stage_loss_2d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/stage_loss_2d/reference)
    set_tests_properties(stage_loss_2d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 37 项：成桥态无应力曲率 benchmark（B5 第二十六次，opt-in）。
    # *EN 单元行 K0Y=/K0Z= 录入无应力曲率，固端力并入 elementEquivalentLoad3D，
    # 端内力自动为 M = EI(κ − κ0)；κ0=0 时逐位不变。opt-in 独立输出
    # <stem>.curvature.txt（example_3d_report --curvature-detail）。
    # 独立 Python 闭式 oracle + 冻结证据 + 硬负例 + 自平衡/刚度不变/J 不变
    # 强断言 + 默认零污染 + FE/BS 源码契约守卫。离线，不跑 ANSYS。
    add_test(NAME unstressed_curvature_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/unstressed_curvature_3d/compare_unstressed_curvature_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/unstressed_curvature_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/unstressed_curvature_3d/reference)
    set_tests_properties(unstressed_curvature_3d_benchmark PROPERTIES
        TIMEOUT 900
        LABELS "benchmark")

    # 第 38 项：双主缆悬索桥 benchmark（双索面 3D opt-in）。
    # 冻结双/单索面两份 .dat；当前 cable_solver.exe 重放断言两平面独立
    # 找形、H 相等、全部索登记（2×(主缆段+吊索)）；单索面输出不得出现
    # susp_plane_count（默认零污染）。源码契约锁住「附加平面索不并入
    # cfg.mainCable 单链」——否则施工阶段 findErectionState 会把两条链
    # 当成一条并在接缝处误报不连续。离线，不跑 ANSYS。
    add_test(NAME dual_plane_suspension_3d_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/dual_plane_suspension_3d/compare_dual_plane_suspension_3d.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}
                     --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/dual_plane_suspension_3d/cases
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmark/dual_plane_suspension_3d/reference)
    set_tests_properties(dual_plane_suspension_3d_benchmark PROPERTIES
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
# L 阶段：球面弧长法平衡路径跟踪探针（stability_path_3d.hpp），越过荷载
# 极限点走完整 snap-through 下降段；与 ANSYS LINK180 + ARCLEN 对标。
# 独立新文件，不改 fem_3d_nl.hpp 主求解器与 NonlinearOptions3D 默认值。
fem_add_executable(snap_probe snap_probe.cpp)
# 第二十一次力学③：材料非线性纤维梁 / 截面本构探针（section.hpp +
# fem_3d_nl.hpp 的 crFiberLocalForce/Tangent），文件驱动 CSV，供
# tests/benchmark/material_nonlinear_3d 做闭式锚点校验。
fem_add_executable(mat_nl_probe mat_nl_probe.cpp)
# 第二十一次力学⑤：Newmark-β 隐式动力时程 / 动力 snap-through 探针
# （dynamic_newmark_3d.hpp），文件驱动 CSV，供 tests/benchmark/
# dynamic_snapthrough_3d 做独立 SDOF RK4 锚点校验。独立新模块，
# 不改 snap_probe / stability_path_3d / fem_3d_nl 的任何既有路径。
fem_add_executable(dyn_probe dyn_probe.cpp)
# N 阶段：把同一个球面弧长法内核用于**共转 CR 梁刚接框架**后屈曲
# （Williams toggle 型，两根多单元 BEAM，端部铰支 / 顶点刚接），
# 与 ANSYS BEAM188 + ARCLEN 端到端对标。力学仍完全复用内核，
# 不改 fem_3d_nl.hpp / stability_path_3d.hpp。
fem_add_executable(frame_probe frame_probe.cpp)
# 第二十三次 B4：3D 分项预应力损失只读探针（fem_3d.hpp + prestress_loss.hpp），
# 把同一组 *PRESTRESS_LOSS 参数的 3D 内核结果暴露成 CSV，供第 33 项
# CTest 与 2D 明细段做同参同值对拍；不新增力学、不改任何既有路径。
fem_add_executable(loss_probe_3d loss_probe_3d.cpp)
# 第二十四次 B4：2D 结构时变徐变/收缩只读探针（fem_2d.hpp + code_loads.hpp），
# 把 JTG 3362 附录 C 连续拟合的 φ0/βH/εcs0 与逐龄期 φ(t)/εcs(t)/E_eff
# 暴露成 CSV，供第 35 项 CTest 独立闭式 oracle 逐键核对；不新增力学、
# 不改任何既有路径，与 σl6 表值链保持口径隔离。
fem_add_executable(creep_probe_2d creep_probe_2d.cpp)

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
