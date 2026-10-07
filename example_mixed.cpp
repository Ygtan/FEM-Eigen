#include "fem_io.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace fem;

static SolverType parseSolver(const std::string& name) {
    std::string s = name;
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (s == "LDLT" || s == "SIMPLICIAL_LDLT") return SolverType::SimplicialLDLT;
    if (s == "LLT" || s == "SIMPLICIAL_LLT") return SolverType::SimplicialLLT;
    if (s == "LU" || s == "SPARSE_LU") return SolverType::SparseLU;
    if (s == "QR" || s == "SPARSE_QR") return SolverType::SparseQR;
    if (s == "CG" || s == "CONJUGATE_GRADIENT") return SolverType::ConjugateGradient;
    if (s == "BICGSTAB" || s == "BICG") return SolverType::BiCGSTAB;
    if (s == "AMGCL") return SolverType::AMGCL;
    if (s == "CHOLMOD") return SolverType::CHOLMOD;
    throw std::runtime_error("未知求解器: " + name);
}

static const char* solverName(SolverType t) {
    switch (t) {
        case SolverType::SimplicialLDLT: return "SimplicialLDLT";
        case SolverType::SimplicialLLT: return "SimplicialLLT";
        case SolverType::SparseLU: return "SparseLU";
        case SolverType::SparseQR: return "SparseQR";
        case SolverType::ConjugateGradient: return "ConjugateGradient";
        case SolverType::BiCGSTAB: return "BiCGSTAB";
        case SolverType::AMGCL: return "AMGCL";
        case SolverType::CHOLMOD: return "CHOLMOD";
    }
    return "Unknown";
}

static ModelResult timedSolve(const Model2D& model, SolverType type, double& ms) {
    const auto t0 = std::chrono::steady_clock::now();
    ModelResult r = solveModel2D(model, type);
    const auto t1 = std::chrono::steady_clock::now();
    ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return r;
}

static int runMixed(int argc, char** argv) {
    setupConsoleUTF8();

    std::string path = "mixed_truss_beam.dat";
    SolverType solverType = SolverType::SimplicialLDLT;
    bool compare = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--solver") {
            if (i + 1 >= argc) throw std::runtime_error("--solver 需要一个参数");
            solverType = parseSolver(argv[++i]);
        } else if (a == "--compare") {
            compare = true;
        } else {
            path = a;
        }
    }

    const Model2D model = readModel2D(path);
    const std::string title = model.title.empty() ? path : model.title;

    if (compare) {
        const SolverType solvers[] = {
            SolverType::SimplicialLDLT, SolverType::CHOLMOD, SolverType::AMGCL};
        std::cout << "求解器对比（含组装+求解）: " << title << "\n";
        std::cout << "节点数 = " << model.nodes.size()
                  << "，单元数 = " << model.elements.size() << "\n";
        double refMs = -1.0;
        for (const SolverType t : solvers) {
            try {
                double ms = 0.0;
                const ModelResult r = timedSolve(model, t, ms);
                double maxDisp = 0.0;
                for (int i = 0; i < static_cast<int>(r.U.size()); ++i)
                    maxDisp = std::max(maxDisp, std::abs(r.U(i)));
                std::cout << "  " << std::left << std::setw(16) << solverName(t)
                          << "  " << std::right << std::setw(12) << ms << " ms";
                if (t == SolverType::SimplicialLDLT) refMs = ms;
                else if (refMs > 0.0)
                    std::cout << "  (" << (ms / refMs) << "x LDLT)";
                std::cout << "  最大位移=" << maxDisp << "\n";
            } catch (const std::exception& e) {
                std::cout << "  " << std::left << std::setw(16) << solverName(t)
                          << "  不可用: " << e.what() << "\n";
            }
        }
        return 0;
    }

    double solveMs = 0.0;
    const ModelResult r = timedSolve(model, solverType, solveMs);
    const std::vector<ElementForce> forces = elementForces2D(model, r.U);
    std::vector<ElementStress> stresses;
    bool hasStress = true;
    try {
        stresses = elementStress2D(model, forces);
    } catch (const std::runtime_error& e) {
        hasStress = false;
        std::cerr << "提示: " << e.what() << "（跳过应力输出）\n";
    }

    std::string outPath = path;
    const std::size_t p = outPath.rfind(".dat");
    if (p != std::string::npos && p + 4 == outPath.size())
        outPath = outPath.substr(0, p);
    outPath += ".result";
    writeResult(outPath, model, r, forces);

    std::cout << std::setprecision(8);
    std::cout << "模型: " << title << "\n";
    std::cout << "节点数 = " << r.nodeIds.size()
              << "，单元数 = " << model.elements.size() << "\n";
    std::cout << "求解器 = " << solverName(solverType)
              << "，耗时 = " << solveMs << " ms（含组装）\n";
    if (model.analysisType == "NONLINEAR") {
        std::cout << "非线性 Newton 迭代 = " << r.nonlinearIterations
                  << "，自由残差 = " << r.nonlinearResidual << "\n";
    }
    std::cout << "结果已写入: " << outPath << "\n\n";

    std::cout << "节点位移:\n";
    std::cout << "节点    u(x)          v(y)          θ(转角)\n";
    for (std::size_t i = 0; i < r.nodeIds.size(); ++i) {
        const int d = 3 * static_cast<int>(i);
        std::cout << std::setw(4) << r.nodeIds[i]
                  << "  " << std::setw(14) << r.U(d + 0)
                  << "  " << std::setw(14) << r.U(d + 1)
                  << "  " << std::setw(14) << r.U(d + 2) << "\n";
    }

    std::cout << "\n单元内力:\n";
    std::cout << "单元  类型   N1(轴力i)     N2(轴力j)    V1(剪力)     M1(弯矩)     V2(剪力)     M2(弯矩)\n";
    for (const auto& f : forces) {
        if (isRigid(f.type)) continue;
        std::cout << std::setw(4) << f.id << "  " << f.type
                  << "  " << std::setw(13) << f.N1
                  << "  " << std::setw(13) << f.N2;
        if (f.type == "BEAM") {
            std::cout << "  " << std::setw(13) << f.V1
                      << "  " << std::setw(13) << f.M1
                      << "  " << std::setw(13) << f.V2
                      << "  " << std::setw(13) << f.M2;
        }
        std::cout << "\n";
    }

    if (hasStress) {
        std::cout << "\n单元应力:\n";
        for (const auto& s : stresses) {
            if (isRigid(s.type)) continue;
            std::cout << std::setw(4) << s.id << "  " << s.type;
            if (isAxial(s.type)) {
                std::cout << "  σ1=" << std::setw(13) << s.s1
                          << "  σ2=" << std::setw(13) << s.s2;
            } else {
                std::cout << "  上翼缘(1)=" << std::setw(13) << s.stop1
                          << " 下翼缘(1)=" << std::setw(13) << s.sbot1
                          << " 上翼缘(2)=" << std::setw(13) << s.stop2
                          << " 下翼缘(2)=" << std::setw(13) << s.sbot2;
            }
            std::cout << "\n";
        }
    }

    std::cout << "\n支座反力（非零项）:\n";
    std::cout << "节点   dof          反力\n";
    for (std::size_t i = 0; i < r.nodeIds.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            const double v = r.R(3 * static_cast<int>(i) + d);
            if (std::abs(v) > 1e-9)
                std::cout << std::setw(4) << r.nodeIds[i]
                          << "   " << (d + 1)
                          << "    " << std::setw(14) << v << "\n";
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    try {
        return runMixed(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
