#include "fem.hpp"

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

static const char* solverName(fem::SolverType t) {
    using fem::SolverType;
    switch (t) {
        case SolverType::SimplicialLDLT: return "SimplicialLDLT";
        case SolverType::SimplicialLLT: return "SimplicialLLT";
        case SolverType::SparseLU: return "SparseLU";
        case SolverType::SparseQR: return "SparseQR";
        case SolverType::ConjugateGradient: return "ConjugateGradient";
        case SolverType::BiCGSTAB: return "BiCGSTAB";
    }
    return "Unknown";
}

int main() {
    using namespace fem;
    setupConsoleUTF8();

    // 构造一个已知解析解的 SPD 系统：1D 杆（3 单元、4 节点）。
    const double E = 200e9, A = 0.01;
    const int dim = 1;
    const int nNodes = 4;
    const int dofsPerNode = dim;
    const int nDofs = nNodes * dofsPerNode;

    std::vector<Eigen::Vector3d> x = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
    std::vector<std::pair<int, int>> conn = {{0, 1}, {1, 2}, {2, 3}};

    std::vector<DenseMat> kes;
    std::vector<std::vector<int>> gdofs;
    for (const auto& [a, b] : conn) {
        TrussElement e{E, A, dim, x[a], x[b]};
        kes.push_back(e.stiffness());
        gdofs.push_back(nodeDofs({a, b}, dofsPerNode));
    }

    SparseMat K = assemble(kes, gdofs, nDofs);
    DenseVec F = DenseVec::Zero(nDofs);
    F(2) = 1000.0; // 节点 2 施加 1000 N
    applyDirichlet(K, F, {0}, {0.0});

    const SolverType types[] = {
        SolverType::SimplicialLDLT, SolverType::SimplicialLLT,
        SolverType::SparseLU, SolverType::SparseQR,
        SolverType::ConjugateGradient, SolverType::BiCGSTAB};

    std::cout << "各求解后端结果（节点2位移，解析解 u2 = 1e-6）:\n";
    bool allPass = true;
    for (const SolverType t : types) {
        const DenseVec U = solve(K, F, t);
        const double err = std::abs(U(2) - 1e-6);
        const bool ok = err < 1e-9;
        allPass = allPass && ok;
        std::cout << "  " << solverName(t) << ": u2 = " << U(2)
                  << ", err = " << err << (ok ? "  [PASS]" : "  [FAIL]") << "\n";
    }
    std::cout << (allPass ? "全部通过" : "存在失败") << "\n";
    return allPass ? 0 : 1;
}
