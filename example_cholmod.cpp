#include "fem.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    using namespace fem;
    setupConsoleUTF8();

    const double E0 = 200e9, A = 0.01;
    const int dim = 1;
    const int nNodes = 51;
    const int dofsPerNode = dim;
    const int nDofs = nNodes;

    // 构造 1D 杆的 K 与 F（给定弹性模量 E）。
    auto makeSystem = [&](double E, SparseMat& K, DenseVec& F) {
        std::vector<Eigen::Vector3d> x(nNodes);
        for (int i = 0; i < nNodes; ++i) x[i] = {double(i), 0.0, 0.0};

        std::vector<DenseMat> kes;
        std::vector<std::vector<int>> gdofs;
        for (int i = 0; i < nNodes - 1; ++i) {
            TrussElement e{E, A, dim, x[i], x[i + 1]};
            kes.push_back(e.stiffness());
            gdofs.push_back(nodeDofs({i, i + 1}, dofsPerNode));
        }

        K = assemble(kes, gdofs, nDofs);
        F = DenseVec::Zero(nDofs);
        F(nNodes / 2) = 1000.0;
        applyDirichlet(K, F, {0}, {0.0});
    };

    SparseMat K, K2;
    DenseVec F, F2;
    makeSystem(E0, K, F);
    makeSystem(2.0 * E0, K2, F2); // 结构不变，仅刚度值翻倍

    // 多右端项：三个不同载荷，每列一个。
    DenseMat B = DenseMat::Zero(nDofs, 3);
    B(nNodes / 4, 0) = 1.0;
    for (int i = 1; i < nNodes; ++i) B(i, 1) = 1.0;
    B(nNodes / 2, 2) = 1.0;

    CholmodSolver solver;
    solver.analyzePattern(K); // 符号分解：只做一次

    solver.factorize(K);
    const DenseVec u1 = solver.solve(F); // 单右端项

    solver.factorize(K2); // 复用符号分解
    const DenseVec u2 = solver.solve(F2);

    solver.factorize(K);
    const DenseMat U = solver.solve(B); // 多右端项

    // 参考解（Eigen LDLT）。
    const DenseVec r1 = solve(K, F, SolverType::SimplicialLDLT);
    const DenseVec r2 = solve(K2, F2, SolverType::SimplicialLDLT);
    DenseMat R(nDofs, 3);
    for (int j = 0; j < 3; ++j)
        R.col(j) = solve(K, B.col(j), SolverType::SimplicialLDLT);

    const double e1 = (u1 - r1).norm();
    const double e2 = (u2 - r2).norm();
    const double e3 = (U - R).norm();

    std::cout << std::setprecision(10);
    std::cout << "单右端项误差       = " << e1 << "\n";
    std::cout << "复用符号分解误差   = " << e2 << "\n";
    std::cout << "多右端项误差       = " << e3 << "\n";
    std::cout << "u2(中点)/u1(中点)  = " << u2(nNodes / 2) / u1(nNodes / 2)
              << "（应为 0.5）\n";

    const bool ok = e1 < 1e-12 && e2 < 1e-12 && e3 < 1e-12;
    std::cout << (ok ? "[PASS]" : "[FAIL]") << "\n";
    return ok ? 0 : 1;
}
