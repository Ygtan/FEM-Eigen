#include "fem.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

int main(int argc, char** argv) {
    using namespace fem;
    setupConsoleUTF8();

    // 二维结构化 CST 网格，每边 N 格；自由度约 2*(N+1)^2。
    int N = 100;
    if (argc > 1) N = std::atoi(argv[1]);
    if (N < 2) N = 2;

    const double E = 1000.0, nu = 0.3, t = 1.0;
    const int dofsPerNode = 2;
    const int nPer = N + 1;
    const int nNodes = nPer * nPer;
    const int nDofs = nNodes * dofsPerNode;

    auto id = [&](int i, int j) { return i * nPer + j; };

    std::vector<Eigen::Vector2d> x(nNodes);
    for (int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            x[id(i, j)] = {double(i) / N, double(j) / N};

    std::vector<std::array<int, 3>> tris;
    tris.reserve(2 * N * N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            const int a = id(i, j), b = id(i + 1, j);
            const int c = id(i + 1, j + 1), d = id(i, j + 1);
            tris.push_back({a, b, c});
            tris.push_back({a, c, d});
        }

    std::vector<DenseMat> kes;
    std::vector<std::vector<int>> gdofs;
    kes.reserve(tris.size());
    gdofs.reserve(tris.size());

    for (const auto& tri : tris) {
        TriangleCST e{E, nu, t, x[tri[0]], x[tri[1]], x[tri[2]]};
        kes.push_back(e.stiffness());
        gdofs.push_back(nodeDofs({tri[0], tri[1], tri[2]}, dofsPerNode));
    }

    SparseMat K = assemble(kes, gdofs, nDofs);
    DenseVec F = DenseVec::Zero(nDofs);

    // 右边缘施加向下剪力。
    for (int j = 0; j <= N; ++j)
        F(2 * id(N, j) + 1) -= 1.0;

    // 左边缘固定 u = v = 0。
    std::vector<int> bc;
    std::vector<double> bcv;
    for (int j = 0; j <= N; ++j) {
        bc.push_back(2 * id(0, j) + 0); bcv.push_back(0.0);
        bc.push_back(2 * id(0, j) + 1); bcv.push_back(0.0);
    }
    applyDirichlet(K, F, bc, bcv);

    std::cout << std::setprecision(6);
    std::cout << "自由度 n = " << nDofs
              << "，CST 单元数 = " << tris.size() << "\n";

    auto run = [&](SolverType type) {
        const auto t0 = std::chrono::steady_clock::now();
        DenseVec U = solve(K, F, type);
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        return std::make_pair(U, ms);
    };

    const auto [Uref, tLDLT] = run(SolverType::SimplicialLDLT);
    std::cout << "SimplicialLDLT : " << tLDLT << " ms\n";

    const auto [Ucho, tCHO] = run(SolverType::CHOLMOD);
    std::cout << "CHOLMOD        : " << tCHO << " ms\n";

    const auto [Uamg, tAMG] = run(SolverType::AMGCL);
    std::cout << "AMGCL(BiCGStab): " << tAMG << " ms\n";

    const double errCho = (Uref - Ucho).norm() / (Uref.norm() + 1e-300);
    const double errAmg = (Uref - Uamg).norm() / (Uref.norm() + 1e-300);
    std::cout << "CHOLMOD 相对误差 = " << errCho << "\n";
    std::cout << "AMGCL  相对误差 = " << errAmg << "\n";
    std::cout << "右上角位移: u=" << Uref(2 * id(N, N))
              << ", v=" << Uref(2 * id(N, N) + 1) << "\n";
    const bool ok = errCho < 1e-6 && errAmg < 1e-6;
    std::cout << (ok ? "[PASS]" : "[FAIL]") << "\n";
    return ok ? 0 : 1;
}
