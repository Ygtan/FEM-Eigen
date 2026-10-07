#include "fem.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    using namespace fem;
    setupConsoleUTF8();

    const double E = 100.0, nu = 0.25, t = 1.0;
    const int dofsPerNode = 2;
    const int nNodes = 4;
    const int nDofs = nNodes * dofsPerNode;

    // 单位正方形，剖分为两个 CST 三角形。
    std::vector<Eigen::Vector2d> x = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    std::vector<std::array<int, 3>> tris = {{0, 1, 3}, {1, 2, 3}};

    // 目标应力 sigma_x = 1（单向拉伸）。
    Eigen::Vector3d sigma(1.0, 0.0, 0.0);

    std::vector<DenseMat> kes;
    std::vector<std::vector<int>> gdofs;
    DenseVec F = DenseVec::Zero(nDofs);

    for (const auto& tri : tris) {
        TriangleCST e{E, nu, t, x[tri[0]], x[tri[1]], x[tri[2]]};
        const DenseMat B = e.B();
        const double A = e.area();
        kes.push_back(e.stiffness());

        const std::vector<int> g = nodeDofs({tri[0], tri[1], tri[2]}, dofsPerNode);
        gdofs.push_back(g);

        // 等效节点载荷：f = t * A * B^T * sigma。
        const DenseVec fe = t * A * B.transpose() * sigma;
        for (int i = 0; i < 6; ++i) F(g[i]) += fe(i);
    }

    SparseMat K = assemble(kes, gdofs, nDofs);

    // 约束刚体位移：固定节点0，并约束节点1的 y 向位移。
    applyDirichlet(K, F, {0, 1, 3}, {0.0, 0.0, 0.0});

    const DenseVec U = solveSPD(K, F);

    std::cout << std::setprecision(10);
    std::cout << "2D 平面应力 CST（单位正方形 patch test）节点位移:\n";
    for (int i = 0; i < nNodes; ++i)
        std::cout << "  节点" << i << ": u=" << U(2 * i)
                  << ", v=" << U(2 * i + 1) << "\n";

    // 期望：u = 0.01*x，v = -0.0025*y。
    double err = 0.0;
    for (int i = 0; i < nNodes; ++i) {
        const double ue = 0.01 * x[i](0);
        const double ve = -0.0025 * x[i](1);
        err += std::hypot(U(2 * i) - ue, U(2 * i + 1) - ve);
    }
    std::cout << "累计位移误差 = " << err << "\n";
    std::cout << (err < 1e-12 ? "[PASS]" : "[FAIL]") << "\n";
    return 0;
}
