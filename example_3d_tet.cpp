#include "fem.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    using namespace fem;
    setupConsoleUTF8();

    const double E = 100.0, nu = 0.25;
    const int dofsPerNode = 3;
    const int nNodes = 4;
    const int nDofs = nNodes * dofsPerNode;

    // 单位四面体（由原点出发的三条单位棱边）。
    std::vector<Eigen::Vector3d> x = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    const std::array<int, 4> tet = {0, 1, 2, 3};

    // 目标应力 sigma_x = 1。
    Eigen::Matrix<double, 6, 1> sigma;
    sigma << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Tetra4 e{E, nu, x[0], x[1], x[2], x[3]};
    const DenseMat B = e.B();
    const double V = e.volume();
    const std::vector<int> g = nodeDofs({0, 1, 2, 3}, dofsPerNode);

    DenseVec F = DenseVec::Zero(nDofs);
    const DenseVec fe = V * B.transpose() * sigma; // 等效节点载荷
    for (int i = 0; i < 12; ++i) F(g[i]) += fe(i);

    SparseMat K = assemble({e.stiffness()}, {g}, nDofs);

    // 约束刚体位移：固定节点0，节点1 的 v/w，节点2 的 w。
    applyDirichlet(K, F, {0, 1, 2, 4, 5, 8}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    const DenseVec U = solveSPD(K, F);

    std::cout << std::setprecision(10);
    std::cout << "3D 四面体（patch test）节点位移:\n";
    for (int i = 0; i < nNodes; ++i)
        std::cout << "  节点" << i << ": u=" << U(3 * i)
                  << ", v=" << U(3 * i + 1) << ", w=" << U(3 * i + 2) << "\n";

    // 期望：u = 0.01*x，v = -0.0025*y，w = -0.0025*z。
    double err = 0.0;
    for (int i = 0; i < nNodes; ++i) {
        const double ue = 0.01 * x[i](0);
        const double ve = -0.0025 * x[i](1);
        const double we = -0.0025 * x[i](2);
        err += std::hypot(U(3 * i) - ue,
                          std::hypot(U(3 * i + 1) - ve, U(3 * i + 2) - we));
    }
    std::cout << "累计位移误差 = " << err << "\n";
    std::cout << (err < 1e-12 ? "[PASS]" : "[FAIL]") << "\n";
    return 0;
}
