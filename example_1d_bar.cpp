#include "fem.hpp"

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    using namespace fem;
    setupConsoleUTF8();

    const double E = 200e9, A = 0.01;
    const int dim = 1;
    const int dofsPerNode = dim;
    const int nNodes = 4;
    const int nDofs = nNodes * dofsPerNode;

    // 节点坐标与单元连接（一维杆，沿 x 轴）。
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

    applyDirichlet(K, F, {0}, {0.0}); // 固定节点 0

    const DenseVec U = solveSPD(K, F);

    std::cout << std::setprecision(10);
    std::cout << "1D 杆（3 单元）节点位移:\n";
    for (int i = 0; i < nNodes; ++i)
        std::cout << "  u" << i << " = " << U(i) << "\n";

    DenseVec expected(nDofs);
    expected << 0.0, 5e-7, 1e-6, 1e-6;
    const double err = (U - expected).norm();
    std::cout << "误差 ||U - 解析解|| = " << err << "\n";
    std::cout << (err < 1e-15 ? "[PASS]" : "[FAIL]") << "\n";
    return 0;
}
