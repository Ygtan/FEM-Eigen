#include "fem.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    using namespace fem;
    setupConsoleUTF8();

    const double E = 200e9, A = 0.01;
    const int dim = 2;
    const int dofsPerNode = dim;
    const int nNodes = 3;
    const int nDofs = nNodes * dofsPerNode;

    // 两杆三角形桁架：节点 0(0,0)、1(1,0) 固定，节点 2(0.5,1) 受向下 1000 N。
    std::vector<Eigen::Vector3d> x = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 1.0, 0.0}};
    std::vector<std::pair<int, int>> conn = {{0, 2}, {1, 2}};

    std::vector<DenseMat> kes;
    std::vector<std::vector<int>> gdofs;
    for (const auto& [a, b] : conn) {
        TrussElement e{E, A, dim, x[a], x[b]};
        kes.push_back(e.stiffness());
        gdofs.push_back(nodeDofs({a, b}, dofsPerNode));
    }

    SparseMat K = assemble(kes, gdofs, nDofs);
    DenseVec F = DenseVec::Zero(nDofs);
    F(2 * 2 + 1) = -1000.0; // 节点 2 的 y 方向

    // 固定节点 0 与节点 1 的全部自由度。
    applyDirichlet(K, F, {0, 1, 2, 3}, {0.0, 0.0, 0.0, 0.0});

    const DenseVec U = solveSPD(K, F);

    std::cout << std::setprecision(10);
    std::cout << "2D 平面桁架（两杆三角形）节点位移:\n";
    for (int i = 0; i < nNodes; ++i)
        std::cout << "  节点" << i << ": u=" << U(2 * i)
                  << ", v=" << U(2 * i + 1) << "\n";

    // 解析解：顶点(节点2) u=0, v=-3.49385621e-7。
    const double err = std::hypot(U(4) - 0.0, U(5) + 3.49385621e-7);
    std::cout << "顶点位移误差 = " << err << "\n";
    std::cout << (err < 1e-14 ? "[PASS]" : "[FAIL]") << "\n";
    return 0;
}
