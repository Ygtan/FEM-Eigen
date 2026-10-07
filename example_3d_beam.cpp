#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_beam.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    auto U = solveModel3D(m);
    // 梁模型为 6 DOF/节点，节点2竖向位移索引 = 6*1 + 1
    std::cout << "node2 uy = " << U[7] << "\n";
    return 0;
}
