#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_stage.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    auto stages = solveModel3DStages(m);
    for (const auto& s : stages) {
        std::cout << "STAGE " << s.stage << " node2 ux=" << s.U[3] << "\n";
    }
    return 0;
}
