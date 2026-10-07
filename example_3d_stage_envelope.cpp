#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_stage.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    auto stages = solveModel3DStages(m);
    auto env = compute3DStageEnvelope(m, stages);
    auto camber = compute3DCamber(m, stages);
    for (const auto& e : env)
        std::cout << "E" << e.id << " Nmin=" << e.Nmin << " Nmax=" << e.Nmax << "\n";
    for (const auto& c : camber)
        std::cout << "N" << c.node << " dof" << c.dof << " camber=" << c.camber << "\n";
    return 0;
}
