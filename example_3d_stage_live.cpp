#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_stage_live.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    for (const auto& spec : m.influences) {
        auto envs = compute3DStageLiveEnvelope(m, spec);
        for (const auto& e : envs)
            std::cout << "STAGE " << e.stage << " N" << spec.node
                      << " dof" << spec.dof << " max=" << e.maxVal << " min=" << e.minVal << "\n";
    }
    return 0;
}
