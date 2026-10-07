#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_influence.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    for (const auto& spec : m.influences) {
        double totalMax=-1e300, totalMin=1e300;
        for (const auto& lane : m.lanes) {
            auto eta = compute3DInfluence(m, spec, lane);
            auto env = envelope3DFromInfluence(lane, eta);
            totalMax = std::max(totalMax, env.first);
            totalMin = std::min(totalMin, env.second);
        }
        std::cout << "N" << spec.node << " dof" << spec.dof
                  << " max=" << totalMax << " min=" << totalMin << "\n";
    }
    return 0;
}
