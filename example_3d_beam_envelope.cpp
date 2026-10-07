#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_beam_stage.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    auto stages = solveModel3DStages(m);
    auto env = compute3DDisplacementEnvelope(m, stages);
    auto beamEnv = compute3DBeamEnvelope(m, stages);
    for (const auto& e : env)
        std::cout << "N" << e.node << " dof" << e.dof
                  << " min=" << e.minVal << " max=" << e.maxVal << "\n";
    for (const auto& e : beamEnv)
        std::cout << "BEAM" << e.id
                  << " Nmin=" << e.Nmin << " Nmax=" << e.Nmax
                  << " MzMin=" << e.MzMin << " MzMax=" << e.MzMax << "\n";
    return 0;
}
