#include "fem_3d.hpp"

#include <fstream>
#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "bridge_3d_grillage.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    auto stages = solveModel3DStages(m);
    auto dispEnv = compute3DDisplacementEnvelope(m, stages);
    auto beamEnv = compute3DBeamEnvelope(m, stages);
    std::ofstream out(path.substr(0,path.rfind(".dat"))+"_summary.txt");
    out << "stage displacement envelope\n";
    for (auto& e : dispEnv)
        out << "N" << e.node << " dof" << e.dof << " min=" << e.minVal << " max=" << e.maxVal << "\n";
    out << "beam force envelope\n";
    for (auto& e : beamEnv)
        out << "E" << e.id << " Nmin=" << e.Nmin << " Nmax=" << e.Nmax
            << " MzMin=" << e.MzMin << " MzMax=" << e.MzMax << "\n";
    if (!m.influences.empty()) {
        auto live = compute3DStageLiveEnvelope(m, m.influences.front());
        out << "stage live envelope\n";
        for (auto& e : live)
            out << "STAGE" << e.stage << " max=" << e.maxVal << " min=" << e.minVal << "\n";
    }
    std::cout << "summary written: " << path.substr(0,path.rfind(".dat"))+"_summary.txt" << "\n";
    return 0;
}
