#include "fem_3d.hpp"

#include <iostream>
#include <fstream>

using namespace fem;

int main(int argc, char** argv) {
    std::string path = "verify_3d_beam_live.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    std::ofstream report(path.substr(0,path.rfind(".dat"))+"_live_envelope.csv");
    report << "elem,component,end,max,min\n";
    for (const auto& ei : m.elementInfluences) {
        for (const auto& lane : m.lanes) {
            auto eta = compute3DBeamAxialInfluence(m, ei.elem, lane);
            auto env = envelope3DFromInfluence(lane, eta);
            std::cout << "E" << ei.elem << " Nmax=" << env.first << " Nmin=" << env.second << "\n";
            report << ei.elem << ",N,1," << env.first << "," << env.second << "\n";
        }
    }
    for (const auto& e : m.elements) {
        if (e.type != "BEAM") continue;
        for (const auto& lane : m.lanes) {
            for (const std::string comp : {"MY","MZ"}) {
                auto eta = compute3DBeamMomentInfluence(m, e.id, comp, 1, lane);
                auto env = envelope3DFromInfluence(lane, eta);
                std::cout << "E" << e.id << " " << comp << "1 max=" << env.first << " min=" << env.second << "\n";
                report << e.id << "," << comp << ",1," << env.first << "," << env.second << "\n";
            }
            for (const std::string comp : {"VY","VZ"}) {
                auto eta = compute3DBeamShearInfluence(m, e.id, comp, 1, lane);
                auto env = envelope3DFromInfluence(lane, eta);
                std::cout << "E" << e.id << " " << comp << "1 max=" << env.first << " min=" << env.second << "\n";
                report << e.id << "," << comp << ",1," << env.first << "," << env.second << "\n";
            }
            if (e.Wtop > 0 && e.Wbot > 0) {
                for (const std::string comp : {"TOP","BOT"}) {
                    auto eta = compute3DBeamStressInfluence(m, e.id, comp, 1, lane);
                    auto env = envelope3DFromInfluence(lane, eta);
                    std::cout << "E" << e.id << " " << comp << "1 max=" << env.first << " min=" << env.second << "\n";
                    report << e.id << "," << comp << ",1," << env.first << "," << env.second << "\n";
                }
            }
        }
    }
    return 0;
}
