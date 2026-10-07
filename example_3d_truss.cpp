#include "fem_3d.hpp"

#include <iostream>

using namespace fem;

int main(int argc, char** argv) {
    Model3D m;
    if (argc > 1) m = readModel3D(argv[1]);
    else {
        m.nodes[1] = Vec3(0, 0, 0);
        m.nodes[2] = Vec3(1, 0, 0);
        Element3D e;
        e.id = 1; e.type="TRUSS"; e.E=200e9; e.A=0.001; e.n1=1; e.n2=2;
        m.elements.push_back(e);
        m.constraints.push_back({1,1}); m.constraints.push_back({1,2}); m.constraints.push_back({1,3});
        m.constraints.push_back({2,2}); m.constraints.push_back({2,3});
        m.loads.push_back({2,1,1000.0});
    }
    const auto U = solveModel3D(m);
    std::cout << "node2 ux = " << U[3] << "\n";
    return 0;
}
