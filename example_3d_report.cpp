#include "fem_3d.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace fem;

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "verify_3d_beam.dat";
        if (argc > 1) path = argv[1];
        Model3D m = readModel3D(path);
        const std::vector<double> U = solveModel3D(m);
        const DofMap3D dm = buildDofMap3D(m);

        std::ofstream out(path.substr(0, path.rfind(".dat")) + ".report.txt");
        out << std::setprecision(8);

        out << "# FEM 3D 计算结果\n";
        out << "节点位移 (m, rad)\n";
        for (const auto& [id, x] : m.nodes) {
            out << id;
            const int n = dm.count(id);
            for (int k = 0; k < n; ++k)
                out << " " << U[dm.dof(id, k)];
            out << "\n";
        }

        out << "单元内力 (N, N*m)\n";
        for (const auto& e : m.elements) {
            out << e.id << " " << e.type;
            if (e.type == "BEAM") {
                const Eigen::VectorXd fl = beam3DLocalEndForces(m, dm, e, U);
                const BeamEndForces fi = beamEndForcesAt(fl, 12, 1);
                out << " N=" << fi.N << " Vy=" << fi.Vy << " Vz=" << fi.Vz
                    << " T=" << fi.T << " My=" << fi.My << " Mz=" << fi.Mz;
            } else if (e.type == "BEAM7") {
                const Eigen::VectorXd fl = beam7LocalEndForces(m, dm, e, U);
                const BeamEndForces fi = beamEndForcesAt(fl, 14, 1);
                out << " N=" << fi.N << " Vy=" << fi.Vy << " Vz=" << fi.Vz
                    << " T=" << fi.T << " My=" << fi.My << " Mz=" << fi.Mz
                    << " B=" << fi.B;
            } else {
                const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
                const double L = std::max(d.norm(), 1e-12);
                const Vec3 c = d / L;
                double dl = 0.0;
                for (int k = 0; k < 3; ++k)
                    dl += (U[dm.dof(e.n2, k)] - U[dm.dof(e.n1, k)]) * c(k);
                // 单元分布荷载/自重下需把 i 端等效节点荷载沿轴向的分量加回，
                // 与 fem_2d.cpp 的 f.N1 = Nrestore + fe_i 保持一致。
                const Eigen::VectorXd fe = elementEquivalentLoad3D(m, e);
                const double fe_i = fe(0) * c(0) + fe(1) * c(1) + fe(2) * c(2);
                out << " N=" << (e.E * e.A / L * dl + fe_i);
            }
            out << "\n";
        }
        std::cout << "report written: " << path.substr(0, path.rfind(".dat")) + ".report.txt" << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
