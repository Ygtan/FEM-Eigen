// 非线性施工阶段的文件接口。
// 读入与线性阶段相同的 .dat，调用 solveModel3DStagesNonlinear。
// 每个阶段的位移是当前结构的全量平衡，不能按线性增量累加。
#include "fem_3d_nl.hpp"
#include "newton_output_3d.hpp"   // [NEWTON]：每阶段内层共转列式 Newton 过程

#include <fstream>
#include <iostream>
#include <string>

using namespace fem;

int main(int argc, char** argv) {
    setupConsoleUTF8();
    if (argc < 2) {
        std::cerr << "用法: example_3d_stage_nl <model.dat>\n";
        return 2;
    }
    const std::string path = argv[1];
    try {
        Model3D m = readModel3D(path);
        NonlinearOptions3D nl;
        nl.loadSteps = 4;
        nl.maxIter = 25;
        const auto stages = solveModel3DStagesNonlinear(m, nl);

        const std::string outPath = path + ".nlstage.txt";
        std::ofstream out(outPath);
        if (!out) throw std::runtime_error("无法写出 " + outPath);
        out << "# total equilibrium; do not sum stage displacements\n";
        out << "NSTAGES " << stages.size() << "\n";
        for (const auto& s : stages) {
            out << "STAGE " << s.id << " converged=" << (s.converged ? 1 : 0)
                << " ndof=" << s.u.size() << "\n";
            out << "[U]\n";
            for (std::size_t i = 0; i < s.dm.nodeIds.size(); ++i) {
                const int node = s.dm.nodeIds[i];
                const int n = std::min(s.dm.dofCount[i], 3);
                for (int k = 0; k < n; ++k) {
                    const int d = s.dm.dofStart[i] + k;
                    if (d < 0 || d >= static_cast<int>(s.u.size())) continue;
                    out << node << " " << (k + 1) << " " << s.u[d] << "\n";
                }
            }
            out << "[CABLE]\n";
            for (const auto& [eid, T] : s.cableForce) {
                const auto it = s.cableL0.find(eid);
                const double L0 = (it == s.cableL0.end()) ? 0.0 : it->second;
                out << eid << " " << T << " " << L0 << "\n";
            }
            out << "[BEAM]\n";
            for (const auto& [eid, ends] : s.beamEnds) {
                const BeamEnd3D* pair[2] = {&ends.first, &ends.second};
                for (int end = 0; end < 2; ++end) {
                    const BeamEnd3D& f = *pair[end];
                    out << eid << " " << (end + 1)
                        << " " << f.N << " " << f.Vy << " " << f.Vz
                        << " " << f.T << " " << f.My << " " << f.Mz << "\n";
                }
            }
            // 内层 Newton 过程：每阶段一段。与 [LENGTH] 的外层迭代不是一层；
            // 本算例的 TENSION= 调索另由 CableLengthResult3D 处理，不在结果
            // 文件里重复堆中间求解。阶段位移是全量平衡，不能跨阶段累加。
            if (!s.newton.steps.empty())
                out << formatNewtonSection3D(s.newton, "[NEWTON]", s.id);
            else
                out << newtonErrorSection3D("[NEWTON]", "no_newton_solution");
            out << "[ENDSTAGE]\n";
        }
        // 没有索时，活载仍是参考几何上的线弹性影响线，与 example_3d_grillage 相同。
        // 有索时用末阶段切线刚度；冲击系数的基频也取这套切线，不用 EA/L。
        bool hasCable = false;
        for (const auto& e : m.elements)
            if (e.type == "CABLE") hasCable = true;
        if (!m.lanes.empty() && hasCable && !stages.empty() && stages.back().liveTangent) {
            auto endName = [](int end) { return end == 1 ? " i" : " j"; };
            auto nz = [](double v) { return v == 0.0 ? 0.0 : v; };
            const auto& back = stages.back();
            out << "live basis\nbasis=tangent\n";
            if (back.liveImpact.automatic) {
                out << "lane impact\n";
                out << "f1_vertical=" << back.liveImpact.f1
                    << " impact_mu=" << back.liveImpact.mu
                    << " impact_factor=" << back.liveImpact.factor << "\n";
            }
            out << "live beam envelope\n";
            for (const auto& e : back.liveBeam)
                out << "E" << e.elem << endName(e.end)
                    << " MyMin=" << nz(e.MyMin) << " MyMax=" << nz(e.MyMax)
                    << " VzMin=" << nz(e.VzMin) << " VzMax=" << nz(e.VzMax) << "\n";
            if (!back.liveDisp.empty()) {
                out << "stage live envelope\n";
                for (const auto& e : back.liveDisp)
                    out << "STAGE" << e.stage << " max=" << e.maxVal
                        << " min=" << e.minVal << "\n";
            }
        }
        if (!m.lanes.empty() && !hasCable) {
            auto endName = [](int end) { return end == 1 ? " i" : " j"; };
            auto nz = [](double v) { return v == 0.0 ? 0.0 : v; };
            LaneImpact3D imp;
            auto liveBeam = compute3DLiveBeamEnvelope(m, &imp);
            if (imp.automatic) {
                out << "lane impact\n";
                out << "f1_vertical=" << imp.f1
                    << " impact_mu=" << imp.mu
                    << " impact_factor=" << imp.factor << "\n";
            }
            out << "live beam envelope\n";
            for (auto& e : liveBeam)
                out << "E" << e.elem << endName(e.end)
                    << " MyMin=" << nz(e.MyMin) << " MyMax=" << nz(e.MyMax)
                    << " VzMin=" << nz(e.VzMin) << " VzMax=" << nz(e.VzMax) << "\n";
            if (!m.influences.empty()) {
                auto live = compute3DStageLiveEnvelope(m, m.influences.front());
                out << "stage live envelope\n";
                for (auto& e : live)
                    out << "STAGE" << e.stage << " max=" << e.maxVal
                        << " min=" << e.minVal << "\n";
            }
        }
        std::cout << "wrote " << outPath << " stages=" << stages.size() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
