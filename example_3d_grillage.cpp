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
    // 六个分量全输出，不只 N/Mz。
    //   梁格里的纵向梁沿 X 走，竖向弯曲的弯矩落在 **My**（绕局部 y，水平轴），
    //   而 Mz 是横向平面内的弯曲。旧版只打印 Mz，于是竖向受弯主梁的弯矩
    //   全是 0——报告上看着"内力包络 222 个单元"，实际主弯矩一个没写。
    for (auto& e : beamEnv) {
        out << "E" << e.id
            << " Nmin=" << e.Nmin << " Nmax=" << e.Nmax
            << " VyMin=" << e.VyMin << " VyMax=" << e.VyMax
            << " VzMin=" << e.VzMin << " VzMax=" << e.VzMax
            << " Tmin=" << e.Tmin << " Tmax=" << e.Tmax
            << " MyMin=" << e.MyMin << " MyMax=" << e.MyMax
            << " MzMin=" << e.MzMin << " MzMax=" << e.MzMax
            << (e.valid ? "" : " invalid=1") << "\n";
    }
    // 成桥态（最后阶段）的恒载效应：运营期荷载组合的基准。
    // 上面的 beam force 是跨阶段包络，含半成品构型的临时内力，
    // 与运营期活载相加没有物理意义。
    auto nz = [](double v) { return v == 0.0 ? 0.0 : v; };
    auto endName = [](int end) { return end == 1 ? " i" : " j"; };
    auto finalBeam = compute3DFinalStageBeamForces(m, stages);
    out << "final stage beam force\n";
    for (auto& e : finalBeam)
        out << "E" << e.elem << endName(e.end)
            << " N=" << nz(e.N) << " Vy=" << nz(e.Vy) << " Vz=" << nz(e.Vz)
            << " T=" << nz(e.T) << " My=" << nz(e.My) << " Mz=" << nz(e.Mz) << "\n";
    auto finalDisp = compute3DFinalStageDisplacement(m, stages);
    out << "final stage displacement\n";
    for (auto& e : finalDisp)
        out << "N" << e.node << " dof" << e.dof << " u=" << nz(e.minVal) << "\n";

    if (!m.influences.empty()) {
        auto live = compute3DStageLiveEnvelope(m, m.influences.front());
        out << "stage live envelope\n";
        for (auto& e : live)
            out << "STAGE" << e.stage << " max=" << e.maxVal << " min=" << e.minVal << "\n";
    }
    // 活载内力包络：荷载组合（1.2S_G + 1.4S_Q）必须有活载产生的 M、V，
    // 而上面那段只有影响点的位移。位移包络算不出承载力组合。
    if (!m.lanes.empty()) {
        auto liveBeam = compute3DLiveBeamEnvelope(m);
        out << "live beam envelope\n";
        for (auto& e : liveBeam)
            out << "E" << e.elem << endName(e.end)
                << " MyMin=" << nz(e.MyMin) << " MyMax=" << nz(e.MyMax)
                << " VzMin=" << nz(e.VzMin) << " VzMax=" << nz(e.VzMax) << "\n";
    }
    std::cout << "summary written: " << path.substr(0,path.rfind(".dat"))+"_summary.txt" << "\n";
    return 0;
}
