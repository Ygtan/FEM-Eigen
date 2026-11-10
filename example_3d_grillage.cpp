#include "fem_3d.hpp"
#include "modal_output_3d.hpp"
#include "buckling_output_3d.hpp"

#include <fstream>
#include <iomanip>
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

    // 各阶段末的索状态。
    //
    // 这一段是"施工阶段能不能调索"的唯一证据，也是与前几段性质不同的地方：
    // 上面的位移/内力都是**增量**，可以逐阶段累加；索力是**状态量**，
    // 累加没有意义（加出来是"张拉动作的总和"）。所以这里逐阶段列，
    // 不提供包络。
    //
    // 两个偏差各有用处：
    //   dev    本阶段末索力 vs 本阶段目标（0 = 本阶段没调这根索）
    //   devDes 本阶段末索力 vs 该索**设计目标**（跨阶段继承）
    // 后者才是施工控制关心的问题："先张好的索被后来的索带跑了多少"。
    // 逐对张拉时它一般不为零——这正是需要二次调索的原因，藏起来就等于
    // 假装一次张拉到位。
    {
        bool any = false;
        for (const auto& sr : stages) if (!sr.cables.empty()) { any = true; break; }
        if (any) {
            out << "stage cable tension\n";
            for (const auto& sr : stages)
                for (const auto& cs : sr.cables)
                    out << "STAGE" << sr.stage << " E" << cs.elem
                        << " T=" << cs.T << " L0=" << cs.L0
                        << " chord=" << cs.chord
                        << " target=" << cs.target << " dev=" << cs.dev
                        << " Tdes=" << cs.Tdes << " devDes=" << cs.devDes << "\n";
        }
    }

    // 施工阶段预应力：每阶段每束的有效预加力（束号 → P·loss）。
    // 前端靠它校核"有效预应力 vs 张拉控制应力"，以及分阶段看哪批束何时张拉。
    {
        bool any = false;
        for (const auto& sr : stages)
            if (!sr.prestressForces.empty()) { any = true; break; }
        if (any) {
            out << "stage prestress force\n";
            for (const auto& sr : stages)
                for (const auto& [id, Pe] : sr.prestressForces)
                    out << "STAGE" << sr.stage << " P" << id << " Pe=" << Pe << "\n";
        }
    }

    // 分项预应力损失（JTG 3362-2018 §6.2）：只有 LOSSREF 路径才输出。
    // 只给 LOSS= 或不写损失时内核根本没算分项，不输出一堆 0 冒充「各项为 0」。
    // 段头必须登记进 Python 的 sections 元组，否则这些 STAGE/P 行会被上一个
    // section 的正则收走，污染旧段。
    {
        bool any = false;
        for (const auto& sr : stages)
            if (!sr.prestressLossDetails.empty()) { any = true; break; }
        if (any) {
            out << "stage prestress loss\n";
            // 只给本段提高精度：σ 与 RATIO 要供报告复算，6 位有效数字在
            // 大张拉力下会让 P·RATIO 与另一段的 Pe 差出几十牛。写完立刻
            // 恢复，既有段（含 stage prestress force）逐位不变。
            const std::streamsize oldPrec = out.precision();
            out << std::setprecision(12);
            for (const auto& sr : stages)
                for (const auto& [id, d, sc] : sr.prestressLossDetails)
                    out << "STAGE" << sr.stage << " P" << id
                        << " SIGMA_CON=" << sc
                        << " SL1=" << d.sigmaL1
                        << " SL2=" << d.sigmaL2
                        << " SL3=" << d.sigmaL3
                        << " SL4=" << d.sigmaL4
                        << " SL5=" << d.sigmaL5
                        << " SL6=" << d.sigmaL6
                        << " SPE=" << d.sigmaPe
                        << " RATIO=" << (sc > 0.0 ? d.sigmaPe / sc : 0.0)
                        << "\n";
            out.precision(oldPrec);
        }
    }

    if (!m.influences.empty()) {
        auto live = compute3DStageLiveEnvelope(m, m.influences.front());
        out << "stage live envelope\n";
        for (auto& e : live)
            out << "STAGE" << e.stage << " max=" << e.maxVal << " min=" << e.minVal << "\n";
    }
    // 活载内力包络：荷载组合（1.2S_G + 1.4S_Q）必须有活载产生的 M、V，
    // 而上面那段只有影响点的位移。位移包络算不出承载力组合。
    if (!m.lanes.empty()) {
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
    }

    // ---- 模态：把内核本来就在算的东西写出来 -------------------------------
    // resolveLaneImpact3D 内部已经对**同一个成桥态构型**调了
    // computeNaturalFrequencies3D(mc, 8)，但只留下一个竖向基频用于 μ。这里
    // 用同样的 K/M 再输出前 8 阶的频率、周期、参与质量与节点振型，报告里
    // 的逐阶数值与 f1_vertical 因而天然同源，不存在"两套频率"。
    //
    // 不静默：无质量 / 无自由自由度 / 特征求解抛异常都要写一段 ERROR=…，
    // 让下游明确知道"这次没有模态结果"，而不是解析不到就当没有这回事。
    {
        const Model3D mc = finalStageModel3D(m);
        const DofMap3D dmM = buildDofMap3D(mc);
        const SparseMat3D Mm = assembleMassMatrix3D(mc, dmM);
        if (Mm.nonZeros() == 0) {
            out << modalErrorSection3D("modal", "no_mass");
        } else {
            try {
                const ModalResult3D modal = computeNaturalFrequencies3D(mc, 8);
                if (modal.frequencies.empty())
                    out << modalErrorSection3D("modal", "no_free_dof");
                else
                    out << formatModalSection3D(dmM, Mm, modal, "modal");
            } catch (const std::exception& e) {
                out << modalErrorSection3D("modal", "solver_exception");
                std::cerr << "[modal] " << e.what() << "\n";
            }
        }
    }

    // ---- 线性屈曲：把内核本来就在算的 λ 写出来 ---------------------------
    // computeBucklingFactors3D 解 (K + λ·Kg)·φ = 0，之前只在
    // example_3d_modal / modal_bench 里打印到控制台，报告拿不到。这里用
    // **成桥态构型 + 成桥态恒载合力**作参考荷载（阶段是增量法，见
    // finalStageBucklingModel3D），输出前 6 阶 λ 与临界荷载 λ·P_ref。
    //
    // 这是线弹性屈曲，不是规范稳定验算：λ 的含义、上限性质和未覆盖项由
    // 报告层写明；这里只保证数值来自内核同一次特征求解，不另算一套。
    // 不静默：无参考荷载 / 无屈曲模态 / 求解抛异常都写 ERROR 段。
    {
        const Model3D mb = finalStageBucklingModel3D(m);
        const DofMap3D dmB = buildDofMap3D(mb);
        const BucklingReference3D ref = bucklingReference3D(mb);
        if (ref.prefN <= 0.0) {
            out << bucklingErrorSection3D("buckling", "no_reference_load");
        } else {
            try {
                // λ 为空不写 ERROR：没有受压杆件时本就不该有屈曲模态，
                // STATUS=no_compression 把「物理上无压杆」与「求解失败」
                // 分开，报告层据此给 info 而不是假的"验算失败"。
                const BucklingResult3D bk = computeBucklingResult3D(mb, 6);
                out << formatBucklingSection3D(mb, dmB, bk, "buckling");
            } catch (const std::exception& e) {
                out << bucklingErrorSection3D("buckling", "solver_exception");
                std::cerr << "[buckling] " << e.what() << "\n";
            }
        }
    }
    std::cout << "summary written: " << path.substr(0,path.rfind(".dat"))+"_summary.txt" << "\n";
    return 0;
}
