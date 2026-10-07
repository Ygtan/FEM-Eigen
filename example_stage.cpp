#include "fem_io.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace fem;

static SolverType parseSolver(const std::string& name) {
    std::string s = name;
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (s == "LDLT" || s == "SIMPLICIAL_LDLT") return SolverType::SimplicialLDLT;
    if (s == "LLT" || s == "SIMPLICIAL_LLT") return SolverType::SimplicialLLT;
    if (s == "LU" || s == "SPARSE_LU") return SolverType::SparseLU;
    if (s == "QR" || s == "SPARSE_QR") return SolverType::SparseQR;
    if (s == "CG" || s == "CONJUGATE_GRADIENT") return SolverType::ConjugateGradient;
    if (s == "BICGSTAB" || s == "BICG") return SolverType::BiCGSTAB;
    if (s == "AMGCL") return SolverType::AMGCL;
    if (s == "CHOLMOD") return SolverType::CHOLMOD;
    throw std::runtime_error("未知求解器: " + name);
}

static std::string baseName(const std::string& path) {
    std::string b = path;
    const std::size_t p = b.rfind(".dat");
    if (p != std::string::npos && p + 4 == b.size()) b = b.substr(0, p);
    return b;
}

static void writeStageResults(const std::string& path,
                              const Model2D& model,
                              const std::vector<StageResult>& stages) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法写入施工阶段结果: " + path);
    out << std::setprecision(12);
    out << "# 施工阶段正装分析结果\n";
    out << "# 模型: " << (model.title.empty() ? "未命名" : model.title) << "\n\n";

    for (const auto& sr : stages) {
        out << "STAGE " << sr.stage << "\n";
        out << "节点位移\n节点\tu\tv\ttheta\n";
        std::map<int, int> index;
        int n = 0;
        for (const auto& [id, c] : model.nodes) index[id] = n++;
        for (const auto& [id, c] : model.nodes) {
            const int d = 3 * index[id];
            out << id << '\t' << sr.U(d) << '\t' << sr.U(d + 1) << '\t' << sr.U(d + 2) << '\n';
        }

        out << "单元内力\n单元\t类型\tN1\tN2\tV1\tM1\tV2\tM2\n";
        for (const auto& f : sr.forces) {
            if (isRigid(f.type)) continue;
            out << f.id << '\t' << f.type << '\t' << f.N1 << '\t' << f.N2;
            if (f.type == "BEAM")
                out << '\t' << f.V1 << '\t' << f.M1 << '\t' << f.V2 << '\t' << f.M2;
            else
                out << "\t0\t0\t0\t0";
            out << '\n';
        }
        out << "\n";
    }
}

static void writeCamberResults(const std::string& path,
                               const std::vector<CamberResult>& cambers) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法写入预拱度结果: " + path);
    out << std::setprecision(12);
    out << "node\tdof\toriginal\tfinal_displacement\tcamber\tinstallation\n";
    for (const auto& c : cambers)
        out << c.node << '\t' << c.dof << '\t' << c.original << '\t'
            << c.finalDisplacement << '\t' << c.camber << '\t' << c.installation << '\n';
}

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "bridge_tied_arch_80m.dat";
        SolverType solver = SolverType::SimplicialLDLT;
        bool backward = false;
        bool iterate = false;
        bool camberIterate = false;
        bool feedback = false;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--solver" && i + 1 < argc) solver = parseSolver(argv[++i]);
            else if (a == "--backward") backward = true;
            else if (a == "--iterate") iterate = true;
            else if (a == "--camber-iterate") camberIterate = true;
            else if (a == "--feedback") feedback = true;
            else path = a;
        }

        const Model2D model = readModel2D(path);
        if (feedback) {
            const auto fc = computeFeedbackCorrection(model, solver);
            std::cout << "反馈修正：当前响应=" << fc.current
                      << " 灵敏度=" << fc.derivative
                      << " 修正量=" << fc.correction << "\n";
        } else if (camberIterate) {
            const auto it = iterateCamber(model, solver);
            const std::string outPath = baseName(path) + "_stages_camber_iterated.txt";
            writeStageResults(outPath, model, it.stages);
            std::cout << "已写入非线性预拱度迭代结果: " << outPath << "\n";
            for (int k = 0; k < it.cambers.size(); ++k)
                std::cout << "节点 " << model.cambers[k].node
                          << " dof " << model.cambers[k].dof
                          << " 预拱度=" << it.cambers(k)
                          << " 最终位移=" << it.finalDisplacements(k) << "\n";
            std::cout << "迭代次数=" << it.iterations << "\n";
        } else if (iterate) {
            const auto it = iterateStagePrestress(model, solver);
            const std::string outPath = baseName(path) + "_stages_iterated.txt";
            writeStageResults(outPath, model, it.stages);
            std::cout << "已写入迭代后施工阶段结果: " << outPath << "\n";
            for (int k = 0; k < it.scales.size(); ++k)
                std::cout << "控制变量 " << (k + 1) << " 缩放系数=" << it.scales(k) << "\n";
            for (int j = 0; j < it.finalValues.size(); ++j)
                std::cout << "目标 " << (j + 1) << " 响应=" << it.finalValues(j) << "\n";
            std::cout << "迭代次数=" << it.iterations << "\n";
        } else {
            const auto stages = backward ? solveModel2DStagesBackward(model, solver)
                                         : solveModel2DStages(model, solver);
            const std::string outPath = baseName(path) + (backward ? "_stages_backward.txt" : "_stages.txt");
            writeStageResults(outPath, model, stages);
            std::cout << "已写入施工阶段结果: " << outPath << "\n";
            for (const auto& sr : stages) {
                double maxAbs = 0.0;
                for (int i = 0; i < sr.U.size(); ++i) maxAbs = std::max(maxAbs, std::abs(sr.U(i)));
                std::cout << "STAGE " << sr.stage << " 最大绝对位移=" << maxAbs << "\n";
            }
            if (!backward) {
                const double res = forwardBackwardStageResidual(model, solver);
                std::cout << "正装-倒拆状态残差范数=" << res << "\n";
            }
            if (!model.cambers.empty()) {
                const auto cambers = computeCamber(model, solver);
                const std::string camberPath = baseName(path) + "_camber.txt";
                writeCamberResults(camberPath, cambers);
                std::cout << "已写入预拱度结果: " << camberPath << "\n";
                for (const auto& c : cambers)
                    std::cout << "节点 " << c.node << " dof " << c.dof
                              << " 预拱度=" << c.camber
                              << " 安装标高=" << c.installation << "\n";
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
