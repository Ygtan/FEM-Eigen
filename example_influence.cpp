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

static void writeInfluenceFile(const std::string& path,
                               const Model2D& model,
                               const Model2D::LaneLoad& lane,
                               const std::vector<Model2D::InfluenceRequest>& reqs,
                               SolverType solver) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法写入影响线文件: " + path);
    out << std::setprecision(12);
    out << "# 影响线（单位荷载沿 " << lane.dir << " 方向）\n";
    out << "# LANE " << lane.id << "  q=" << lane.q << "  P=" << lane.P << "\n";

    for (const auto& req : reqs) {
        const InfluenceCurve c = computeInfluenceCurve(model, lane, req, solver);
        out << "\n# " << req.kind << " id=" << req.id;
        if (req.kind == "DISP" || req.kind == "REACTION")
            out << " node=" << req.node << " dof=" << req.dof;
        else
            out << " elem=" << req.elem << " comp=" << req.component << " end=" << req.end;
        out << "\n# max=" << c.vmax << " @s=" << c.smax
            << "  min=" << c.vmin << " @s=" << c.smin << "\n";
        out << "s\teta\n";
        for (std::size_t k = 0; k < c.s.size(); ++k)
            out << c.s[k] << '\t' << c.eta[k] << '\n';
    }
}

static void writeEnvelopeFile(const std::string& path,
                              const Model2D& model,
                              const Model2D::LaneLoad& lane,
                              SolverType solver) {
    const auto specs = defaultEnvelopeSpecs(model);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法写入包络文件: " + path);
    out << std::setprecision(12);
    out << "# 车道荷载内力/应力包络\n";
    out << "# LANE " << lane.id << "  dir=" << lane.dir
        << "  q=" << lane.q << "  P=" << lane.P << "\n";
    out << "label\tmax\tmin\tsmax\tsmin\n";

    for (const auto& sp : specs) {
        EnvelopeValue ev;
        if (lane.method == "NONLINEAR") {
            ev = nonlinearEnvelopeFromLane(model, lane, sp.req, solver);
        } else {
            const InfluenceCurve c = computeInfluenceCurve(model, lane, sp.req, solver);
            ev = envelopeFromInfluence(model, lane, c);
        }
        out << sp.label << '\t' << ev.maxVal << '\t' << ev.minVal
            << '\t' << ev.smax << '\t' << ev.smin << '\n';
    }
}

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "mixed_truss_beam.dat";
        SolverType solver = SolverType::SimplicialLDLT;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--solver" && i + 1 < argc) {
                solver = parseSolver(argv[++i]);
            } else {
                path = a;
            }
        }

        const Model2D model = readModel2D(path);
        if (model.laneLoads.empty())
            throw std::runtime_error("输入文件没有定义 *LANE_LOAD");
        if (model.influences.empty())
            std::cout << "未定义 *INFLUENCE，仅输出默认内力/应力包络。\n";

        const std::string base = baseName(path);
        for (const auto& lane : model.laneLoads) {
            if (!model.influences.empty()) {
                const std::string fname = base + "_influence_lane" + std::to_string(lane.id) + ".txt";
                writeInfluenceFile(fname, model, lane, model.influences, solver);
                std::cout << "已写入影响线: " << fname << "\n";
            }

            const std::string fname = base + "_envelope_lane" + std::to_string(lane.id) + ".txt";
            writeEnvelopeFile(fname, model, lane, solver);
            std::cout << "已写入包络: " << fname << "\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
