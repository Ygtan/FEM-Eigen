// ============================================================
// beam_nl_probe.cpp —— 大转动梁柱效应（共转列式）内核探针 CLI
//
// 目的：把 fem_3d_nl.hpp 里**已实现但只在 test_3d_nl.cpp 硬编码断言**
// 的 CR 梁非线性求解能力，暴露成文件驱动、机器可读的 CSV，交给
// tests/benchmark/beam_nl_ansys_3d/ 与 ANSYS BEAM188 + NLGEOM 做端到端
// 对标。这里不重写任何力学：模型构造、约束、荷载与
// test_3d_nl.cpp::testCRCantileverCircle / testCRBeamColumn 完全一致，
// 平衡求解直接调 solveModel3DNonlinear()（共转 + 解析几何刚度）。
//
// 两个工况族（都是桥梁"梁柱效应 / 大转动"的经典解析标尺）：
//   CIRCLE : 悬臂梁端弯矩 M = lambda * 2*pi*E*Iz/L
//            纯弯曲 elastica -> 圆弧闭式 theta=M*L/(E*Iz),
//            x_tip = R*sin(theta), y_tip = R*(1-cos(theta)), R=L/theta.
//            lambda=1 时梁卷成整圆、尖端回到根部、转角 2*pi。
//   PDELTA : 简支梁轴压 P = lambda*Pcr + 跨中横向力 Q
//            Timoshenko 梁柱精确解
//            y_ex = Q/(2*P*k) * (tan(k*L/2) - k*L/2), k=sqrt(P/(E*Iz));
//            几何刚度缺失时放大系数会退化到 1。
//
// 输入（纯 ASCII 参数文件，由 BridgeStudio/ansys_bridge_beamnl.py 同源导出）：
//     G <gid> E=... G=... L=... A=... IY=... IZ=... J=... NELEM=...
//            LOADSTEPS=... MAXITER=... TOL=... Q=...
//     C <cid> <gid> <CIRCLE|PDELTA> <lambda>
//
// 输出 CSV（无表头，首列为记录类型）：
//     N,<cid>,<node>,<ux>,<uy>,<uz>,<rotx>,<roty>,<rotz>
//     S,<cid>,<key>,<value>
//   S 行给出同源参数、解析标尺、FEM 派生量与收敛信息，便于
//   compare 脚本同时做「ANSYS 黄金值」与「FEM 重跑逐位一致」两层校验。
//
// 用法：beam_nl_probe.exe <params.beam> <out.csv>
// ============================================================
#include "fem_3d.hpp"
#include "fem_3d_nl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace fem;

namespace {

constexpr double PI = 3.14159265358979323846;

struct Group {
    std::string id;
    double E = 0.0, G = 0.0, L = 0.0;
    double A = 0.0, IY = 0.0, IZ = 0.0, J = 0.0;
    double Q = 0.0, TOL = 1e-10;
    int NELEM = 0, LOADSTEPS = 0, MAXITER = 0;
};

struct Case {
    std::string id, gid, kind;
    double lambda = 0.0;
};

std::string trim(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::map<std::string, double> parseKV(const std::string& s) {
    std::map<std::string, double> out;
    std::istringstream ss(s);
    std::string tok;
    while (ss >> tok) {
        const std::size_t p = tok.find('=');
        if (p == std::string::npos) continue;
        out[tok.substr(0, p)] = std::stod(tok.substr(p + 1));
    }
    return out;
}

double get(const std::map<std::string, double>& kv, const std::string& k,
           double def = 0.0) {
    const auto it = kv.find(k);
    return it == kv.end() ? def : it->second;
}

struct CaseResult {
    std::vector<std::array<double, 6>> nodal;
    std::vector<std::pair<std::string, double>> summary;
};

CaseResult solveOne(const Group& g, const Case& c) {
    Model3D m;
    const double dx = g.L / static_cast<double>(g.NELEM);
    for (int i = 0; i <= g.NELEM; ++i)
        m.nodes[i + 1] = Vec3(i * dx, 0.0, 0.0);
    for (int i = 0; i < g.NELEM; ++i) {
        Element3D e;
        e.id = i + 1;
        e.type = "BEAM";
        e.E = g.E; e.A = g.A; e.Iy = g.IY; e.Iz = g.IZ; e.J = g.J; e.G = g.G;
        e.n1 = i + 1; e.n2 = i + 2;
        m.elements.push_back(e);
    }
    const int tip = g.NELEM + 1;
    const int mid = g.NELEM / 2 + 1;

    double M = 0.0, P = 0.0, theta = 0.0, Pcr = 0.0;
    double xtip_ex = 0.0, ytip_ex = 0.0, rotex = 0.0;
    double yex = 0.0, ylin = 0.0;

    if (c.kind == "CIRCLE") {
        for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0, 0.0});
        M = c.lambda * 2.0 * PI * g.E * g.IZ / g.L;
        m.loads.emplace_back(tip, 6, M);           // 尖端绕 z 的死弯矩
        theta = M * g.L / (g.E * g.IZ);
        const double R = (std::abs(theta) > 1e-14) ? g.L / theta : 1e300;
        xtip_ex = R * std::sin(theta);
        ytip_ex = R * (1.0 - std::cos(theta));
        rotex = theta;
    } else if (c.kind == "PDELTA") {
        for (int d = 1; d <= 4; ++d) m.constraints.push_back({1, d, 0.0, 0.0});
        m.constraints.push_back({tip, 2, 0.0, 0.0});
        m.constraints.push_back({tip, 3, 0.0, 0.0});
        Pcr = PI * PI * g.E * g.IZ / (g.L * g.L);
        P = c.lambda * Pcr;
        m.loads.emplace_back(tip, 1, -P);
        m.loads.emplace_back(mid, 2, -g.Q);
        ylin = g.Q * g.L * g.L * g.L / (48.0 * g.E * g.IZ);
        if (std::abs(P) > 1e-12) {
            const double k = std::sqrt(P / (g.E * g.IZ));
            yex = (g.Q / (2.0 * P * k)) *
                  (std::tan(k * g.L / 2.0) - k * g.L / 2.0);
        } else {
            yex = ylin;
        }
    } else {
        throw std::runtime_error("unknown case kind: " + c.kind);
    }

    NonlinearOptions3D opt;
    opt.loadSteps = g.LOADSTEPS;
    opt.maxIter = g.MAXITER;
    opt.tolForce = g.TOL;
    const NonlinearResult3D r = solveModel3DNonlinear(m, opt);
    if (!r.converged)
        throw std::runtime_error("nonlinear solve did not converge: " + c.id);

    const DofMap3D dm = buildDofMap3D(m);
    CaseResult out;
    out.nodal.assign(g.NELEM + 1, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    for (int i = 1; i <= g.NELEM + 1; ++i) {
        for (int k = 0; k < 6; ++k) {
            const int d = dm.safeDof(i, k);
            out.nodal[i - 1][k] = (d >= 0 && d < static_cast<int>(r.u.size()))
                                      ? r.u[d] : 0.0;
        }
    }
    const auto U = [&](int node, int k) {
        const int d = dm.safeDof(node, k);
        return (d >= 0 && d < static_cast<int>(r.u.size())) ? r.u[d] : 0.0;
    };
    const double xtip_fem = g.L + U(tip, 0);
    const double ytip_fem = U(tip, 1);
    const double rot_fem = U(tip, 5);
    const double ymid_fem = -U(mid, 1);

    out.summary = {
        {"E", g.E}, {"G", g.G}, {"L", g.L}, {"A", g.A},
        {"IY", g.IY}, {"IZ", g.IZ}, {"J", g.J}, {"NELEM", static_cast<double>(g.NELEM)},
        {"LOADSTEPS", static_cast<double>(g.LOADSTEPS)},
        {"MAXITER", static_cast<double>(g.MAXITER)}, {"TOL", g.TOL},
        {"Q", g.Q}, {"lambda", c.lambda},
        {"M", M}, {"P", P}, {"Pcr", Pcr}, {"theta", theta},
        {"xtip_ex", xtip_ex}, {"ytip_ex", ytip_ex}, {"rot_ex", rotex},
        {"y_pd_ex", yex}, {"y_lin", ylin},
        {"xtip_fem", xtip_fem}, {"ytip_fem", ytip_fem},
        {"rot_fem", rot_fem}, {"ymid_fem", ymid_fem},
        {"amp_fem", (c.kind == "PDELTA" && std::abs(ylin) > 1e-30)
                        ? ymid_fem / ylin : 0.0},
        {"converged", 1.0},
        {"iterations", static_cast<double>(r.totalIterations)},
        {"steps", static_cast<double>(r.steps)},
    };
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: beam_nl_probe.exe <params.beam> <out.csv>\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open params: " << argv[1] << "\n";
        return 2;
    }
    std::map<std::string, Group> groups;
    std::vector<Case> cases;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        std::istringstream ss(t);
        std::string tag;
        ss >> tag;
        if (tag == "G") {
            Group g;
            ss >> g.id;
            std::string rest;
            std::getline(ss, rest);
            const auto kv = parseKV(rest);
            g.E = get(kv, "E");       g.G = get(kv, "G");
            g.L = get(kv, "L");       g.A = get(kv, "A");
            g.IY = get(kv, "IY");     g.IZ = get(kv, "IZ");
            g.J = get(kv, "J");       g.Q = get(kv, "Q");
            g.TOL = get(kv, "TOL", 1e-10);
            g.NELEM = static_cast<int>(get(kv, "NELEM"));
            g.LOADSTEPS = static_cast<int>(get(kv, "LOADSTEPS"));
            g.MAXITER = static_cast<int>(get(kv, "MAXITER"));
            if (!(g.E > 0.0 && g.L > 0.0 && g.NELEM > 0 &&
                  g.LOADSTEPS > 0 && g.MAXITER > 0))
                throw std::runtime_error("invalid group: " + g.id);
            groups[g.id] = g;
        } else if (tag == "C") {
            Case c;
            if (!(ss >> c.id >> c.gid >> c.kind >> c.lambda))
                throw std::runtime_error("invalid case line: " + t);
            cases.push_back(c);
        }
    }
    if (cases.empty()) {
        std::cerr << "no cases in params file\n";
        return 2;
    }

    std::ofstream out(argv[2]);
    if (!out) {
        std::cerr << "cannot write: " << argv[2] << "\n";
        return 2;
    }
    out << std::setprecision(14);
    for (const auto& c : cases) {
        const auto it = groups.find(c.gid);
        if (it == groups.end())
            throw std::runtime_error("unknown group: " + c.gid);
        const CaseResult r = solveOne(it->second, c);
        for (std::size_t i = 0; i < r.nodal.size(); ++i) {
            out << "N," << c.id << "," << (i + 1);
            for (double v : r.nodal[i]) out << "," << v;
            out << "\n";
        }
        for (const auto& kv : r.summary)
            out << "S," << c.id << "," << kv.first << "," << kv.second << "\n";
    }
    return 0;
}
