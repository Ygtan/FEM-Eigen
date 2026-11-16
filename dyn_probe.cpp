// dyn_probe.cpp — 动力 snap-through / Newmark-β 隐式时程探针（文件驱动 CLI）
//
// 第二十一次力学⑤：独立 opt-in 模块 dynamic_newmark_3d.hpp 的入口。
// 不修改 snap_probe.cpp / stability_path_3d.hpp / fem_3d_nl.hpp 的任何既有路径。
//
// 输入 *.dyn（纯 ASCII）：
//   G <gid> E=.. A=.. B=.. H=.. PREF=.. DENS=.. [ADDMASS=..]
//   C <cid> <gid> DYN=1 DT=.. TEND=.. [BETA=.. GAMMA=.. MASSMODE=..
//          CMASS=.. CSTIFF=.. LOADTYPE=.. TOL=.. MAXITER=.. SNAP=.. STRAIN=.. TAG=..]
//   DYN=0 的 case 只输出 S,dyn,0，不跑时程（负例 / 默认零污染守卫用）。
//
// 物理模型：von Mises 两杆浅桁架，两端固定，顶点 Z 约束、Y 向荷载 −PREF。
// 独立解析锚点（与 crElementInternal 轴向单元同本构 N=EA(L−L0)/L0）：
//   L(v)=hypot(B,H−v),  P(v)=2EA(L0−L)(H−v)/(L0·L)
//   一致质量顶点平动等效质量 m=2·(dens·A+addmass)·L0/3
//   集中质量（行和）顶点等效质量 m=(dens·A+addmass)·L0
//
// 输出 CSV（无表头，首列记录类型）：
//   T,<cid>,<step>,<t>,<lambda>,<v>,<vel>,<acc>,<residual>,<iter>,
//       <kinetic>,<strain>,<external_work>
//   S,<cid>,<key>,<value>
#include "dynamic_newmark_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Group {
    double E = 2.1e11, A = 1.0e-3, B = 2.0, H = 0.5;
    double PREF = 1.0e6, DENS = 7850.0, ADDMASS = 0.0;
};

struct Case {
    std::string id, gid, tag;
    int DYN = 0;
    double DT = 1.0e-4, TEND = 0.05;
    double BETA = 0.25, GAMMA = 0.5;
    int MASSMODE = 0, LOADTYPE = 0, STRAIN = 0;
    double CMASS = 0.0, CSTIFF = 0.0, SNAP = 0.0;
    double TOL = 1.0e-7;
    int MAXITER = 30;
};

std::vector<std::string> split(const std::string& s) {
    std::istringstream is(s);
    std::vector<std::string> t;
    std::string x;
    while (is >> x) t.push_back(x);
    return t;
}

std::map<std::string, double> parseKV(const std::vector<std::string>& tok, int from) {
    std::map<std::string, double> kv;
    for (std::size_t i = from; i < tok.size(); ++i) {
        const auto pos = tok[i].find('=');
        if (pos == std::string::npos) continue;
        try {
            kv[tok[i].substr(0, pos)] = std::stod(tok[i].substr(pos + 1));
        } catch (const std::exception&) {
            // 非数值键（TAG=...）忽略
        }
    }
    return kv;
}

double get(const std::map<std::string, double>& kv, const std::string& k, double d) {
    auto it = kv.find(k);
    return it == kv.end() ? d : it->second;
}

double pOfV(double v, const Group& g, double L0) {
    const double L = std::hypot(g.B, g.H - v);
    return 2.0 * g.E * g.A * (L0 - L) * (g.H - v) / (L0 * L);
}

double analyticLimitV(const Group& g, double L0) {
    double lo = 0.0, hi = g.H;
    for (int i = 0; i < 200; ++i) {
        const double m1 = lo + (hi - lo) * 0.3819660112501051;
        const double m2 = lo + (hi - lo) * 0.6180339887498949;
        if (pOfV(m1, g, L0) < pOfV(m2, g, L0)) lo = m1;
        else hi = m2;
    }
    return 0.5 * (lo + hi);
}

fem::Model3D makeModel(const Group& g, double L0) {
    fem::Model3D m;
    m.nodes[1] = fem::Vec3(-g.B, 0.0, 0.0);
    m.nodes[2] = fem::Vec3(0.0, g.H, 0.0);
    m.nodes[3] = fem::Vec3(g.B, 0.0, 0.0);
    for (int i = 0; i < 2; ++i) {
        fem::Element3D e;
        e.id = i + 1;
        e.type = "TRUSS";
        e.E = g.E; e.A = g.A; e.L0 = L0;
        e.dens = g.DENS; e.addmass = g.ADDMASS;
        e.n1 = i + 1; e.n2 = i + 2;
        m.elements.push_back(e);
    }
    for (int d = 1; d <= 3; ++d) {
        m.constraints.push_back({1, d, 0.0, 0.0});
        m.constraints.push_back({3, d, 0.0, 0.0});
    }
    m.constraints.push_back({2, 3, 0.0, 0.0});
    m.loads.emplace_back(2, 2, -g.PREF);
    return m;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: dyn_probe <input.dyn> <output.csv>\n");
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::fprintf(stderr, "cannot open input: %s\n", argv[1]);
        return 2;
    }
    std::map<std::string, Group> groups;
    std::vector<Case> cases;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto t = split(line);
        if (t.empty() || t[0][0] == '#' || t[0][0] == '!') continue;
        if (t[0] == "G" && t.size() >= 2) {
            Group g;
            const auto kv = parseKV(t, 2);
            g.E = get(kv, "E", g.E); g.A = get(kv, "A", g.A);
            g.B = get(kv, "B", g.B); g.H = get(kv, "H", g.H);
            g.PREF = get(kv, "PREF", g.PREF); g.DENS = get(kv, "DENS", g.DENS);
            g.ADDMASS = get(kv, "ADDMASS", g.ADDMASS);
            groups[t[1]] = g;
        } else if (t[0] == "C" && t.size() >= 3) {
            Case c;
            c.id = t[1]; c.gid = t[2];
            const auto kv = parseKV(t, 3);
            c.DYN = static_cast<int>(get(kv, "DYN", c.DYN));
            c.DT = get(kv, "DT", c.DT); c.TEND = get(kv, "TEND", c.TEND);
            c.BETA = get(kv, "BETA", c.BETA); c.GAMMA = get(kv, "GAMMA", c.GAMMA);
            c.MASSMODE = static_cast<int>(get(kv, "MASSMODE", c.MASSMODE));
            c.LOADTYPE = static_cast<int>(get(kv, "LOADTYPE", c.LOADTYPE));
            c.STRAIN = static_cast<int>(get(kv, "STRAIN", c.STRAIN));
            c.CMASS = get(kv, "CMASS", c.CMASS);
            c.CSTIFF = get(kv, "CSTIFF", c.CSTIFF);
            c.SNAP = get(kv, "SNAP", c.SNAP);
            c.TOL = get(kv, "TOL", c.TOL);
            c.MAXITER = static_cast<int>(get(kv, "MAXITER", c.MAXITER));
            for (const auto& x : t)
                if (x.rfind("TAG=", 0) == 0) c.tag = x.substr(4);
            cases.push_back(c);
        }
    }
    if (cases.empty()) {
        std::fprintf(stderr, "no cases in %s\n", argv[1]);
        return 2;
    }

    std::ofstream out(argv[2]);
    out << std::setprecision(14);
    auto S = [&](const std::string& cid, const char* k, double v) {
        out << "S," << cid << "," << k << "," << v << "\n";
    };

    for (const auto& c : cases) {
        auto gi = groups.find(c.gid);
        if (gi == groups.end()) {
            std::fprintf(stderr, "unknown group %s\n", c.gid.c_str());
            return 2;
        }
        const Group& g = gi->second;
        const double L0 = std::hypot(g.B, g.H);
        const double vcr = analyticLimitV(g, L0);
        const double pcr = pOfV(vcr, g, L0);
        const double rhoA = g.DENS * g.A + g.ADDMASS;

        S(c.id, "dyn", static_cast<double>(c.DYN));
        S(c.id, "E", g.E); S(c.id, "A", g.A);
        S(c.id, "B", g.B); S(c.id, "H", g.H); S(c.id, "L0", L0);
        S(c.id, "PREF", g.PREF); S(c.id, "DENS", g.DENS);
        S(c.id, "ADDMASS", g.ADDMASS);
        S(c.id, "analytic_v_cr", vcr); S(c.id, "analytic_P_cr", pcr);
        S(c.id, "analytic_lambda_cr", pcr / g.PREF);
        S(c.id, "rhoA", rhoA);
        S(c.id, "m_consistent", 2.0 * rhoA * L0 / 3.0);
        S(c.id, "m_lumped", rhoA * L0);
        if (c.DYN == 0) continue;   // 零污染：只输出静态锚点，无运行时键

        fem::Model3D m = makeModel(g, L0);
        fem::NewmarkOptions3D opt;
        opt.dt = c.DT; opt.tEnd = c.TEND;
        opt.beta = c.BETA; opt.gamma = c.GAMMA;
        opt.massMode = c.MASSMODE;
        opt.cMass = c.CMASS; opt.cStiff = c.CSTIFF;
        opt.loadType = c.LOADTYPE;
        opt.tolForce = c.TOL; opt.maxIter = c.MAXITER;
        opt.ctrlNode = 2; opt.ctrlDof = 2;
        opt.snapThreshold = c.SNAP;
        opt.strainMeasure = c.STRAIN;

        fem::NewmarkResult3D r;
        bool threw = false;
        try {
            r = fem::solveNewmarkDynamic3D(m, opt);
        } catch (const std::exception& e) {
            threw = true;
            std::fprintf(stderr, "case %s threw: %s\n", c.id.c_str(), e.what());
        }
        S(c.id, "exception", threw ? 1.0 : 0.0);
        if (threw) {
            S(c.id, "converged", 0.0);
            continue;
        }

        double vMax = 0.0, tVmax = 0.0, lamMax = 0.0, kinMax = 0.0;
        for (const auto& p : r.points) {
            const double vv = -p.ctrlDisp;   // 向下为正
            out << "T," << c.id << "," << p.step << "," << p.t << ","
                << p.lambda << "," << vv << "," << p.ctrlVel << ","
                << p.ctrlAcc << "," << p.residual << "," << p.iterations
                << "," << p.kinetic << "," << p.strain << ","
                << p.externalWork << "\n";
            if (vv > vMax) { vMax = vv; tVmax = p.t; }
            lamMax = std::max(lamMax, p.lambda);
            kinMax = std::max(kinMax, p.kinetic);
        }
        S(c.id, "DT", c.DT); S(c.id, "TEND", c.TEND);
        S(c.id, "BETA", c.BETA); S(c.id, "GAMMA", c.GAMMA);
        S(c.id, "MASSMODE", static_cast<double>(c.MASSMODE));
        S(c.id, "CMASS", c.CMASS); S(c.id, "CSTIFF", c.CSTIFF);
        S(c.id, "LOADTYPE", static_cast<double>(c.LOADTYPE));
        S(c.id, "SNAP", c.SNAP); S(c.id, "TOL", c.TOL);
        S(c.id, "MAXITER", static_cast<double>(c.MAXITER));
        S(c.id, "STRAIN", static_cast<double>(c.STRAIN));
        S(c.id, "npoints", static_cast<double>(r.points.size()));
        S(c.id, "steps", static_cast<double>(r.steps));
        S(c.id, "converged", r.converged ? 1.0 : 0.0);
        S(c.id, "error_empty", r.error.empty() ? 1.0 : 0.0);
        S(c.id, "mass_total", r.massTotal);
        S(c.id, "v_max", vMax); S(c.id, "t_at_v_max", tVmax);
        S(c.id, "lambda_max", lamMax);
        S(c.id, "v_final", r.points.empty() ? 0.0 : -r.points.back().ctrlDisp);
        S(c.id, "vel_final", r.points.empty() ? 0.0 : r.points.back().ctrlVel);
        S(c.id, "acc_final", r.points.empty() ? 0.0 : r.points.back().ctrlAcc);
        S(c.id, "kinetic_max", kinMax);
        S(c.id, "strain_final", r.points.empty() ? 0.0 : r.points.back().strain);
        S(c.id, "work_final", r.points.empty() ? 0.0 : r.points.back().externalWork);
        S(c.id, "snap_time", r.snapTime);
        S(c.id, "snap_disp", r.snapDisp);
        S(c.id, "snap_lambda", r.snapLambda);
    }
    return 0;
}
