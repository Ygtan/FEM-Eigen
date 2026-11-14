// frame_probe.cpp — CR 梁刚接框架后屈曲 / snap-through 球面弧长法探针
//
// 与 snap_probe.cpp 同构的文件驱动 CLI，但把单元从「两杆 TRUSS 桁架」
// 升级为「两根多单元 BEAM 刚接框架」（Williams toggle frame 型）：
//
//   左支座 (-B,0) -- NELEM 个 BEAM -- 顶点 (0,H) -- NELEM 个 BEAM -- 右支座 (B,0)
//
//   端部铰支：UX/UY/UZ/ROTX/ROTY 约束，释放 ROTZ
//   顶点与中间节点：UZ/ROTX/ROTY 约束（保持 XY 平面内），UX/UY/ROTZ 自由
//   参考荷载：顶点 -Y，幅值 PREF；控制自由度：顶点 UY，终止 |v| >= VMAX
//
// 力学完全复用 fem_3d_nl.hpp 的共转 CR 梁（随动三元组 + 解析一致切线
// crBeamTangent12）与 stability_path_3d.hpp 的球面弧长法
// （lambda 为未知量 + Crisfield 约束）。本文件不重写任何力学，
// 也不改 fem_3d_nl.hpp / stability_path_3d.hpp。
//
// 截面刻意取 IY = IZ 的对称截面：beamRefTriad 用 global Z 作参考 ez，
// 而 ANSYS BEAM188 无 K 节点时的缺省局部轴未必同向；IY=IZ 时任何绕
// 杆轴的局部轴旋转都不改变刚度，从根上消除两套求解器的局部轴口径差。
// 平面内弯曲绕局部 z，G 只影响扭转（ROTX/ROTY 已约束），不参与本问题。
//
// 输入 *.frame（纯 ASCII，与 .snap / .beam 同构）：
//   G <gid> E=.. G=.. A=.. IY=.. IZ=.. J=.. B=.. H=.. NELEM=.. PREF=..
//          TOL=.. MAXSTEPS=.. MAXITER=.. VMAX=.. DSMAX=.. DSMIN=..
//   C <cid> <gid> DS=.. PSI=.. TAG=..
//
// 输出 CSV（无表头，首列记录类型）：
//   P,<cid>,<step>,<lambda>,<v>,<residual>,<iterations>,<ds>,<pivot>
//   S,<cid>,<key>,<value>
//   E,<cid>,<error_code>                （仅失败时）
//
// 用法：frame_probe.exe <input.frame> <output.csv>
#include "stability_path_3d.hpp"

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
    double E = 2.1e11;
    double G = 8.1e13;       // 大 G：ANSYS 侧用于退化 Timoshenko 剪切柔度
    double A = 1.0e-3;
    double IY = 1.0e-7;
    double IZ = 1.0e-7;
    double J = 2.0e-7;
    double B = 2.0;
    double H = 0.5;
    double PREF = 1.0e6;
    double VMAX = 1.0;
    double TOL = 1.0e-9;
    double DSMAX = 0.02;
    double DSMIN = 1.0e-8;
    int NELEM = 8;
    int MAXSTEPS = 2000;
    int MAXITER = 40;
};

struct Case {
    std::string id, gid, tag;
    double DS = 0.005;
    double PSI = 1.0;
};

std::map<std::string, double> parseKV(const std::vector<std::string>& tok, int from) {
    std::map<std::string, double> kv;
    for (std::size_t i = from; i < tok.size(); ++i) {
        const auto pos = tok[i].find('=');
        if (pos == std::string::npos) continue;
        try {
            kv[tok[i].substr(0, pos)] = std::stod(tok[i].substr(pos + 1));
        } catch (const std::exception&) {
            // 非数值键（TAG=...）只作标签
        }
    }
    return kv;
}

double get(const std::map<std::string, double>& kv, const std::string& k, double d) {
    auto it = kv.find(k);
    return it == kv.end() ? d : it->second;
}

std::vector<std::string> split(const std::string& s) {
    std::istringstream is(s);
    std::vector<std::string> t;
    std::string x;
    while (is >> x) t.push_back(x);
    return t;
}

// 构造 Williams toggle 型刚接框架；返回顶点节点号。
int makeFrame(const Group& g, fem::Model3D& m) {
    const int n = g.NELEM;
    const int apex = n + 1;
    const int right = 2 * n + 1;
    const fem::Vec3 pl(-g.B, 0.0, 0.0);
    const fem::Vec3 pa(0.0, g.H, 0.0);
    const fem::Vec3 pr(g.B, 0.0, 0.0);

    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        m.nodes[1 + i] = pl + t * (pa - pl);
    }
    for (int i = 1; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        m.nodes[apex + i] = pa + t * (pr - pa);
    }

    int eid = 0;
    for (int i = 1; i <= 2 * n; ++i) {
        fem::Element3D e;
        e.id = ++eid;
        e.type = "BEAM";
        e.E = g.E;
        e.G = g.G;
        e.A = g.A;
        e.Iy = g.IY;
        e.Iz = g.IZ;
        e.J = g.J;
        e.n1 = i;
        e.n2 = i + 1;
        m.elements.push_back(e);
    }

    // 全部节点保持 XY 平面内（约束 UZ / ROTX / ROTY）
    for (int nd = 1; nd <= right; ++nd) {
        for (int d : {3, 4, 5}) m.constraints.push_back({nd, d, 0.0, 0.0});
    }
    // 两端铰支：约束平动，释放绕 Z 的转角
    for (int nd : {1, right}) {
        m.constraints.push_back({nd, 1, 0.0, 0.0});
        m.constraints.push_back({nd, 2, 0.0, 0.0});
    }
    m.loads.emplace_back(apex, 2, -g.PREF);
    return apex;
}

// 三点二次 Lagrange 插值（与 BridgeStudio/ansys_bridge_snap.py::interp_lambda 同口径）
double interpLambda(const std::vector<std::pair<double, double>>& curve, double vq) {
    if (curve.empty()) return 0.0;
    if (curve.size() < 3) {
        for (std::size_t i = 0; i + 1 < curve.size(); ++i) {
            const double a = curve[i].first, b = curve[i + 1].first;
            if (a <= vq && vq <= b && b > a) {
                const double t = (vq - a) / (b - a);
                return curve[i].second + t * (curve[i + 1].second - curve[i].second);
            }
        }
        return 0.0;
    }
    std::size_t i = 0;
    for (std::size_t k = 1; k < curve.size(); ++k)
        if (std::abs(curve[k].first - vq) < std::abs(curve[i].first - vq)) i = k;
    const std::size_t j = std::min(std::max<std::size_t>(i, 1) - 1, curve.size() - 3);
    const double xs[3] = {curve[j].first, curve[j + 1].first, curve[j + 2].first};
    const double ys[3] = {curve[j].second, curve[j + 1].second, curve[j + 2].second};
    for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 3; ++b)
            if (std::abs(xs[a] - xs[b]) < 1e-14) return 0.0;
    double s = 0.0;
    for (int a = 0; a < 3; ++a) {
        double t = ys[a];
        for (int b = 0; b < 3; ++b)
            if (a != b) t *= (vq - xs[b]) / (xs[a] - xs[b]);
        s += t;
    }
    return s;
}

// 相邻三点二次插值求离散路径极值，消除弧长步采样误差
void quadExtremum(const std::vector<std::pair<double, double>>& curve, int i,
                  double& vStar, double& lStar) {
    const int n = static_cast<int>(curve.size());
    if (i <= 0 || i + 1 >= n) {
        vStar = curve[i].first; lStar = curve[i].second; return;
    }
    const double x0 = curve[i - 1].first, x1 = curve[i].first, x2 = curve[i + 1].first;
    const double y0 = curve[i - 1].second, y1 = curve[i].second, y2 = curve[i + 1].second;
    if (std::abs(x1 - x0) < 1e-14 || std::abs(x2 - x1) < 1e-14) {
        vStar = x1; lStar = y1; return;
    }
    const double s1 = (y1 - y0) / (x1 - x0);
    const double s2 = (y2 - y1) / (x2 - x1);
    const double d2 = (s2 - s1) / (x2 - x0);
    if (std::abs(d2) < 1e-14) { vStar = x1; lStar = y1; return; }
    vStar = 0.5 * (x0 + x1) - s1 / (2.0 * d2);
    lStar = y0 + s1 * (vStar - x0) + d2 * (vStar - x0) * (vStar - x1);
}

// FEM 内核在参考构形上的**精确线性切线刚度**：单位增量下解 K0*du = Fref，
// 由顶点柔度反算 k0 = PREF / v。这不是独立锚点（它复用内核的单元切线），
// 只作为 FEM 侧的基准量，供 benchmark 与 BridgeStudio 侧独立 Python
// 线性刚架刚度矩阵严格对拍。
double linearSlope(const fem::Model3D& m, int apex, double pref) {
    const fem::DofMap3D dm = fem::buildDofMap3D(m);
    const int nd = dm.nDofs;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nd, nd);
    Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
    Eigen::Matrix<double, 14, 1> fe;
    Eigen::Matrix<double, 14, 14> Ke;
    for (const auto& e : m.elements) {
        const fem::CRElemData ed = fem::crMakeElemData(m, e, dm);
        fem::crElementInternal(ed, ue, fe, true, &Ke);
        for (int r = 0; r < ed.n; ++r)
            for (int c = 0; c < ed.n; ++c)
                K(ed.gdofs[r], ed.gdofs[c]) += Ke(r, c);
    }
    Eigen::VectorXd F = Eigen::VectorXd::Zero(nd);
    for (const auto& ld : m.loads) {
        const int d = dm.safeDof(std::get<0>(ld), std::get<1>(ld) - 1);
        if (d >= 0) F(d) += std::get<2>(ld);
    }
    std::vector<int> cdofs;
    for (const auto& c : m.constraints) {
        if (c.k > 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0 && std::find(cdofs.begin(), cdofs.end(), d) == cdofs.end())
            cdofs.push_back(d);
    }
    for (int d : cdofs) {
        K.row(d).setZero();
        K.col(d).setZero();
        K(d, d) = 1.0;
        F(d) = 0.0;
    }
    const Eigen::VectorXd du = Eigen::PartialPivLU<Eigen::MatrixXd>(K).solve(F);
    const int ctrl = dm.safeDof(apex, 1);
    if (ctrl < 0) return 0.0;
    const double v = -du(ctrl);       // Fref = -PREF，故 v > 0
    return (v != 0.0) ? pref / v : 0.0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: frame_probe <input.frame> <output.csv>\n");
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
            g.E = get(kv, "E", g.E);         g.G = get(kv, "G", g.G);
            g.A = get(kv, "A", g.A);         g.IY = get(kv, "IY", g.IY);
            g.IZ = get(kv, "IZ", g.IZ);      g.J = get(kv, "J", g.J);
            g.B = get(kv, "B", g.B);         g.H = get(kv, "H", g.H);
            g.PREF = get(kv, "PREF", g.PREF);
            g.VMAX = get(kv, "VMAX", g.VMAX);
            g.TOL = get(kv, "TOL", g.TOL);
            g.DSMAX = get(kv, "DSMAX", g.DSMAX);
            g.DSMIN = get(kv, "DSMIN", g.DSMIN);
            g.NELEM = static_cast<int>(get(kv, "NELEM", g.NELEM));
            g.MAXSTEPS = static_cast<int>(get(kv, "MAXSTEPS", g.MAXSTEPS));
            g.MAXITER = static_cast<int>(get(kv, "MAXITER", g.MAXITER));
            groups[t[1]] = g;
        } else if (t[0] == "C" && t.size() >= 3) {
            Case c;
            c.id = t[1];
            c.gid = t[2];
            const auto kv = parseKV(t, 3);
            c.DS = get(kv, "DS", c.DS);
            c.PSI = get(kv, "PSI", c.PSI);
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

    for (const auto& c : cases) {
        auto gi = groups.find(c.gid);
        if (gi == groups.end()) {
            std::fprintf(stderr, "unknown group %s\n", c.gid.c_str());
            return 2;
        }
        const Group& g = gi->second;
        fem::Model3D m;
        const int apex = makeFrame(g, m);

        fem::ArcLengthOptions3D opt;
        opt.ds0 = c.DS;
        opt.dsMin = g.DSMIN;
        opt.dsMax = g.DSMAX;
        opt.psi = c.PSI;
        opt.tolForce = g.TOL;
        opt.maxSteps = g.MAXSTEPS;
        opt.maxIter = g.MAXITER;
        opt.ctrlNode = apex;
        opt.ctrlDof = 2;
        opt.ctrlTarget = g.VMAX;

        fem::ArcLengthResult3D r;
        try {
            r = fem::solveArcLengthPath3D(m, opt);
        } catch (const std::exception& e) {
            out << "E," << c.id << ",solver_exception\n";
            out << "S," << c.id << ",converged,0\n";
            std::fprintf(stderr, "case %s threw: %s\n", c.id.c_str(), e.what());
            continue;
        }
        if (!r.error.empty()) out << "E," << c.id << "," << r.error << "\n";

        double worstRes = 0.0;
        std::vector<std::pair<double, double>> curve;
        std::vector<double> lams;
        for (const auto& p : r.points) {
            const double v = -p.ctrlDisp;
            out << "P," << c.id << "," << p.step << ","
                << p.lambda << "," << v << "," << p.residual << ","
                << p.iterations << "," << p.ds << "," << p.pivotRatio << "\n";
            curve.emplace_back(v, p.lambda);
            lams.push_back(p.lambda);
            worstRes = std::max(worstRes, p.residual);
        }

        // 极限点 = **第一个** lambda 局部极大。刚接框架翻越后弯曲刚度
        // 强烈再硬化，全局最大通常落在终点，不能当极限荷载。
        int iPeak = 0, iMin = 0, iGMax = 0;
        double lamPeak = 0.0, lamMin = 0.0, lamGMax = 0.0;
        double k0Secant = 0.0;
        if (!curve.empty()) {
            lamPeak = lamMin = lamGMax = lams[0];
            for (int i = 0; i < static_cast<int>(lams.size()); ++i) {
                if (lams[i] > lamGMax) { lamGMax = lams[i]; iGMax = i; }
                if (lams[i] < lamMin) { lamMin = lams[i]; iMin = i; }
            }
            bool found = false;
            for (int i = 1; i + 1 < static_cast<int>(lams.size()); ++i) {
                if (lams[i] > lams[i - 1] && lams[i] >= lams[i + 1]) {
                    lamPeak = lams[i]; iPeak = i; found = true; break;
                }
            }
            if (!found) { lamPeak = lamGMax; iPeak = iGMax; }
            if (curve.size() >= 2 && curve[1].first > curve[0].first)
                k0Secant = g.PREF * (lams[1] - lams[0]) /
                           (curve[1].first - curve[0].first);
        }
        const double vAtPeak = curve.empty() ? 0.0 : curve[iPeak].first;
        double vPeakFit = vAtPeak, lamPeakFit = lamPeak;
        double vMinFit = curve.empty() ? 0.0 : curve[iMin].first;
        double lamMinFit = lamMin;
        if (!curve.empty()) {
            quadExtremum(curve, iPeak, vPeakFit, lamPeakFit);
            quadExtremum(curve, iMin, vMinFit, lamMinFit);
        }

        auto S = [&](const char* k, double v) {
            out << "S," << c.id << "," << k << "," << v << "\n";
        };
        S("E", g.E);              S("G", g.G);
        S("A", g.A);              S("IY", g.IY);
        S("IZ", g.IZ);            S("J", g.J);
        S("B", g.B);              S("H", g.H);
        S("NELEM", static_cast<double>(g.NELEM));
        S("L0", std::hypot(g.B, g.H));
        S("PREF", g.PREF);        S("VMAX", g.VMAX);
        S("DS0", c.DS);           S("PSI", c.PSI);
        S("TOL", g.TOL);          S("DSMAX", g.DSMAX);
        S("DSMIN", g.DSMIN);
        S("MAXSTEPS", static_cast<double>(g.MAXSTEPS));
        S("MAXITER", static_cast<double>(g.MAXITER));
        S("apex_node", static_cast<double>(apex));
        S("k0_fem", linearSlope(m, apex, g.PREF));
        S("k0_secant", k0Secant);
        S("npoints", static_cast<double>(r.points.size()));
        S("steps", static_cast<double>(r.steps));
        S("converged", r.converged ? 1.0 : 0.0);
        S("error_empty", r.error.empty() ? 1.0 : 0.0);
        S("lambda_final", r.lambda);
        S("v_final", r.points.empty() ? 0.0 : -r.points.back().ctrlDisp);
        S("lambda_max", lamPeak);
        S("lambda_global_max", lamGMax);
        S("v_at_lambda_max", vAtPeak);
        S("lambda_max_fit", lamPeakFit);
        S("v_at_lambda_max_fit", vPeakFit);
        S("lambda_min", lamMin);
        S("lambda_min_fit", lamMinFit);
        S("v_at_lambda_min_fit", vMinFit);
        S("lambda_at_01H", interpLambda(curve, 0.1 * g.H));
        S("lambda_at_02H", interpLambda(curve, 0.2 * g.H));
        S("lambda_at_H", interpLambda(curve, g.H));
        S("lambda_at_2H", interpLambda(curve, 2.0 * g.H));
        S("max_residual", worstRes);
        S("max_residual_rel", g.PREF > 0.0 ? worstRes / g.PREF : 0.0);
    }
    return 0;
}
