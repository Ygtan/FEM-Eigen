// snap_probe.cpp — 后屈曲 / snap-through 球面弧长法对标探针（文件驱动 CLI）
//
// 输入 *.snap（纯 ASCII，与 .sag / .beam 同构）：
//   G <gid> E=.. A=.. B=.. H=.. PREF=.. TOL=.. MAXSTEPS=.. MAXITER=.. VMAX=..
//   C <cid> <gid> DS=.. PSI=.. TAG=..
//
// 物理模型：von Mises 两杆浅桁架
//   节点 1(-B,0,0) — 2(0,H,0) — 3(B,0,0)，两端固定、顶点 Z 约束；
//   两根 TRUSS，参考荷载 Fref 施加在顶点 -Y 方向，幅值 PREF。
// 独立解析锚点（与 fem_3d_nl.hpp 轴向单元同本构 N=EA(L-L0)/L0）：
//   L(v)=sqrt(B^2+(H-v)^2),  P(v) = -2*N*(H-v)/L
// 极值点 dP/dv=0 由黄金分割求得，P(H)=0，P(2H)=0。
//
// 输出 CSV（无表头，首列记录类型）：
//   P,<cid>,<step>,<lambda>,<v>,<residual>,<iterations>,<ds>,<pivot>
//   S,<cid>,<key>,<value>
//
// 由 stability_path_3d.hpp 求解；不修改任何既有求解器或输入格式。
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
    double E = 2.1e11, A = 1.0e-3, B = 2.0, H = 0.5;
    double PREF = 2.0e6;
    double DSMAX = 0.2, DSMIN = 1.0e-7, TOL = 1.0e-7;
    double VMAX = 1.0;
    int MAXSTEPS = 400, MAXITER = 30;
};

struct Case {
    std::string id, gid;
    double DS = 0.02, PSI = 1.0;
    int STRAIN = 0;              // 0=工程应变（内核口径） 1=对数应变（ANSYS KEYOPT(2)=1）
    // 第二十一次力学①（opt-in，默认 0/0/0 → 既有 .snap 逐位不变）：
    int CTRLMODE = 0;            // 1=位移控制（控制 DOF 位移增量锁定）
    int LINESEARCH = 0;          // 1=Armijo 回溯线搜索
    double CTRLSTEP = 0.0;       // 每步控制位移增量（>0 生效，否则回退 DS）
    int CTRLNODE = 2;            // 控制节点（负例可指向固定端触发 singular_control_dof）
    int CTRLDOF = 2;             // 控制自由度 1~6
    // 第二十一次力学②（opt-in，默认 0/0/0 → 既有 .snap 逐位不变）：
    int IMPERF = 0;              // 1=施加线性屈曲模态初始缺陷
    int IMPMODE = 0;             // 屈曲模态序号（0=最低阶）
    double IMPAMP = 0.0;         // 缺陷幅值（长度，加到节点坐标）
    // 第二十一次力学④（opt-in，默认全 0 → 既有 .snap 逐位不变）：
    int BIFURC = 0;              // 1=启用分岔检测与分支切换
    double BIFPIV = 0.05;        // pivotRatio 触发阈值
    int BIFMODE = 1;             // 扰动模态序号（1=反对称）
    double BIFAMP = 0.0;         // 扰动幅值；0=只检测不切换（负例）
    std::string tag;
};

std::map<std::string, double> parseKV(const std::vector<std::string>& tok, int from) {
    std::map<std::string, double> kv;
    for (std::size_t i = from; i < tok.size(); ++i) {
        const auto pos = tok[i].find('=');
        if (pos == std::string::npos) continue;
        try {
            kv[tok[i].substr(0, pos)] = std::stod(tok[i].substr(pos + 1));
        } catch (const std::exception&) {
            // 非数值键（TAG=... 等）只作标签，不进参数表
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

double pOfV(double v, const Group& g, double L0, int strain) {
    const double L = std::hypot(g.B, g.H - v);
    const double N = (strain == 1)
        ? g.E * g.A * std::log(L / L0)
        : g.E * g.A * (L - L0) / L0;
    return -2.0 * N * (g.H - v) / L;
}

double analyticLimitV(const Group& g, double L0, int strain) {
    double lo = 0.0, hi = g.H;
    for (int i = 0; i < 200; ++i) {
        const double m1 = lo + (hi - lo) * 0.3819660112501051;
        const double m2 = lo + (hi - lo) * 0.6180339887498949;
        if (pOfV(m1, g, L0, strain) < pOfV(m2, g, L0, strain)) lo = m1;
        else hi = m2;
    }
    return 0.5 * (lo + hi);
}

fem::Model3D makeTruss(const Group& g, double L0) {
    fem::Model3D m;
    m.nodes[1] = fem::Vec3(-g.B, 0.0, 0.0);
    m.nodes[2] = fem::Vec3(0.0, g.H, 0.0);
    m.nodes[3] = fem::Vec3(g.B, 0.0, 0.0);
    for (int i = 0; i < 2; ++i) {
        fem::Element3D e;
        e.id = i + 1;
        e.type = "TRUSS";
        e.E = g.E;
        e.A = g.A;
        e.L0 = L0;
        e.n1 = i + 1;
        e.n2 = i + 2;
        m.elements.push_back(e);
    }
    for (int d = 1; d <= 3; ++d) {
        m.constraints.push_back({1, d, 0.0, 0.0});
        m.constraints.push_back({3, d, 0.0, 0.0});
    }
    m.constraints.push_back({2, 3, 0.0, 0.0});   // 顶点 Z 约束，保持平面内
    m.loads.emplace_back(2, 2, -g.PREF);
    return m;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: snap_probe <input.snap> <output.csv>\n");
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
            g.E = get(kv, "E", g.E);
            g.A = get(kv, "A", g.A);
            g.B = get(kv, "B", g.B);
            g.H = get(kv, "H", g.H);
            g.PREF = get(kv, "PREF", g.PREF);
            g.DSMAX = get(kv, "DSMAX", g.DSMAX);
            g.DSMIN = get(kv, "DSMIN", g.DSMIN);
            g.TOL = get(kv, "TOL", g.TOL);
            g.VMAX = get(kv, "VMAX", g.VMAX);
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
            c.STRAIN = static_cast<int>(get(kv, "STRAIN", c.STRAIN));
            c.CTRLMODE = static_cast<int>(get(kv, "CTRLMODE", c.CTRLMODE));
            c.LINESEARCH = static_cast<int>(get(kv, "LINESEARCH", c.LINESEARCH));
            c.CTRLSTEP = get(kv, "CTRLSTEP", c.CTRLSTEP);
            c.CTRLNODE = static_cast<int>(get(kv, "CTRLNODE", c.CTRLNODE));
            c.CTRLDOF = static_cast<int>(get(kv, "CTRLDOF", c.CTRLDOF));
            c.IMPERF = static_cast<int>(get(kv, "IMPERF", c.IMPERF));
            c.IMPMODE = static_cast<int>(get(kv, "IMPMODE", c.IMPMODE));
            c.IMPAMP = get(kv, "IMPAMP", c.IMPAMP);
            c.BIFURC = static_cast<int>(get(kv, "BIFURC", c.BIFURC));
            c.BIFPIV = get(kv, "BIFPIV", c.BIFPIV);
            c.BIFMODE = static_cast<int>(get(kv, "BIFMODE", c.BIFMODE));
            c.BIFAMP = get(kv, "BIFAMP", c.BIFAMP);
            for (const auto& x : t) {
                if (x.rfind("TAG=", 0) == 0) c.tag = x.substr(4);
            }
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
        const double L0 = std::hypot(g.B, g.H);
        fem::Model3D m = makeTruss(g, L0);

        // 第二十一次力学②：opt-in 初始缺陷（线性屈曲模态按幅值加到节点坐标）。
        // IMPERF=0 时整段跳过，模型 / 路径 / 输出与既有 .snap 逐位一致。
        double impLambda = 0.0, impDx2 = 0.0, impDy2 = 0.0, impMaxPhi = 0.0;
        int impOk = 0;
        if (c.IMPERF == 1) {
            const fem::BucklingModeShape3D bk =
                fem::computeBucklingModeShape3D(m, c.IMPMODE);
            if (bk.error.empty()) {
                const fem::DofMap3D dm = fem::buildDofMap3D(m);
                const fem::Vec3 p2before = m.nodes.at(2);
                // 无应力缺陷：同步更新 TRUSS 的 L0，缺陷构形在 λ=0 平衡。
                if (fem::applyImperfectionToNodes(m, dm, bk.phi, c.IMPAMP, true)) {
                    impOk = 1;
                    impLambda = bk.lambda;
                    impDx2 = m.nodes.at(2)(0) - p2before(0);
                    impDy2 = m.nodes.at(2)(1) - p2before(1);
                    for (double x : bk.phi) impMaxPhi = std::max(impMaxPhi, std::abs(x));
                }
            }
        }

        fem::ArcLengthOptions3D opt;
        opt.ds0 = c.DS;
        opt.dsMin = g.DSMIN;
        opt.dsMax = g.DSMAX;
        opt.psi = c.PSI;
        opt.tolForce = g.TOL;
        opt.maxSteps = g.MAXSTEPS;
        opt.maxIter = g.MAXITER;
        opt.ctrlNode = 2;
        opt.ctrlDof = 2;
        opt.ctrlTarget = g.VMAX;
        opt.strainMeasure = c.STRAIN;
        opt.controlMode = c.CTRLMODE;
        opt.lineSearch = c.LINESEARCH;
        opt.ctrlDispStep = c.CTRLSTEP;
        opt.ctrlNode = c.CTRLNODE;
        opt.ctrlDof = c.CTRLDOF;
        // 第二十一次力学④：opt-in 分岔检测与分支切换。
        opt.bifurcate = c.BIFURC;
        opt.bifPivot = c.BIFPIV;
        opt.bifMode = c.BIFMODE;
        opt.bifAmp = c.BIFAMP;

        fem::ArcLengthResult3D r;
        try {
            r = fem::solveArcLengthPath3D(m, opt);
        } catch (const std::exception& e) {
            out << "S," << c.id << ",exception," << 1.0 << "\n";
            out << "S," << c.id << ",converged,0\n";
            std::fprintf(stderr, "case %s threw: %s\n", c.id.c_str(), e.what());
            continue;
        }

        const double vcr = analyticLimitV(g, L0, c.STRAIN);
        const double pcr = pOfV(vcr, g, L0, c.STRAIN);
        double lamMax = 0.0, lamMin = 0.0, vAtLamMax = 0.0, worst = 0.0;
        int iMax = 0, iMin = 0;
        std::vector<std::pair<double, double>> curve;   // (v, lambda)
        for (const auto& p : r.points) {
            const double v = -p.ctrlDisp;
            out << "P," << c.id << "," << p.step << ","
                << p.lambda << "," << v << "," << p.residual << ","
                << p.iterations << "," << p.ds << "," << p.pivotRatio << "\n";
            // 力学④：分支路径的模态投影逐点输出（仅 opt-in，默认无 M 行）。
            if (c.BIFURC == 1)
                out << "M," << c.id << "," << p.step << "," << p.bifModalAmp << "\n";
            curve.emplace_back(v, p.lambda);
            const int idx = static_cast<int>(curve.size()) - 1;
            if (p.lambda > lamMax) { lamMax = p.lambda; vAtLamMax = v; iMax = idx; }
            if (p.lambda < lamMin) { lamMin = p.lambda; iMin = idx; }
            if (v >= 0.0 && v <= 2.0 * g.H)
                worst = std::max(worst, std::abs(p.lambda * g.PREF
                                                 - pOfV(v, g, L0, c.STRAIN)));
        }

        // 离散平衡点未必恰好落在极值点上：用相邻三点二次插值求极值，
        // 消除弧长步采样误差（路径上 v 单调递增，插值良态）。
        auto quadExtremum = [&](int i, double& vStar, double& lStar) {
            const int n = static_cast<int>(curve.size());
            if (i <= 0 || i + 1 >= n) {
                vStar = curve[i].first; lStar = curve[i].second; return;
            }
            const double x0 = curve[i - 1].first, x1 = curve[i].first, x2 = curve[i + 1].first;
            const double y0 = curve[i - 1].second, y1 = curve[i].second, y2 = curve[i + 1].second;
            const double s1 = (y1 - y0) / (x1 - x0);
            const double s2 = (y2 - y1) / (x2 - x1);
            const double d2 = (s2 - s1) / (x2 - x0);
            if (std::abs(d2) < 1.0e-14 || std::abs(x1 - x0) < 1e-14) {
                vStar = x1; lStar = y1; return;
            }
            vStar = 0.5 * (x0 + x1) - s1 / (2.0 * d2);
            lStar = y0 + s1 * (vStar - x0) + d2 * (vStar - x0) * (vStar - x1);
        };
        auto interpLambda = [&](double vq) {
            for (std::size_t i = 0; i + 1 < curve.size(); ++i) {
                const double a = curve[i].first, b = curve[i + 1].first;
                if (a <= vq && vq <= b && b > a) {
                    const double t = (vq - a) / (b - a);
                    return curve[i].second + t * (curve[i + 1].second - curve[i].second);
                }
            }
            return 0.0;
        };
        double vMaxFit = vAtLamMax, lamMaxFit = lamMax;
        double vMinFit = 0.0, lamMinFit = lamMin;
        if (!curve.empty()) {
            quadExtremum(iMax, vMaxFit, lamMaxFit);
            quadExtremum(iMin, vMinFit, lamMinFit);
        }
        const double lamAtH = interpLambda(g.H);
        const double lamAt2H = interpLambda(2.0 * g.H);

        auto S = [&](const char* k, double v) { out << "S," << c.id << "," << k << "," << v << "\n"; };
        S("E", g.E);
        S("A", g.A);
        S("B", g.B);
        S("H", g.H);
        S("L0", L0);
        S("PREF", g.PREF);
        S("DS0", c.DS);
        S("PSI", c.PSI);
        S("STRAIN", static_cast<double>(c.STRAIN));
        S("TOL", g.TOL);
        S("VMAX", g.VMAX);
        S("analytic_v_cr", vcr);
        S("analytic_P_cr", pcr);
        S("analytic_lambda_cr", pcr / g.PREF);
        S("analytic_P_at_H", pOfV(g.H, g, L0, c.STRAIN));
        S("analytic_P_at_2H", pOfV(2.0 * g.H, g, L0, c.STRAIN));
        S("analytic_P_cr_eng", pOfV(analyticLimitV(g, L0, 0), g, L0, 0));
        S("analytic_P_cr_log", pOfV(analyticLimitV(g, L0, 1), g, L0, 1));
        S("npoints", static_cast<double>(r.points.size()));
        S("steps", static_cast<double>(r.steps));
        S("converged", r.converged ? 1.0 : 0.0);
        S("error_empty", r.error.empty() ? 1.0 : 0.0);
        S("lambda_final", r.lambda);
        S("v_final", r.points.empty() ? 0.0 : -r.points.back().ctrlDisp);
        S("lambda_max", lamMax);
        S("lambda_min", lamMin);
        S("v_at_lambda_max", vAtLamMax);
        S("lambda_max_fit", lamMaxFit);
        S("v_at_lambda_max_fit", vMaxFit);
        S("lambda_min_fit", lamMinFit);
        S("v_at_lambda_min_fit", vMinFit);
        S("lambda_at_H", lamAtH);
        S("lambda_at_2H", lamAt2H);
        S("worst_dP", worst);
        S("worst_dP_rel", pcr > 0.0 ? worst / pcr : 0.0);
        // 第二十一次力学①：位移控制 / 线搜索诊断，只在 opt-in 时输出，
        // CTRLMODE=0 且 LINESEARCH=0 的既有 .snap 输出逐位不变。
        if (c.CTRLMODE == 1 || c.LINESEARCH == 1) {
            int totalCuts = 0, maxCuts = 0, stepsCut = 0;
            double minAlpha = 1.0;
            for (const auto& p : r.points) {
                totalCuts += p.lineCuts;
                maxCuts = std::max(maxCuts, p.lineCuts);
                if (p.lineCuts > 0) ++stepsCut;
                minAlpha = std::min(minAlpha, p.alpha);
            }
            S("ctrlmode", static_cast<double>(c.CTRLMODE));
            S("linesearch", static_cast<double>(c.LINESEARCH));
            S("ctrlstep", c.CTRLSTEP > 0.0 ? c.CTRLSTEP : c.DS);
            S("ls_total_cuts", static_cast<double>(totalCuts));
            S("ls_max_cuts", static_cast<double>(maxCuts));
            S("ls_steps_cut", static_cast<double>(stepsCut));
            S("ls_min_alpha", minAlpha);
        }
        // 第二十一次力学②：初始缺陷诊断，只在 IMPERF=1 时输出，
        // IMPERF=0 的既有 .snap 输出逐位不变。
        if (c.IMPERF == 1) {
            S("imperf", static_cast<double>(c.IMPERF));
            S("impmode", static_cast<double>(c.IMPMODE));
            S("impamp", c.IMPAMP);
            S("imp_ok", static_cast<double>(impOk));
            S("imp_lambda", impLambda);          // 线性屈曲因子（与 λ_max 不同源）
            S("imp_dx2", impDx2);                // 缺陷引起的顶点 X 坐标变化
            S("imp_dy2", impDy2);                // 缺陷引起的顶点 Y 坐标变化
            S("imp_maxphi", impMaxPhi);          // 模态最大分量（归一化 =1）
        }
        // 第二十一次力学④：分岔检测 / 分支切换诊断，只在 BIFURC=1 时输出，
        // BIFURC=0 的既有 .snap 输出逐位不变。
        if (c.BIFURC == 1) {
            S("bifurc", 1.0);
            S("bifpiv", c.BIFPIV);
            S("bifmode", static_cast<double>(c.BIFMODE));
            S("bifamp", c.BIFAMP);
            S("bif_detected", static_cast<double>(r.bifDetected));
            S("bif_switched", static_cast<double>(r.bifSwitched));
            S("bif_step", static_cast<double>(r.bifStep));
            S("bif_lambda", r.bifLambda);
            S("bif_ctrl_disp", r.bifCtrlDisp);
            S("bif_eig_ratio", r.bifMinEigRatio);
            S("bif_pivot", r.bifPivot);
            double ampBefore = 0.0, ampAfter = 0.0, ampAtSwitch = 0.0;
            for (const auto& p : r.points) {
                if (p.step < r.bifStep) ampBefore = std::max(ampBefore, std::abs(p.bifModalAmp));
                else if (p.step > r.bifStep) ampAfter = std::max(ampAfter, std::abs(p.bifModalAmp));
                else ampAtSwitch = p.bifModalAmp;
            }
            S("bif_max_amp_before", ampBefore);
            S("bif_amp_at_switch", ampAtSwitch);
            S("bif_max_amp_after", ampAfter);
            S("bif_amp_final", r.points.empty() ? 0.0 : r.points.back().bifModalAmp);
        }
    }
    return 0;
}
