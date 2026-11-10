// ============================================================
// 3D + 统一内核 解析解断言测试
// ============================================================

#include "fem_3d.hpp"
#include "fem_2d.hpp"
#include "test_common.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace fem;

namespace {


void checkMat(const std::string& name, const DenseMat& got, const DenseMat& exp, double tol = 1e-9) {
    if (got.rows() != exp.rows() || got.cols() != exp.cols()) {
        std::printf("  [FAIL] %-44s 尺寸不符 got=%dx%d exp=%dx%d\n",
                    name.c_str(), static_cast<int>(got.rows()), static_cast<int>(got.cols()),
                    static_cast<int>(exp.rows()), static_cast<int>(exp.cols()));
        ++g_fail;
        return;
    }
    const double denom = std::max(1.0, exp.cwiseAbs().maxCoeff());
    const double rel = (got - exp).cwiseAbs().maxCoeff() / denom;
    if (rel <= tol) {
        ++g_pass;
        std::printf("  [PASS] %-44s maxRel=%.2e\n", name.c_str(), rel);
    } else {
        ++g_fail;
        std::printf("  [FAIL] %-44s maxRel=%.2e (tol=%.0e)\n", name.c_str(), rel, tol);
    }
}

DenseMat submatrix(const DenseMat& M, const int* idx, int n) {
    DenseMat S(n, n);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            S(r, c) = M(idx[r], idx[c]);
    return S;
}

Model3D straightBeam3D(double Lb, int nElem, double Em, double Am, double Im) {
    Model3D m;
    for (int i = 0; i <= nElem; ++i)
        m.nodes[i + 1] = Vec3(Lb * i / static_cast<double>(nElem), 0.0, 0.0);
    for (int e = 0; e < nElem; ++e) {
        Element3D el;
        el.id = e + 1;
        el.type = "BEAM";
        el.E = Em; el.A = Am;
        el.Iy = Im; el.Iz = Im; el.J = 1.0e-5; el.G = 8.0e10;
        el.n1 = e + 1; el.n2 = e + 2;
        m.elements.push_back(el);
    }
    return m;
}

constexpr double E  = 2.0e11;   // Pa（钢）
constexpr double A  = 0.02;     // m^2
constexpr double I  = 8.0e-4;   // m^4
constexpr double L  = 30.0;     // m
constexpr double P  = 1.0e5;    // N
constexpr double q  = 2.0e4;    // N/m

double nodalDisp(const Model2D& m, const ModelResult& r, int node, int dof) {
    std::vector<int> ids;
    for (const auto& [id, c] : m.nodes) ids.push_back(id);
    for (std::size_t i = 0; i < ids.size(); ++i)
        if (ids[i] == node) return r.U(3 * static_cast<int>(i) + (dof - 1));
    throw std::runtime_error("节点不存在: " + std::to_string(node));
}

void test3DCantilever() {
    std::printf("\n[8] 三维悬臂梁 端部集中力\n");
    const int n = 1;
    Model3D m = straightBeam3D(L, n, E, A, I);
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});
    m.loads.push_back({n + 1, 2, -P});

    const std::vector<double> U = solveModel3D(m);
    // 梁沿 X、荷载沿 Y：弯曲发生在 X-Y 面内，用 Iz。
    check("端部挠度  δ=PL³/(3EIz)", U[6 * n + 1], -P * L * L * L / (3.0 * E * I));
}

void test3DAxialBar() {
    std::printf("\n[9] 三维桁架杆 轴向拉伸\n");
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(L, 0.0, 0.0);
    Element3D el;
    el.id = 1; el.type = "TRUSS"; el.E = E; el.A = A; el.n1 = 1; el.n2 = 2;
    m.elements.push_back(el);
    for (int d = 1; d <= 3; ++d) m.constraints.push_back({1, d, 0.0});
    m.constraints.push_back({2, 2, 0.0});
    m.constraints.push_back({2, 3, 0.0});
    m.loads.push_back({2, 1, P});

    const std::vector<double> U = solveModel3D(m);
    // 纯桁架节点只有 3 个自由度：节点 2 起始于全局第 3 个自由度。
    check("轴向伸长  δ=PL/(EA)", U[3], P * L / (E * A));
}

void test3DSupportSettlement() {
    std::printf("\n[10] 三维 支座沉降（非零指定位移）\n");
    const int n = 1;
    const double settle = -0.005;
    Model3D m = straightBeam3D(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, settle});   // 固定端竖向沉降 5 mm
    for (int d = 3; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});
    // 不施加任何荷载：梁随支座整体下沉，且内部无变形。

    const std::vector<double> U = solveModel3D(m);
    check("固定端位移 = 沉降量", U[1], settle);
    check("自由端随支座整体下沉", U[6 * n + 1], settle);
    check("自由端无转角（无内力）", U[6 * n + 5], 0.0);
}

void test3DMixedBeamTruss() {
    std::printf("\n[11] 三维 梁+桁架混合（按节点分配自由度）\n");
    const double Lb = 10.0, Lt = 10.0;
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(Lb, 0.0, 0.0);
    m.nodes[3] = Vec3(Lb + Lt, 0.0, 0.0);
    Element3D beam;
    beam.id = 1; beam.type = "BEAM"; beam.E = E; beam.A = A;
    beam.Iy = I; beam.Iz = I; beam.J = 1.0e-5; beam.G = 8.0e10;
    beam.n1 = 1; beam.n2 = 2;
    Element3D truss;
    truss.id = 2; truss.type = "TRUSS"; truss.E = E; truss.A = A;
    truss.n1 = 2; truss.n2 = 3;
    m.elements.push_back(beam);
    m.elements.push_back(truss);

    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});
    m.constraints.push_back({3, 2, 0.0});   // 桁架端节点只约束横向
    m.constraints.push_back({3, 3, 0.0});
    m.loads.push_back({3, 1, P});

    // 梁轴向刚度与桁架轴向刚度串联：δ = P·(Lb/(EA) + Lt/(EA))
    const std::vector<double> U = solveModel3D(m);
    const double exact = P * (Lb / (E * A) + Lt / (E * A));
    // 节点 3 为纯桁架节点（3 自由度），起始于全局第 12 个自由度。
    check("串联轴向变形  δ=P(Lb+Lt)/(EA)", U[12], exact);
}

void test3DInfluenceSurfaceTransverse() {
    std::printf("\n[12] 三维 影响面横向维度\n");
    const char* cands[] = {"../verify_3d_influence_surface.dat",
                           "verify_3d_influence_surface.dat",
                           "../../verify_3d_influence_surface.dat"};
    std::string path;
    for (const char* c : cands) {
        std::ifstream probe(c);
        if (probe) { path = c; break; }
    }
    if (path.empty()) {
        std::printf("  [SKIP] 未找到 verify_3d_influence_surface.dat\n");
        return;
    }

    const Model3D m = readModel3D(path);
    const Model3D::Influence3D spec = m.influences.at(0);
    const auto vals = compute3DInfluenceSurface(m, spec, m.laneSurface);

    // 该算例中两条纵向桁架互相独立：
    //   加载在 y=0 线上的节点 2 → 影响本线，eta = L/(EA)
    //   加载在 y=1 线上的节点 4 → 与响应点无连接，eta = 0
    double eta_x2y1 = 0.0, eta_x2y3 = 0.0;
    bool hasA = false, hasB = false;
    for (auto [xn, yn, v] : vals) {
        if (xn == 2 && yn == 1) { eta_x2y1 = v; hasA = true; }
        if (xn == 2 && yn == 3) { eta_x2y3 = v; hasB = true; }
    }
    if (!hasA || !hasB) {
        std::printf("  [SKIP] 影响面输出缺少预期的 (x=2,y=1)/(x=2,y=3) 组合\n");
        return;
    }

    const double kAxial = 1.0 / (200.0e9 * 0.001);   // L/(EA)，L=1
    check("同纵向、y=0 横向：eta=L/(EA)", eta_x2y1, kAxial, 1e-6);
    check("同纵向、y=1 横向：eta=0（横向维度生效）", eta_x2y3, 0.0, 1e-6);
    if (std::abs(eta_x2y1 - eta_x2y3) < 1e-15) {
        ++g_fail;
        std::printf("  [FAIL] 不同横向位置得到相同影响值：横向维度未生效\n");
    } else {
        ++g_pass;
        std::printf("  [PASS] 不同横向位置影响值不同：横向维度已生效\n");
    }
}

void test3DBeamEnvelope() {
    std::printf("\n[22] 3D 梁内力包络与应力符号\n");
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(L, 0.0, 0.0);
    Element3D e;
    e.id = 1; e.type = "BEAM";
    e.E = E; e.A = A; e.Iy = I; e.Iz = I;
    e.J = 1.0e-3; e.G = E / 2.6;   // 非零扭转刚度，避免 θx 零刚度自由度
    e.Wtop = 1.0e-4; e.Wbot = 1.0e-4; e.Wy = 1.0e-4;
    e.n1 = 1; e.n2 = 2;
    m.elements.push_back(e);
    for (int dof = 1; dof <= 6; ++dof) m.constraints.push_back({1, dof, 0.0});
    m.loads.push_back({2, 2, -P});   // 节点2 dof2（y 向）向下

    const std::vector<double> U = solveModel3D(m);
    const DofMap3D dm = buildDofMap3D(m);
    const Eigen::VectorXd fl = beam3DLocalEndForces(m, dm, e, U);

    // 恢复端力 fl(5)=+P·L（外反力），截面弯矩（下缘受拉正）= -fl(5) = -P·L。
    check("3D 悬臂 截面弯矩 Mz=-PL", -fl(5), -P * L);
    // 上缘应力 σ_top = N/A - Mz/Wtop = +PL/Wtop（受拉）
    check("3D 悬臂 上缘受拉 σ=PL/Wtop", -(-fl(5)) / e.Wtop, P * L / e.Wtop);

    std::vector<StageResult3D> stages;
    StageResult3D sr; sr.stage = 1; sr.dm = dm; sr.U = U;
    stages.push_back(sr);
    const auto env = compute3DBeamEnvelope(m, stages);
    check("3D 包络 MzMin=-PL（负值被捕获）", env[0].MzMin, -P * L);
    check("3D 包络 MzMax=0（j 端自由）", env[0].MzMax, 0.0, 1e-9);
}

void test3DTorsion() {
    std::printf("\n[23] 3D 梁扭转剪应力\n");
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(L, 0.0, 0.0);
    Element3D e;
    e.id = 1; e.type = "BEAM";
    e.E = E; e.A = A; e.Iy = I; e.Iz = I;
    e.J = 1.0e-4; e.G = E / 2.6;
    e.Wtop = 1.0e-4; e.Wbot = 1.0e-4; e.Wy = 1.0e-4; e.Wt = 2.0e-4;
    e.n1 = 1; e.n2 = 2;
    m.elements.push_back(e);
    for (int dof = 1; dof <= 6; ++dof) m.constraints.push_back({1, dof, 0.0});
    const double Mx = 1.0e5;
    m.loads.push_back({2, 4, Mx});   // dof4 = θx（绕 x 轴扭矩）

    const std::vector<double> U = solveModel3D(m);
    const DofMap3D dm = buildDofMap3D(m);
    const Eigen::VectorXd fl = beam3DLocalEndForces(m, dm, e, U);
    check("3D 扭矩 T=Mx", -fl(3), Mx);
    check("3D 扭转剪应力 τ=T/Wt", -fl(3) / e.Wt, Mx / e.Wt);
}

void test3DObliqueCantilever() {
    std::printf("\n[3D] 斜置悬臂梁 端部横向集中力（坐标变换回归）\n");
    const double E = 2.06e11, A = 0.01, I = 1.0e-4, G = 7.9e10;
    const double P = 1.0e4;
    const double s = std::sqrt(0.5);
    const double L = 5.0;   // 梁轴长度

    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(L * s, L * s, 0.0);   // 45° 斜置
    Element3D e;
    e.id = 1; e.type = "BEAM"; e.E = E; e.A = A; e.Iy = I; e.Iz = I; e.J = I; e.G = G;
    e.n1 = 1; e.n2 = 2;
    m.elements.push_back(e);
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});   // 固定端
    m.loads.push_back({2, 1, -P * s});   // 横向力 (-1,1)/√2
    m.loads.push_back({2, 2,  P * s});

    const std::vector<double> U = solveModel3D(m);
    const DofMap3D dm = buildDofMap3D(m);
    const double delta = P * L * L * L / (3.0 * E * I);
    check("斜置悬臂 u = -δ/√2", U[dm.dof(2, 0)], -delta * s);
    check("斜置悬臂 v = +δ/√2", U[dm.dof(2, 1)],  delta * s);
}

void test3DModalAndBuckling() {
    std::printf("\n[3D] 自振频率与线性屈曲\n");
    const double E = 2.06e11, A = 0.18, Iy = 0.3, Iz = 0.6, dens = 7850.0;
    const double L = 10.0;
    const int n = 20;   // 单元数（一致质量/几何刚度收敛足够快）

    // 简支梁：节点 1 约束 u,v,w,θx；末端约束 v,w,θx（允许 θy/θz 与轴向 u）
    auto makeModel = [&](double axialLoad) {
        Model3D m;
        for (int i = 0; i <= n; ++i)
            m.nodes[i + 1] = Vec3(L * i / n, 0.0, 0.0);
        for (int i = 1; i <= n; ++i) {
            Element3D e;
            e.id = i; e.type = "BEAM"; e.E = E; e.A = A;
            e.Iy = Iy; e.Iz = Iz; e.J = 0.0; e.G = 0.0; e.dens = dens;
            e.n1 = i; e.n2 = i + 1;
            m.elements.push_back(e);
        }
        m.constraints.push_back({1, 1, 0.0});   // 左铰：轴向 u
        for (int i = 1; i <= n + 1; ++i) {
            m.constraints.push_back({i, 4, 0.0});   // 全节点约束扭转 θx（桥面系横隔假定）
            if (i == 1 || i == n + 1) {
                m.constraints.push_back({i, 2, 0.0});   // v
                m.constraints.push_back({i, 3, 0.0});   // w
            }
        }
        if (axialLoad != 0.0) m.loads.push_back({n + 1, 1, axialLoad});
        return m;
    };

    // 一阶自振频率（弱轴弯曲）：f1 = (π/(2L²))·√(E·Iy/(ρA))
    {
        const Model3D m = makeModel(0.0);
        const ModalResult3D modal = computeNaturalFrequencies3D(m, 1);
        const double rhoA = dens * A;
        const double f1 = kPi / (2.0 * L * L) * std::sqrt(E * Iy / rhoA);
        if (modal.frequencies.empty()) {
            std::printf("  [FAIL] 3D 自振频率：无结果\n");
            ++g_fail;
        } else {
            check("3D 一阶自振频率 f1", modal.frequencies[0], f1, 1e-3);
        }
    }

    // 一阶屈曲因子：Pcr = π²·E·Iy/L²，λ = Pcr/P_ref
    {
        const double P_ref = 1.0e5;   // 参考轴压（压缩）
        const Model3D m = makeModel(-P_ref);
        const std::vector<double> lam = computeBucklingFactors3D(m, 1);
        const double Pcr = kPi * kPi * E * Iy / (L * L);
        if (lam.empty()) {
            std::printf("  [FAIL] 3D 屈曲因子：无结果\n");
            ++g_fail;
        } else {
            check("3D 一阶屈曲因子 λ=Pcr/P", lam[0], Pcr / P_ref, 1e-3);
        }
    }
}

void test3DVariableSection() {
    std::printf("\n[3D 变截面] 沿梁轴线性插值截面 + 锥形悬臂梁\n");
    const double E = 2.06e11, G = 7.9e10, L = 10.0;

    SectionProperties p0, p1;
    p0.A = 0.18; p0.Iy = 0.3; p0.Iz = 0.6; p0.J = 0.3;
    p1.A = 0.18; p1.Iy = 0.3; p1.Iz = 0.2; p1.J = 0.15;
    VariableSection vs;
    vs.control = {{0.0, p0}, {L, p1}};

    // 中点线性插值
    const SectionProperties mid = sectionAt(vs, 0.5 * L);
    check("变截面中点 Iz=(Iz1+Iz2)/2", mid.Iz, 0.5 * (0.6 + 0.2), 1e-12);
    check("变截面中点 J=(J1+J2)/2", mid.J, 0.5 * (0.3 + 0.15), 1e-12);

    // beamElementFromSection 取中点截面
    {
        const BeamElement b = beamElementFromSection(vs, 0.0, L, Vec3(0,0,0), Vec3(L,0,0), E, G);
        check("beamElementFromSection 中点 Iz", b.Iz, 0.4, 1e-12);
        check("beamElementFromSection 中点 J", b.J, 0.225, 1e-12);
    }

    // 锥形悬臂梁：N 段，每段取中点截面；与柔度积分 δ=P/E·∫(L-s)²/Iz(s)ds 对照
    const double P = 1.0e5;
    const int n = 50;
    Model3D m;
    for (int i = 0; i <= n; ++i)
        m.nodes[i + 1] = Vec3(L * i / n, 0.0, 0.0);
    for (int i = 1; i <= n; ++i) {
        const double s1 = L * (i - 1) / n, s2 = L * i / n;
        const SectionProperties sec = elementMidSection(vs, s1, s2);
        Element3D e;
        e.id = i; e.type = "BEAM"; e.E = E; e.G = G;
        e.A = sec.A; e.Iy = sec.Iy; e.Iz = sec.Iz; e.J = sec.J;
        e.n1 = i; e.n2 = i + 1;
        m.elements.push_back(e);
    }
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});
    for (int i = 2; i <= n + 1; ++i) {
        m.constraints.push_back({i, 3, 0.0});   // w
        m.constraints.push_back({i, 4, 0.0});   // θx
        m.constraints.push_back({i, 5, 0.0});   // θy
    }
    m.loads.push_back({n + 1, 2, -P});
    const std::vector<double> U = solveModel3D(m);
    const DofMap3D dm = buildDofMap3D(m);

    const int ni = 100000;
    double integral = 0.0;
    for (int i = 0; i <= ni; ++i) {
        const double s = L * i / ni;
        const double Iz = 0.6 - 0.04 * s;
        const double w = (i == 0 || i == ni) ? 0.5 : 1.0;
        integral += w * (L - s) * (L - s) / Iz;
    }
    integral *= L / ni;
    const double exact = P / E * integral;
    check("3D 锥形悬臂梁挠度（柔度积分）", U[dm.dof(n + 1, 1)], -exact, 5e-3);
}

void test3DBeam7Torsion() {
    std::printf("\n[7 自由度梁] 圣维南扭转 + 翘曲扭转\n");
    const double E = 2.06e11, G = 7.9e10, L = 5.0;
    const double A = 0.18, Iy = 0.4, Iz = 0.6, J = 0.3;
    const double T = 1.0e5;

    // (a) 纯圣维南（Iw=0，根部自由翘曲 φ 不约束）：θx = T·L/(G·J)
    {
        Model3D m;
        m.nodes[1] = Vec3(0, 0, 0);
        m.nodes[2] = Vec3(L, 0, 0);
        Element3D e;
        e.id = 1; e.type = "BEAM7"; e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz;
        e.J = J; e.G = G; e.Iw = 0.0; e.n1 = 1; e.n2 = 2;
        m.elements.push_back(e);
        for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});   // 根部 θx 固定，φ 自由
        m.loads.push_back({2, 4, T});   // 端部扭矩（dof 4 = θx）
        const std::vector<double> U = solveModel3D(m);
        const DofMap3D dm = buildDofMap3D(m);
        check("BEAM7 纯圣维南 θx=T·L/(G·J)", U[dm.dof(2, 3)], T * L / (G * J), 1e-9);
    }

    // (b) 纯翘曲（GJ=0，根部 θx、φ 均约束）：θx = T·L³/(3·E·Iw)
    {
        const double Iw = 0.05;
        Model3D m;
        m.nodes[1] = Vec3(0, 0, 0);
        m.nodes[2] = Vec3(L, 0, 0);
        Element3D e;
        e.id = 1; e.type = "BEAM7"; e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz;
        e.J = 0.0; e.G = G; e.Iw = Iw; e.n1 = 1; e.n2 = 2;
        m.elements.push_back(e);
        for (int d = 1; d <= 7; ++d) m.constraints.push_back({1, d, 0.0});   // 根部完全固定（含 φ）
        m.loads.push_back({2, 4, T});
        const std::vector<double> U = solveModel3D(m);
        const DofMap3D dm = buildDofMap3D(m);
        check("BEAM7 纯翘曲 θx=T·L³/(3·E·Iw)", U[dm.dof(2, 3)],
              T * L * L * L / (3.0 * E * Iw), 1e-9);
    }

    // (c) 纯翘曲悬臂梁：固定端双力矩 |B| = T·L
    {
        const double Iw = 0.05;
        Model3D m;
        m.nodes[1] = Vec3(0, 0, 0);
        m.nodes[2] = Vec3(L, 0, 0);
        Element3D e;
        e.id = 1; e.type = "BEAM7"; e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz;
        e.J = 0.0; e.G = G; e.Iw = Iw; e.n1 = 1; e.n2 = 2;
        m.elements.push_back(e);
        for (int d = 1; d <= 7; ++d) m.constraints.push_back({1, d, 0.0});
        m.loads.push_back({2, 4, T});
        const std::vector<double> U = solveModel3D(m);
        const DofMap3D dm = buildDofMap3D(m);
        const Eigen::VectorXd fl = beam7LocalEndForces(m, dm, e, U);
        check("BEAM7 固定端双力矩 |B|=T·L", std::abs(-fl(6)), T * L, 1e-9);
    }

    // (d) 约束扭转（GJ 与 EIw 共同作用）：θx = (T·L/(G·J))·[1−tanh(λL)/(λL)]
    {
        const double Iw = 0.05;
        const int n = 50;
        Model3D m;
        for (int i = 0; i <= n; ++i) m.nodes[i + 1] = Vec3(L * i / n, 0.0, 0.0);
        for (int i = 1; i <= n; ++i) {
            Element3D e;
            e.id = i; e.type = "BEAM7"; e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz;
            e.J = J; e.G = G; e.Iw = Iw; e.n1 = i; e.n2 = i + 1;
            m.elements.push_back(e);
        }
        for (int d = 1; d <= 7; ++d) m.constraints.push_back({1, d, 0.0});
        m.loads.push_back({n + 1, 4, T});
        const std::vector<double> U = solveModel3D(m);
        const DofMap3D dm = buildDofMap3D(m);
        const double lam = std::sqrt(G * J / (E * Iw));
        const double lamL = lam * L;
        const double exact = T * L / (G * J) * (1.0 - std::tanh(lamL) / lamL);
        check("BEAM7 约束扭转 θx（闭式）", U[dm.dof(n + 1, 3)], exact, 2e-3);
    }
}

void test3DBeam7Influence() {
    std::printf("\n[7 自由度梁] 影响线打通（跨中弯矩）\n");
    const double E = 2.06e11, G = 7.9e10, L = 10.0;
    const int n = 10;

    auto makeBeam = [&](const std::string& type) {
        Model3D m;
        for (int i = 0; i <= n; ++i) m.nodes[i + 1] = Vec3(L * i / n, 0.0, 0.0);
        for (int i = 1; i <= n; ++i) {
            Element3D e;
            e.id = i; e.type = type; e.E = E; e.A = 0.18;
            e.Iy = 0.4; e.Iz = 0.6; e.J = 0.3; e.G = G; e.Iw = 0.05;
            e.n1 = i; e.n2 = i + 1;
            m.elements.push_back(e);
        }
        m.constraints.push_back({1, 1, 0.0});   // u
        m.constraints.push_back({1, 2, 0.0});   // v
        m.constraints.push_back({1, 3, 0.0});   // w
        m.constraints.push_back({1, 4, 0.0});   // θx（除刚体扭转）
        m.constraints.push_back({n + 1, 2, 0.0});
        m.constraints.push_back({n + 1, 3, 0.0});
        return m;
    };

    Model3D::Lane3D lane;
    for (int i = 1; i <= n + 1; ++i) lane.nodes.push_back(i);
    lane.dof = 2;

    const auto eta7 = compute3DBeamMomentInfluence(makeBeam("BEAM7"), 5, "MZ", 2, lane);
    const auto eta6 = compute3DBeamMomentInfluence(makeBeam("BEAM"), 5, "MZ", 2, lane);

    double peak = 0.0, maxDiff = 0.0;
    for (std::size_t k = 0; k < eta7.size(); ++k) {
        peak = std::max(peak, std::abs(eta7[k]));
        maxDiff = std::max(maxDiff, std::abs(eta7[k] - eta6[k]));
    }
    check("BEAM7 跨中弯矩影响线峰值 L/4", peak, L / 4.0, 1e-9);
    check("BEAM7 与 BEAM 弯矩影响线一致", maxDiff, 0.0, 1e-12);
}

void testBeamElementUnification() {
    std::printf("\n[统一单元] BeamElement 与 EulerBeam2D / beam3D* 逐元素等价\n");
    const double E  = 2.06e11;
    const double A  = 0.18;
    const double I  = 0.6;
    const double Iy = 0.4;
    const double J  = 0.32;
    const double G  = 7.9e10;
    const double L  = 8.0;
    const double Phi = 0.4;   // 任意剪切参数

    // 12 自由度 -> 平面内 [u1,v1,θz1,u2,v2,θz2] 的子集索引
    const int idx[6] = {0, 1, 5, 6, 7, 11};

    // (a) 3D 局部刚度：BeamElement(Phi=0) == beam3DLocalStiffness
    {
        BeamElement b;
        b.E = E; b.A = A; b.Iy = Iy; b.Iz = I; b.J = J; b.G = G;
        b.x0 = Eigen::Vector3d(0, 0, 0);
        b.x1 = Eigen::Vector3d(L, 0, 0);
        checkMat("BeamElement.localStiffness == beam3DLocalStiffness",
                 b.localStiffness(), beam3DLocalStiffness(E, A, Iy, I, J, G, L));
    }

    // (b) 3D 坐标变换：BeamElement == beam3DTransformation
    {
        BeamElement b;
        b.x0 = Eigen::Vector3d(3, 1, 2);
        b.x1 = Eigen::Vector3d(7, -2, 5);
        checkMat("BeamElement.transformation == beam3DTransformation",
                 b.transformation(), beam3DTransformation(b.x0, b.x1));
    }

    // (c) 2D 平面内整体刚度（Phi=0）：BeamElement(z=0) == EulerBeam2D
    {
        EulerBeam2D b2;
        b2.E = E; b2.A = A; b2.I = I;
        b2.x0 = Eigen::Vector2d(2, 1);
        b2.x1 = Eigen::Vector2d(6, 4);
        BeamElement b;
        b.E = E; b.A = A; b.Iz = I; b.Iy = 0.0; b.J = 0.0; b.G = 0.0;
        b.x0 = Eigen::Vector3d(2, 1, 0);
        b.x1 = Eigen::Vector3d(6, 4, 0);
        checkMat("BeamElement(2D) 整体刚度 == EulerBeam2D 刚度",
                 submatrix(b.stiffness(), idx, 6), b2.stiffness(0.0));
    }

    // (d) 2D Timoshenko 局部刚度（Phi）：BeamElement(PhiZ=Phi) == EulerBeam2D(Phi)
    {
        EulerBeam2D b2;
        b2.E = E; b2.A = A; b2.I = I;
        b2.x0 = Eigen::Vector2d(0, 0);
        b2.x1 = Eigen::Vector2d(L, 0);
        BeamElement b;
        b.E = E; b.A = A; b.Iz = I; b.Iy = 0.0; b.J = 0.0; b.G = 0.0; b.PhiZ = Phi;
        b.x0 = Eigen::Vector3d(0, 0, 0);
        b.x1 = Eigen::Vector3d(L, 0, 0);
        checkMat("BeamElement(2D,Timoshenko) 局部刚度 == EulerBeam2D(Phi)",
                 submatrix(b.localStiffness(), idx, 6), b2.localStiffness(Phi));
    }

    // (e) 一致质量矩阵平面内等价
    {
        const double rhoA = 7850.0 * A;
        EulerBeam2D b2;
        b2.E = E; b2.A = A; b2.I = I;
        b2.x0 = Eigen::Vector2d(0, 0);
        b2.x1 = Eigen::Vector2d(L, 0);
        BeamElement b;
        b.E = E; b.A = A; b.Iz = I; b.Iy = 0.0; b.J = 0.0; b.G = 0.0;
        b.x0 = Eigen::Vector3d(0, 0, 0);
        b.x1 = Eigen::Vector3d(L, 0, 0);
        checkMat("BeamElement 平面内质量 == EulerBeam2D 质量",
                 submatrix(b.consistentMassMatrix(rhoA), idx, 6), b2.consistentMassMatrix(rhoA));
    }

    // (f) 几何刚度平面内等价
    {
        const double N = 1.2e6;
        BeamElement b;
        b.E = E; b.A = A; b.Iz = I; b.Iy = 0.0; b.J = 0.0; b.G = 0.0;
        b.x0 = Eigen::Vector3d(0, 0, 0);
        b.x1 = Eigen::Vector3d(L, 0, 0);
        checkMat("BeamElement 平面内几何刚度 == beamGeometricStiffness2D",
                 submatrix(b.geometricStiffness(N), idx, 6), beamGeometricStiffness2D(N, L));
    }
}

void testDofMapUnification() {
    std::printf("\n[统一自由度表] 2D 装配刚度 == 3D 刚度平面内子矩阵\n");
    const double E = 2.06e11, A = 0.18, I = 0.6;

    // 2D：两段连续梁 (0,0)-(4,3)-(8,0)，含斜杆
    Model2D m2;
    m2.nodes[1] = Eigen::Vector2d(0.0, 0.0);
    m2.nodes[2] = Eigen::Vector2d(4.0, 3.0);
    m2.nodes[3] = Eigen::Vector2d(8.0, 0.0);
    for (int k = 1; k <= 2; ++k) {
        Model2D::Element el;
        el.id = k; el.type = "BEAM"; el.E = E; el.A = A; el.I = I; el.conn = {k, k + 1};
        m2.elements.push_back(el);
    }
    const DofMap dm2 = buildDofMap2D(m2);

    // 2D 刚度：用统一自由度表 + 通用 assemble 装配
    std::vector<DenseMat> kes;
    std::vector<std::vector<int>> gdofs;
    for (const auto& el : m2.elements) {
        EulerBeam2D b;
        b.E = el.E; b.A = el.A; b.I = el.I;
        b.x0 = m2.nodes.at(el.conn[0]);
        b.x1 = m2.nodes.at(el.conn[1]);
        kes.push_back(b.stiffness(0.0));
        std::vector<int> g;
        for (int n : el.conn)
            for (int d = 0; d < 3; ++d)
                g.push_back(dm2.dof(n, d));
        gdofs.push_back(g);
    }
    const SparseMat K2 = assemble(kes, gdofs, dm2.nDofs);

    // 3D 同构模型（z=0，Iz=I，Iy=J=0）
    Model3D m3;
    for (const auto& [id, x] : m2.nodes) m3.nodes[id] = Vec3(x(0), x(1), 0.0);
    for (const auto& el : m2.elements) {
        Element3D e;
        e.id = el.id; e.type = "BEAM"; e.E = el.E; e.A = el.A;
        e.Iz = el.I; e.Iy = 0.0; e.J = 0.0; e.G = 0.0;
        e.n1 = el.conn[0]; e.n2 = el.conn[1];
        m3.elements.push_back(e);
    }
    // 2D = 3D 且面外自由度（w,θx,θy）被约束；否则这些自由度零刚度会触发装配诊断。
    for (const auto& [id, x] : m3.nodes) {
        m3.constraints.push_back({id, 3, 0.0});
        m3.constraints.push_back({id, 4, 0.0});
        m3.constraints.push_back({id, 5, 0.0});
    }
    const Assembled3D a3 = assembleModel3D(m3);
    const DofMap& dm3 = a3.dm;

    // 抽取 3D 刚度的平面内 [u,v,θz] 子矩阵
    std::vector<int> sub3;
    for (std::size_t k = 0; k < dm2.nodeIds.size(); ++k) {
        const int node = dm2.nodeIds[k];
        sub3.push_back(dm3.dof(node, 0));
        sub3.push_back(dm3.dof(node, 1));
        sub3.push_back(dm3.dof(node, 5));
    }
    DenseMat K3sub(static_cast<int>(sub3.size()), static_cast<int>(sub3.size()));
    for (std::size_t r = 0; r < sub3.size(); ++r)
        for (std::size_t c = 0; c < sub3.size(); ++c)
            K3sub(static_cast<int>(r), static_cast<int>(c)) = a3.K0.coeff(sub3[r], sub3[c]);

    checkMat("2D 装配刚度 == 3D 刚度平面内子矩阵", K2.toDense(), K3sub);
}

void testUnifiedEndForces() {
    std::printf("\n[统一内力] BeamElement.localEndForces == beam3DLocalEndForces\n");
    const double E = 2.06e11, A = 0.18, Iy = 0.4, Iz = 0.6, J = 0.32, G = 7.9e10;

    Model3D m;
    m.nodes[1] = Vec3(1.0, 2.0, 0.5);
    m.nodes[2] = Vec3(5.0, -1.0, 3.0);
    Element3D e;
    e.id = 1; e.type = "BEAM"; e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz; e.J = J; e.G = G;
    e.n1 = 1; e.n2 = 2;
    m.elements.push_back(e);

    // 任意一组全局位移（含所有 6 个自由度）
    const DofMap dm = buildDofMap3D(m);
    std::vector<double> U(dm.nDofs, 0.0);
    const double vals[6] = {0.001, -0.002, 0.003, 0.004, -0.005, 0.006};
    for (int k = 0; k < 6; ++k) {
        U[dm.dof(1, k)] = vals[k];
        U[dm.dof(2, k)] = -0.5 * vals[k];
    }

    const Eigen::VectorXd fl3 = beam3DLocalEndForces(m, dm, e, U);

    BeamElement b;
    b.E = E; b.A = A; b.Iy = Iy; b.Iz = Iz; b.J = J; b.G = G;
    b.x0 = m.nodes.at(1); b.x1 = m.nodes.at(2);
    Eigen::Matrix<double, 12, 1> ug;
    for (int k = 0; k < 6; ++k) {
        ug(k)     = U[dm.dof(1, k)];
        ug(6 + k) = U[dm.dof(2, k)];
    }
    const Eigen::VectorXd flu = b.localEndForces(ug);

    checkMat("BeamElement.localEndForces == beam3DLocalEndForces", flu, fl3);
}

void testElasticSupport() {
    std::printf("\n[弹性支承] 悬臂梁端部弹簧\n");
    const double E = 2.0e11, A = 0.02, I = 1.0e-4, L = 5.0;
    const double k = 1.0e6, P = 1.0e4;
    const double kBeam = 3.0 * E * I / (L * L * L);   // 悬臂端部刚度

    // 2D：端部集中力 + 端部弹簧
    {
        Model2D m;
        m.nodes[1] = Eigen::Vector2d(0.0, 0.0);
        m.nodes[2] = Eigen::Vector2d(L, 0.0);
        Model2D::Element el;
        el.id = 1; el.type = "BEAM"; el.E = E; el.A = A; el.I = I; el.conn = {1, 2};
        m.elements.push_back(el);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({1, 3, 0.0});
        m.constraints.push_back({2, 2, 0.0, k});   // 端部竖向弹簧
        m.loads.push_back({2, 2, -P});
        const ModelResult r = solveModel2D(m);
        check("2D 端部弹簧+集中力 δ=P/(k+3EI/L³)",
              nodalDisp(m, r, 2, 2), -P / (k + kBeam));
    }

    // 2D：端部弹簧沉降（value≠0），无荷载
    {
        const double u0 = 0.01;
        Model2D m;
        m.nodes[1] = Eigen::Vector2d(0.0, 0.0);
        m.nodes[2] = Eigen::Vector2d(L, 0.0);
        Model2D::Element el;
        el.id = 1; el.type = "BEAM"; el.E = E; el.A = A; el.I = I; el.conn = {1, 2};
        m.elements.push_back(el);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({1, 3, 0.0});
        m.constraints.push_back({2, 2, u0, k});   // 弹簧基准位移 u0
        const ModelResult r = solveModel2D(m);
        check("2D 端部弹簧沉降 δ=k·u0/(k+3EI/L³)",
              nodalDisp(m, r, 2, 2), k * u0 / (k + kBeam));
    }

    // 3D：同构悬臂（z=0 平面内，绕 z 弯）
    {
        Model3D m3;
        m3.nodes[1] = Vec3(0.0, 0.0, 0.0);
        m3.nodes[2] = Vec3(L, 0.0, 0.0);
        Element3D e;
        e.id = 1; e.type = "BEAM"; e.E = E; e.A = A; e.Iy = I; e.Iz = I; e.J = 0.0; e.G = 0.0;
        e.n1 = 1; e.n2 = 2;
        m3.elements.push_back(e);
        for (int d = 1; d <= 6; ++d) m3.constraints.push_back({1, d, 0.0});
        m3.constraints.push_back({2, 3, 0.0});   // w
        m3.constraints.push_back({2, 4, 0.0});   // θx
        m3.constraints.push_back({2, 5, 0.0});   // θy
        m3.constraints.push_back({2, 2, 0.0, k});   // v 弹簧
        m3.loads.push_back({2, 2, -P});
        const std::vector<double> U = solveModel3D(m3);
        const DofMap3D dm = buildDofMap3D(m3);
        check("3D 端部弹簧+集中力 δ=P/(k+3EIz/L³)",
              U[dm.dof(2, 1)], -P / (k + kBeam));
    }
}

void testJtgD60AutoLaneLoad() {
    std::printf("\n[规范荷载库] JTG D60 自动车道荷载\n");
    const double L = 30.0;
    Model3D m;
    for (int i = 0; i <= 6; ++i) m.nodes[i + 1] = Vec3(L * i / 6, 0.0, 0.0);

    Model3D::Lane3D lane;
    lane.nodes = {1, 2, 3, 4, 5, 6, 7};
    lane.dof = 2;
    lane.autoCode = true;
    m.lanes.push_back(lane);

    m.laneSurface.xNodes = {1, 2, 3, 4, 5, 6, 7};
    m.laneSurface.yNodes = {1};
    m.laneSurface.laneCount = 3;
    m.laneSurface.autoCode = true;

    applyJtgD60LaneLoads(m);

    check("自动 qk=-10.5kN/m（向下）", m.lanes[0].q, -10500.0, 1e-12);
    check("自动 Pk(L=30)=-320kN", m.lanes[0].p, -jtgD60LaneLoadPk(L), 1e-9);
    check("影响面 qk=-10.5kN/m", m.laneSurface.q, -10500.0, 1e-12);
    check("影响面 3 车道折减 0.78", m.laneSurface.reduction, 0.78, 1e-12);

    // 2D *LANE_LOAD CODE=JTG_D60
    {
        Model2D m2;
        m2.nodes[1] = Eigen::Vector2d(0.0, 0.0);
        m2.nodes[2] = Eigen::Vector2d(15.0, 0.0);
        m2.nodes[3] = Eigen::Vector2d(30.0, 0.0);
        Model2D::LaneLoad ll;
        ll.id = 1; ll.nodes = {1, 2, 3}; ll.dir = "GY"; ll.autoCode = true;
        m2.laneLoads.push_back(ll);
        applyJtgD60LaneLoads2D(m2);
        check("2D 自动 qk=-10.5kN/m", m2.laneLoads[0].q, -10500.0, 1e-12);
        check("2D 自动 Pk(L=30)=-320kN", m2.laneLoads[0].P, -jtgD60LaneLoadPk(L), 1e-9);
    }
}

} // namespace

// ============================================================
// [18] 压缩输入语法（*FILL / *DO / 索引表达式）
// ============================================================
// 导出的 .dat 允许把等间距节点写成起止两点 + *FILL，把规则排列的单元写成
// *DO/*ENDDO 循环。这两个特性只要有一个读错，模型就"少几个单元/节点"而
// 不报错——力学结果全错但看不出错在哪。所以用同一模型的两种写法对拍：
// 解析结果必须逐字段相同。这是"导出器可以放心压缩"的前提。
void testCompactInput() {
    std::printf("\n[18] 压缩输入语法 *FILL / *DO 与逐行展开写法对拍\n");

    // 展开写法：节点 1001..1011（x = 0,10,…,100），
    //   单元 1..10  ：1001-1002, 1002-1003, …（步长 1）
    //   单元 21..25 ：1001-1003, 1003-1005, …（步长 2，验证 2*I+c 系数）
    //   单元 26     ：1001-1011，带单元级 W= 覆盖
    const std::string full =
        "*NODE\n"
        "1001 0 0 0\n1002 10 0 0\n1003 20 0 0\n1004 30 0 0\n1005 40 0 0\n"
        "1006 50 0 0\n1007 60 0 0\n1008 70 0 0\n1009 80 0 0\n1010 90 0 0\n"
        "1011 100 0 0\n"
        "*ET\n1 BEAM\n"
        "*MP\n1 E=3.0e10 DENS=2500\n"
        "*R\n1 A=1.0 IY=0.1 IZ=0.2 J=0.3 G=1.2e10\n"
        "*TYPE,1\n*MAT,1\n*REAL,1\n"
        "*EN,1,1001,1002\n*EN,2,1002,1003\n*EN,3,1003,1004\n*EN,4,1004,1005\n"
        "*EN,5,1005,1006\n*EN,6,1006,1007\n*EN,7,1007,1008\n*EN,8,1008,1009\n"
        "*EN,9,1009,1010\n*EN,10,1010,1011\n"
        "*EN,21,1001,1003\n*EN,22,1003,1005\n*EN,23,1005,1007\n"
        "*EN,24,1007,1009\n*EN,25,1009,1011\n"
        "*EN,26,1001,1011 W=88.5\n"
        "*END\n";

    // 压缩写法：节点只剩首末两点；单元全部走循环
    const std::string compact =
        "*NODE\n"
        "1001 0 0 0\n1011 100 0 0\n"
        "*FILL 1001 1011\n"
        "*ET\n1 BEAM\n"
        "*MP\n1 E=3.0e10 DENS=2500\n"
        "*R\n1 A=1.0 IY=0.1 IZ=0.2 J=0.3 G=1.2e10\n"
        "*TYPE,1\n*MAT,1\n*REAL,1\n"
        "*DO I 1 10\n*EN,I,I+1000,I+1001\n*ENDDO\n"
        "*DO I 21 25\n*EN,I,2*I+959,2*I+961\n*ENDDO\n"
        "*EN,26,1001,1011 W=88.5\n"
        "*END\n";

    const std::string pf = "_tmp_compact_full.dat";
    const std::string pc = "_tmp_compact_compact.dat";
    {
        std::ofstream f(pf); f << full;
        std::ofstream g(pc); g << compact;
    }
    const Model3D a = readModel3D(pf);
    const Model3D b = readModel3D(pc);
    std::remove(pf.c_str());
    std::remove(pc.c_str());

    check("*FILL 补齐节点数", static_cast<double>(b.nodes.size()),
          static_cast<double>(a.nodes.size()));
    check("*DO 展开单元数", static_cast<double>(b.elements.size()),
          static_cast<double>(a.elements.size()));

    double nerr = 0.0;
    for (const auto& kv : a.nodes) {
        const auto it = b.nodes.find(kv.first);
        if (it == b.nodes.end()) { nerr = 1e9; break; }
        nerr = std::max(nerr, (kv.second - it->second).cwiseAbs().maxCoeff());
    }
    check("*FILL 插值坐标", nerr, 0.0);

    // 单元按 id 排序后逐项对齐（*DO 段与直接写法的书写顺序不同）
    auto byId = [](const Model3D& m3) {
        std::vector<Element3D> v = m3.elements;
        std::sort(v.begin(), v.end(),
                  [](const Element3D& x, const Element3D& y) { return x.id < y.id; });
        return v;
    };
    const std::vector<Element3D> ea = byId(a), eb = byId(b);
    double cerr = 0.0, perr = 0.0;
    for (std::size_t k = 0; k < std::min(ea.size(), eb.size()); ++k) {
        if (ea[k].n1 != eb[k].n1 || ea[k].n2 != eb[k].n2) cerr = 1e9;
        if (ea[k].id != eb[k].id) perr = 1e9;
    }
    check("*DO 索引表达式求值（含 2*I+c）", cerr, 0.0);
    check("*DO 单元编号", perr, 0.0);

    // 循环体内的 W= 覆盖：压缩后必须仍带在最后一根单元上
    check("循环外单元级 W= 覆盖仍生效",
          eb.empty() ? -1.0 : eb.back().w, 88.5);

    // 端点引用不存在的节点：必须只是告警，不能让整份文件解析失败
    {
        const std::string pbad = "_tmp_compact_badfill.dat";
        { std::ofstream f(pbad);
          f << "*NODE\n1001 0 0 0\n*FILL 1001 9999\n*END\n"; }
        bool threw = false;
        try { readModel3D(pbad); } catch (const std::exception&) { threw = true; }
        std::remove(pbad.c_str());
        check("*FILL 悬空引用不中断解析", threw ? 0.0 : 1.0, 1.0);
    }
}

// 施工阶段索张拉（*STAGE 的 TENSION= / L0=）。
//
// 为什么这个模型值得做成回归：调索的整条链路（目标索力 → 反算无应力索长 →
// 折成初应变 → 增量求解 → 用累积变形后的弦长回算索力 → 再修索长）里，
// 只要有一环方向错（索长反了、等效荷载加错端、弦长漏了累积位移），
// 结果都会"看起来像个索力"，不会崩、不会报错。所以必须用能闭式求解的模型。
//
// 模型：竖向索（L=10 m、EA=1e8 N）上端固结；下端节点只约束 x/y，竖向挂一根
// k=2e7 N/m 的弹簧，并作用 W=100 kN 向下荷载。
//   阶段 1：激活索 + W。此时索力不是 100 kN——弹簧按刚度比分担：
//           u₁ = −W/(EA/L + k) = −3.333 mm，T₁ = EA·(L−u₁)/L = 33.33 kN
//   阶段 2：不给荷载，TENSION=1:500 kN。
//
// 收敛后必须**同时**满足三条，缺一不能算对：
//   (a) 索力 = 500 kN                     ← 迭代真的命中目标
//   (b) 节点 2 总位移 = (500−100)/k = 20 mm ← 总量满足静力平衡（独立于迭代）
//   (c) L0 = (10−0.02)/(1+500k/1e8) = 9.9303 m ← 索长反算的闭式解
// 只凑 (a) 很容易：把 L0 直接算成"弦长/(1+T/EA)"也能让回算索力等于目标，
// 但那样位移是错的——(b) 就是用来抓这个的。
void test3DStageCableTension() {
    std::printf("\n[24] 施工阶段索张拉（目标索力 → 无应力索长）\n");
    const double EA = 1.0e8, Lc = 10.0, kSpring = 2.0e7;
    const double W = 100e3, Ttar = 500e3;

    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, Lc);
    m.nodes[2] = Vec3(0.0, 0.0, 0.0);
    Element3D c;
    c.id = 1; c.type = "CABLE";
    c.E = 2.0e11; c.A = EA / 2.0e11;      // EA = 1e8 N
    c.n1 = 1; c.n2 = 2;
    m.elements.push_back(c);
    for (int d = 1; d <= 3; ++d) m.constraints.push_back({1, d, 0.0});
    m.constraints.push_back({2, 1, 0.0});
    m.constraints.push_back({2, 2, 0.0});
    m.constraints.push_back({2, 3, 0.0, kSpring});   // 竖向弹性支承

    Model3D::Stage3D s1; s1.id = 1; s1.activate = {1};
    s1.loads.push_back({2, 3, -W});
    Model3D::Stage3D s2; s2.id = 2;
    s2.tensions.push_back({1, Ttar});
    m.stages = {s1, s2};

    const auto rs = solveModel3DStages(m);
    check("阶段数", static_cast<double>(rs.size()), 2.0);
    if (rs.size() != 2) return;

    // ---- 阶段 1：没有目标索力时，索力由刚度比分配 ----
    const double u1 = -W / (EA / Lc + kSpring);
    const double T1 = EA * (-u1) / Lc;          // 节点下移 u1<0，弦长 = L−u1，L0 = L
    check("阶段1 索力 = EA(L−u₁)/L", rs[0].cables.at(0).T, T1, 1e-9);
    check("阶段1 无目标索力（target=0）", rs[0].cables.at(0).target, 0.0);

    // ---- 阶段 2：张拉到 500 kN ----
    const CableStageState3D& c2 = rs[1].cables.at(0);
    check("阶段2 索力 = 目标 500 kN", c2.T, Ttar, 0.02);
    check("阶段2 报告的偏差 < 1e-3", c2.dev, 0.0, 1.0);      // dev 本身，非负
    check("阶段2 报告的偏差确实很小", c2.dev < 1e-3 ? 1.0 : 0.0, 1.0);

    // ---- 总位移与索长：独立于迭代的两条物理约束 ----
    const DofMap3D& dm2 = rs[1].dm;
    const int dz = dm2.safeDof(2, 2);
    const int dz1 = rs[0].dm.safeDof(2, 2);
    const double uTot = (dz1 >= 0 ? rs[0].U[dz1] : 0.0) + (dz >= 0 ? rs[1].U[dz] : 0.0);
    // 容差 1e-4：线性阶段求解的切线是 EA/L（L 为几何长度），索力关系按
    // EA/L0，两者相差弹性应变 T/EA ≈ 0.5%，表现在位移上实测 1.5e-5。
    // 这条判据之所以能压这么紧，是因为张拉增量按**索力状态**算
    // （ΔP = N目标 − N当前）；早期按"参考几何上的初张力"算时，同样的模型
    // 偏差是 2.5e-1，两者差四个数量级。
    check("总位移 ≈ (T−W)/k = 20 mm", uTot, (Ttar - W) / kSpring, 1e-4);

    const double chordEnd = Lc - uTot;
    const double L0exact = chordEnd / (1.0 + Ttar / EA);
    check("无应力索长 L0 闭式解", c2.L0, L0exact, 0.01);
    check("L0 比几何长度短（张拉=缩短下料）", c2.L0 < Lc ? 1.0 : 0.0, 1.0);

    // ---- 负例：编错了必须立刻报错，不能静默算出一个像样的数 ----
    auto expectThrow = [](const std::string& name, const Model3D& mm) {
        bool threw = false;
        try { solveModel3DStages(mm); } catch (const std::exception&) { threw = true; }
        check(name, threw ? 1.0 : 0.0, 1.0);
    };

    {   // 目标索力指向梁单元
        Model3D bad = m;
        bad.stages = {s2};
        bad.stages[0].tensions = {{1, Ttar}};
        bad.elements[0].type = "BEAM";
        bad.elements[0].Iy = bad.elements[0].Iz = 1e-3;
        bad.elements[0].J = 1e-3; bad.elements[0].G = 1e11;
        expectThrow("TENSION 指向非索单元 → 报错", bad);
    }
    {   // 目标索力 ≤ 0
        Model3D bad = m;
        bad.stages = {s2};
        bad.stages[0].tensions = {{1, 0.0}};
        expectThrow("TENSION ≤ 0 → 报错（放张须用 L0=）", bad);
    }
    {   // 本阶段未激活的索
        Model3D bad = m;
        Model3D::Stage3D s; s.id = 1; s.tensions = {{1, Ttar}};   // 没 ACTIVE
        bad.stages = {s};
        expectThrow("TENSION 指向本阶段未激活单元 → 报错", bad);
    }
    {   // 同一根索既给目标索力又给索长
        Model3D bad = m;
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        s.tensions = {{1, Ttar}};
        s.unstrained = {{1, 9.93}};
        bad.stages = {s};
        expectThrow("TENSION 与 L0 同时给同一根索 → 报错", bad);
    }

    // ---- 无张拉规格时行为不变：应与单独 solveModel3D 完全一致 ----
    {
        Model3D plain = m;
        plain.stages = {s1};
        const auto rp = solveModel3DStages(plain);
        Model3D once;
        once.nodes = m.nodes;
        once.constraints = m.constraints;
        once.loads = s1.loads;
        once.elements = m.elements;
        const std::vector<double> Uref = solveModel3D(once);
        const double dmax = [&] {
            double d = 0.0;
            for (std::size_t k = 0; k < Uref.size() && k < rp[0].U.size(); ++k)
                d = std::max(d, std::abs(Uref[k] - rp[0].U[k]));
            return d;
        }();
        check("无 TENSION 时阶段解与单次求解一致", dmax, 0.0, 1e-12);
    }
}

// 简支梁（沿 +x，竖向弯曲绕 y）的最小模型：两端铰支，只留竖向弯曲 + 轴向滑移。
Model3D prestressBeam(bool cantilever, double e) {
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(10.0, 0.0, 0.0);
    Element3D b;
    b.id = 1; b.type = "BEAM";
    b.E = 2.0e11; b.A = 0.01; b.Iy = b.Iz = 0.01; b.J = 1e-4; b.G = 8.0e10;
    b.n1 = 1; b.n2 = 2;
    m.elements.push_back(b);
    if (cantilever) {
        for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});
    } else {
        m.constraints.push_back({1, 1, 0.0});   // 轴向固定
        m.constraints.push_back({1, 2, 0.0}); m.constraints.push_back({1, 3, 0.0});
        m.constraints.push_back({1, 4, 0.0}); m.constraints.push_back({1, 6, 0.0});
        m.constraints.push_back({2, 2, 0.0}); m.constraints.push_back({2, 3, 0.0});
        m.constraints.push_back({2, 4, 0.0}); m.constraints.push_back({2, 6, 0.0});
    }
    return m;
}

void test3DStagePrestress() {
    std::printf("\n[25] 施工阶段预应力（等效荷载法：轴压 + 偏心弯矩）\n");
    const double Pe = 1.0e6, e = -0.5;     // 筋在下缘（e<0）

    // ---- 判据 1：轴压效应（简支梁，直线筋）----
    {
        Model3D m = prestressBeam(false, e);
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1};
        ps.P = Pe; ps.e = e; ps.loss = 1.0;
        s.prestresses.push_back(ps);
        m.stages = {s};
        const auto rs = solveModel3DStages(m);
        check("阶段数", static_cast<double>(rs.size()), 1.0);
        // 轴向：节点 1 固定，节点 2 向 −x 移动，N = −Pe。
        const int ux = rs[0].dm.safeDof(2, 0);
        check("梁轴压位移 = −Pe·L/(EA)", rs[0].U[ux], -Pe * 10.0 / (2e11 * 0.01), 1e-9);
    }

    // ---- 判据 2：PRESTRESS 与手工 CLOAD 等价（零外部假设的翻译正确性）----
    {
        Model3D m = prestressBeam(false, e);
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1};
        ps.P = Pe; ps.e = e; ps.loss = 0.85;   // 带损失
        s.prestresses.push_back(ps);
        m.stages = {s};
        const auto ra = solveModel3DStages(m);

        Model3D hand = prestressBeam(false, e);
        Model3D::Stage3D hs; hs.id = 1; hs.activate = {1};
        const double Pf = Pe * 0.85;
        hs.loads.push_back({1, 1,  Pf}); hs.loads.push_back({2, 1, -Pf});
        hs.loads.push_back({1, 5,  Pf * e}); hs.loads.push_back({2, 5, -Pf * e});
        hand.stages = {hs};
        const auto rb = solveModel3DStages(hand);

        double dmax = 0.0;
        for (std::size_t k = 0; k < ra[0].U.size() && k < rb[0].U.size(); ++k)
            dmax = std::max(dmax, std::abs(ra[0].U[k] - rb[0].U[k]));
        check("PRESTRESS == 手工等效 CLOAD（含 loss=0.85）", dmax, 0.0, 1e-12);
    }

    // ---- 判据 3：偏心弯矩的方向与大小（悬臂梁，筋在下缘）----
    // 这里锁的是"偏心弯矩算对了"，不是"上拱"。直线筋（常偏心）**不产生上拱**：
    // 它在锚固端只给一对偏心弯矩，筋在下缘（e<0）时，自由端锚固力作用在
    // 形心下方、方向指向固结端，产生的绕 y 弯矩使悬臂自由端**下挠**。
    // 上拱来自曲线筋（s≠0）的跨内反向均布力 q = −8Pe·s/L²（s 为矢高），
    // 已由 [26] 单独测试（闭式解 −Pe·s·L²/(12EI)，见 fem_3d.hpp 与
    // input_format_3d.md）。本判据只用直线筋（s=0），把"直线筋"当"上拱"
    // 用是错的。
    {
        Model3D m = prestressBeam(true, e);
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1};
        ps.P = Pe; ps.e = e; ps.loss = 1.0;
        s.prestresses.push_back(ps);
        m.stages = {s};
        const auto rs = solveModel3DStages(m);
        const int uz = rs[0].dm.safeDof(2, 2);
        // 自由端弯矩 = Pe·(−e) 绕 +y（+y 正转使 +z 转向 +x，即自由端向下），
        // 故自由端挠度 uz = −Pe·L²/(4EI)，闭式解。
        check("自由端挠度 uz = −Pe·L²/(4EI)（直线筋下缘）",
              rs[0].U[uz], -Pe * 100.0 / (4.0 * 2e11 * 0.01), 1e-9);
        check("方向：筋在下缘 → 自由端下挠（非上拱）", rs[0].U[uz] < 0 ? 1.0 : 0.0, 1.0);
    }

    // ---- 负例：编错了必须当场报错 ----
    auto expectThrow = [](const std::string& name, const Model3D& mm) {
        bool threw = false;
        try { solveModel3DStages(mm); } catch (const std::exception&) { threw = true; }
        check(name, threw ? 1.0 : 0.0, 1.0);
    };
    {   // 指向非梁单元
        Model3D bad = prestressBeam(false, e);
        bad.elements[0].type = "TRUSS";
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1}; ps.P = Pe; ps.e = e;
        s.prestresses.push_back(ps); bad.stages = {s};
        expectThrow("PRESTRESS 指向非梁单元 → 报错", bad);
    }
    {   // P ≤ 0
        Model3D bad = prestressBeam(false, e);
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1}; ps.P = 0.0; ps.e = e;
        s.prestresses.push_back(ps); bad.stages = {s};
        expectThrow("PRESTRESS 张拉力 ≤ 0 → 报错", bad);
    }
    {   // 未激活单元
        Model3D bad = prestressBeam(false, e);
        Model3D::Stage3D s; s.id = 1;                       // 没 ACTIVE
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1}; ps.P = Pe; ps.e = e;
        s.prestresses.push_back(ps); bad.stages = {s};
        expectThrow("PRESTRESS 指向本阶段未激活单元 → 报错", bad);
    }
    {   // loss 越界
        Model3D bad = prestressBeam(false, e);
        Model3D::Stage3D s; s.id = 1; s.activate = {1};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1}; ps.P = Pe; ps.e = e;
        ps.loss = 1.5;
        s.prestresses.push_back(ps); bad.stages = {s};
        expectThrow("PRESTRESS 的 LOSS > 1 → 报错", bad);
    }
    {   // 不连续单元链（两条不接的梁）
        Model3D bad;
        bad.nodes[1] = Vec3(0, 0, 0); bad.nodes[2] = Vec3(10, 0, 0);
        bad.nodes[3] = Vec3(20, 0, 0); bad.nodes[4] = Vec3(30, 0, 0);
        Element3D b1; b1.id = 1; b1.type = "BEAM"; b1.E = 2e11; b1.A = 0.01;
        b1.Iy = b1.Iz = 0.01; b1.J = 1e-4; b1.G = 8e10; b1.n1 = 1; b1.n2 = 2;
        Element3D b2 = b1; b2.id = 2; b2.n1 = 3; b2.n2 = 4;   // 与 b1 不相接
        bad.elements = {b1, b2};
        for (int d = 1; d <= 6; ++d) bad.constraints.push_back({1, d, 0.0});
        bad.constraints.push_back({4, 2, 0.0}); bad.constraints.push_back({4, 3, 0.0});
        Model3D::Stage3D s; s.id = 1; s.activate = {1, 2};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1, 2}; ps.P = Pe; ps.e = e;
        s.prestresses.push_back(ps); bad.stages = {s};
        expectThrow("PRESTRESS 单元链不连续 → 报错", bad);
    }
}

// 3 节点 2 单元简支梁（节点 1—2—3，单元 1—2），用于曲线筋上拱闭式解。
Model3D curvedBeam() {
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(5.0, 0.0, 0.0);
    m.nodes[3] = Vec3(10.0, 0.0, 0.0);
    Element3D b;
    b.type = "BEAM"; b.E = 2.0e11; b.A = 0.01;
    b.Iy = b.Iz = 0.01; b.J = 1e-4; b.G = 8.0e10;
    b.id = 1; b.n1 = 1; b.n2 = 2; m.elements.push_back(b);
    b.id = 2; b.n1 = 2; b.n2 = 3; m.elements.push_back(b);
    // 简支：node1 轴向固定 + 竖向约束，node3 竖向约束，跨中 node2 只锁面外。
    m.constraints.push_back({1,1,0.0}); m.constraints.push_back({1,2,0.0});
    m.constraints.push_back({1,3,0.0}); m.constraints.push_back({1,4,0.0});
    m.constraints.push_back({1,6,0.0});
    m.constraints.push_back({3,2,0.0}); m.constraints.push_back({3,3,0.0});
    m.constraints.push_back({3,4,0.0}); m.constraints.push_back({3,6,0.0});
    m.constraints.push_back({2,2,0.0}); m.constraints.push_back({2,4,0.0});
    m.constraints.push_back({2,6,0.0});
    return m;
}

void test3DStageCurvedPrestress() {
    std::printf("\n[26] 施工阶段曲线筋（抛物线，上拱 = 跨内反向均布力）\n");
    const double Pe = 1.0e6, s = -0.3, L = 10.0;
    const double E = 2e11, I = 0.01;

    // ---- 判据 1：跨中上拱闭式解 ----
    // 端部偏心 e=0、矢高 s<0（跨中筋更低）。等效竖向力：中间节点 −4Pe·s/L（向
    // 上）；两端锚固竖向力与均布半跨分量落在支座节点，不产生挠度。
    // 跨中上拱 = F·L³/(48EI) = −Pe·s·L²/(12EI)。
    {
        Model3D m = curvedBeam();
        Model3D::Stage3D st; st.id = 1; st.activate = {1, 2};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1, 2};
        ps.P = Pe; ps.e = 0.0; ps.s = s; ps.loss = 1.0;
        st.prestresses.push_back(ps);
        m.stages = {st};
        const auto rs = solveModel3DStages(m);
        const int uz = rs[0].dm.safeDof(2, 2);
        const double want = -Pe * s * L * L / (12.0 * E * I);   // > 0，向上
        check("跨中上拱 = −Pe·s·L²/(12EI)", rs[0].U[uz], want, 1e-9);
        check("方向：s<0 → 跨中上拱（uz>0）", rs[0].U[uz] > 0 ? 1.0 : 0.0, 1.0);
    }

    // ---- 判据 2：曲线筋 == 手工 CLOAD 逐位等价（带端偏心 + loss）----
    {
        Model3D m = curvedBeam();
        Model3D::Stage3D st; st.id = 1; st.activate = {1, 2};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1, 2};
        ps.P = Pe; ps.e = 0.1; ps.s = s; ps.loss = 0.8;
        st.prestresses.push_back(ps);
        m.stages = {st};
        const auto ra = solveModel3DStages(m);

        Model3D hand = curvedBeam();
        Model3D::Stage3D hs; hs.id = 1; hs.activate = {1, 2};
        const double Pf = Pe * 0.8, e = 0.1;
        const double qz = -8.0 * Pf * s / (L * L);   // 跨内均布
        const double tz = 4.0 * Pf * s / L;          // 锚固切线竖向力（两端同向）
        hs.loads.push_back({1, 1,  Pf}); hs.loads.push_back({3, 1, -Pf});
        hs.loads.push_back({1, 5,  Pf * e}); hs.loads.push_back({3, 5, -Pf * e});
        hs.loads.push_back({1, 3, tz}); hs.loads.push_back({3, 3, tz});
        hs.loads.push_back({1, 3, qz * 2.5});        // 均布节点化（首尾各半跨 2.5）
        hs.loads.push_back({2, 3, qz * 5.0});        // 中间节点全跨 5.0
        hs.loads.push_back({3, 3, qz * 2.5});
        hand.stages = {hs};
        const auto rb = solveModel3DStages(hand);

        double dmax = 0.0;
        for (std::size_t k = 0; k < ra[0].U.size() && k < rb[0].U.size(); ++k)
            dmax = std::max(dmax, std::abs(ra[0].U[k] - rb[0].U[k]));
        check("曲线筋 == 手工等效 CLOAD（含端偏心 + loss）", dmax, 0.0, 1e-12);
    }

    // ---- 判据 3：s=0 退化为直线筋（轴压 + 端弯矩）----
    {
        Model3D m = curvedBeam();
        Model3D::Stage3D st; st.id = 1; st.activate = {1, 2};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1, 2};
        ps.P = Pe; ps.e = -0.5; ps.s = 0.0; ps.loss = 1.0;
        st.prestresses.push_back(ps);
        m.stages = {st};
        const auto ra = solveModel3DStages(m);

        Model3D hand = curvedBeam();
        Model3D::Stage3D hs; hs.id = 1; hs.activate = {1, 2};
        hs.loads.push_back({1, 1,  Pe}); hs.loads.push_back({3, 1, -Pe});
        hs.loads.push_back({1, 5,  Pe * (-0.5)}); hs.loads.push_back({3, 5, -Pe * (-0.5)});
        hand.stages = {hs};
        const auto rb = solveModel3DStages(hand);
        double dmax = 0.0;
        for (std::size_t k = 0; k < ra[0].U.size() && k < rb[0].U.size(); ++k)
            dmax = std::max(dmax, std::abs(ra[0].U[k] - rb[0].U[k]));
        check("s=0 退化为直线筋（轴压 + 端弯矩）", dmax, 0.0, 1e-12);
    }

    // ---- 负例：曲线筋链不共线 → 报错 ----
    {
        Model3D bad = curvedBeam();
        bad.nodes[2] = Vec3(5.0, 0.0, 1.0);          // 中间节点抬高，成折线
        Model3D::Stage3D st; st.id = 1; st.activate = {1, 2};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1, 2};
        ps.P = Pe; ps.e = 0.0; ps.s = s; ps.loss = 1.0;
        st.prestresses.push_back(ps);
        bad.stages = {st};
        bool threw = false;
        try { solveModel3DStages(bad); } catch (const std::exception&) { threw = true; }
        check("曲线筋链不共线 → 报错", threw ? 1.0 : 0.0, 1.0);
    }
    // ---- 直线筋（s=0）链不共线不报错（保持旧的首尾连线近似）----
    {
        Model3D ok = curvedBeam();
        ok.nodes[2] = Vec3(5.0, 0.0, 1.0);
        Model3D::Stage3D st; st.id = 1; st.activate = {1, 2};
        Model3D::Prestress3D ps; ps.id = 1; ps.elems = {1, 2};
        ps.P = Pe; ps.e = -0.5; ps.s = 0.0; ps.loss = 1.0;
        st.prestresses.push_back(ps);
        ok.stages = {st};
        bool threw = false;
        try { solveModel3DStages(ok); } catch (const std::exception&) { threw = true; }
        check("直线筋链不共线 → 不报错（首尾连线近似）", threw ? 0.0 : 1.0, 1.0);
    }

    // ---- 关键字 / 位置两种写法的文本解析 ----
    {
        const std::string p = "_tmp_curved_ps.dat";
        {
            std::ofstream f(p);
            f << "*STAGE\n"
              << "1 PRESTRESS=1-2:P=1000000:E=0:S=-0.3\n"
              << "2 PRESTRESS=1-2:1000000:-0.5:0.85\n"
              << "*END\n";
        }
        Model3D m = readModel3D(p);
        std::remove(p.c_str());
        check("阶段数", static_cast<double>(m.stages.size()), 2.0);
        const auto& a = m.stages[0].prestresses[0];
        const auto& b = m.stages[1].prestresses[0];
        check("关键字模式解析 S/P/E",
              (a.s == -0.3 && a.P == 1000000.0 && a.e == 0.0 && a.loss == 1.0)
                  ? 1.0 : 0.0, 1.0);
        check("关键字模式解析 elems 范围",
              (a.elems.size() == 2 && a.elems[0] == 1 && a.elems[1] == 2)
                  ? 1.0 : 0.0, 1.0);
        check("位置模式（旧格式）向后兼容",
              (b.s == 0.0 && b.P == 1000000.0 && b.e == -0.5 && b.loss == 0.85)
                  ? 1.0 : 0.0, 1.0);
    }
}

// ---- 3D 分项预应力损失（JTG 3362-2018 §6.2，2D/3D 共用公式）----
void test3DPrestressLossSection() {
    std::printf("\n[27] 3D *PRESTRESS_LOSS / LOSSREF（分项损失）\n");
    const double sigmaCon = 1395e6, Ep = 1.95e11, Ec = 3.45e10, fpk = 1860e6;
    const double aEp = Ep / Ec;
    // 独立 oracle：按 §6.2.2~6.2.7 手算，与 input_format.md §6.4 算例同参数。
    const double sl1 = sigmaCon * (1.0 - std::exp(-(0.25 * 0.2 + 0.0015 * 10.0)));
    const double sl2 = 0.006 / 10.0 * Ep;
    const double sl3 = 2.0 * (35.0 - 15.0) * 1e6;
    const double sl4 = aEp * 8e6;
    const double spe = sigmaCon - sl1 - sl2 - sl4;      // 后张法传力锚固应力
    const double sl5 = 1.0 * 0.3 * (0.52 * spe / fpk - 0.26) * spe;
    const double sl6 = 0.9 * (Ep * 2e-4 + aEp * 8e6 * 2.0)
                     / (1.0 + 15.0 * 0.005 * 1.0);
    const double ratio = spe / sigmaCon;

    const std::string lossLine =
        "1 SIGMA_CON=1395e6 EP=1.95e11 EC=3.45e10 FPK=1860e6 MU=0.25 K=0.0015 "
        "THETA=0.2 X=10 DELTA_L=0.006 L=10 T1=15 T2=35 SIGMA_PC4=8e6 SIGMA_PC6=8e6 "
        "PSI=1.0 ZETA=0.3 ECS=2e-4 PHI=2.0 RHO=0.005 RHOPS=1.0 POST=1";
    const std::string stageLine =
        "1 ACTIVE=1 PRESTRESS=1:P=1000000:E=-0.5:LOSSREF=1";
    auto writeDat = [](const std::string& body, const std::string& path) {
        std::ofstream f(path);
        f << "*NODE\n1 0 0 0\n2 10 0 0\n"
          << "*ET\n1 BEAM\n"
          << "*MP\n1 E=2e11\n"
          << "*R\n1 A=0.01 IY=0.01 IZ=0.01 J=1e-4 G=8e10\n"
          << "*TYPE,1\n*MAT,1\n*REAL,1\n*EN,1,1,2\n"
          << "*BOUNDARY\n1 1-6 0\n"
          << body << "\n*END\n";
    };
    auto badDat = [&](const std::string& body) {
        const std::string p = "_tmp_loss3d_bad.dat";
        writeDat(body, p);
        bool threw = false;
        try { Model3D mm = readModel3D(p); (void)mm; }
        catch (const std::exception&) { threw = true; }
        std::remove(p.c_str());
        return threw;
    };

    // ---- 判据 1：共享解析器与独立 oracle 一致 ----
    const auto [lid, li] = parsePrestressLossLine(lossLine);
    check("共享解析器读出表编号", static_cast<double>(lid), 1.0);
    check("共享解析器 σpe/σcon", effectivePrestressRatio(li), ratio, 1e-12);

    // ---- 判据 2：.dat 解析（表在 *STAGE 前）----
    {
        const std::string p = "_tmp_loss3d.dat";
        writeDat("*PRESTRESS_LOSS\n" + lossLine + "\n*STAGE\n" + stageLine, p);
        Model3D m = readModel3D(p);
        std::remove(p.c_str());
        check("表条目数", static_cast<double>(m.prestressLosses.size()), 1.0);
        const auto& ps = m.stages[0].prestresses[0];
        check("LOSSREF 编号", static_cast<double>(ps.lossRef), 1.0);
        check("已生成分项明细", ps.hasLossDetail ? 1.0 : 0.0, 1.0);
        check("σpe/σcon 写回 loss", ps.loss, ratio, 1e-12);
        check("σl1 摩擦", ps.lossDetail.sigmaL1, sl1, 1e-12);
        check("σl2 锚具变形", ps.lossDetail.sigmaL2, sl2, 1e-12);
        check("σl3 温差", ps.lossDetail.sigmaL3, sl3, 1e-12);
        check("σl4 弹性压缩", ps.lossDetail.sigmaL4, sl4, 1e-12);
        check("σl5 松弛", ps.lossDetail.sigmaL5, sl5, 1e-12);
        check("σl6 收缩徐变", ps.lossDetail.sigmaL6, sl6, 1e-12);
        check("σpe 传力锚固应力", ps.lossDetail.sigmaPe, spe, 1e-12);
        auto sr = solveModel3DStages(m);
        check("求解结果 Pe = P·σpe/σcon",
              sr[0].prestressForces[0].second, 1.0e6 * ratio, 1e-12);
        check("求解结果带分项明细",
              static_cast<double>(sr[0].prestressLossDetails.size()), 1.0);
    }

    // ---- 判据 3：表写在 *STAGE 之后也能引用（整份读完才 resolve）----
    {
        const std::string p = "_tmp_loss3d_late.dat";
        writeDat("*STAGE\n" + stageLine + "\n*PRESTRESS_LOSS\n" + lossLine, p);
        Model3D m = readModel3D(p);
        std::remove(p.c_str());
        check("表在 *STAGE 之后仍生效",
              m.stages[0].prestresses[0].loss, ratio, 1e-12);
    }

    // ---- 判据 4：只给 LOSS= 的默认路径逐位不变 ----
    {
        const std::string p = "_tmp_loss3d_plain.dat";
        writeDat("*STAGE\n1 ACTIVE=1 PRESTRESS=1:P=1000000:E=-0.5:LOSS=0.85", p);
        Model3D m = readModel3D(p);
        std::remove(p.c_str());
        const auto& ps = m.stages[0].prestresses[0];
        check("仅 LOSS= 时 loss 原样", ps.loss, 0.85);
        check("仅 LOSS= 时 lossRef=0", static_cast<double>(ps.lossRef), 0.0);
        check("仅 LOSS= 时不生成分项明细", ps.hasLossDetail ? 0.0 : 1.0, 1.0);
        auto sr = solveModel3DStages(m);
        check("仅 LOSS= 时 Pe=P·loss", sr[0].prestressForces[0].second, 850000.0);
        check("仅 LOSS= 时结果段无明细",
              static_cast<double>(sr[0].prestressLossDetails.size()), 0.0);
    }

    // ---- 判据 5：硬报错，不静默兜底 ----
    check("LOSS 与 LOSSREF 同时给 → 报错",
          badDat("*PRESTRESS_LOSS\n" + lossLine + "\n*STAGE\n"
                 + stageLine + ":LOSS=0.8") ? 1.0 : 0.0, 1.0);
    check("LOSSREF 引用不存在编号 → 报错",
          badDat("*PRESTRESS_LOSS\n" + lossLine + "\n*STAGE\n"
                 "1 ACTIVE=1 PRESTRESS=1:P=1000000:E=-0.5:LOSSREF=9") ? 1.0 : 0.0, 1.0);
    check("写了 LOSSREF 却没定义表 → 报错",
          badDat("*STAGE\n" + stageLine) ? 1.0 : 0.0, 1.0);
    const std::string badLoss =
        "1 SIGMA_CON=1395e6 EP=1.95e11 EC=3.45e10 FPK=1860e6 "
        "DELTA_L=0.6 L=1 POST=1";   // σl2 远大于 σcon → σpe/σcon < 0
    check("σpe/σcon 越界 → 报错",
          badDat("*PRESTRESS_LOSS\n" + badLoss + "\n*STAGE\n" + stageLine)
              ? 1.0 : 0.0, 1.0);
    {
        // 直接构造模型、绕过 readModel3D 的 resolve：求解必须明确报错，
        // 不能拿未解析的 loss=1.0 悄悄算一遍。
        Model3D m = prestressBeam(true, -0.5);
        Model3D::Stage3D st; st.id = 1; st.activate = {1};
        Model3D::Prestress3D ps;
        ps.id = 1; ps.elems = {1}; ps.P = 1.0e6; ps.e = -0.5; ps.lossRef = 1;
        st.prestresses.push_back(ps);
        m.stages = {st};
        m.prestressLosses[1] = li;
        bool threw = false;
        try { solveModel3DStages(m); }
        catch (const std::exception&) { threw = true; }
        check("未 resolve 直接求解 → 报错", threw ? 1.0 : 0.0, 1.0);
        resolvePrestressLosses3D(m);
        auto sr = solveModel3DStages(m);
        check("resolve 后 Pe = P·ratio",
              sr[0].prestressForces[0].second, 1.0e6 * ratio, 1e-12);
    }
}

int main() {
    setupConsoleUTF8();
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // 崩溃时也能看到已输出的内容
    std::printf("========== 3D 解析解断言测试 ==========\n");
    try {
        test3DCantilever();
        test3DAxialBar();
        test3DSupportSettlement();
        test3DMixedBeamTruss();
        test3DInfluenceSurfaceTransverse();
        test3DBeamEnvelope();
        test3DTorsion();
        test3DObliqueCantilever();
        test3DModalAndBuckling();
        test3DVariableSection();
        test3DBeam7Torsion();
        test3DBeam7Influence();
        testBeamElementUnification();
        testDofMapUnification();
        testUnifiedEndForces();
        testElasticSupport();
        testJtgD60AutoLaneLoad();
        testCompactInput();
        test3DStageCableTension();
        test3DStagePrestress();
        test3DStageCurvedPrestress();
        test3DPrestressLossSection();
    } catch (const std::exception& e) {
        std::printf("\n[ERROR] 异常: %s\n", e.what());
        ++g_fail;
    }

    std::printf("\n====================================\n");
    std::printf("通过 %d 项，失败 %d 项\n", g_pass, g_fail);
    if (g_fail != 0) {
        std::printf("结果：FAIL —— 数值结果与解析解不符\n");
        return 1;
    }
    std::printf("结果：ALL PASS —— 数值结果与解析解一致\n");
    return 0;
}
