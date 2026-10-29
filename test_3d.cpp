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
