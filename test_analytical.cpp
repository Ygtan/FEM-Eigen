// ============================================================
// 解析解断言测试
//
// 目的：run_tests.ps1 原先只校验"程序退出码 == 0"，任何数值回归
// （梁刚度系数写错、符号翻转、坐标变换出错）都不会被发现。这里用
// 材料力学闭式解做断言，使数值错误能被自动检出。
//
// 说明：欧拉-伯努利梁单元对"节点荷载"和"一致荷载向量"给出的节点
// 位移/转角/杆端力都是精确的，因此下面的比较应当达到机器精度量级。
// 判据取相对误差 1e-8 —— 远小于任何实际建模错误会引入的偏差。
// ============================================================

#include "fem_io.hpp"
#include "fem_3d.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace fem;

namespace {

int g_pass = 0;
int g_fail = 0;

void check(const std::string& name, double got, double exact, double tol = 1e-8) {
    const double denom = std::max(1.0, std::abs(exact));
    const double rel = std::abs(got - exact) / denom;
    if (rel <= tol) {
        ++g_pass;
        std::printf("  [PASS] %-44s got=%14.8g  exact=%14.8g  rel=%.2e\n",
                    name.c_str(), got, exact, rel);
    } else {
        ++g_fail;
        std::printf("  [FAIL] %-44s got=%14.8g  exact=%14.8g  rel=%.2e (tol=%.0e)\n",
                    name.c_str(), got, exact, rel, tol);
    }
}

// 等分直梁：nElem 个 BEAM 单元，沿 +X 布置，节点号 1..nElem+1。
Model2D straightBeam(double L, int nElem, double E, double A, double I) {
    Model2D m;
    m.title = "straight beam";
    for (int i = 0; i <= nElem; ++i)
        m.nodes[i + 1] = Eigen::Vector2d(L * i / static_cast<double>(nElem), 0.0);
    for (int e = 0; e < nElem; ++e) {
        Model2D::Element el;
        el.id = e + 1;
        el.type = "BEAM";
        el.E = E; el.A = A; el.I = I;
        el.conn = {e + 1, e + 2};
        m.elements.push_back(el);
    }
    return m;
}

double nodalDisp(const Model2D& m, const ModelResult& r, int node, int dof) {
    std::vector<int> ids;
    for (const auto& [id, c] : m.nodes) ids.push_back(id);
    for (std::size_t i = 0; i < ids.size(); ++i)
        if (ids[i] == node) return r.U(3 * static_cast<int>(i) + (dof - 1));
    throw std::runtime_error("节点不存在: " + std::to_string(node));
}

double nodalReaction(const Model2D& m, const ModelResult& r, int node, int dof) {
    std::vector<int> ids;
    for (const auto& [id, c] : m.nodes) ids.push_back(id);
    for (std::size_t i = 0; i < ids.size(); ++i)
        if (ids[i] == node) return r.R(3 * static_cast<int>(i) + (dof - 1));
    throw std::runtime_error("节点不存在: " + std::to_string(node));
}

constexpr double E  = 2.0e11;   // Pa（钢）
constexpr double A  = 0.02;     // m^2
constexpr double I  = 8.0e-4;   // m^4
constexpr double L  = 30.0;     // m
constexpr double P  = 1.0e5;    // N
constexpr double q  = 2.0e4;    // N/m

// ---- 1. 简支梁跨中集中力：δ = P L³ / (48 E I) ----
void testSimplySupportedPointLoad() {
    std::printf("\n[1] 简支梁 跨中集中力 P\n");
    const int n = 2;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});      // 左铰：u
    m.constraints.push_back({1, 2, 0.0});      // 左铰：v
    m.constraints.push_back({n + 1, 2, 0.0});  // 右滚：v
    const int mid = n / 2 + 1;                 // 节点 2
    m.loads.push_back({mid, 2, -P});

    const ModelResult r = solveModel2D(m);
    const double exact = -P * L * L * L / (48.0 * E * I);
    check("跨中挠度  δ=PL³/(48EI)", nodalDisp(m, r, mid, 2), exact);

    // 支座反力各 P/2
    check("左支座反力  R=P/2", nodalReaction(m, r, 1, 2), P / 2.0);
    check("右支座反力  R=P/2", nodalReaction(m, r, n + 1, 2), P / 2.0);
}

// ---- 2. 悬臂梁端部集中力：δ = P L³ / (3 E I)，θ = P L² / (2 E I) ----
void testCantileverPointLoad() {
    std::printf("\n[2] 悬臂梁 端部集中力 P\n");
    const int n = 1;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({1, 3, 0.0});   // 固定端
    m.loads.push_back({n + 1, 2, -P});

    const ModelResult r = solveModel2D(m);
    check("端部挠度  δ=PL³/(3EI)", nodalDisp(m, r, n + 1, 2), -P * L * L * L / (3.0 * E * I));
    check("端部转角  θ=PL²/(2EI)", nodalDisp(m, r, n + 1, 3), -P * L * L / (2.0 * E * I));
    check("固定端反力  R=P",       nodalReaction(m, r, 1, 2), P);
}

// ---- 3. 简支梁均布荷载：δ = 5 q L⁴ / (384 E I)，M中 = q L² / 8 ----
void testSimplySupportedUDL() {
    std::printf("\n[3] 简支梁 均布荷载 q\n");
    const int n = 4;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({n + 1, 2, 0.0});
    for (int e = 1; e <= n; ++e)
        m.elemLoads.push_back({e, "GY", -q, -q});

    const ModelResult r = solveModel2D(m);
    const int mid = n / 2 + 1;
    check("跨中挠度  δ=5qL⁴/(384EI)", nodalDisp(m, r, mid, 2),
          -5.0 * q * L * L * L * L / (384.0 * E * I));
    check("支座反力  R=qL/2",         nodalReaction(m, r, 1, 2), q * L / 2.0);

    // 跨中弯矩（取绝对值比较，避免符号约定纠缠）
    double mMax = 0.0;
    for (const auto& f : elementForces2D(m, r.U))
        mMax = std::max(mMax, std::max(std::abs(f.M1), std::abs(f.M2)));
    check("跨中弯矩  |M|=qL²/8", mMax, q * L * L / 8.0);
}

// ---- 4. 悬臂梁均布荷载：δ = q L⁴ / (8 E I) ----
void testCantileverUDL() {
    std::printf("\n[4] 悬臂梁 均布荷载 q\n");
    const int n = 4;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({1, 3, 0.0});
    for (int e = 1; e <= n; ++e)
        m.elemLoads.push_back({e, "GY", -q, -q});

    const ModelResult r = solveModel2D(m);
    check("端部挠度  δ=qL⁴/(8EI)", nodalDisp(m, r, n + 1, 2),
          -q * L * L * L * L / (8.0 * E * I));
    check("固定端弯矩  |M|=qL²/2",
          std::abs(elementForces2D(m, r.U).front().M1), q * L * L / 2.0);
}

// ---- 5. 梁端力符号约定（P2-6） ----
// 注释约定：M1>0 表示"下侧受拉"；N1>0 表示"受拉"。
// 这些约定一旦写反，输出"看起来合理"但物理意义全错，因此用解析解锁定。
void testBeamSignConvention() {
    std::printf("\n[5] 梁端力符号约定\n");

    // 5a. 悬臂梁端部向下集中力：固定端为"上侧受拉"→ M1 应为负。
    {
        Model2D m = straightBeam(L, 1, E, A, I);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({1, 3, 0.0});
        m.loads.push_back({2, 2, -P});
        const ModelResult r = solveModel2D(m);
        const auto forces = elementForces2D(m, r.U);
        const auto& f = forces.front();
        // 固定端弯矩 = -P*L（上侧受拉，按约定应为负）
        check("悬臂固定端 M1 = -PL（上侧受拉为负）", f.M1, -P * L);
        // 剪力约定：V1/V2 均为"截面剪力向上作用于单元"为正。
        // 悬臂常剪力 = P，因此两端均为 +P。
        check("悬臂左端剪力 V1 = P",  f.V1, P);
        check("悬臂右端剪力 V2 = P",  f.V2, P);
        // 轴力为零
        check("悬臂轴力 N1 = 0", f.N1, 0.0, 1e-12);
    }

    // 5b. 简支梁跨中向下集中力：跨中为"下侧受拉"→ M 应为正。
    {
        const int n = 2;
        Model2D m = straightBeam(L, n, E, A, I);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({n + 1, 2, 0.0});
        m.loads.push_back({2, 2, -P});            // 节点 2 为跨中
        const ModelResult r = solveModel2D(m);
        const auto forces = elementForces2D(m, r.U);
        // 单元 1（左半跨）：节点 2 端弯矩 M2 = +PL/4（下侧受拉）
        check("简支左半跨 M2 = +PL/4（下侧受拉为正）", forces[0].M2, P * L / 4.0);
        // 单元 2（右半跨）：节点 2 端弯矩 M1 = +PL/4
        check("简支右半跨 M1 = +PL/4（下侧受拉为正）", forces[1].M1, P * L / 4.0);
    }

    // 5c. 轴力杆受拉：N1 = N2 = +F（拉力为正）。
    {
        Model2D m;
        m.nodes[1] = Eigen::Vector2d(0.0, 0.0);
        m.nodes[2] = Eigen::Vector2d(L, 0.0);
        Model2D::Element el;
        el.id = 1; el.type = "TRUSS"; el.E = E; el.A = A; el.conn = {1, 2};
        m.elements.push_back(el);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({2, 2, 0.0});
        m.loads.push_back({2, 1, P});   // 节点 2 受向右拉力
        const ModelResult r = solveModel2D(m);
        const auto forces = elementForces2D(m, r.U);
        check("桁架拉杆 N1 = +F（受拉为正）", forces.front().N1, P);
        check("桁架拉杆 N2 = +F（受拉为正）", forces.front().N2, P);
    }

    // 5d. 轴力杆受压：N1 = N2 = -F（压力为负）。
    {
        Model2D m;
        m.nodes[1] = Eigen::Vector2d(0.0, 0.0);
        m.nodes[2] = Eigen::Vector2d(L, 0.0);
        Model2D::Element el;
        el.id = 1; el.type = "TRUSS"; el.E = E; el.A = A; el.conn = {1, 2};
        m.elements.push_back(el);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({2, 2, 0.0});
        m.loads.push_back({2, 1, -P});  // 节点 2 受向左压力
        const ModelResult r = solveModel2D(m);
        const auto forces = elementForces2D(m, r.U);
        check("桁架压杆 N1 = -F（受压为负）", forces.front().N1, -P);
        check("桁架压杆 N2 = -F（受压为负）", forces.front().N2, -P);
    }
}

// ---- 6. 轴力杆：δ = P L / (E A) ----
void testAxialBar() {
    std::printf("\n[6] 桁架杆 轴向拉伸\n");
    Model2D m;
    m.nodes[1] = Eigen::Vector2d(0.0, 0.0);
    m.nodes[2] = Eigen::Vector2d(L, 0.0);
    Model2D::Element el;
    el.id = 1; el.type = "TRUSS"; el.E = E; el.A = A; el.conn = {1, 2};
    m.elements.push_back(el);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({2, 2, 0.0});
    m.loads.push_back({2, 1, P});

    const ModelResult r = solveModel2D(m);
    check("轴向伸长  δ=PL/(EA)", nodalDisp(m, r, 2, 1), P * L / (E * A));

    const auto forces = elementForces2D(m, r.U);
    check("轴力  N=P", forces.front().N1, P);
}

// ---- 6. 两跨连续梁（超静定）：中支座反力 = 5qL/8（每跨 q、跨度 L）----
// 用于验证组装与约束处理在多跨体系下仍然正确。
void testContinuousBeam() {
    std::printf("\n[7] 两跨连续梁 均布荷载（超静定）\n");
    const int nPerSpan = 4;
    const int n = 2 * nPerSpan;
    Model2D m = straightBeam(2.0 * L, n, E, A, I);
    const int left = 1, midNode = nPerSpan + 1, right = n + 1;
    m.constraints.push_back({left, 1, 0.0});
    m.constraints.push_back({left, 2, 0.0});
    m.constraints.push_back({midNode, 2, 0.0});
    m.constraints.push_back({right, 2, 0.0});
    for (int e = 1; e <= n; ++e)
        m.elemLoads.push_back({e, "GY", -q, -q});

    const ModelResult r = solveModel2D(m);
    // 两跨等跨等载：中支座反力 1.25 qL，端支座 0.375 qL，合计 2qL
    check("中支座反力  R=1.25qL", nodalReaction(m, r, midNode, 2), 1.25 * q * L, 1e-6);
    check("端支座反力  R=0.375qL", nodalReaction(m, r, left, 2), 0.375 * q * L, 1e-6);
    const double total = nodalReaction(m, r, left, 2) + nodalReaction(m, r, midNode, 2)
                       + nodalReaction(m, r, right, 2);
    check("反力合计  ΣR=2qL（整体平衡）", total, 2.0 * q * L, 1e-6);
}

// ==================== 三维部分 ====================

// 沿 +X 的等分直梁（BEAM）。
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

// ---- 7. 三维悬臂梁端部集中力：δ = P L³ / (3 E Iz) ----
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

// ---- 8. 三维桁架杆：δ = P L / (E A) ----
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

// ---- 9. 支座沉降（非零指定位移）----
// 旧实现的 Model3D 约束只有 (node, dof)，且求解时硬编码 F(d)=0，
// 三维模型完全无法施加支座沉降。
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

// ---- 10. 梁+桁架混合：按节点分配自由度 ----
// 旧实现"存在一个 BEAM 就给所有节点 6 个自由度"，纯桁架节点的转动
// 自由度恒为零刚度，总体矩阵奇异 —— 这个算例在旧实现下会直接报错。
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

// ---- 11. 影响面横向维度 ----
// 旧实现的内层循环不使用 yn，同一条一维影响线被复制 N 遍。
// 这里验证不同横向位置会得到不同的影响值。
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

// ---- 13. 徐变 AEMM：可重入 + 逐阶段时效 + 多 spec 不连乘 ----
//
// 物理判据：加载龄期 0、徐变系数 φ 的梁，长期挠度 = 弹性挠度 × (1 + χφ)。
// 这是 AEMM 的直接推论，也是最容易验错的一项。
void testCreepAEMM() {
    std::printf("\n[14] 徐变 AEMM（可重入 / 逐阶段 / 多 spec）\n");

    const int n = 2;
    const double chi = 0.8;
    const double phi = 1.5;

    // 基准：无徐变的弹性解
    Model2D base = straightBeam(L, n, E, A, I);
    base.constraints.push_back({1, 1, 0.0});
    base.constraints.push_back({1, 2, 0.0});
    base.constraints.push_back({n + 1, 2, 0.0});
    base.loads.push_back({2, 2, -P});
    const double d0 = std::abs(nodalDisp(base, solveModel2D(base), 2, 2));

    auto makeCreepModel = [&]() {
        Model2D m = straightBeam(L, n, E, A, I);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({n + 1, 2, 0.0});
        m.loads.push_back({2, 2, -P});
        Model2D::CreepSpec cs;
        cs.all = true;
        cs.phi = phi;
        cs.chi = chi;
        m.creepSpecs = {cs};
        snapshotBaseModulus(m);
        return m;
    };

    // (a) 徐变放大倍数 (1 + χφ)
    {
        Model2D m = makeCreepModel();
        processTimeEffects(m);
        const double d = std::abs(nodalDisp(m, solveModel2D(m), 2, 2));
        check("徐变长期挠度 = 弹性挠度 × (1+χφ)", d, d0 * (1.0 + chi * phi), 1e-10);
    }

    // (b) 可重入：再处理一次，结果必须完全不变
    //     （旧实现 `e.E *= factor` 会把模量折减两次）
    {
        Model2D m = makeCreepModel();
        processTimeEffects(m);
        const double E1 = m.elements[0].E;
        processTimeEffects(m);
        processTimeEffects(m);
        const double E3 = m.elements[0].E;
        check("processTimeEffects 幂等（调用 3 次模量不变）", E3, E1, 1e-14);
        check("折减后模量 = E0/(1+χφ)", E1, E / (1.0 + chi * phi), 1e-12);
    }

    // (c) 多 spec 不连乘：两条相同 spec 折减一次，而不是平方
    {
        Model2D m = makeCreepModel();
        m.creepSpecs.push_back(m.creepSpecs.front());
        processTimeEffects(m);
        check("两条相同 *CREEP 仅折减一次", m.elements[0].E,
              E / (1.0 + chi * phi), 1e-12);
    }

    // (d) 逐阶段 AEMM：Δφ = φ(t_now) - φ(t_load)
    //     晚加载的单元折减少（Δφ 小），早加载的折减多。
    {
        const double tLoad = 100.0;
        const double tNow  = 1000.0;
        Model2D m = makeCreepModel();
        m.creepSpecs.front().phi = 0.0;   // 强制走龄期路径
        std::map<int, double> loadAge;
        for (const auto& e : m.elements) loadAge[e.id] = tLoad;
        processTimeEffects(m, &loadAge, tNow);

        const double phiNow = creepCoefficientAt(m, tNow);
        const double phiT0  = creepCoefficientAt(m, tLoad);
        const double dphi   = std::max(0.0, phiNow - phiT0);
        check("逐阶段 AEMM：E_eff = E0/(1+χΔφ)", m.elements[0].E,
              E / (1.0 + chi * dphi), 1e-12);

        // 晚加载（Δφ 小）的有效模量必须高于早加载（Δφ 大）
        std::map<int, double> early;
        for (const auto& e : m.elements) early[e.id] = 0.0;
        Model2D m2 = makeCreepModel();
        m2.creepSpecs.front().phi = 0.0;
        processTimeEffects(m2, &early, tNow);
        if (m.elements[0].E > m2.elements[0].E && dphi > 0.0) {
            ++g_pass;
            std::printf("  [PASS] %-44s late=%.6g > early=%.6g\n",
                        "晚加载单元模量高于早加载单元",
                        m.elements[0].E, m2.elements[0].E);
        } else {
            ++g_fail;
            std::printf("  [FAIL] %-44s late=%.6g early=%.6g dphi=%g\n",
                        "晚加载单元模量高于早加载单元",
                        m.elements[0].E, m2.elements[0].E, dphi);
        }
    }

    // (e) 多阶段：阶段 2 的有效模量必须是 E0/(1+χ·(φ(t2)-φ(t0)))，
    //     而不是旧实现的逐级累乘 E0/((1+χφ1)(1+χφ2))。
    //     判据用简支梁跨中挠度 δ = PL³/(48·Eeff·I)。
    {
        Model2D m = straightBeam(L, 2, E, A, I);
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({3, 2, 0.0});
        // 时变徐变只在给了时变模型时才启用；CURVE 且无数据点时会自动关闭。
        m.creepModel = "CEB90";
        m.creepPhi0 = 2.0; m.creepT0 = 7.0;
        m.creepBetaH = 700.0; m.creepPower = 0.3;

        Model2D::Stage s1;
        s1.id = 1;
        s1.time = 100.0;
        s1.activate = {1, 2};          // 第 1 阶段安装，加载龄期 = 0
        Model2D::Stage s2;
        s2.id = 2;
        s2.time = 1000.0;              // 第 2 阶段才加载
        s2.loads.push_back({2, 2, -P});
        m.stages = {s1, s2};

        const auto res = solveModel2DStages(m);
        // 节点 2 是 0-based 索引 1，dof 2 → 3*1 + 1 = 4
        const double d2 = std::abs(res.back().U(4));

        const double phi1 = creepCoefficientAt(m, 100.0);
        const double phi2 = creepCoefficientAt(m, 1000.0);
        const double EeffNew = E / (1.0 + 0.8 * (phi2 - 0.0));            // AEMM
        const double EeffOld = E / ((1.0 + 0.8 * phi1) * (1.0 + 0.8 * phi2)); // 旧累乘
        const double dNew = P * L * L * L / (48.0 * EeffNew * I);
        const double dOld = P * L * L * L / (48.0 * EeffOld * I);

        check("两阶段徐变 阶段2 挠度（AEMM）", d2, dNew, 1e-8);
        const double gap = std::abs(d2 - dOld) / dOld;
        if (gap > 1e-2) {
            ++g_pass;
            std::printf("  [PASS] %-44s 与旧累乘值差 %.1f%%\n",
                        "已不是旧的逐级累乘（修复生效）", gap * 100.0);
        } else {
            ++g_fail;
            std::printf("  [FAIL] %-44s 差仅 %.3g%%，疑似仍在累乘\n",
                        "已不是旧的逐级累乘（修复生效）", gap * 100.0);
        }
    }
}

// ---- 14. 非线性：两杆桁架（von Mises）大位移解析解 + 收敛保障 ----
//
// 几何：两杆长 L0，下端铰支于 (±b, 0)，顶点在 (0, h)。顶点受竖向荷载 P 下沉 v。
// 变形后杆长 L = √(b² + (h-v)²)，应变 ε = (L-L0)/L0，轴力 N = EA·ε。
// 顶点竖向平衡：P = 2·N·(h-v)/L
// 这是几何非线性的经典闭式解，用来锁住 Newton 内核。
void testNonlinearTwoBarTruss() {
    std::printf("\n[15] 非线性 两杆桁架（von Mises）+ 收敛保障\n");

    const double b = 2.0;      // 半跨
    const double h = 0.5;      // 初始矢高（浅桁架 → 强几何非线性）
    const double Abar = 1.0e-3;
    const double L0 = std::sqrt(b * b + h * h);

    Model2D m;
    m.analysisType = "NONLINEAR";
    m.nodes[1] = Eigen::Vector2d(-b, 0.0);
    m.nodes[2] = Eigen::Vector2d(0.0, h);
    m.nodes[3] = Eigen::Vector2d(b, 0.0);
    for (int e = 0; e < 2; ++e) {
        Model2D::Element el;
        el.id = e + 1;
        el.type = "TRUSS";
        el.E = E; el.A = Abar; el.I = 0.0;
        el.conn = {e + 1, e + 2};
        m.elements.push_back(el);
    }
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({3, 1, 0.0});
    m.constraints.push_back({3, 2, 0.0});

    // 解析解：给定顶点下沉 v，反算所需荷载 P(v)。
    // 顶点下移 → 杆长 L 变短 → 杆受压（N<0）→ 竖向分力向上抵抗外荷载，
    // 故 P = -2·N·(h-v)/L，N = EA·(L-L0)/L0。
    auto loadForDisp = [&](double v) {
        const double L = std::sqrt(b * b + (h - v) * (h - v));
        const double N = E * Abar * (L - L0) / L0;
        return -2.0 * N * (h - v) / L;
    };

    // 取一个位于"上升段"（远低于失稳临界点）的荷载，二分求解析位移。
    const double Pdown = 8.0e5;
    const double hiBound = 0.12;   // P(v) 在 [0, 0.12] 上单调递增
    auto solveExact = [&]() {
        double lo = 0.0, hi = hiBound;
        for (int i = 0; i < 200; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (loadForDisp(mid) < Pdown) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    };
    const double vExact = solveExact();

    // 线性解作为对照：确认这个算例确实处在非线性显著区（否则测不出东西）
    const double kLin = 2.0 * (E * Abar / L0) * (h / L0) * (h / L0);
    const double vLin = Pdown / kLin;
    std::printf("  参考：线性解 v=%.6g，非线性解 v=%.6g（相差 %.1f%%）\n",
                vLin, vExact, 100.0 * (vExact - vLin) / vLin);
    if (std::abs(vExact - vLin) / vLin < 0.05) {
        ++g_fail;
        std::printf("  [FAIL] 算例非线性效应不足 5%%，测不出非线性内核的问题\n");
    }

    // (a) 不分級 + 线搜索：与解析解对比
    m.loads.push_back({2, 2, -Pdown});
    const ModelResult r1 = solveModel2D(m);
    const double v1 = -nodalDisp(m, r1, 2, 2);
    check("两杆桁架 顶点下沉（不分級）", v1, vExact, 1e-6);

    // (b) 荷载分级 8 步：结果必须与不分級一致
    {
        Model2D ms = m;
        ms.nonlinearLoadSteps = 8;
        const ModelResult r2 = solveModel2D(ms);
        const double v2 = -nodalDisp(ms, r2, 2, 2);
        check("两杆桁架 顶点下沉（8 级加载）", v2, vExact, 1e-6);
        check("分级与不分级结果一致", v2, v1, 1e-7);
    }

    // (c) 关掉线搜索：正常收敛工况下结果仍应一致
    {
        Model2D ms = m;
        ms.nonlinearLineSearch = false;
        const ModelResult r3 = solveModel2D(ms);
        check("关闭线搜索 结果仍一致",
              -nodalDisp(ms, r3, 2, 2), vExact, 1e-6);
    }

    // (d) 不收敛时返回"部分成功结果"而不是直接抛异常
    {
        Model2D ms = m;
        ms.nonlinearThrowOnFailure = false;
        ms.nonlinearTol = 1e-30;          // 容差低于机器精度 → 必然达不到
        ms.nonlinearMaxIter = 5;
        ms.nonlinearLoadSteps = 1;
        NonlinearDiag diag;
        const ModelResult r4 = solveModel2D(ms, SolverType::SimplicialLDLT,
                                            nullptr, nullptr, &diag);
        if (!r4.nonlinearConverged) {
            ++g_pass;
            std::printf("  [PASS] %-44s 已施加%.1f%%荷载\n",
                        "不收敛时返回部分结果（不抛异常）", diag.loadRatio * 100.0);
        } else {
            ++g_fail;
            std::printf("  [FAIL] %-44s 竟然判定为收敛\n", "不收敛时返回部分结果（不抛异常）");
        }
        if (!diag.failureReason.empty() && diag.maxResidualNode >= 0 &&
            !diag.residualHistory.empty()) {
            ++g_pass;
            std::printf("  [PASS] %-44s 最大残差 节点%d dof%d\n",
                        "诊断信息含原因/残差自由度/残差历史",
                        diag.maxResidualNode, diag.maxResidualDof);
        } else {
            ++g_fail;
            std::printf("  [FAIL] %-44s\n", "诊断信息不完整");
        }
    }

    // (e) 抛异常时，异常文本必须带诊断信息（而不是一句"未收敛"）
    {
        Model2D ms = m;
        ms.nonlinearTol = 1e-30;          // 同上，强制不收敛
        ms.nonlinearMaxIter = 5;
        ms.nonlinearLoadSteps = 1;
        bool threw = false;
        std::string msg;
        try {
            solveModel2D(ms);
        } catch (const std::exception& ex) {
            threw = true;
            msg = ex.what();
        }
        if (threw && msg.find("残差最大自由度") != std::string::npos &&
            msg.find("残差历史") != std::string::npos) {
            ++g_pass;
            std::printf("  [PASS] %-44s\n", "未收敛异常带诊断信息");
        } else {
            ++g_fail;
            std::printf("  [FAIL] %-44s threw=%d msg=%.120s\n",
                        "未收敛异常带诊断信息", static_cast<int>(threw), msg.c_str());
        }
    }
}

// ---- 15. 求解器复用：符号分解只做一次 ----
//
// P0-1/P1-1 已经把"重复装配+重复分解"从热路径上拿掉了，但真正省掉的是
// **数值**分解那一半；符号分解（AMD 排序 + 消去树）只有结构变化时才需要。
// 这里用 analyze/factorize 计数把这件事钉死，防止以后改回 compute()。
void testSolverReuse() {
    std::printf("\n[16] 求解器复用（analyze 一次 / factorize 多次）\n");

    const int n = 4;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({n + 1, 2, 0.0});
    m.loads.push_back({3, 2, -P});

    // 同一结构、不同数值（模拟徐变折减 / 阶段推进 / Newton 迭代）
    const Assembled2D a1 = assembleModel2DLinear(m);

    Model2D m2 = m;
    for (auto& e : m2.elements) e.E *= 0.5;      // 数值变了，结构没变
    const Assembled2D a2 = assembleModel2DLinear(m2);

    ReusableSolver rs(SolverType::SimplicialLDLT);
    rs.factorizeReuse(a1.K);
    const DenseVec u1 = rs.solve(a1.F);
    rs.factorizeReuse(a2.K);
    const DenseVec u2 = rs.solve(a2.F);

    if (rs.analyzeCount() == 1 && rs.factorizeCount() == 2) {
        ++g_pass;
        std::printf("  [PASS] %-44s analyze=%d factorize=%d\n",
                    "结构不变：符号分解只做一次",
                    rs.analyzeCount(), rs.factorizeCount());
    } else {
        ++g_fail;
        std::printf("  [FAIL] %-44s analyze=%d factorize=%d（期望 1 / 2）\n",
                    "结构不变：符号分解只做一次",
                    rs.analyzeCount(), rs.factorizeCount());
    }

    // 结果与一次性求解一致
    const DenseVec ref1 = solve(a1.K, a1.F, SolverType::SimplicialLDLT);
    const DenseVec ref2 = solve(a2.K, a2.F, SolverType::SimplicialLDLT);
    check("复用求解与一次性求解一致（工况1）", (u1 - ref1).norm(), 0.0, 1e-12);
    check("复用求解与一次性求解一致（工况2）", (u2 - ref2).norm(), 0.0, 1e-12);

    // 结构变化时必须重新 analyze
    Model2D m3 = m;
    Model2D::Element extra;
    extra.id = 100;
    extra.type = "BEAM";
    extra.E = E; extra.A = A; extra.I = I;
    extra.conn = {1, 3};                       // 新增单元 → 非零结构改变
    m3.elements.push_back(extra);
    const Assembled2D a3 = assembleModel2DLinear(m3);
    rs.factorizeReuse(a3.K);
    if (rs.analyzeCount() == 2 && rs.factorizeCount() == 3) {
        ++g_pass;
        std::printf("  [PASS] %-44s analyze=%d factorize=%d\n",
                    "结构变化：自动重新符号分解",
                    rs.analyzeCount(), rs.factorizeCount());
    } else {
        ++g_fail;
        std::printf("  [FAIL] %-44s analyze=%d factorize=%d（期望 2 / 3）\n",
                    "结构变化：自动重新符号分解",
                    rs.analyzeCount(), rs.factorizeCount());
    }

    // 多右端项
    DenseMat B(a1.K.rows(), 3);
    B.col(0) = a1.F;
    B.col(1) = a1.F * 2.0;
    B.col(2) = DenseVec::Zero(a1.K.rows());
    ReusableSolver rs2(SolverType::SimplicialLDLT);
    rs2.factorizeReuse(a1.K);
    const DenseMat X = rs2.solveMulti(B);
    check("多右端项 第1列 = 单右端解", (X.col(0) - u1).norm(), 0.0, 1e-12);
    check("多右端项 第2列 = 2×(单右端解)", (X.col(1) - 2.0 * u1).norm(), 0.0, 1e-12);
    check("多右端项 第3列 = 零右端 → 零解", X.col(2).norm(), 0.0, 1e-12);
}

// ---- 17. 梁温度荷载：均匀升温 / 温度梯度 ----
// 温度等效荷载此前只靠"与旧基线逐字节一致"兜底，从未用闭式解验证过符号与数值。
// 这里补两个最基础、最无歧义的解析解，锁定均匀升温轴力与温度梯度反拱。
void testBeamTemperature() {
    std::printf("\n[17] 梁温度荷载\n");
    const double alpha = 1.2e-5;  // 热膨胀系数
    const double h = 0.5;         // 梁高
    const double dT = 30.0;       // 温升 / 温差

    // 17a. 两端固定梁均匀升温：轴向被约束 → 轴力 N = -E·A·α·ΔT（受压）
    {
        Model2D m = straightBeam(L, 1, E, A, I);
        m.elements[0].alpha = alpha;
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({1, 3, 0.0});
        m.constraints.push_back({2, 1, 0.0});
        m.constraints.push_back({2, 2, 0.0});
        m.constraints.push_back({2, 3, 0.0});
        m.uniformTemp = dT;

        const ModelResult r = solveModel2D(m);
        const auto f = elementForces2D(m, r.U).front();
        check("两端固定 均匀升温 N=-EAαΔT", f.N1, -E * A * alpha * dT);
        check("两端固定 均匀升温 N2=N1", f.N2, -E * A * alpha * dT);
    }

    // 17b. 简支梁温度梯度：端弯矩应为 0，跨中反拱 |δ|=α·ΔT·L²/(8h)
    {
        const int n = 4;
        Model2D m = straightBeam(L, n, E, A, I);
        for (auto& e : m.elements) { e.alpha = alpha; e.h = h; }
        m.constraints.push_back({1, 1, 0.0});
        m.constraints.push_back({1, 2, 0.0});
        m.constraints.push_back({n + 1, 2, 0.0});
        for (int e = 1; e <= n; ++e)
            m.elemTemps.push_back({e, "GRADIENT", dT});

        const ModelResult r = solveModel2D(m);
        const auto forces = elementForces2D(m, r.U);
        const double MT = E * I * alpha * dT / h;   // 温度梯度等效弯矩
        double mMax = 0.0;
        for (const auto& f : forces)
            mMax = std::max(mMax, std::max(std::abs(f.M1), std::abs(f.M2)));
        check("温度梯度 端弯矩≈0", mMax / MT, 0.0, 1e-8);

        const int mid = n / 2 + 1;
        check("温度梯度 跨中反拱 |δ|=αΔTL²/(8h)",
              std::abs(nodalDisp(m, r, mid, 2)),
              alpha * dT * L * L / (8.0 * h));
    }
}

// ---- 18. 均布荷载梁端力符号 ----
// 补锁均布荷载下的 V1/V2/M1/M2 符号（此前只在集中力工况下锁定过）。
void testUDLSignConvention() {
    std::printf("\n[18] 均布荷载梁端力符号\n");
    const int n = 4;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({n + 1, 2, 0.0});
    for (int e = 1; e <= n; ++e)
        m.elemLoads.push_back({e, "GY", -q, -q});

    const ModelResult r = solveModel2D(m);
    const auto forces = elementForces2D(m, r.U);
    const double Mscale = q * L * L / 8.0;
    // 左端单元（单元1）的节点1端 = 左支座：M=0，V=+qL/2
    check("左支座 M1=0", forces.front().M1 / Mscale, 0.0, 1e-8);
    check("左支座 V1=+qL/2", forces.front().V1, q * L / 2.0);
    // 右端单元（单元n）的节点2端 = 右支座：M=0，V=-qL/2
    check("右支座 M2=0", forces.back().M2 / Mscale, 0.0, 1e-8);
    check("右支座 V2=-qL/2", forces.back().V2, -q * L / 2.0);
}

// ---- 19. 梁 P-Δ / 几何刚度（屈曲放大）----
// 悬臂柱：轴压 P + 端部横力 Q。轴压会软化弯曲刚度，使挠度按
// 1/(1-P/Pcr) 放大（Pcr = 悬臂欧拉临界力 π²EI/(4L²)）。
// 线性解无此效应，非线性解（含几何刚度）应有放大。
void testBeamPDelta() {
    std::printf("\n[19] 梁 P-Δ 几何刚度\n");
    const int n = 10;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({1, 3, 0.0});

    const double pi = 3.14159265358979323846;
    const double Pcr = pi * pi * E * I / (4.0 * L * L);
    const double P = 0.25 * Pcr;
    const double Q = 1.0e3;

    Model2D mlin = m;
    mlin.loads.push_back({n + 1, 1, -P});
    mlin.loads.push_back({n + 1, 2, -Q});
    mlin.analysisType = "LINEAR";
    const double dLin = nodalDisp(mlin, solveModel2D(mlin), n + 1, 2);

    Model2D mnl = m;
    mnl.loads.push_back({n + 1, 1, -P});
    mnl.loads.push_back({n + 1, 2, -Q});
    mnl.analysisType = "NONLINEAR";
    mnl.nonlinearLoadSteps = 4;
    const double dNl = nodalDisp(mnl, solveModel2D(mnl), n + 1, 2);

    check("挠度放大系数 ≈ 1/(1-P/Pcr)", dNl / dLin, 1.0 / (1.0 - P / Pcr), 0.03);
    if (std::abs(dNl) > std::abs(dLin)) {
        ++g_pass;
        std::printf("  [PASS] P-Δ 使挠度增大（受压软化）  dNl/dLin=%.4f\n", dNl / dLin);
    } else {
        ++g_fail;
        std::printf("  [FAIL] P-Δ 未增大挠度  dLin=%g dNl=%g\n", dLin, dNl);
    }
}

// ---- 20. 梁端部释放（铰）----
// 一端固定、一端铰支（moment release）的固支-铰支梁（propped cantilever）。
// 均布荷载 q 下解析解：固定端弯矩 |M1|=qL²/8，铰端 M2=0。
void testBeamEndRelease() {
    std::printf("\n[20] 梁端部释放（铰）\n");
    Model2D m = straightBeam(L, 1, E, A, I);
    m.elements[0].releases.push_back(5);   // 释放 j 端 θ（局部自由度 5 = 铰）
    m.constraints.push_back({1, 1, 0.0});  // 固定端 u
    m.constraints.push_back({1, 2, 0.0});  // 固定端 v
    m.constraints.push_back({1, 3, 0.0});  // 固定端 θ
    m.constraints.push_back({2, 1, 0.0});  // 铰支端 u
    m.constraints.push_back({2, 2, 0.0});  // 铰支端 v
    for (int e = 1; e <= 1; ++e)
        m.elemLoads.push_back({e, "GY", -q, -q});

    const ModelResult r = solveModel2D(m);
    const auto f = elementForces2D(m, r.U).front();
    check("固支-铰支梁 固定端 |M1|=qL²/8", std::abs(f.M1), q * L * L / 8.0);
    check("固支-铰支梁 铰端 M2=0", f.M2, 0.0, 1e-8);
}

// ---- 21. Timoshenko 梁（剪切变形）----
// 深梁（h/L=0.4，剪切显著）：简支梁跨中集中力挠度
//   δ = PL³/(48EI) + PL/(4G·As)，G=E/(2(1+ν))，As=K·A。
void testBeamTimoshenko() {
    std::printf("\n[21] Timoshenko 梁（剪切变形）\n");
    const double Ld = 1.0;           // 短跨深梁
    const double b = 0.1, h = 0.4;   // 矩形截面
    const double Ad = b * h;         // A
    const double Id = b * h * h * h / 12.0; // I
    const double nu = 0.3, K = 5.0 / 6.0;
    const double Pd = 1.0e6;         // 集中力

    const int n = 4;
    Model2D m = straightBeam(Ld, n, E, Ad, Id);
    for (auto& e : m.elements) { e.nu = nu; e.K = K; }
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({n + 1, 2, 0.0});
    const int mid = n / 2 + 1;
    m.loads.push_back({mid, 2, -Pd});

    const ModelResult r = solveModel2D(m);
    const double G = E / (2.0 * (1.0 + nu));
    const double As = K * Ad;
    const double exact = -Pd * Ld * Ld * Ld / (48.0 * E * Id) - Pd * Ld / (4.0 * G * As);
    check("深梁跨中挠度（含剪切）", nodalDisp(m, r, mid, 2), exact);
}

// ---- 22. 3D 梁内力包络与应力符号 ----
// 悬臂梁（固定 i 端、自由 j 端）+ 向下的端部集中力：
//   固定端截面弯矩（下缘受拉正）应为 -P·L（上缘受拉），自由端为 0。
//   该负值必须被包络正确捕获（验证"0 初值吞负值"已修），且上缘应力为 +PL/Wtop。
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

// ---- 23. 3D 梁扭转剪应力 ----
// 悬臂梁端部施加绕 x 轴扭矩 Mx：梁内扭矩恒为 Mx，扭转剪应力 τ=T/Wt=Mx/Wt。
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

// ---- 24. 自振频率（质量矩阵 + 广义特征值）----
// 简支梁基频 f1 = π/(2L²)·√(EI/(ρA))（欧拉-伯努利梁闭式解）。
void testNaturalFrequency() {
    std::printf("\n[24] 自振频率（质量矩阵 + 特征值）\n");
    const double dens = 7850.0;
    const int n = 16;
    Model2D m = straightBeam(L, n, E, A, I);
    for (auto& e : m.elements) e.dens = dens;
    m.constraints.push_back({1, 1, 0.0});
    m.constraints.push_back({1, 2, 0.0});
    m.constraints.push_back({n + 1, 2, 0.0});

    const ModalResult2D mr = computeNaturalFrequencies2D(m, 3);
    const double rhoA = dens * A;
    const double f1 = kPi / (2.0 * L * L) * std::sqrt(E * I / rhoA);
    check("简支梁基频 f1=π/(2L²)√(EI/ρA)", mr.frequencies[0], f1, 1e-3);
    // 二阶频率应为 4 倍基频（简支梁 f_n = n²·f1）
    check("简支梁二阶频率 f2=4·f1", mr.frequencies[1], 4.0 * f1, 1e-3);
}

// ---- 25. 线性屈曲（欧拉临界力）----
// 两端铰支柱轴压：欧拉临界力 Pcr = π²EI/L²。
void testBuckling() {
    std::printf("\n[25] 线性屈曲（欧拉临界力）\n");
    const int n = 16;
    Model2D m = straightBeam(L, n, E, A, I);
    m.constraints.push_back({1, 1, 0.0});       // 左端 u=0
    m.constraints.push_back({1, 2, 0.0});       // 左端 v=0
    m.constraints.push_back({n + 1, 2, 0.0});   // 右端 v=0（滚轴）
    m.loads.push_back({n + 1, 1, -1.0});        // 参考轴压 P_ref=1

    const auto bf = computeBucklingFactors2D(m, 1);
    const double Pcr = kPi * kPi * E * I / (L * L);
    check("简支柱屈曲 Pcr=π²EI/L²", bf[0], Pcr, 1e-3);
}

// ---- 26. 冲击系数（按基频）----
void testImpactCoefficient() {
    std::printf("\n[26] 冲击系数（按基频）\n");
    check("冲击系数 f=5Hz", impactCoefficientFromFrequency(5.0),
          0.1767 * std::log(5.0) - 0.0157, 1e-12);
    check("冲击系数 f=1Hz（下限 0.05）", impactCoefficientFromFrequency(1.0), 0.05, 1e-12);
    check("冲击系数 f=20Hz（上限 0.45）", impactCoefficientFromFrequency(20.0), 0.45, 1e-12);
}

} // namespace

int main() {
    setupConsoleUTF8();
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // 崩溃时也能看到已输出的内容
    std::printf("========== 解析解断言测试 ==========\n");

    try {
        testSimplySupportedPointLoad();
        testCantileverPointLoad();
        testSimplySupportedUDL();
        testCantileverUDL();
        testBeamSignConvention();
        testAxialBar();
        testContinuousBeam();
        test3DCantilever();
        test3DAxialBar();
        test3DSupportSettlement();
        test3DMixedBeamTruss();
        test3DInfluenceSurfaceTransverse();
        testCreepAEMM();
        testNonlinearTwoBarTruss();
        testSolverReuse();
        testBeamTemperature();
        testUDLSignConvention();
        testBeamPDelta();
        testBeamEndRelease();
        testBeamTimoshenko();
        test3DBeamEnvelope();
        test3DTorsion();
        testNaturalFrequency();
        testBuckling();
        testImpactCoefficient();
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
