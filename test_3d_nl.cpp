// ============================================================
// 3D 几何非线性（共转列式 CR）解析解断言测试
//
// 这些用例是 P3 的"地基"：CR 的 B 矩阵与几何刚度一旦写错，
// 静力精度、成桥态索力、屈曲临界荷载会一起错，而且错得不明显。
// 因此这里全部对标**闭式解**，而不是"跑通不崩"。
// ============================================================

#include "fem_3d_nl.hpp"
#include "test_common.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace fem;

namespace {

constexpr double PI = 3.14159265358979323846;

// 单根梁的最小模型（供单元级测试用）
Model3D oneBeam(const Vec3& p1, const Vec3& p2, double E, double A,
                double Iy, double Iz, double J, double G) {
    Model3D m;
    m.nodes[1] = p1;
    m.nodes[2] = p2;
    Element3D e;
    e.id = 1; e.type = "BEAM";
    e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz; e.J = J; e.G = G;
    e.n1 = 1; e.n2 = 2;
    m.elements.push_back(e);
    return m;
}

// ---------------------------------------------------------------
// 1. CR 切线刚度在参考态必须等于线性刚度（验证 B 矩阵推导）
//    取斜梁：一旦坐标变换或 B 的符号写反，这里立刻暴露。
// ---------------------------------------------------------------
void testCRTangentEqualsLinear() {
    std::printf("\n[CR 切线刚度 = 线性刚度（斜梁）]\n");
    const Vec3 p1(0.0, 0.0, 0.0), p2(3.0, 4.0, 5.0);
    const double E = 2.1e11, A = 0.02, Iy = 3e-4, Iz = 2e-4, J = 5e-4, G = 8.1e10;
    Model3D m = oneBeam(p1, p2, E, A, Iy, Iz, J, G);
    const DofMap3D dm = buildDofMap3D(m);
    const CRElemData ed = crMakeElemData(m, m.elements[0], dm);

    Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
    Eigen::Matrix<double, 14, 1> fe;
    Eigen::Matrix<double, 14, 14> Ke;
    crElementInternal(ed, ue, fe, true, &Ke);

    const double L = (p2 - p1).norm();
    const Eigen::Matrix<double, 12, 12> Kl =
        beam3DLocalStiffness(E, A, Iy, Iz, J, G, L);
    const Eigen::Matrix<double, 12, 12> T = beam3DTransformation(p1, p2);
    const Eigen::Matrix<double, 12, 12> Klin = T.transpose() * Kl * T;

    const Eigen::Matrix<double, 12, 12> dif = Ke.topLeftCorner<12, 12>() - Klin;
    const double scale = Klin.cwiseAbs().maxCoeff();
    check("CR 切线(u=0) 与线性刚度最大相对差",
          dif.cwiseAbs().maxCoeff() / scale, 0.0, 1e-8);
    // 参考态内力必须为零（以 EA 归一，消除刚度量级带来的量纲影响）
    check("CR 参考态内力范数/EA", fe.head<12>().norm() / (E * A), 0.0, 1e-12);
}

// ---------------------------------------------------------------
// 2. 刚体运动（平动 + 任意大转动）下内力恒为零
//    这是 CR 的立身之本：刚体运动必须从位移场里被精确剔除。
// ---------------------------------------------------------------
void testCRRigidBody() {
    std::printf("\n[CR 刚体运动不变性]\n");
    const Vec3 p1(0.0, 0.0, 0.0), p2(2.5, 1.0, -0.5);
    const double E = 2.1e11, A = 0.02, Iy = 3e-4, Iz = 2e-4, J = 5e-4, G = 8.1e10;
    Model3D m = oneBeam(p1, p2, E, A, Iy, Iz, J, G);
    const DofMap3D dm = buildDofMap3D(m);
    const CRElemData ed = crMakeElemData(m, m.elements[0], dm);

    const Vec3 psi(0.7, -0.4, 0.9);          // 约 66° 的有限转动
    const Eigen::Matrix3d R = expMapRot(psi);
    const Vec3 t(0.3, -0.2, 0.15);           // 附加刚体平动

    Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
    // 刚体运动：u_i = (R - I) p_i + t，两节点转角向量同为 ψ
    const Vec3 u1 = (R - Eigen::Matrix3d::Identity()) * p1 + t;
    const Vec3 u2 = (R - Eigen::Matrix3d::Identity()) * p2 + t;
    ue.segment<3>(0) = u1;
    ue.segment<3>(3) = psi;
    ue.segment<3>(6) = u2;
    ue.segment<3>(9) = psi;

    Eigen::Matrix<double, 14, 1> fe;
    crElementInternal(ed, ue, fe, false);
    const double fscale = E * A;              // 内力自然量级
    check("刚体平动+有限转动：内力范数/EA", fe.head<12>().norm() / fscale, 0.0, 1e-8);
}

// ---------------------------------------------------------------
// 3. 纯扭转刚体运动：绕自身轴任意转动，内力为零
// ---------------------------------------------------------------
void testCRRigidTorsion() {
    std::printf("\n[CR 刚体扭转不变性]\n");
    const Vec3 p1(0.0, 0.0, 0.0), p2(0.0, 0.0, 4.0);   // 竖直梁，轴 = z
    const double E = 2.1e11, A = 0.02, Iy = 3e-4, Iz = 2e-4, J = 5e-4, G = 8.1e10;
    Model3D m = oneBeam(p1, p2, E, A, Iy, Iz, J, G);
    const DofMap3D dm = buildDofMap3D(m);
    const CRElemData ed = crMakeElemData(m, m.elements[0], dm);

    Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
    const double alpha = 1.3;                 // 74.5° 自转
    const Vec3 psi(0.0, 0.0, alpha);
    ue.segment<3>(3) = psi;                   // 节点 1 转角
    ue.segment<3>(9) = psi;                   // 节点 2 转角（与节点 1 相同 → 刚体）

    Eigen::Matrix<double, 14, 1> fe;
    crElementInternal(ed, ue, fe, false);
    check("刚体扭转：内力范数/EA", fe.head<12>().norm() / (E * A), 0.0, 1e-8);
}

// ---------------------------------------------------------------
// 4. 轴压梁柱 P-Δ：与 Timoshenko 精确解对比
//    y_max = (Q/(2Pk))·(tan(kL/2) - kL/2)， k = √(P/EI)
//    这是桥梁"梁柱效应"最直接的解析标尺。
// ---------------------------------------------------------------
void testCRBeamColumn() {
    std::printf("\n[CR 梁柱效应 P-Δ vs Timoshenko 精确解]\n");
    const double Lx = 10.0;
    // A 取大值：Timoshenko 梁柱解假定轴线不可伸长，轴向压缩使梁变短会让
    // Pcr 略微提高（本内核含此效应，故更精确）。实测残差 = 0.85/n² + C(A)，
    // C ∝ 1/A；取 A=0.2 时 C≈1.4e-4，可让"精确解"名副其实。
    const double E = 2.1e11, A = 0.2;
    const double Iz = 2e-4;        // 面内弯曲（x-y 平面）
    const double Iy = 5e-2;        // 面外弯曲刚度放大，避免面外屈曲混入
    const double J = 1e-3, G = 8.1e10;
    const int n = 32;              // 离散误差 ≈ 0.85/n² ≈ 8e-4

    Model3D m;
    for (int i = 0; i <= n; ++i) m.nodes[i + 1] = Vec3(i * Lx / n, 0.0, 0.0);
    for (int i = 0; i < n; ++i) {
        Element3D e;
        e.id = i + 1; e.type = "BEAM";
        e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz; e.J = J; e.G = G;
        e.n1 = i + 1; e.n2 = i + 2;
        m.elements.push_back(e);
    }
    // 支座：左端铰（含 θx 抗扭约束），右端滚动（uy/uz 约束，ux 自由以施加轴压）
    for (int d = 1; d <= 4; ++d) m.constraints.push_back({1, d, 0.0, 0.0});   // ux,uy,uz,θx
    m.constraints.push_back({n + 1, 2, 0.0, 0.0});
    m.constraints.push_back({n + 1, 3, 0.0, 0.0});

    const double Pcr = PI * PI * E * Iz / (Lx * Lx);
    const double P = 0.5 * Pcr;
    const double Q = 2.0e4;
    m.loads.push_back({n + 1, 1, -P});                 // 轴压（右端向 -x 推）
    const int mid = n / 2 + 1;
    m.loads.push_back({mid, 2, -Q});                   // 跨中横向力

    NonlinearOptions3D opt;
    opt.loadSteps = 6;
    opt.maxIter = 30;
    opt.tolForce = 1e-10;
    const NonlinearResult3D r = solveModel3DNonlinear(m, opt);

    const DofMap3D dm = buildDofMap3D(m);
    const double wy = r.u[dm.dof(mid, 1)];             // 跨中 y 位移（dof 2 → index 1）

    const double k = std::sqrt(P / (E * Iz));
    const double yex = (Q / (2.0 * P * k)) * (std::tan(k * Lx / 2.0) - k * Lx / 2.0);
    const double ylin = Q * Lx * Lx * Lx / (48.0 * E * Iz);

    check("梁柱：跨中挠度（含 P-Δ 放大）", -wy, yex, 2e-3);
    check("梁柱：放大系数（相对线性解）", -wy / ylin, yex / ylin, 2e-3);
    // 几何刚度若缺失，放大系数会退化到 1.0
    check("梁柱：放大系数显著大于 1", -wy / ylin, yex / ylin, 2e-3);

    // 轴力与外荷载平衡（偏差为单元倾斜的二阶几何效应）
    double sumN = 0.0;
    for (std::size_t i = 0; i < r.axial.size(); ++i) sumN += r.axial[i];
    check("梁柱：单元平均轴力 = -P", sumN / n, -P, 1e-3);
}

// ---------------------------------------------------------------
// 5. 悬臂梁端部弯矩卷成整圆（经典 CR benchmark）
//    M = 2πEI/L 时，梁应闭合成整圆：尖端回到根部，尖端转角 = 2π。
// ---------------------------------------------------------------
void testCRCantileverCircle() {
    std::printf("\n[CR 悬臂梁端弯矩 2πEI/L 卷成整圆]\n");
    const double Lx = 10.0;
    const double E = 2.1e11, A = 0.02, Iy = 5e-2, Iz = 2e-4, J = 1e-2, G = 8.1e10;
    const int n = 8;

    Model3D m;
    for (int i = 0; i <= n; ++i) m.nodes[i + 1] = Vec3(i * Lx / n, 0.0, 0.0);
    for (int i = 0; i < n; ++i) {
        Element3D e;
        e.id = i + 1; e.type = "BEAM";
        e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz; e.J = J; e.G = G;
        e.n1 = i + 1; e.n2 = i + 2;
        m.elements.push_back(e);
    }
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0, 0.0});   // 固定端

    const double M = 2.0 * PI * E * Iz / Lx;
    m.loads.push_back({n + 1, 6, M});                  // 尖端绕 z 的弯矩

    NonlinearOptions3D opt;
    opt.loadSteps = 40;
    opt.maxIter = 40;
    opt.tolForce = 1e-8;
    const NonlinearResult3D r = solveModel3DNonlinear(m, opt);

    const DofMap3D dm = buildDofMap3D(m);
    const int tip = n + 1;
    const Vec3 utip(r.u[dm.dof(tip, 0)], r.u[dm.dof(tip, 1)], r.u[dm.dof(tip, 2)]);
    const double rot = r.u[dm.dof(tip, 5)];
    const Vec3 pos = Vec3(Lx, 0.0, 0.0) + utip;        // 尖端变形后位置

    check("卷圆：尖端回到根部（位置范数）", pos.norm(), 0.0, 5e-3);
    check("卷圆：尖端 x 位移 = -L", utip.x(), -Lx, 5e-3);
    check("卷圆：尖端 y 位移 = 0", utip.y(), 0.0, 5e-3);
    check("卷圆：尖端转角 = 2π", std::abs(rot), 2.0 * PI, 5e-3);
    // 内力矩应与外加弯矩平衡
    double mom = 0.0;
    for (std::size_t i = 0; i < r.localForces.size(); ++i)
        mom = r.localForces[i][11];                     // 单元 j 端绕 z 的局部弯矩
    check("卷圆：单元内力矩 = 2πEI/L", std::abs(mom), M, 1e-2);
}

// ---------------------------------------------------------------
// 6. 小荷载退化：CR 结果必须收敛到线性解
// ---------------------------------------------------------------
void testCRSmallLoadMatchesLinear() {
    std::printf("\n[CR 小荷载退化到线性解]\n");
    const double Lx = 6.0;
    const double E = 2.1e11, A = 0.02, Iy = 5e-2, Iz = 3e-4, J = 1e-2, G = 8.1e10;
    const int n = 6;

    Model3D m;
    for (int i = 0; i <= n; ++i) m.nodes[i + 1] = Vec3(i * Lx / n, 0.0, 0.0);
    for (int i = 0; i < n; ++i) {
        Element3D e;
        e.id = i + 1; e.type = "BEAM";
        e.E = E; e.A = A; e.Iy = Iy; e.Iz = Iz; e.J = J; e.G = G;
        e.n1 = i + 1; e.n2 = i + 2;
        m.elements.push_back(e);
    }
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0, 0.0});
    const double Pt = 5.0e3;                            // 小集中力
    m.loads.push_back({n + 1, 2, -Pt});

    NonlinearOptions3D opt;
    opt.loadSteps = 3;
    opt.maxIter = 20;
    opt.tolForce = 1e-10;
    const NonlinearResult3D r = solveModel3DNonlinear(m, opt);
    const std::vector<double> ulin = solveModel3D(m);

    const DofMap3D dm = buildDofMap3D(m);
    const int tip = n + 1;
    const double wtip = r.u[dm.dof(tip, 1)];
    const double wlin = ulin[dm.dof(tip, 1)];
    const double wth = -Pt * Lx * Lx * Lx / (3.0 * E * Iz);

    check("小荷载：CR 尖端挠度 vs 线性内核", wtip, wlin, 1e-5);
    check("小荷载：CR 尖端挠度 vs 悬臂理论解", wtip, wth, 1e-4);
}

// ---------------------------------------------------------------
// 7. 索单元：无应力长度 L0 与"只受拉"行为
//    （成桥态索长迭代的基础：给定 L0 才能得到正确的成桥索力）
// ---------------------------------------------------------------
void testCableUnstrainedLength() {
    std::printf("\n[索单元：无应力长度与只受拉]\n");
    const double E = 1.95e11, A = 1.0e-3, Lg = 20.0;

    auto run = [&](const std::string& type, double L0) {
        Model3D m;
        m.nodes[1] = Vec3(0.0, 0.0, 0.0);
        m.nodes[2] = Vec3(Lg, 0.0, 0.0);
        Element3D e;
        e.id = 1; e.type = type; e.E = E; e.A = A; e.L0 = L0;
        e.n1 = 1; e.n2 = 2;
        m.elements.push_back(e);
        for (int d = 1; d <= 3; ++d) {
            m.constraints.push_back({1, d, 0.0, 0.0});
            m.constraints.push_back({2, d, 0.0, 0.0});
        }
        NonlinearOptions3D opt;
        opt.loadSteps = 1; opt.maxIter = 10; opt.tolForce = 1e-10;
        return solveModel3DNonlinear(m, opt).axial[0];
    };

    // L0 = 0.9L：索被预拉（无应力长度短于几何长度 → 张紧）
    const double Nten = run("CABLE", 0.9 * Lg);
    check("索：L0=0.9L 张紧轴力 = EA·(L/L0-1)", Nten, E * A * (1.0 / 0.9 - 1.0), 1e-6);
    // L0 = 1.1L：松弛 → 只受拉单元不传力
    const double Nslack = run("CABLE", 1.1 * Lg);
    check("索：L0=1.1L 松弛轴力 = 0", Nslack, 0.0, 1e-9);
    // 对照：TRUSS 同样 L0 会受压（说明只受拉逻辑确实是 CABLE 专属）
    const double Ntruss = run("TRUSS", 1.1 * Lg);
    check("桁架：L0=1.1L 轴力为压（对照）", Ntruss, E * A * (1.0 - 1.1) / 1.1, 1e-6);
}

}  // namespace

int main() {
    setupConsoleUTF8();
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("========== 3D 几何非线性（共转列式）解析解断言测试 ==========\n");
    try {
        testCRTangentEqualsLinear();
        testCRRigidBody();
        testCRRigidTorsion();
        testCRBeamColumn();
        testCRCantileverCircle();
        testCRSmallLoadMatchesLinear();
        testCableUnstrainedLength();
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
