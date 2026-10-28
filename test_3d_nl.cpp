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
#include <ctime>
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

// ---------------------------------------------------------------
// 8. 缆索垂度：与**精确悬链线**独立对照（不是用近似验近似）
//    不可伸长时：弧长 s = (2H/w)·sinh(w·Lh/(2H))
//    可伸长时：弧长 s = L0·(1 + H/(EA))，再代入上式
// ---------------------------------------------------------------
void testCableSagVsExactCatenary() {
    std::printf("\n[缆索垂度 vs 精确悬链线]\n");
    // 精确悬链线：由弧长 s 与水平投影 Lh 反解水平力 H（两端等高时弦 = Lh）
    auto exactH = [](double s, double Lh, double w) {
        // g(H) = (2H/w)·sinh(w·Lh/(2H)) - s = 0，H ∈ (0, ∞)
        double H = w * Lh * Lh / (8.0 * std::sqrt(3.0 * Lh * std::max(s - Lh, 1e-9) / 8.0));
        if (!(H > 0.0) || !std::isfinite(H)) H = w * Lh;
        for (int it = 0; it < 200; ++it) {
            const double u = w * Lh / (2.0 * H);
            const double g = (2.0 * H / w) * std::sinh(u) - s;
            const double gp = (2.0 / w) * std::sinh(u) - (Lh / H) * std::cosh(u);
            const double step = g / gp;
            H -= step;
            if (H < 1e-6) H = 1e-6;
            if (std::abs(step) <= 1e-14 * std::max(1.0, H)) break;
        }
        return H;
    };

    const double Lh = 100.0, w = 100.0;         // 水平跨 100m，重度 100 N/m
    // (a) 近似不可伸长：EA 极大 → 弧长≈L0
    {
        const double L0 = 100.001;              // 弧长比弦长长 1mm → 垂度 ≈ 0.19m
        const double EA = 1e15;
        const double T = cableTensionFromChord(Lh, L0, EA, w, Lh);
        const double Hex = exactH(L0, Lh, w);
        check("索：垂度索力 vs 精确悬链线（不可伸长）", T, Hex, 2e-3);
        // 与抛物线手算式 H = wLh²/(8f)、f = √(3Lh(s-Lh)/8) 对照
        const double f0 = std::sqrt(3.0 * Lh * (L0 - Lh) / 8.0);
        check("索：垂度索力 vs 抛物线公式", T, w * Lh * Lh / (8.0 * f0), 1e-2);
    }
    // (b) 可伸长：EA 有限，弧长 = L0(1 + H/EA)，与代入拉伸弧长的精确悬链线对照
    {
        const double L0 = 100.02, EA = 2.0e8;
        const double T = cableTensionFromChord(Lh, L0, EA, w, Lh);
        // 精确弹性悬链线：g(H) = (2H/w)·sinh(wLh/(2H)) - L0(1+H/EA) = 0
        auto gElas = [&](double hh) {
            return (2.0 * hh / w) * std::sinh(w * Lh / (2.0 * hh))
                   - L0 * (1.0 + hh / EA);
        };
        double lo = 1.0e3, hi = 1.0e7;            // g(lo) > 0 > g(hi)
        for (int i = 0; i < 200; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (gElas(mid) > 0.0) lo = mid; else hi = mid;
        }
        const double H = 0.5 * (lo + hi);
        check("索：垂度索力 vs 精确悬链线（可伸长）", T, H, 5e-3);
    }
    // (c) Ernst 等效模量 = 力-伸长关系的切线（保证内力与切线共轭）
    {
        const double L0 = 100.01, EA = 2.0e8, Ac = 1.0e-3, E = EA / Ac;
        const double T = cableTensionFromChord(Lh, L0, EA, w, Lh);
        // 固定水平投影 Lh、只变弦长 Ln 求偏导（这正是 Ernst 式的定义）
        // 步长取 1e-6：T(Ln) 曲率极大，截断误差 ~36·h²（相对），
        // 与舍入误差 ~1e-11/h 平衡点在 h≈7e-7，故 1e-6 是最佳量级。
        const double h = 1e-6;
        const double dTdL = (cableTensionFromChord(Lh + h, L0, EA, w, Lh)
                             - cableTensionFromChord(Lh - h, L0, EA, w, Lh)) / (2.0 * h);
        const double Eeq = cableErnstModulus(E, Ac, w, Lh, T, L0);
        check("索：Ernst 模量 = 力-伸长切线", Eeq * Ac / L0, dTdL, 1e-7);
        // 垂度必须让索"变软"：Eeq < E
        check("索：Ernst 模量小于 E（垂度软化）", Eeq < E ? 1.0 : 0.0, 1.0, 0.0);
        std::printf("        Eeq=%.6g Pa < E=%.6g Pa（垂度使索软化 %.1f%%）\n",
                    Eeq, E, 100.0 * (1.0 - Eeq / E));
    }
}

// ---------------------------------------------------------------
// 9. 成桥态：无应力索长迭代（双索，各索目标索力不同）
//    两根竖索 + 悬臂主梁（一维体系，收敛后有闭式解）：
//      梁端竖向平衡：索力 T1+T2 向上，恒载 P 与索自重 wi·L0i/2 向下
//        v = (T1 + T2 - P - Σwi·L0i/2)/kb，  kb = 3EI/L³（悬臂梁端刚度）
//      竖索无水平投影 → 无垂度：Ti = EA·(Hi - v - L0i)/L0i
//    联立即得各索 L0i 与 v 的闭式解，逐项与迭代结果对照。
// ---------------------------------------------------------------
void testCableLengthIteration() {
    std::printf("\n[成桥态：无应力索长迭代（双索）]\n");
    const double Eb = 3.45e10, Ab = 8.0, Ib = 6.0;       // 主梁：A=8m²，I=6m^4
    const double Ec = 1.95e11, Ac = 5.0e-3, wc = 46.0;   // 拉索（竖直 → 无垂度）
    const double Lb = 30.0;                              // 悬臂 30m
    const double H1 = 40.0, H2 = 35.0;                   // 两索索高
    const double T1 = 1.5e6, T2 = 1.0e6;                 // 目标成桥索力
    const double P = 2.45e6;                             // 梁端恒载 2450 kN

    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);        // 梁根（固结）
    m.nodes[2] = Vec3(Lb, 0.0, 0.0);         // 梁端 = 两索下锚点
    m.nodes[3] = Vec3(Lb, 0.0, H1);          // 索 1 上锚点（固结）
    m.nodes[4] = Vec3(Lb, 0.0, H2);          // 索 2 上锚点（固结）
    {
        Element3D g; g.id = 1; g.type = "BEAM";
        g.E = Eb; g.A = Ab; g.Iy = Ib; g.Iz = Ib; g.J = 12.0; g.G = 1.4e10;
        g.n1 = 1; g.n2 = 2; m.elements.push_back(g);
        Element3D c1; c1.id = 10; c1.type = "CABLE";
        c1.E = Ec; c1.A = Ac; c1.w = wc; c1.n1 = 2; c1.n2 = 3;
        m.elements.push_back(c1);
        Element3D c2; c2.id = 11; c2.type = "CABLE";
        c2.E = Ec; c2.A = Ac; c2.w = wc; c2.n1 = 2; c2.n2 = 4;
        m.elements.push_back(c2);
    }
    for (int d = 1; d <= 6; ++d) {
        m.constraints.push_back({1, d, 0.0, 0.0});
        m.constraints.push_back({3, d, 0.0, 0.0});
        m.constraints.push_back({4, d, 0.0, 0.0});
    }
    m.loads.push_back({2, 3, -P});

    // ---- 闭式解（索自重是 L0 的一阶小量，代入后直接解出 v）----
    const double kb = 3.0 * Eb * Ib / (Lb * Lb * Lb);     // 悬臂梁端刚度
    const double EAc = Ec * Ac;
    const double a1 = 1.0 + T1 / EAc, a2 = 1.0 + T2 / EAc;
    const double vth = (T1 + T2 - P - wc * H1 / (2.0 * a1) - wc * H2 / (2.0 * a2))
                       / (kb + wc / (2.0 * a1) + wc / (2.0 * a2));
    const double L01th = (H1 - vth) / a1;
    const double L02th = (H2 - vth) / a2;

    CableLengthOptions3D co;
    co.maxIter = 40;
    co.tolTension = 1e-8;
    co.nl.loadSteps = 4;
    co.nl.maxIter = 30;
    co.nl.tolForce = 1e-11;
    const CableLengthResult3D cl =
        iterateCableUnstrainedLengths(m, {{10, T1}, {11, T2}}, co);

    check("索长迭代：收敛", cl.converged ? 1.0 : 0.0, 1.0, 0.0);

    const DofMap3D dm = buildDofMap3D(m);
    const double vz = cl.u[dm.dof(2, 2)];                 // 梁端竖向位移
    const double Tt[2] = {T1, T2};
    const double Hh[2] = {H1, H2};
    const double L0th[2] = {L01th, L02th};
    for (int i = 0; i < 2; ++i) {
        check("索长迭代：成桥索力 = 目标索力", cl.tension[i], Tt[i], 1e-7);
        check("索长迭代：无应力索长 vs 闭式解", cl.L0[i], L0th[i], 1e-8);
        check("索长迭代：弦长 = 索高 - 梁端上挠", cl.chord[i], Hh[i] - vth, 1e-8);
        // 本构自检：竖索 T = EA·(Ln - L0)/L0
        check("索长迭代：竖索本构 T = EA(Ln-L0)/L0",
              EAc * (cl.chord[i] - cl.L0[i]) / cl.L0[i], Tt[i], 1e-7);
        // 无应力索长必须短于弦长（索被张拉）
        check("索长迭代：无应力索长 < 弦长（已张拉）",
              cl.L0[i] < cl.chord[i] ? 1.0 : 0.0, 1.0, 0.0);
    }
    check("索长迭代：梁端竖向位移 vs 闭式解", vz, vth, 1e-7);
    std::printf("  %d 轮（首轮含 1 次差分探测），索力误差=%.2e\n"
                "  索1：L0=%.5f m，弦长=%.5f m；索2：L0=%.5f m，弦长=%.5f m\n"
                "  梁端上挠=%.4f mm（预抬高迭代的一阶估计）\n",
                cl.iterations, cl.maxTensionErr,
                cl.L0[0], cl.chord[0], cl.L0[1], cl.chord[1], 1.0e3 * vz);
}

// ---------------------------------------------------------------
// 10. 刚性支承连续梁法：自动定目标成桥索力
//     两跨等跨连续梁（跨长 L），各跨跨中一个集中力 P，三处竖索（两端 + 中支
//     点）。把索换成刚性竖向支承后即标准的"两跨连续梁"，反力有教科书闭式解
//     （三弯矩方程 / 力法，两种推导互相印证）：
//        M_B = -3PL/16，  R_B = 11P/8，  R_A = R_C = 5P/16
//     竖索的竖向分量就是索力本身，故 T = R。
//     再叠加一条严格的整体平衡：ΣT = 2P（梁竖向无其他支承）。
// ---------------------------------------------------------------
namespace {
// 三索 + 两跨连续梁（每跨 2 个单元），L 为单跨跨长，H 为索高
Model3D makeThreeCableBeam(double L, double H, double Eb, double A, double I,
                           double Ec, double Ac, double P) {
    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(0.5 * L, 0.0, 0.0);
    m.nodes[3] = Vec3(L, 0.0, 0.0);
    m.nodes[4] = Vec3(1.5 * L, 0.0, 0.0);
    m.nodes[5] = Vec3(2.0 * L, 0.0, 0.0);
    m.nodes[6] = Vec3(0.0, 0.0, H);
    m.nodes[7] = Vec3(L, 0.0, H);
    m.nodes[8] = Vec3(2.0 * L, 0.0, H);
    for (int i = 1; i <= 4; ++i) {
        Element3D g; g.id = i; g.type = "BEAM";
        g.E = Eb; g.A = A; g.Iy = I; g.Iz = I; g.J = 2.0 * I; g.G = 1.4e10;
        g.n1 = i; g.n2 = i + 1; m.elements.push_back(g);
    }
    const int cn[3] = {1, 3, 5};
    for (int k = 0; k < 3; ++k) {
        Element3D c; c.id = 10 + k; c.type = "CABLE";
        c.E = Ec; c.A = Ac; c.w = 0.0;      // 竖索：无垂度、无自重
        c.n1 = cn[k]; c.n2 = 6 + k;
        m.elements.push_back(c);
    }
    // 约束：A 处水平定位 + 防扭转；面外（y、θz）全约束；塔顶锚点固结
    m.constraints.push_back({1, 1, 0.0, 0.0});
    m.constraints.push_back({1, 4, 0.0, 0.0});
    for (int n = 1; n <= 5; ++n) {
        m.constraints.push_back({n, 2, 0.0, 0.0});
        m.constraints.push_back({n, 6, 0.0, 0.0});
    }
    for (int n = 6; n <= 8; ++n)
        for (int d = 1; d <= 6; ++d) m.constraints.push_back({n, d, 0.0, 0.0});
    m.loads.push_back({2, 3, -P});
    m.loads.push_back({4, 3, -P});
    return m;
}
}  // namespace

void testRigidSupportCableTension() {
    std::printf("\n[刚性支承连续梁法：自动定目标索力]\n");
    const double L = 30.0, H = 40.0, P = 1.0e6;
    const Model3D m = makeThreeCableBeam(L, H, 3.45e10, 6.0, 4.0, 1.95e11, 5.0e-3, P);

    const std::vector<double> T =
        cableTensionsRigidSupport(m, {{10, 1}, {11, 3}, {12, 5}});

    // 闭式解（三弯矩方程）：R_A = R_C = 5P/16，R_M = 11P/8
    check("刚性支承法：中支点索力 = 11P/8", T[1], 11.0 * P / 8.0, 1e-9);
    check("刚性支承法：端支点索力 = 5P/16", T[0], 5.0 * P / 16.0, 1e-9);
    check("刚性支承法：另一端索力 = 5P/16", T[2], 5.0 * P / 16.0, 1e-9);
    check("刚性支承法：索力之和 = 总恒载 2P", T[0] + T[1] + T[2], 2.0 * P, 1e-9);
    std::printf("  T端=%.6g kN，T中=%.6g kN，ΣT=%.6g kN（恒载 %.6g kN）\n",
                T[0] / 1e3, T[1] / 1e3, (T[0] + T[1] + T[2]) / 1e3, 2.0 * P / 1e3);
}

// ---------------------------------------------------------------
// 11. 成桥态闭环：单索悬臂梁（有完整闭式解）
//     梁根固结、梁端竖向恒载 P、竖索到上方固定点，目标成桥索力 T*。
//       · 下料索长：L0 = Hc/(1 + T*/EA)（竖索无垂度，闭式）
//       · 梁端平衡：T* + kb·z_m - P = 0，kb = 3EI/L³
//         ⇒ 预抬高 z_m = (P - T*)/kb   （索力偏大 → 制造时把梁端放低）
//     收敛后成桥位置必须回到设计线形：弦长 = Hc、索力 = T*。
// ---------------------------------------------------------------
void testCompletionStateSingleCable() {
    std::printf("\n[成桥态闭环：单索悬臂梁 vs 闭式解]\n");
    const double Eb = 3.45e10, Ab = 6.0, Ib = 6.0;
    const double Ec = 1.95e11, Ac = 5.0e-3;
    const double Lb = 30.0, Hc = 40.0;
    const double P = 2.4e6, Tt = 2.5e6;

    Model3D m;
    m.nodes[1] = Vec3(0.0, 0.0, 0.0);
    m.nodes[2] = Vec3(Lb, 0.0, 0.0);
    m.nodes[3] = Vec3(Lb, 0.0, Hc);
    {
        Element3D g; g.id = 1; g.type = "BEAM";
        g.E = Eb; g.A = Ab; g.Iy = Ib; g.Iz = Ib; g.J = 12.0; g.G = 1.4e10;
        g.n1 = 1; g.n2 = 2; m.elements.push_back(g);
        Element3D c; c.id = 10; c.type = "CABLE";
        c.E = Ec; c.A = Ac; c.w = 0.0;
        c.n1 = 2; c.n2 = 3; m.elements.push_back(c);
    }
    for (int d = 1; d <= 6; ++d) {
        m.constraints.push_back({1, d, 0.0, 0.0});
        m.constraints.push_back({3, d, 0.0, 0.0});
    }
    m.loads.push_back({2, 3, -P});

    CompletionStateOptions3D co;
    co.maxIter = 20;
    co.tolClosure = 1e-10;
    co.nl.loadSteps = 2;
    co.nl.maxIter = 40;
    co.nl.tolForce = 1e-11;
    const CompletionStateResult3D cs = solveCompletionState3D(m, {{10, 2, Tt}}, co);

    check("成桥态：收敛（几何闭合）", cs.converged ? 1.0 : 0.0, 1.0, 0.0);

    const DofMap3D dm = buildDofMap3D(m);
    const int dz = dm.dof(2, 2);                 // 梁端竖向自由度
    const double kb = 3.0 * Eb * Ib / (Lb * Lb * Lb);
    const double zm = (P - Tt) / kb;             // 闭式预抬高

    check("成桥态：成桥索力 = 目标索力", cs.tension[0], Tt, 1e-8);
    check("成桥态：成桥弦长 = 设计弦长（几何闭合）", cs.chord[0], Hc, 1e-9);
    check("成桥态：下料索长 = Hc/(1+T*/EA)", cs.L0[0], Hc / (1.0 + Tt / (Ec * Ac)), 1e-12);
    check("成桥态：预抬高 vs 闭式解 (P-T*)/kb", cs.camber[dz], zm, 1e-6);
    // 位移与预抬高严格反号：制造线形 + 位移 = 设计线形
    check("成桥态：预抬高 + 位移 = 0", cs.camber[dz] + cs.u[dz], 0.0, 1e-10);
    check("成桥态：下料索长 < 设计弦长（已张拉）",
          cs.L0[0] < cs.Ltarget[0] ? 1.0 : 0.0, 1.0, 0.0);
    std::printf("  %d 轮几何迭代；L0=%.5f m，弦长=%.5f m，"
                "梁端预抬高=%.4f mm（闭式 %.4f mm）\n",
                cs.iterations, cs.L0[0], cs.chord[0],
                1.0e3 * cs.camber[dz], 1.0e3 * zm);
}

// ---------------------------------------------------------------
// 12. 成桥态闭环：三索两跨连续梁（索力自动 + 预抬高）
//     索力由刚性支承连续梁法自动给定，再解制造线形。独立检验两条：
//       · 几何闭合（成桥位置 = 设计线形，弦长回到设计值）
//       · 整体平衡：梁竖向无其他支承，故 ΣT = 2P 必须在成桥态仍成立
// ---------------------------------------------------------------
void testCompletionStateThreeCable() {
    std::printf("\n[成桥态闭环：三索两跨连续梁（自动定索力）]\n");
    const double L = 30.0, H = 40.0, P = 1.0e6;
    Model3D m = makeThreeCableBeam(L, H, 3.45e10, 6.0, 4.0, 1.95e11, 5.0e-3, P);

    CompletionStateOptions3D co;
    co.maxIter = 20;
    co.tolClosure = 1e-10;
    co.nl.loadSteps = 2;
    co.nl.maxIter = 40;
    co.nl.tolForce = 1e-11;
    // tension = 0 → 由刚性支承连续梁法自动定
    const CompletionStateResult3D cs =
        solveCompletionState3D(m, {{10, 1, 0.0}, {11, 3, 0.0}, {12, 5, 0.0}}, co);

    check("三索成桥态：收敛", cs.converged ? 1.0 : 0.0, 1.0, 0.0);
    check("三索成桥态：索力 = 刚性支承法目标值",
          cs.tension[1], 11.0 * P / 8.0, 1e-7);
    check("三索成桥态：整体平衡 ΣT = 2P",
          cs.tension[0] + cs.tension[1] + cs.tension[2], 2.0 * P, 1e-7);
    for (int k = 0; k < 3; ++k)
        check("三索成桥态：弦长回到设计值（索 " + std::to_string(k + 1) + "）",
              cs.chord[k], cs.Ltarget[k], 1e-9);
    std::printf("  %d 轮；T=(%.4g, %.4g, %.4g) kN，闭合误差=%.2e m\n",
                cs.iterations, cs.tension[0] / 1e3, cs.tension[1] / 1e3,
                cs.tension[2] / 1e3, cs.maxClosureErr);
}

// ---------------------------------------------------------------
// 13. 索影响矩阵：解析式 J = M·D 必须等于"逐根扰动 + 完整重算"的黄金标准。
//     这条是用以确保"提速不损失正确性"的硬约束：P3 提速把建立 J 的代价从
//     nc 次非线性求解降为 1 次分解 + nc 次回代，二者必须给出同一个 J。
//     三索模型索间耦合最强（中索一动，两根边索都跟着变），是理想试金石。
// ---------------------------------------------------------------
void testCableInfluenceMatrix() {
    std::printf("\n[索影响矩阵：解析 J = M·D vs 逐根扰动重算]\n");
    const double L = 30.0, H = 40.0, P = 1.0e6;
    Model3D m = makeThreeCableBeam(L, H, 3.45e10, 6.0, 4.0, 1.95e11, 5.0e-3, P);
    const std::vector<int> ids = {10, 11, 12};
    const std::vector<double> Ttar = {5.0 * P / 16.0, 11.0 * P / 8.0, 5.0 * P / 16.0};
    for (std::size_t k = 0; k < 3; ++k) {
        for (auto& e : m.elements) {
            if (e.id != ids[k]) continue;
            const Vec3 d3 = m.nodes.at(e.n2) - m.nodes.at(e.n1);
            e.L0 = cableUnstrainedLengthForTension(d3.norm(), Ttar[k], e.E * e.A, e.w,
                                                   std::sqrt(d3.x() * d3.x() + d3.y() * d3.y()));
        }
    }

    NonlinearOptions3D nl;
    nl.loadSteps = 2; nl.maxIter = 40; nl.tolForce = 1e-11;
    NLTangent3D tan;
    const NonlinearResult3D r0 = solveModel3DNonlinear(m, nl, &tan);
    std::vector<double> T0(3, 0.0);
    for (int k = 0; k < 3; ++k) {
        std::size_t i = 0;
        while (r0.elemIds[i] != ids[k]) ++i;
        T0[k] = r0.axial[i];
    }

    const DofMap3D dm = buildDofMap3D(m);
    const int nd = dm.nDofs;
    std::vector<char> isC(nd, 0), isT(nd, 0);
    for (const auto& c : m.constraints) {
        if (c.k > 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0) isC[d] = 1;
    }
    for (std::size_t i = 0; i < dm.nodeIds.size(); ++i)
        for (int a = 0; a < std::min(dm.dofCount[i], 3); ++a) isT[dm.dofStart[i] + a] = 1;
    for (int d = 0; d < nd; ++d) if (isC[d]) isT[d] = 0;

    const CableSensitivity3D s = cableSensitivity3D(m, dm, ids, r0.u, T0, isT, tan);
    const Eigen::MatrixXd J = s.M * s.D.asDiagonal();

    // 黄金标准：逐根扰动 L0_j，完整重算，取 (T_i(L0_j+δ) - T_i)/δ
    const double relp = 1e-4;
    double maxRel = 0.0;
    for (int j = 0; j < 3; ++j) {
        Model3D mj = m;
        double dL = 0.0;
        for (auto& e : mj.elements)
            if (e.id == ids[j]) { dL = -relp * e.L0; e.L0 += dL; break; }
        const NonlinearResult3D rj = solveModel3DNonlinear(mj, nl, nullptr);
        for (int i = 0; i < 3; ++i) {
            std::size_t q = 0;
            while (rj.elemIds[q] != ids[i]) ++q;
            const double Jnum = (rj.axial[q] - T0[i]) / dL;
            const double den = std::max(std::abs(Jnum), std::abs(J(i, j)));
            maxRel = std::max(maxRel, std::abs(Jnum - J(i, j)) / std::max(den, 1e-12));
        }
    }
    // 一阶影响矩阵：扰动量 1e-4 下的二阶残差应在 1e-3 以内
    check("影响矩阵：解析 J = M·D vs 扰动重算（最大相对差）", maxRel, 0.0, 2e-3);
    // 对角元必须是 dT/dL0 的量级且为负（索做短 → 索力增大）
    check("影响矩阵：对角元为负（L0 缩短 → 索力增大）",
          (J(0, 0) < 0.0 && J(1, 1) < 0.0 && J(2, 2) < 0.0) ? 1.0 : 0.0, 1.0, 0.0);
    // 非对角元非零 = 索间耦合确实存在（否则退化成对角割线，收敛会变慢）
    check("影响矩阵：非对角元非零（索间耦合）",
          std::abs(J(0, 1)) > 1e-8 * std::abs(J(1, 1)) ? 1.0 : 0.0, 1.0, 0.0);
    std::printf("  J 对角=(%.4g, %.4g, %.4g) N/m；最大相对差=%.2e\n",
                J(0, 0), J(1, 1), J(2, 2), maxRel);
}

// ---------------------------------------------------------------
// 10. 解析切线 vs 差分切线（黄金基线）
//
// 解析式（P3 提速）必须与"逐自由度中心差分"逐元素一致，否则 Newton 的
// 二次收敛就没了——而且这种错误不会体现在任何静力结果里，只体现在迭代
// 次数上，极难察觉。因此把它钉成断言，三种单元各比一次：
//   · BEAM(12)   —— 材料项 + 几何项（含平分三元组的自旋 ω）
//   · BEAM7(14)  —— 14 维材料项（翘曲—扭转耦合不能靠嵌入代替）
//   · CABLE(6)   —— 含垂度，N = N(Ln, Lh) 的两个偏导
// 比较状态取**有限变形**（弯曲 + 扭转 + 轴向），此时几何刚度非零。
// ---------------------------------------------------------------
namespace {

// 解析 vs 差分：返回最大相对差，同时回填两侧最大绝对量级
double tangentGap(const CRElemData& ed, const Eigen::Matrix<double, 14, 1>& ue,
                  int n, double* scaleOut = nullptr) {
    Eigen::Matrix<double, 14, 1> fe;
    Eigen::Matrix<double, 14, 14> Ka, Kf;
    crElementInternal(ed, ue, fe, true, &Ka, nullptr, nullptr, false);
    crElementInternal(ed, ue, fe, true, &Kf, nullptr, nullptr, true);   // 差分
    double md = 0.0, sc = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            md = std::max(md, std::abs(Ka(i, j) - Kf(i, j)));
            sc = std::max(sc, std::abs(Kf(i, j)));
        }
    if (scaleOut) *scaleOut = sc;
    return md / std::max(sc, 1e-30);
}

}  // namespace

void testAnalyticTangentVsFD() {
    std::printf("\n[解析切线 vs 差分切线（黄金基线）]\n");
    const double E = 2.1e11, A = 0.02, Iy = 3e-4, Iz = 2e-4, J = 5e-4, G = 8.1e10;
    const Vec3 p1(0.0, 0.0, 0.0), p2(3.0, 4.0, 5.0);

    // 有限变形的位移场（弯曲 5e-2 rad、扭转 2e-3 rad、轴向 1.2e-2 m）
    Eigen::Matrix<double, 14, 1> u = Eigen::Matrix<double, 14, 1>::Zero();
    u[0] = 1.2e-2; u[1] = 3.0e-3; u[2] = -2.0e-3;
    u[3] = 2.0e-3; u[4] = 5.0e-2; u[5] = -3.0e-2;
    u[6] = -8.0e-3; u[7] = 5.0e-3; u[8] = 4.0e-3;
    u[9] = -3.0e-3; u[10] = 5.5e-2; u[11] = 2.5e-2;

    // (a) 参考态：切线必须等于线性刚度，解析式与差分都必须命中
    {
        Model3D m = oneBeam(p1, p2, E, A, Iy, Iz, J, G);
        const DofMap3D dm = buildDofMap3D(m);
        const CRElemData ed = crMakeElemData(m, m.elements[0], dm);
        Eigen::Matrix<double, 14, 1> z = Eigen::Matrix<double, 14, 1>::Zero();
        check("解析切线：参考态 vs 差分", tangentGap(ed, z, 12), 0.0, 1e-8);
    }

    // (b) BEAM(12) 有限变形
    {
        Model3D m = oneBeam(p1, p2, E, A, Iy, Iz, J, G);
        const DofMap3D dm = buildDofMap3D(m);
        const CRElemData ed = crMakeElemData(m, m.elements[0], dm);
        double sc = 0.0;
        const double gap = tangentGap(ed, u, 12, &sc);
        check("解析切线：BEAM(12) 大变形 vs 差分", gap, 0.0, 1e-8);
    }

    // (c) BEAM7(14) 有限变形（含翘曲自由度）
    {
        Model3D m = oneBeam(p1, p2, E, A, Iy, Iz, J, G);
        m.elements[0].type = "BEAM7";
        m.elements[0].Iw = 1.0e-8;
        const DofMap3D dm = buildDofMap3D(m);
        const CRElemData ed = crMakeElemData(m, m.elements[0], dm);
        Eigen::Matrix<double, 14, 1> u7 = Eigen::Matrix<double, 14, 1>::Zero();
        for (int i = 0; i < 6; ++i) { u7[i] = u[i]; u7[7 + i] = u[6 + i]; }
        u7[6] = 1.0e-3;
        u7[13] = -1.5e-3;
        check("解析切线：BEAM7(14) 大变形 vs 差分", tangentGap(ed, u7, 14), 0.0, 1e-8);
    }

    // (d) 含垂度的索（N = N(Ln, Lh)，切线含 k2 = ∂N/∂Lh）
    {
        Model3D m;
        m.nodes[1] = Vec3(0.0, 0.0, 0.0);
        m.nodes[2] = Vec3(100.0, 0.0, 40.0);
        Element3D e;
        e.id = 1; e.type = "CABLE"; e.E = 1.95e11; e.A = 6.0e-3; e.w = 50.0;
        e.L0 = 105.0; e.n1 = 1; e.n2 = 2;
        m.elements.push_back(e);
        const DofMap3D dm = buildDofMap3D(m);
        const CRElemData ed = crMakeElemData(m, m.elements[0], dm);
        Eigen::Matrix<double, 14, 1> uc = Eigen::Matrix<double, 14, 1>::Zero();
        uc[0] = 1e-2; uc[1] = -2e-2; uc[2] = 1e-2;
        uc[3] = -1e-3; uc[4] = 2e-3; uc[5] = 1e-3;
        check("解析切线：垂度索 vs 差分", tangentGap(ed, uc, 6), 0.0, 1e-8);
    }
}

// ---------------------------------------------------------------
// 11. 解析切线的提速效果（自检：解析路径不得慢于差分路径）
//     这里是粗粒度计时，只用于防止"优化变慢"的回归；
//     正式对标见 bench/nl_bench。
// ---------------------------------------------------------------
void testAnalyticTangentSpeed() {
    std::printf("\n[解析切线 vs 差分切线（单元级耗时）]\n");
    const double E = 2.1e11, A = 0.02, Iy = 3e-4, Iz = 2e-4, J = 5e-4, G = 8.1e10;
    Model3D m = oneBeam(Vec3(0, 0, 0), Vec3(3, 4, 5), E, A, Iy, Iz, J, G);
    const DofMap3D dm = buildDofMap3D(m);
    const CRElemData ed = crMakeElemData(m, m.elements[0], dm);
    Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
    ue[0] = 1.2e-2; ue[4] = 5.0e-2; ue[10] = 5.5e-2;

    const int N = 20000;
    double sink = 0.0;
    auto timeIt = [&](bool byFD) {
        Eigen::Matrix<double, 14, 1> fe;
        Eigen::Matrix<double, 14, 14> Ke;
        const clock_t t0 = std::clock();
        for (int k = 0; k < N; ++k) {
            crElementInternal(ed, ue, fe, true, &Ke, nullptr, nullptr, byFD);
            sink += Ke(0, 0);
        }
        return static_cast<double>(std::clock() - t0) / CLOCKS_PER_SEC / N * 1e6;
    };
    const double tA = timeIt(false);
    const double tF = timeIt(true);
    if (sink == 1.2345) std::printf("");
    std::printf("  解析 %.2f µs / 差分 %.2f µs → %.1fx\n", tA, tF, tF / tA);
    // 门限放宽（CI 机器噪声大），只要求解析路径不慢于差分的一半
    check("解析切线：单元级耗时 ≤ 差分的一半", tA / tF, 0.0, 0.5);
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
        testCableSagVsExactCatenary();
        testCableLengthIteration();
        testRigidSupportCableTension();
        testCompletionStateSingleCable();
        testCompletionStateThreeCable();
        testAnalyticTangentVsFD();
        testAnalyticTangentSpeed();
        testCableInfluenceMatrix();
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
