// ce_3d_benchmark.cpp —— 3D 约束方程（索鞍顶推口径 C）内核验收
//
// 独立于悬索找形，直接构造模型，覆盖线性与非线性两条求解路径：
//   1. 线性 solveModel3D：三节点级联 CE 命中 Δ、缆鞍固结、固结方向锁定；
//   2. 非线性 solveModel3DNonlinear：同样命中 Δ（守住常数项只进内力侧、
//      不重复进 Fext 的 2b 回归）；
//   3. 无效 dof（0/7）与纯轴向节点转动 dof 被拒绝；
//   4. 重复系数合并正确。
//
// 模型映射（与 suspension_3d 顶推一致）：
//   2=塔顶结构 T（连 BEAM），3=主缆端 C（连 TRUSS），5=索鞍 S（孤立，靠 CE）
//   CE1 UX(S)−UX(T)=Δ；CE2 UX(C)−UX(S)=0；CE3~6 UY/UZ 固结。
//
// 退出码 0 = 全通过；非 0 = 失败（CTest 直接判定）。
#include "fem_3d.hpp"
#include "fem_3d_nl.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace fem;

namespace {

int g_fail = 0;

void expect(const std::string& name, bool ok, const std::string& detail = "") {
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  " << detail;
    std::cout << "\n";
    if (!ok) ++g_fail;
}

bool almostEq(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

// 构造三节点索鞍顶推模型；jack_on=false 时退化为纯梁+桁架（无 CE）。
Model3D saddleModel(double Delta, bool jack_on) {
    Model3D m;
    m.nodes[1] = Vec3(0, 0, 0);
    m.nodes[2] = Vec3(0, 0, 10);
    m.nodes[3] = Vec3(0, 0, 10);
    m.nodes[4] = Vec3(3, 1, 3);   // y 分量≠0，使节点3的 UY 有桁架刚度（jack_on=false 亦自洽）

    Element3D eb; eb.id = 1; eb.type = "BEAM"; eb.n1 = 1; eb.n2 = 2;
    eb.E = 3.45e10; eb.A = 12.0; eb.Iy = 16.0; eb.Iz = 16.0;
    eb.J = 25.0; eb.G = 1.4e10;
    m.elements.push_back(eb);
    Element3D ec; ec.id = 2; ec.type = "TRUSS"; ec.n1 = 3; ec.n2 = 4;
    ec.E = 1.95e11; ec.A = 0.06;
    m.elements.push_back(ec);

    for (int d = 1; d <= 6; ++d) {
        m.constraints.push_back({1, d, 0.0, 0.0});
        m.constraints.push_back({4, d, 0.0, 0.0});
    }

    if (jack_on) {
        m.nodes[5] = Vec3(0, 0, 10);   // 索鞍 S：孤立节点，仅顶推时存在，靠 CE 罚刚度
        auto add = [&](int id, double b, int n1, int d1, double c1,
                       int n2, int d2, double c2) {
            Model3D::ConstraintEquation ce;
            ce.id = id; ce.constant = b;
            ce.terms.push_back({n1, d1, c1});
            ce.terms.push_back({n2, d2, c2});
            m.ces.push_back(ce);
        };
        add(1, Delta, 5, 1, 1.0, 2, 1, -1.0);   // UX(S)-UX(T)=Δ
        add(2, 0.0, 3, 1, 1.0, 5, 1, -1.0);     // UX(C)-UX(S)=0
        add(3, 0.0, 5, 2, 1.0, 2, 2, -1.0);     // UY(S)-UY(T)=0
        add(4, 0.0, 5, 3, 1.0, 2, 3, -1.0);     // UZ(S)-UZ(T)=0
        add(5, 0.0, 3, 2, 1.0, 2, 2, -1.0);     // UY(C)-UY(T)=0
        add(6, 0.0, 3, 3, 1.0, 2, 3, -1.0);     // UZ(C)-UZ(T)=0
        m.cePenalty = 1e12;
    }
    return m;
}

void checkPath(const std::string& tag, bool nonlinear, double Delta) {
    Model3D m = saddleModel(Delta, true);
    std::vector<double> U;
    if (nonlinear) {
        NonlinearOptions3D opt;
        opt.loadSteps = 1;
        opt.maxIter = 50;
        opt.tolForce = 1e-10;
        opt.tolDisp = 1e-12;
        NonlinearResult3D r = solveModel3DNonlinear(m, opt);
        expect(tag + " 非线性收敛", r.converged,
               "iter=" + std::to_string(r.totalIterations));
        U = r.u;
    } else {
        U = solveModel3D(m);
    }
    const DofMap3D dm = buildDofMap3D(m);
    const double uxT = U[dm.safeDof(2, 0)];
    const double uxC = U[dm.safeDof(3, 0)];
    const double uxS = U[dm.safeDof(5, 0)];
    const double uyC = U[dm.safeDof(3, 1)], uyT = U[dm.safeDof(2, 1)];
    const double uzC = U[dm.safeDof(3, 2)], uzT = U[dm.safeDof(2, 2)];

    // 罚刚度 1e12 / 结构刚度 ~4e10 → 相对误差 ~1e-4；50mm 下绝对 ~1e-5。
    expect(tag + " UX(S)-UX(T)=Δ", almostEq(uxS - uxT, Delta, 2e-5),
           "got=" + std::to_string(uxS - uxT) + " target=" + std::to_string(Delta));
    expect(tag + " UX(C)-UX(S)=0（缆鞍固结）", almostEq(uxC - uxS, 0.0, 2e-5),
           "got=" + std::to_string(uxC - uxS));
    expect(tag + " 主缆端相对塔顶滑移 Δ", almostEq(uxC - uxT, Delta, 3e-5),
           "got=" + std::to_string(uxC - uxT));
    expect(tag + " UY/UZ 固结成立",
           almostEq(uyC - uyT, 0.0, 1e-5) && almostEq(uzC - uzT, 0.0, 5e-5),
           "uy=" + std::to_string(uyC - uyT) + " uz=" + std::to_string(uzC - uzT));
}

}  // namespace

int main() {
    std::cout << "=== 3D CE 约束方程 benchmark（索鞍顶推口径 C）===\n";

    const double Delta = 0.05;   // 50 mm

    std::cout << "\n[1] 线性 solveModel3D\n";
    checkPath("linear", false, Delta);

    std::cout << "\n[2] 非线性 solveModel3DNonlinear（常数项不重复）\n";
    checkPath("nonlinear", true, Delta);

    std::cout << "\n[3] 无效 dof 被拒绝\n";
    {
        Model3D m = saddleModel(Delta, false);
        Model3D::ConstraintEquation ce;
        ce.id = 1; ce.constant = 0.0;
        ce.terms.push_back({2, 7, 1.0});   // dof 7 越界
        m.ces.push_back(ce);
        bool threw = false;
        try { (void)solveModel3D(m); } catch (const std::exception&) { threw = true; }
        expect("dof=7 被拒绝", threw);
    }
    {
        Model3D m = saddleModel(Delta, false);
        Model3D::ConstraintEquation ce;
        ce.id = 1; ce.constant = 0.0;
        ce.terms.push_back({3, 6, 1.0});   // 节点3 连 TRUSS，无转动 dof 6
        m.ces.push_back(ce);
        bool threw = false;
        try { (void)solveModel3D(m); } catch (const std::exception&) { threw = true; }
        expect("纯轴向节点转动 dof 被拒绝", threw);
    }

    std::cout << "\n[4] 重复系数合并\n";
    {
        Model3D m = saddleModel(0.0, false);
        // UX(3)-UX(1)=0.001，写成三条同 dof 项合并
        Model3D::ConstraintEquation ce;
        ce.id = 1; ce.constant = 0.001;
        ce.terms.push_back({3, 1, 0.5});
        ce.terms.push_back({3, 1, 0.3});
        ce.terms.push_back({3, 1, 0.2});
        ce.terms.push_back({1, 1, -1.0});
        m.ces.push_back(ce);
        m.cePenalty = 1e12;
        const auto U = solveModel3D(m);
        const DofMap3D dm = buildDofMap3D(m);
        const double diff = U[dm.safeDof(3, 0)] - U[dm.safeDof(1, 0)];
        expect("重复系数合并后 UX(3)-UX(1)=0.001", almostEq(diff, 0.001, 2e-5),
               "got=" + std::to_string(diff));
    }

    std::cout << "\n[5] 空 CE 零污染\n";
    {
        // 用稳定的单悬臂梁验证：无 CE 时可解且逐位确定（不用 saddleModel，
        // 后者在 jack_on=false 时是两个不相连的独立子结构，本身奇异）。
        auto cantilever = []() {
            Model3D m;
            m.nodes[1] = Vec3(0, 0, 0);
            m.nodes[2] = Vec3(3, 1, 10);
            Element3D eb; eb.id = 1; eb.type = "BEAM"; eb.n1 = 1; eb.n2 = 2;
            eb.E = 3.45e10; eb.A = 12.0; eb.Iy = 16.0; eb.Iz = 16.0;
            eb.J = 25.0; eb.G = 1.4e10;
            m.elements.push_back(eb);
            for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0, 0.0});
            return m;
        };
        Model3D a = cantilever();
        Model3D b = cantilever();
        const auto Ua = solveModel3D(a);
        const auto Ub = solveModel3D(b);
        bool same = Ua.size() == Ub.size();
        for (std::size_t i = 0; same && i < Ua.size(); ++i)
            same = almostEq(Ua[i], Ub[i], 1e-14);
        expect("空 CE 模型可解且确定", same);
    }

    std::cout << "\n=== " << (g_fail == 0 ? "CE_3D_BENCHMARK_PASS" : "CE_3D_BENCHMARK_FAIL")
              << " (fail=" << g_fail << ") ===\n";
    return g_fail == 0 ? 0 : 1;
}
