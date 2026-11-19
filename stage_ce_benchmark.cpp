// stage_ce_benchmark.cpp —— 分级顶推 CE 注入施工阶段的内核验收（第 41 项 CTest）
//
// 独立于悬索找形，直接构造三节点索鞍顶推模型 + 三个施工阶段，覆盖：
//   1. 线性 solveModel3DStages：CE= 写本阶段增量 δ，逐阶段累计滑移命中 Δ_i；
//   2. 非线性 solveModel3DStagesNonlinear：全量平衡下内核把 δ 累加为 Δ_i；
//   3. 无 CE= 覆盖：stage 沿用 *CE 原常量（单次顶推等价）；
//   4. CE= 覆盖不存在的 id 显式拒绝；
//   5. m.ces 为空时 stage 可解且确定（默认零污染）；
//   6. CE 引用孤立节点 S 不被硬约束（靠命中间接证明）。
//
// 口径（与 索鞍顶推CE耦合规格.md §8 一致）：
//   *STAGE 的 CE= 写增量 δ_i = Δ_i − Δ_(i−1)；线性直接覆盖 constant=δ_i，
//   非线性把历阶段 δ 累加为累计 Δ_i 后覆盖 constant。
//
// 退出码 0 = 全通过；非 0 = 失败（CTest 直接判定）。
#include "fem_3d.hpp"
#include "fem_3d_nl.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

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

// 构造三节点索鞍顶推模型；withStages 时按 deltas 生成 N 个阶段，每阶段
// CE= 写 deltas[i]（本阶段增量）。baseDelta 为 *CE 的 CE1 原常量。
Model3D stageJackModel(double baseDelta, bool withStages,
                       const std::vector<double>& deltas) {
    Model3D m;
    m.nodes[1] = Vec3(0, 0, 0);
    m.nodes[2] = Vec3(0, 0, 10);      // T 塔顶（连 BEAM）
    m.nodes[3] = Vec3(0, 0, 10);      // C 主缆端（连 TRUSS）
    m.nodes[4] = Vec3(3, 1, 3);       // 固定端（y≠0 使节点3的 UY 有刚度）
    m.nodes[5] = Vec3(0, 0, 10);      // S 索鞍：孤立，仅靠 CE 罚刚度获得自由度

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

    auto add = [&](int id, double b, int n1, int d1, double c1,
                   int n2, int d2, double c2) {
        Model3D::ConstraintEquation ce;
        ce.id = id; ce.constant = b;
        ce.terms.push_back({n1, d1, c1});
        ce.terms.push_back({n2, d2, c2});
        m.ces.push_back(ce);
    };
    add(1, baseDelta, 5, 1, 1.0, 2, 1, -1.0);   // CE1 UX(S)−UX(T)=baseDelta
    add(2, 0.0, 3, 1, 1.0, 5, 1, -1.0);         // CE2 UX(C)−UX(S)=0
    add(3, 0.0, 5, 2, 1.0, 2, 2, -1.0);         // CE3 UY(S)−UY(T)=0
    add(4, 0.0, 5, 3, 1.0, 2, 3, -1.0);         // CE4 UZ(S)−UZ(T)=0
    add(5, 0.0, 3, 2, 1.0, 2, 2, -1.0);         // CE5 UY(C)−UY(T)=0
    add(6, 0.0, 3, 3, 1.0, 2, 3, -1.0);         // CE6 UZ(C)−UZ(T)=0
    m.cePenalty = 1e12;

    if (withStages) {
        for (std::size_t i = 0; i < deltas.size(); ++i) {
            Model3D::Stage3D st;
            st.id = static_cast<int>(i) + 1;
            st.activate = {1, 2};
            st.ceConstants.push_back({1, deltas[i]});
            m.stages.push_back(st);
        }
    }
    return m;
}

void checkLinear() {
    const std::vector<double> d = {0.01, 0.01, 0.01};
    Model3D m = stageJackModel(0.0, true, d);
    const auto out = solveModel3DStages(m);
    expect("线性 stage 数=3", out.size() == 3,
           "got=" + std::to_string(out.size()));

    // StageResult3D::U 是本阶段增量；按节点累加各阶段位移得累计滑移。
    std::map<int, std::array<double, 6>> uTotal;
    for (std::size_t s = 0; s < out.size(); ++s) {
        const auto& sr = out[s];
        for (const auto& [nid, x] : m.nodes) {
            (void)x;
            for (int k = 0; k < sr.dm.count(nid); ++k) {
                const int dd = sr.dm.safeDof(nid, k);
                if (dd >= 0 && dd < static_cast<int>(sr.U.size()))
                    uTotal[nid][static_cast<std::size_t>(k)] += sr.U[static_cast<std::size_t>(dd)];
            }
        }
        const double slip = uTotal[3][0] - uTotal[2][0];
        const double target = 0.01 * static_cast<double>(s + 1);
        expect("线性 第 " + std::to_string(s + 1) + " 阶段累计滑移=Δ",
               almostEq(slip, target, 3e-5),
               "got=" + std::to_string(slip) + " target=" + std::to_string(target));
    }
}

void checkNonlinear() {
    const std::vector<double> d = {0.01, 0.01, 0.01};
    Model3D m = stageJackModel(0.0, true, d);
    NonlinearOptions3D opt;
    opt.loadSteps = 1;
    opt.maxIter = 50;
    opt.tolForce = 1e-10;
    opt.tolDisp = 1e-12;
    const auto out = solveModel3DStagesNonlinear(m, opt);
    expect("非线性 stage 数=3", out.size() == 3,
           "got=" + std::to_string(out.size()));

    // NLStageResult3D::u 是全量位移（每阶段独立全量平衡）——累计滑移直接读出。
    for (std::size_t s = 0; s < out.size(); ++s) {
        const auto& sr = out[s];
        const double uxC = sr.u[static_cast<std::size_t>(sr.dm.safeDof(3, 0))];
        const double uxT = sr.u[static_cast<std::size_t>(sr.dm.safeDof(2, 0))];
        const double slip = uxC - uxT;
        const double target = 0.01 * static_cast<double>(s + 1);
        expect("非线性 第 " + std::to_string(s + 1) + " 阶段累计滑移=Δ",
               almostEq(slip, target, 3e-5) && sr.converged,
               "got=" + std::to_string(slip) + " target=" + std::to_string(target)
               + " converged=" + std::to_string(sr.converged));
    }
}

void checkNoOverride() {
    // 有 stage、但不写 CE=：应沿用 *CE 原常量 baseDelta（单次顶推等价）。
    const double base = 0.02;
    Model3D m = stageJackModel(base, true, {});
    for (int i = 0; i < 2; ++i) {
        Model3D::Stage3D st;
        st.id = i + 1;
        st.activate = {1, 2};
        m.stages.push_back(st);
    }
    NonlinearOptions3D opt;
    opt.loadSteps = 1; opt.maxIter = 50; opt.tolForce = 1e-10; opt.tolDisp = 1e-12;
    const auto out = solveModel3DStagesNonlinear(m, opt);
    expect("无 CE= 覆盖 stage 数=2", out.size() == 2);
    if (out.size() == 2) {
        const auto& sr = out.back();
        const double slip = sr.u[static_cast<std::size_t>(sr.dm.safeDof(3, 0))]
                          - sr.u[static_cast<std::size_t>(sr.dm.safeDof(2, 0))];
        expect("无 CE= 覆盖：stage 沿用 *CE 原常量",
               almostEq(slip, base, 3e-5),
               "got=" + std::to_string(slip) + " target=" + std::to_string(base));
    }
}

void checkReject() {
    Model3D m = stageJackModel(0.0, true, {0.01});
    m.stages[0].ceConstants[0].first = 999;   // 不存在的 CE id
    bool threw = false;
    try { (void)solveModel3DStages(m); } catch (const std::exception&) { threw = true; }
    expect("线性：CE= 引用不存在的 id 被拒绝", threw);

    Model3D n = stageJackModel(0.0, true, {0.01});
    n.stages[0].ceConstants[0].first = 999;
    threw = false;
    NonlinearOptions3D opt; opt.loadSteps = 1;
    try { (void)solveModel3DStagesNonlinear(n, opt); } catch (const std::exception&) { threw = true; }
    expect("非线性：CE= 引用不存在的 id 被拒绝", threw);
}

void checkNoCe() {
    // 默认零污染：无 CE 时 stage 路径应可解且逐位确定。用稳定单悬臂梁
    // （索鞍模型无 CE 时节点3的 UY/UZ 零刚度，本身不适合做此验证）。
    Model3D m;
    m.nodes[1] = Vec3(0, 0, 0);
    m.nodes[2] = Vec3(1, 0, 0);
    Element3D e; e.id = 1; e.type = "BEAM"; e.n1 = 1; e.n2 = 2;
    e.E = 3.45e10; e.A = 0.01; e.Iy = 1e-4; e.Iz = 1e-4;
    e.J = 1e-4; e.G = 1.4e10;
    m.elements.push_back(e);
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0, 0.0});
    m.loads.push_back({2, 2, -1000.0});
    Model3D::Stage3D st; st.id = 1; st.activate = {1};
    m.stages.push_back(st);

    const auto out = solveModel3DStages(m);
    expect("无 CE stage 可解", out.size() == 1 && !out[0].U.empty());
    const auto again = solveModel3DStages(m);
    bool same = again.size() == out.size();
    if (same)
        for (std::size_t i = 0; i < out.size() && same; ++i) {
            if (out[i].U.size() != again[i].U.size()) { same = false; break; }
            for (std::size_t k = 0; k < out[i].U.size(); ++k)
                if (out[i].U[k] != again[i].U[k]) { same = false; break; }
        }
    expect("无 CE stage 逐位确定", same);
}

}  // namespace

int main() {
    std::cout << "=== 分级顶推 CE 施工阶段内核验收（第 41 项）===\n";
    checkLinear();
    std::cout << "\n[非线性 stage 路径]\n";
    checkNonlinear();
    std::cout << "\n[无 CE= 覆盖]\n";
    checkNoOverride();
    std::cout << "\n[拒绝路径]\n";
    checkReject();
    std::cout << "\n[默认零污染]\n";
    checkNoCe();

    std::cout << "\n结果: " << (g_fail == 0 ? "全部通过" : std::to_string(g_fail) + " 项失败")
              << "\n";
    return g_fail == 0 ? 0 : 1;
}
