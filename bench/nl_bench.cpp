// 非线性内核基准：斜拉桥成桥态的规模—耗时曲线。
//
// 跑法：./nl_bench.exe [半跨索数] [主梁单元数/跨]
// 输出：模型规模、单次非线性求解耗时、索长迭代耗时、成桥态总耗时，
//       以及"复用符号分解"开关的对照（证明提速确实来自该优化而非噪声）。
//
// 这座桥是参数化生成的典型双塔斜拉桥：主梁连续梁 + 两座桥塔 + 扇形双索面
// （按单索面等效），索数 = 2·ncHalf。用它能回答"100 根索的实桥成桥态要多久"。

#include "fem_3d_nl.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>

using namespace fem;

namespace {

double nowMs() {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

// 双塔三跨斜拉桥：主跨 Lm，边跨 Ls，塔高 Ht，每塔 ncHalf 对索
Model3D buildCableStayed(int ncHalf, int nGirderPerSpan, double Lm, double Ls,
                         double Ht, double qGirder) {
    Model3D m;
    const double Eb = 3.45e10, Ab = 8.0, Ib = 4.0;
    const double Ec = 1.95e11, Ac = 6.0e-3, wc = 50.0;
    const double Et = 3.45e10, At = 12.0, It = 100.0;

    // 主梁：x = 0 .. 2Ls+Lm，塔在 x = Ls 与 x = Ls+Lm
    const double xT1 = Ls, xT2 = Ls + Lm;
    const double xEnd = 2.0 * Ls + Lm;
    const int nSpan = 3;
    int nid = 0;
    std::vector<double> xs;
    for (int s = 0; s < nSpan; ++s) {
        const double x0 = (s == 0) ? 0.0 : ((s == 1) ? Ls : Ls + Lm);
        const double x1 = (s == 0) ? Ls : ((s == 1) ? Ls + Lm : xEnd);
        for (int i = 0; i < nGirderPerSpan; ++i)
            xs.push_back(x0 + (x1 - x0) * static_cast<double>(i) / nGirderPerSpan);
    }
    xs.push_back(xEnd);
    std::vector<int> gNodes;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        m.nodes[++nid] = Vec3(xs[i], 0.0, 0.0);
        gNodes.push_back(nid);
    }
    int eid = 0;
    for (std::size_t i = 0; i + 1 < gNodes.size(); ++i) {
        Element3D g; g.id = ++eid; g.type = "BEAM";
        g.E = Eb; g.A = Ab; g.Iy = Ib; g.Iz = 4.0 * Ib; g.J = 8.0; g.G = 1.4e10;
        g.n1 = gNodes[i]; g.n2 = gNodes[i + 1];
        m.elements.push_back(g);
        // 主梁恒载：均布 q 的一致荷载向量（平动 + 弯矩），按单元分摊到两端
        const double L = xs[i + 1] - xs[i];
        m.loads.push_back({gNodes[i], 3, -qGirder * L / 2.0});
        m.loads.push_back({gNodes[i + 1], 3, -qGirder * L / 2.0});
    }

    // 桥塔：塔底（与主梁同位置，另建节点）→ 塔顶
    auto buildTower = [&](double xT, int& topNode) {
        const int nSeg = 8;
        std::vector<int> tNodes;
        for (int i = 0; i <= nSeg; ++i) {
            m.nodes[++nid] = Vec3(xT, 0.0, Ht * static_cast<double>(i) / nSeg);
            tNodes.push_back(nid);
        }
        for (int i = 0; i + 1 < static_cast<int>(tNodes.size()); ++i) {
            Element3D t; t.id = ++eid; t.type = "BEAM";
            t.E = Et; t.A = At; t.Iy = It; t.Iz = It; t.J = 40.0; t.G = 1.4e10;
            t.n1 = tNodes[i]; t.n2 = tNodes[i + 1];
            m.elements.push_back(t);
        }
        topNode = tNodes.back();
        // 塔底固结
        for (int d = 1; d <= 6; ++d) m.constraints.push_back({tNodes[0], d, 0.0, 0.0});
    };
    int top1 = 0, top2 = 0;
    buildTower(xT1, top1);
    buildTower(xT2, top2);

    // 拉索：从塔顶向两侧主梁节点扇形展开
    std::vector<int> cableIds;
    auto addCables = [&](int topNode, double xT, double dir) {
        // 取离塔最近的 ncHalf 个主梁节点（跨内）
        std::vector<int> cand;
        for (int n : gNodes) {
            const double x = m.nodes.at(n).x();
            if (dir * (x - xT) > 1e-9) cand.push_back(n);
        }
        std::sort(cand.begin(), cand.end(), [&](int a, int b) {
            return std::abs(m.nodes.at(a).x() - xT) < std::abs(m.nodes.at(b).x() - xT);
        });
        std::vector<int> pick;
        const int step = std::max(1, static_cast<int>(cand.size()) / ncHalf);
        // 跳过紧邻塔支座的主梁节点：刚性支承连续梁法在"支座与索锚点相距仅
        // 一个梁段"时会在两者间剧烈分配反力（实测 min/max 相差 17 倍）。
        const int skip = 1;
        for (int k = 0; k < ncHalf && (skip + k * step) < static_cast<int>(cand.size()); ++k)
            pick.push_back(cand[skip + k * step]);
        for (int n : pick) {
            Element3D c; c.id = ++eid; c.type = "CABLE";
            c.E = Ec; c.A = Ac; c.w = wc;
            c.n1 = topNode; c.n2 = n;
            m.elements.push_back(c);
            cableIds.push_back(eid);
        }
    };
    addCables(top1, xT1, -1.0);
    addCables(top1, xT1, +1.0);
    addCables(top2, xT2, -1.0);
    addCables(top2, xT2, +1.0);

    // 支座：两端 + 两塔处竖向约束；一端纵向约束；全桥横向约束
    auto support = [&](int n, bool fixX) {
        m.constraints.push_back({n, 3, 0.0, 0.0});
        m.constraints.push_back({n, 2, 0.0, 0.0});
        if (fixX) m.constraints.push_back({n, 1, 0.0, 0.0});
    };
    support(gNodes.front(), true);
    support(gNodes.back(), false);
    for (int n : gNodes) {
        if (std::abs(m.nodes.at(n).x() - xT1) < 1e-9 ||
            std::abs(m.nodes.at(n).x() - xT2) < 1e-9) support(n, false);
        m.constraints.push_back({n, 4, 0.0, 0.0});   // 扭转定位
    }
    return m;
}

}  // namespace

int main(int argc, char** argv) {
    const int ncHalf = (argc > 1) ? std::atoi(argv[1]) : 12;
    const int nG = (argc > 2) ? std::atoi(argv[2]) : 12;
    const double Lm = 200.0, Ls = 80.0, Ht = 60.0;
    const double q = 8.0 * 26.0e3;      // 主梁 A=8m² × 26 kN/m³

    std::printf("========== 非线性内核基准：双塔三跨斜拉桥 ==========\n");
    const Model3D m0 = buildCableStayed(ncHalf, nG, Lm, Ls, Ht, q);
    const DofMap3D dm = buildDofMap3D(m0);
    int nc = 0;
    for (const auto& e : m0.elements) if (e.type == "CABLE") ++nc;
    std::printf("主跨 %.0f m / 边跨 %.0f m / 塔高 %.0f m；单元 %d 个（索 %d 根），"
                "节点 %d，自由度 %d\n",
                Lm, Ls, Ht, static_cast<int>(m0.elements.size()), nc,
                static_cast<int>(m0.nodes.size()), dm.nDofs);

    // ---- 1. 单次非线性求解：复用分解 vs 每次完整分解 ----
    auto timeSolve = [&](bool reuse, SolverType st) {
        Model3D m = m0;
        NonlinearOptions3D opt;
        opt.loadSteps = 3; opt.maxIter = 30; opt.tolForce = 1e-8;
        opt.reuseFactorization = reuse;
        opt.solver = st;
        const double t0 = nowMs();
        const NonlinearResult3D r = solveModel3DNonlinear(m, opt);
        const double t1 = nowMs();
        return std::make_pair(t1 - t0, r.totalIterations);
    };
    {
        timeSolve(true, SolverType::SimplicialLDLT);   // 预热（让 pattern 稳定）
        const auto a = timeSolve(false, SolverType::SimplicialLDLT);
        const auto b = timeSolve(true, SolverType::SimplicialLDLT);
        const auto c = timeSolve(true, SolverType::CHOLMOD);
        std::printf("\n[单次非线性求解] %d 次 Newton 迭代\n", a.second);
        std::printf("  Eigen LDLT 每次完整分解: %8.1f ms\n", a.first);
        std::printf("  Eigen LDLT 复用符号分解: %8.1f ms  →  %.2fx\n", b.first, a.first / b.first);
        std::printf("  CHOLMOD    复用符号分解: %8.1f ms  →  %.2fx%s\n", c.first,
                    a.first / c.first, c.second != a.second ? "（迭代数不同，仅供参考）" : "");
    }

    // ---- 2. 索长迭代（影响矩阵 J = M·D）----
    {
        Model3D m = m0;
        std::vector<CableTarget3D> tg;
        for (const auto& e : m.elements)
            if (e.type == "CABLE") tg.push_back({e.id, 3.0e6});
        CableLengthOptions3D co;
        co.maxIter = 15; co.tolTension = 1e-7;
        co.nl.loadSteps = 3; co.nl.maxIter = 30; co.nl.tolForce = 1e-8;
        const double t0 = nowMs();
        const CableLengthResult3D cl = iterateCableUnstrainedLengths(m, tg, co);
        const double t1 = nowMs();
        std::printf("\n[%d 根索的无应力索长迭代]\n", nc);
        std::printf("  %d 轮，耗时 %.1f ms，索力误差 %.2e\n",
                    cl.iterations, t1 - t0, cl.maxTensionErr);
        std::printf("  （影响矩阵 J=M·D：1 次分解 + %d 次回代，而非 %d 次非线性求解）\n",
                    nc, nc);
    }

    // ---- 3. 成桥态闭环（目标索力自动 + 下料长度 + 制造线形）----
    {
        Model3D m = m0;
        std::vector<CompletionCable3D> cc;
        std::vector<CableAnchor3D> cc2;
        for (const auto& e : m.elements) {
            if (e.type != "CABLE") continue;
            // 梁侧锚点 = 主梁那一端（z 较小者）
            const int anchor = (m.nodes.at(e.n1).z() < m.nodes.at(e.n2).z()) ? e.n1 : e.n2;
            cc.push_back({e.id, anchor, 0.0});
            cc2.push_back({e.id, anchor});
        }
        if (argc > 3) {
            const std::vector<double> T = cableTensionsRigidSupport(m, cc2);
            double tmin = 1e30, tmax = -1e30, sum = 0.0;
            for (double t : T) { tmin = std::min(tmin, t); tmax = std::max(tmax, t); sum += t; }
            std::printf("  诊断：刚性支承法索力 min=%.4g max=%.4g sum=%.4g kN\n",
                        tmin / 1e3, tmax / 1e3, sum / 1e3);
        }
        CompletionStateOptions3D co;
        co.maxIter = 12; co.tolClosure = 1e-8; co.verbose = (argc > 3);
        co.nl.loadSteps = 3; co.nl.maxIter = 30; co.nl.tolForce = 1e-8;
        const double t0 = nowMs();
        const CompletionStateResult3D cs = solveCompletionState3D(m, cc, co);
        const double t1 = nowMs();
        std::printf("\n[成桥态闭环：自动定索力 + 下料长度 + 制造线形]\n");
        std::printf("  %d 轮几何迭代，耗时 %.1f ms；闭合误差 %.2e m，收敛=%d\n",
                    cs.iterations, t1 - t0, cs.maxClosureErr,
                    cs.converged ? 1 : 0);
    }

    std::printf("\n====================================\n");
    return 0;
}
