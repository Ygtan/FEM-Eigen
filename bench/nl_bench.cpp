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

    // ---- 0. 单元级微基准：内力 / 解析切线 / 差分切线 ----
    // 这是定位"耗时随单元数线性增长"的关键实验：差分切线每单元要 2n 次内力
    // 评估，而内力评估本身含 4 个超越函数，于是切线吃掉 96% 的单元时间。
    {
        std::vector<CRElemData> eds;
        for (const auto& e : m0.elements) eds.push_back(crMakeElemData(m0, e, dm));
        int ib = -1, ic = -1;
        for (std::size_t i = 0; i < eds.size(); ++i) {
            if (eds[i].n == 12 && ib < 0) ib = static_cast<int>(i);
            if (eds[i].n == 6 && ic < 0) ic = static_cast<int>(i);
        }
        const int N = 20000;
        double sink = 0.0;
        // mode: 0 = 只算内力, 1 = 解析切线, 2 = 差分切线
        auto timeElem = [&](const CRElemData& ed, int mode) {
            Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
            ue[0] = 2e-4; ue[2] = 1e-3; ue[4] = 3e-3; ue[7] = -2e-4; ue[10] = -1e-4;
            Eigen::Matrix<double, 14, 1> fe;
            Eigen::Matrix<double, 14, 14> Ke;
            double best = 1e30;
            for (int rep = 0; rep < 3; ++rep) {          // 取最小，避开调度噪声
                const double t0 = nowMs();
                for (int k = 0; k < N; ++k) {
                    crElementInternal(ed, ue, fe, mode != 0,
                                      mode != 0 ? &Ke : nullptr, nullptr, nullptr,
                                      mode == 2);
                    sink += fe[0] + (mode != 0 ? Ke(0, 0) + Ke(ed.n - 1, ed.n - 1) : 0.0);
                }
                best = std::min(best, (nowMs() - t0) * 1e3 / N);
            }
            return best;                                 // µs / 次
        };
        if (ib >= 0) {
            const double tf = timeElem(eds[ib], 0);
            const double ta = timeElem(eds[ib], 1);
            const double td = timeElem(eds[ib], 2);
            std::printf("\n[单元级微基准 BEAM(12)]\n");
            std::printf("  只算内力:        %7.3f µs\n", tf);
            std::printf("  内力+解析切线:   %7.3f µs  ← 默认\n", ta);
            std::printf("  内力+差分切线:   %7.3f µs  （每单元 2n=24 次内力评估）\n", td);
            std::printf("  解析切线加速:    %.1fx（切线本身 %.1fx）\n",
                        td / ta, (td - tf) / std::max(ta - tf, 1e-9));
        }
        if (ic >= 0) {
            const double tf = timeElem(eds[ic], 0);
            const double ta = timeElem(eds[ic], 1);
            const double td = timeElem(eds[ic], 2);
            std::printf("[单元级微基准 CABLE(6)] 内力 %.3f / 解析 %.3f / 差分 %.3f µs  → %.1fx\n",
                        tf, ta, td, td / ta);
        }
        if (sink == 1.2345) std::printf("");   // 防优化消除
    }

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
        // 目标索力用最小弯曲能量法的解，而不是"所有索统一 3 MN"——后者在扇形
        // 索面里物理上不自洽（外索与内索的合理索力差可达十倍），80 索时会让
        // 非线性求解不收敛（实测）。
        std::vector<CableAnchor3D> ca;
        for (const auto& e : m.elements) {
            if (e.type != "CABLE") continue;
            const int anchor = (m.nodes.at(e.n1).z() < m.nodes.at(e.n2).z()) ? e.n1 : e.n2;
            ca.push_back({e.id, anchor});
        }
        const std::vector<double> Tref = cableTensionsMinBendingEnergy(m, ca, CableForceOptions3D{});
        std::vector<CableTarget3D> tg;
        for (std::size_t k = 0; k < ca.size(); ++k) tg.push_back({ca[k].elem, Tref[k]});
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
        // ---- 两种自动定索力方法对照 ----
        // 判据两条：① 沿索向竖向分量之和必须等于总恒载（全局平衡，与算法无关）；
        //          ② 索力分布的 max/min（刚性支承法在锚点紧邻支座时会畸变）。
        double Wtot = 0.0;
        for (const auto& [node, dof, val] : m.loads) if (dof == 3) Wtot -= val;
        auto report = [&](const char* tag, const std::vector<double>& T) {
            double tmin = 1e30, tmax = -1e30, sum = 0.0, vsum = 0.0;
            for (std::size_t k = 0; k < T.size(); ++k) {
                tmin = std::min(tmin, T[k]);
                tmax = std::max(tmax, T[k]);
                sum += T[k];
                const Element3D* e = nullptr;
                for (const auto& q : m.elements) if (q.id == cc[k].elem) { e = &q; break; }
                const int A = cc[k].anchorNode;
                const int B = (e->n1 == A) ? e->n2 : e->n1;
                const Vec3 ev = (m.nodes.at(B) - m.nodes.at(A)).normalized();
                vsum += T[k] * ev.z();
            }
            std::printf("  %-16s T∈[%.4g, %.4g] kN，ΣT=%.4g kN，max/min=%.1f\n",
                        tag, tmin / 1e3, tmax / 1e3, sum / 1e3, tmax / std::max(tmin, 1e-9));
            std::printf("  %-16s 竖向分量 Σ=%.4g kN（总恒载 %.4g kN，偏差 %.1f%%）\n",
                        "", vsum / 1e3, Wtot / 1e3,
                        100.0 * (vsum - Wtot) / std::max(Wtot, 1e-9));
        };
        report("刚性支承连续梁法", cableTensionsRigidSupport(m, cc2));
        report("最小弯曲能量法", cableTensionsMinBendingEnergy(m, cc2, CableForceOptions3D{}));

        CompletionStateOptions3D co;
        co.maxIter = 20; co.tolClosure = 1e-8; co.verbose = (argc > 3);
        co.nl.loadSteps = 3; co.nl.maxIter = 30; co.nl.tolForce = 1e-8;
        const double t0 = nowMs();
        const CompletionStateResult3D cs = solveCompletionState3D(m, cc, co);
        const double t1 = nowMs();
        std::printf("\n[成桥态闭环：自动定索力 + 下料长度 + 制造线形]\n");
        std::printf("  %d 轮几何迭代，耗时 %.1f ms；闭合误差 %.2e m，收敛=%d\n",
                    cs.iterations, t1 - t0, cs.maxClosureErr, cs.converged ? 1 : 0);
        std::printf("  成桥索力误差 %.2e（相对目标），最大预抬高 %.1f mm\n",
                    cs.maxTensionErr,
                    1e3 * [&] {
                        double mx = 0.0;
                        for (std::size_t i = 0; i < cs.camber.size(); ++i) mx = std::max(mx, std::abs(cs.camber[i]));
                        return mx;
                    }());
    }

    std::printf("\n====================================\n");
    return 0;
}
