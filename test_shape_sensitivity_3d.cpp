// ============================================================
// P3-3b：完整形状敏感度 K_shape 的独立验收 oracle
//
// 规格见 USER_MANUAL.md §9.8：
//  1) crShapeResidual3D 必须与 solveModel3DNonlinear 的 Fext/Fint 口径同源；
//  2) G_full·v 必须与「扰动参考构形 + 完整重求平衡」的 du/dX·v 对拍；
//  3) 完整 / 低秩 / 纯不动点三条路径逐轮误差对比，最终制造线形与索力
//     必须一致到 tolClosure；
//  4) 默认 LowRank 分支不得因本次改动而改变。
//
// 这里全部是**独立 oracle**：不复用 computeShapeSensitivity3D 的公式，
// 数值对拍靠重跑 solveModel3DNonlinear，误差判据靠实际闭合残差。
// ============================================================

#include "fem_3d_nl.hpp"
#include "test_common.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

using namespace fem;

namespace {

// 三索两跨连续梁：竖向支承全部由索提供，梁上无索节点（节点 2/4）也存在，
// 正好用来证明「扰动非锚点的参考坐标同样会引起内力变化」——这正是低秩 G
// 漏掉的那一项。
Model3D makeModel(double L, double H, double Eb, double A, double I,
                  double Ec, double Ac, double P, bool fixEnd = true) {
    Model3D m;
    for (int i = 1; i <= 5; ++i) m.nodes[i] = Vec3(0.5 * L * (i - 1), 0.0, 0.0);
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
        c.E = Ec; c.A = Ac; c.w = 0.0;
        c.n1 = cn[k]; c.n2 = 6 + k;
        m.elements.push_back(c);
    }
    m.constraints.push_back({1, 1, 0.0, 0.0});
    m.constraints.push_back({1, 4, 0.0, 0.0});
    // 端部固结（x / z / θy）：打断几何闭合的两个刚体规范自由度。
    // ① 整段梁参考构形 + 位移反向等量平移；② 绕 node1 的刚体转动
    //    （ΔX 在梁节点上线性递增、转角自由度由求解器自行取 θ）。
    // 不打断时当前构形与索力唯一，但制造线形 ΔX 本身是一族解，
    // (I+G) 有零空间，不同路径会停在同一族的不同点上。
    // fixEnd=false 用于 testGaugeFreedom 复现这个现象。
    if (fixEnd) {
        m.constraints.push_back({1, 3, 0.0, 0.0});
        m.constraints.push_back({1, 5, 0.0, 0.0});
    }
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

// 与 iterateErectionGeometry 里的 isT 构造逐条一致。
std::vector<char> translationDofs(const Model3D& m, const DofMap3D& dm) {
    const int nd = dm.nDofs;
    std::vector<char> isC(nd, 0), isSpring(nd, 0), isT(nd, 0);
    for (const auto& c : m.constraints) {
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        if (c.k > 0.0) isSpring[d] = 1; else isC[d] = 1;
    }
    for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
        const int cnt = std::min(dm.dofCount[i], 3);
        for (int a = 0; a < cnt; ++a) isT[dm.dofStart[i] + a] = 1;
    }
    for (int d = 0; d < nd; ++d) if (isC[d] || isSpring[d]) isT[d] = 0;
    return isT;
}

struct Locator {
    const DofMap3D& dm;
    bool operator()(int dof, int& nid, int& comp) const {
        for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
            if (dof >= dm.dofStart[i] && dof < dm.dofStart[i] + dm.dofCount[i]) {
                nid = dm.nodeIds[i];
                comp = dof - dm.dofStart[i];
                return true;
            }
        }
        return false;
    }
};

// dof 属于的节点是否挂了索单元（低秩 B 在这些 dof 上有非零行）。
std::vector<char> cableAnchorNodes(const Model3D& m) {
    std::vector<char> anchor(9, 0);
    for (const auto& e : m.elements)
        if (e.type == "CABLE") { anchor[e.n1] = 1; anchor[e.n2] = 1; }
    return anchor;
}

// 三根竖索按下料长度给值（刚性支承连续梁法目标索力）。
void setTargetCableL0(Model3D& m, double P) {
    const std::vector<int> ids = {10, 11, 12};
    const double Ttar[3] = {5.0 * P / 16.0, 11.0 * P / 8.0, 5.0 * P / 16.0};
    for (int k = 0; k < 3; ++k) {
        for (auto& e : m.elements) {
            if (e.id != ids[k]) continue;
            const Vec3 d3 = m.nodes.at(e.n2) - m.nodes.at(e.n1);
            const double Lh = std::sqrt(d3.x() * d3.x() + d3.y() * d3.y());
            e.L0 = cableUnstrainedLengthForTension(d3.norm(), Ttar[k],
                                                   e.E * e.A, e.w, Lh);
        }
    }
}

ErectionGeometryOptions3D pathOptions() {
    ErectionGeometryOptions3D o;
    o.maxIter = 60;
    o.tolClosure = 1e-12;
    o.nl.loadSteps = 2;
    o.nl.maxIter = 40;
    o.nl.tolForce = 1e-11;
    o.nl.tolDisp = 1e-12;
    return o;
}

// ---- oracle 1：残差装配口径 + K_shape 低秩之外的项确实存在 ----
void testResidualAndShapeTerm(const Model3D& m, const DofMap3D& dm,
                              const std::vector<char>& isT,
                              const NLTangent3D& tan,
                              const NonlinearResult3D& r0,
                              ShapeSensitivityResult3D& ss,
                              double P) {
    std::printf("\n[oracle 1] 残差口径同源 + K_shape 确实含索之外的项\n");
    const int nd = dm.nDofs;
    std::vector<CRElemData> eds;
    for (const auto& e : m.elements) eds.push_back(crMakeElemData(m, e, dm));
    Eigen::VectorXd u = Eigen::VectorXd::Zero(nd);
    for (int d = 0; d < nd; ++d) u(d) = r0.u[d];
    const Eigen::VectorXd R = crShapeResidual3D(m, dm, eds, u);

    // 收敛态：自由自由度上的残差必须接近 0（口径若漏项，这里立刻非零）。
    std::vector<char> isC(nd, 0);
    for (const auto& c : m.constraints) {
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0 && c.k <= 0.0) isC[d] = 1;
    }
    double rn = 0.0;
    for (int d = 0; d < nd; ++d) if (!isC[d]) rn += R(d) * R(d);
    check("残差口径：收敛态自由自由度 |R| / 2P", std::sqrt(rn) / (2.0 * P), 0.0, 1e-6);

    ShapeSensitivityOptions3D so;
    so.fdStep = 1e-5;
    ss = computeShapeSensitivity3D(m, dm, r0.u, isT, tan, so);
    const int nt = static_cast<int>(ss.tdofs.size());
    check("K_shape：列数 = isT 自由度数", static_cast<double>(nt),
          static_cast<double>(std::count(isT.begin(), isT.end(), char(1))), 0.0);
    check("K_shape：装配次数 = 2·nt", static_cast<double>(ss.assemblies),
          2.0 * nt, 0.0);
    check("K_shape：|K_shape|max > 0", ss.kShapeMax > 0.0 ? 1.0 : 0.0, 1.0, 0.0);
    // 只清硬约束行；自由转角行必须保留（否则 G_full 仍是近似的）。
    // isC 已在收敛态残差检查时构造，直接复用。
    double cRows = 0.0, rotRows = 0.0;
    for (int d = 0; d < nd; ++d) {
        const double v = ss.Kshape.row(d).cwiseAbs().maxCoeff();
        if (isC[d]) cRows = std::max(cRows, v);
        else if (!isT[d]) rotRows = std::max(rotRows, v);
    }
    check("K_shape：硬约束行已清零", cRows, 0.0, 1e-12);
    check("K_shape：自由转角行非零（平衡右端项保留）",
          rotRows > 0.0 ? 1.0 : 0.0, 1.0, 0.0);
    std::printf("  自由转角行最大 |K_shape| = %.4g N·m/m\n", rotRows);

    // 决定性证据：非索锚点节点的参考坐标扰动，K_shape 列非零。
    // 低秩 B 在这些自由度上整列为 0，所以这一项正是 P3-3b 漏掉的部分。
    Locator loc{dm};
    const std::vector<char> anchor = cableAnchorNodes(m);
    double nonAnchorCol = 0.0;
    int nonAnchorCount = 0;
    for (int j = 0; j < nt; ++j) {
        int nid = -1, c = -1;
        if (!loc(ss.tdofs[j], nid, c)) continue;
        if (anchor[nid]) continue;
        ++nonAnchorCount;
        nonAnchorCol = std::max(nonAnchorCol, ss.Kshape.col(j).cwiseAbs().maxCoeff());
    }
    check("K_shape：存在非索锚点的可动平动 dof",
          nonAnchorCount > 0 ? 1.0 : 0.0, 1.0, 0.0);
    check("K_shape：非索锚点列非零（梁参考形状项）",
          nonAnchorCol > 0.0 ? 1.0 : 0.0, 1.0, 0.0);
    std::printf("  非锚点列最大 |K_shape| = %.4g N/m（低秩 G 里这些列恒为 0）\n",
                nonAnchorCol);

    // 与低秩 G 的差有多大（取证用，不设「必须 >20%」的拍脑袋限值）。
    const std::vector<int> ids = {10, 11, 12};
    std::vector<double> T0(3, 0.0);
    for (int k = 0; k < 3; ++k) {
        std::size_t q = 0;
        while (q < r0.elemIds.size() && r0.elemIds[q] != ids[k]) ++q;
        T0[k] = r0.axial[q];
    }
    const CableSensitivity3D s = cableSensitivity3D(m, dm, ids, r0.u, T0, isT, tan);
    const Eigen::MatrixXd GlowNd = -s.W * s.Kc.asDiagonal() * s.B.transpose();
    Eigen::MatrixXd Glow(nt, nt);
    for (int i = 0; i < nt; ++i)
        for (int j = 0; j < nt; ++j) Glow(i, j) = GlowNd(ss.tdofs[i], ss.tdofs[j]);
    const double relDiff = (ss.Gfull - Glow).norm() /
                           std::max(ss.Gfull.norm(), 1e-30);
    check("G_full 与低秩 G 不同（近似误差真实存在）",
          relDiff > 1e-8 ? 1.0 : 0.0, 1.0, 0.0);
    check("G_full 与低秩 G 同量级（不是数值噪声）",
          ss.Gfull.norm() / std::max(Glow.norm(), 1e-30), 1.0, 5.0);
    std::printf("  ||G_full||=%.4g  ||G_low||=%.4g  相对差=%.4g\n",
                ss.Gfull.norm(), Glow.norm(), relDiff);

    // 条件数取证：端部固结后 (I+G_full) 不能还有近零奇异值，否则
    // 制造线形 ΔX 仍不唯一（见 testGaugeFreedom）。
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        Eigen::MatrixXd::Identity(nt, nt) + ss.Gfull);
    const double smin = svd.singularValues().tail<1>()(0);
    const double smax = svd.singularValues()(0);
    check("I+G_full 无近零奇异值（无规范自由度）",
          smin > 1e-4 ? 1.0 : 0.0, 1.0, 0.0);
    std::printf("  σ(I+G_full) ∈ [%.4g, %.4g]，条件数=%.4g\n",
                smin, smax, smax / std::max(smin, 1e-30));
}

// ---- oracle 2：G_full·v 与数值 du/dX·v 逐分量对拍 ----
void testDirectionalAgainstFD(const Model3D& m, const DofMap3D& dm,
                              const std::vector<char>& isT,
                              const NonlinearResult3D& r0,
                              const ShapeSensitivityResult3D& ss,
                              const NonlinearOptions3D& nl) {
    std::printf("\n[oracle 2] G_full·v vs 中心差分 du/dX·v（独立重求平衡）\n");
    const int nt = static_cast<int>(ss.tdofs.size());
    Locator loc{dm};

    // 方向向量：非平凡组合，避免只用单位基向量掩盖耦合项错误。
    Eigen::VectorXd v(nt);
    for (int j = 0; j < nt; ++j)
        v(j) = ((j % 2) ? -1.0 : 1.0) * (1.0 + 0.17 * j);
    v /= v.cwiseAbs().maxCoeff();

    const Eigen::VectorXd ana = ss.Gfull * v;    // nt 维（isT 行）
    auto modelAt = [&](double h) {
        Model3D mm = m;
        for (int j = 0; j < nt; ++j) {
            int nid = -1, c = -1;
            if (!loc(ss.tdofs[j], nid, c) || c > 2) continue;
            mm.nodes[nid][c] += h * v(j);
        }
        return mm;
    };

    // 三组步长：确认差分随 h 稳定（不是某个步长碰巧命中）。
    const double hs[3] = {1e-2, 1e-3, 1e-4};
    for (double h : hs) {
        const NonlinearResult3D rp = solveModel3DNonlinear(modelAt(+h), nl, nullptr);
        const NonlinearResult3D rm = solveModel3DNonlinear(modelAt(-h), nl, nullptr);
        Eigen::VectorXd num(nt);
        for (int i = 0; i < nt; ++i)
            num(i) = (rp.u[ss.tdofs[i]] - rm.u[ss.tdofs[i]]) / (2.0 * h);
        const double l2 = (num - ana).norm() / std::max(ana.norm(), 1e-30);
        const double dot = num.dot(ana) / std::max(num.norm() * ana.norm(), 1e-30);
        check("方向导数 h=" + std::to_string(h).substr(0, 6) + " 相对 L2 差",
              l2, 0.0, 2e-2);
        check("方向导数 h=" + std::to_string(h).substr(0, 6) + " 方向一致",
              dot, 1.0, 2e-2);
        std::printf("  h=%.1e：||num||=%.6g  ||ana||=%.6g  L2差=%.3e  相关=%.8f\n",
                    h, num.norm(), ana.norm(), l2, dot);
    }
    // 换一个方向再做一次：取 G_full 响应最大的那一列，避免挑到弱响应自由度
    // 让「相对误差」被小分母放大。
    int jbest = 0;
    for (int j = 1; j < nt; ++j)
        if (ss.Gfull.col(j).norm() > ss.Gfull.col(jbest).norm()) jbest = j;
    Eigen::VectorXd v2 = Eigen::VectorXd::Zero(nt);
    v2(jbest) = 1.0;
    const Eigen::VectorXd ana2 = ss.Gfull * v2;
    const double h = 1e-3;
    auto modelAt2 = [&](double hh) {
        Model3D mm = m;
        for (int j = 0; j < nt; ++j) {
            int nid = -1, c = -1;
            if (!loc(ss.tdofs[j], nid, c) || c > 2) continue;
            mm.nodes[nid][c] += hh * v2(j);
        }
        return mm;
    };
    const NonlinearResult3D rp = solveModel3DNonlinear(modelAt2(+h), nl, nullptr);
    const NonlinearResult3D rm = solveModel3DNonlinear(modelAt2(-h), nl, nullptr);
    Eigen::VectorXd num2(nt);
    for (int i = 0; i < nt; ++i)
        num2(i) = (rp.u[ss.tdofs[i]] - rm.u[ss.tdofs[i]]) / (2.0 * h);
    check("单自由度方向逐分量 L2 差",
          (num2 - ana2).norm() / std::max(ana2.norm(), 1e-30), 0.0, 2e-2);
    (void)isT;
}

// ---- oracle 3/4：三条路径逐轮误差、终值一致、默认不变 ----
void testThreePaths(const Model3D& m) {
    std::printf("\n[oracle 3/4] 完整 / 低秩 / 纯不动点三路径对比\n");
    const std::vector<int> ids = {10, 11, 12};
    ErectionGeometryOptions3D base;
    base.maxIter = 60;
    // I+G 条件数偏大时，闭合误差 ε 只把 camber 误差压到 ε/σ_min；
    // 这里收紧到 1e-12，让三条路径真正收敛到同一个制造线形，
    // 而不是各自在 tolClosure=1e-9 处提前停下（实测 camber 仍差 2e-5 m）。
    base.tolClosure = 1e-12;
    base.nl.loadSteps = 2;
    base.nl.maxIter = 40;
    base.nl.tolForce = 1e-11;
    base.nl.tolDisp = 1e-12;

    ErectionGeometryResult3D low, full, fp;
    {
        ErectionGeometryOptions3D o = base;
        o.sensitivity = SensitivityMode3D::LowRank;
        low = iterateErectionGeometry(m, ids, o);
    }
    {
        ErectionGeometryOptions3D o = base;
        o.sensitivity = SensitivityMode3D::Full;
        full = iterateErectionGeometry(m, ids, o);
    }
    {
        ErectionGeometryOptions3D o = base;
        o.sensitivity = SensitivityMode3D::FixedPoint;
        // 纯不动点收敛率 = |G| 的谱半径，本算例约 0.933（远慢于 Newton），
        // 需要更多轮才能到同一个 tolClosure；这是路径特性，不是失败。
        o.maxIter = 500;
        fp = iterateErectionGeometry(m, ids, o);
    }

    check("三路径：默认枚举 = LowRank",
          ErectionGeometryOptions3D{}.sensitivity == SensitivityMode3D::LowRank
              ? 1.0 : 0.0, 1.0, 0.0);
    check("三路径：LowRank 收敛", low.converged ? 1.0 : 0.0, 1.0, 0.0);
    check("三路径：Full 收敛", full.converged ? 1.0 : 0.0, 1.0, 0.0);
    check("三路径：FixedPoint 收敛", fp.converged ? 1.0 : 0.0, 1.0, 0.0);
    check("三路径：同一起点（首轮误差相同）",
          std::abs(low.trace.steps.front().err - full.trace.steps.front().err), 0.0, 1e-14);
    std::printf("  起点误差 LowRank=%.4e  Full=%.4e  FixedPoint=%.4e m\n",
                low.trace.steps.front().err, full.trace.steps.front().err,
                fp.trace.steps.front().err);

    auto fmtErr = [](const ErectionGeometryResult3D& r, int i) {
        char b[32];
        if (i >= static_cast<int>(r.trace.steps.size())) std::snprintf(b, 32, "-");
        else std::snprintf(b, 32, "%.4e", r.trace.steps[i].err);
        return std::string(b);
    };
    std::printf("  逐轮闭合误差（m）：\n");
    const int n = static_cast<int>(std::max({low.trace.steps.size(),
                                             full.trace.steps.size(),
                                             fp.trace.steps.size()}));
    const int nshow = std::min(n, 12);
    for (int i = 0; i < nshow; ++i)
        std::printf("    轮 %2d  LowRank=%s  Full=%s  FixedPoint=%s\n", i + 1,
                    fmtErr(low, i).c_str(), fmtErr(full, i).c_str(),
                    fmtErr(fp, i).c_str());
    if (n > nshow) std::printf("    ...（共 %d 轮，纯不动点尾部已省略）\n", n);
    std::printf("  轮数：LowRank=%d  Full=%d  FixedPoint=%d\n",
                low.iterations, full.iterations, fp.iterations);

    // 开关确实接了线：三条路径的第二轮误差不能全同。
    auto secondErr = [](const ErectionGeometryResult3D& r) {
        return r.trace.steps.size() > 1 ? r.trace.steps[1].err : -1.0;
    };
    check("三路径：Full 与 LowRank 第 2 轮误差不同（开关生效）",
          std::abs(secondErr(full) - secondErr(low)) > 1e-12 ? 1.0 : 0.0, 1.0, 0.0);
    // FixedPoint 的 δ = -Fv：首轮步长不应被 (I+G)⁻¹ 放大，故 stepMax <= err。
    // 这个判据不依赖「两条路径恰好重合」的运气。
    const double stepFP = fp.trace.steps.front().stepMax;
    const double stepLow = low.trace.steps.front().stepMax;
    const double err0 = low.trace.steps.front().err;
    std::printf("  首轮步长：LowRank=%.6e  FixedPoint=%.6e  首轮误差=%.6e\n",
                stepLow, stepFP, err0);
    check("三路径：FixedPoint 首轮步长 = |Fv|（G=0 分支生效）",
          stepFP / err0, 1.0, 1e-6);
    check("三路径：LowRank 首轮步长被 (I+G)⁻¹ 放大（区别于 G=0）",
          stepLow > stepFP * (1.0 + 1e-4) ? 1.0 : 0.0, 1.0, 0.0);

    // 每条路径自身逐轮误差非增（线搜索按实际误差裁决）。
    auto monotone = [](const ErectionGeometryResult3D& r) {
        for (std::size_t i = 1; i < r.trace.steps.size(); ++i)
            if (r.trace.steps[i].err > r.trace.steps[i - 1].err * (1.0 + 1e-9))
                return 0.0;
        return 1.0;
    };
    check("三路径：LowRank 逐轮误差非增", monotone(low), 1.0, 0.0);
    check("三路径：Full 逐轮误差非增", monotone(full), 1.0, 0.0);
    check("三路径：FixedPoint 逐轮误差非增", monotone(fp), 1.0, 0.0);

    // 终值一致：制造线形与索力（规格第 3 条）。
    auto maxDiff = [](const std::vector<double>& a, const std::vector<double>& b) {
        double v = 0.0;
        for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
            v = std::max(v, std::abs(a[i] - b[i]));
        return v;
    };
    const double camberLF = maxDiff(low.camber, full.camber);
    const double camberLP = maxDiff(low.camber, fp.camber);
    // 取证：制造线形若出现差异，定位到具体自由度，看是否属于「闭合方程
    // 规范自由度」（I+G 在该方向零空间、ΔX+u 对它不敏感）。
    {
        const DofMap3D dmg = buildDofMap3D(m);
        Locator loc{dmg};
        std::vector<std::pair<double, int>> diffs;
        for (std::size_t i = 0; i < low.camber.size(); ++i)
            diffs.push_back({std::abs(low.camber[i] - full.camber[i]), (int)i});
        std::sort(diffs.begin(), diffs.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });
        std::printf("  制造线形差异最大自由度（node/dof, LowRank, Full, "
                    "LowRank.u, Full.u, ΔX+u 两条路径）：\n");
        for (int q = 0; q < std::min(5, (int)diffs.size()); ++q) {
            int nid = -1, c = -1;
            loc(diffs[q].second, nid, c);
            const int d = diffs[q].second;
            std::printf("    node=%d comp=%d  Δ=%.4e  camber=(%.6f, %.6f)  "
                        "u=(%.6f, %.6f)  Fv=(%.3e, %.3e)\n",
                        nid, c, diffs[q].first,
                        low.camber[d], full.camber[d], low.u[d], full.u[d],
                        low.camber[d] + low.u[d], full.camber[d] + full.u[d]);
        }
    }
    check("终值：Full 与 LowRank 制造线形一致 (m)", camberLF, 0.0, 1e-6);
    check("终值：FixedPoint 与 LowRank 制造线形一致 (m)", camberLP, 0.0, 1e-6);
    double tenLF = 0.0, tenLP = 0.0;
    for (std::size_t k = 0; k < low.tension.size(); ++k) {
        const double den = std::max(std::abs(low.tension[k]), 1.0);
        tenLF = std::max(tenLF, std::abs(full.tension[k] - low.tension[k]) / den);
        tenLP = std::max(tenLP, std::abs(fp.tension[k] - low.tension[k]) / den);
    }
    check("终值：Full 与 LowRank 索力一致 (相对)", tenLF, 0.0, 1e-6);
    check("终值：FixedPoint 与 LowRank 索力一致 (相对)", tenLP, 0.0, 1e-6);
    check("终值：三条路径闭合误差均 <= tolClosure",
          std::max({low.maxClosureErr, full.maxClosureErr, fp.maxClosureErr}),
          0.0, base.tolClosure);
    std::printf("  终值差：camber(L-F)=%.3e m  camber(L-P)=%.3e m  "
                "T(L-F)=%.3e  T(L-P)=%.3e\n",
                camberLF, camberLP, tenLF, tenLP);

    // 默认路径回归保护：显式 LowRank 与默认构造逐位相同。
    ErectionGeometryOptions3D def;
    def.maxIter = base.maxIter;
    def.tolClosure = base.tolClosure;
    def.nl = base.nl;
    const ErectionGeometryResult3D dflt = iterateErectionGeometry(m, ids, def);
    check("默认路径：与显式 LowRank 制造线形逐位相同",
          maxDiff(dflt.camber, low.camber), 0.0, 0.0);
    check("默认路径：与显式 LowRank 轮数相同",
          static_cast<double>(dflt.iterations), static_cast<double>(low.iterations), 0.0);
    check("默认路径：与显式 LowRank 闭合误差逐位相同",
          std::abs(dflt.maxClosureErr - low.maxClosureErr), 0.0, 0.0);
}

// ---------------------------------------------------------------
// 补充取证：几何闭合的刚体平移规范自由度（不是 bug，是方程本身的性质）
//
// 若整段梁的所有节点都没有竖向/平面内转动约束，把参考构形整体平移 δ、
// 位移取 u = -δ，则当前构形、CR 梁内力、索的当前弦长与索力全都不变。
// 于是 Fv = ΔX + u 对这一方向不敏感，(I+G) 有零空间：
//   · 成桥构形与索力唯一；
//   · 但制造线形 ΔX 是一族解，不同 Newton 路径会停在不同代表点。
// 生产算例里塔/墩的固结约束天然打断它；测试模型必须显式固结，
// 否则「三路径制造线形一致」这条 oracle 本身不成立。
// ---------------------------------------------------------------
void testGaugeFreedom(double L, double H, double P) {
    std::printf("\n[补充取证] 无端部固结时的刚体平移规范自由度\n");
    Model3D m = makeModel(L, H, 3.45e10, 6.0, 4.0, 1.95e11, 5.0e-3, P, false);
    setTargetCableL0(m, P);
    const std::vector<int> ids = {10, 11, 12};

    ErectionGeometryOptions3D oLow = pathOptions();
    oLow.sensitivity = SensitivityMode3D::LowRank;
    oLow.maxIter = 40;
    oLow.tolClosure = 1e-9;
    const ErectionGeometryResult3D low = iterateErectionGeometry(m, ids, oLow);
    ErectionGeometryOptions3D oFull = pathOptions();
    oFull.sensitivity = SensitivityMode3D::Full;
    oFull.maxIter = 40;
    oFull.tolClosure = 1e-9;
    const ErectionGeometryResult3D full = iterateErectionGeometry(m, ids, oFull);

    check("规范自由度：LowRank 收敛", low.converged ? 1.0 : 0.0, 1.0, 0.0);
    check("规范自由度：Full 收敛", full.converged ? 1.0 : 0.0, 1.0, 0.0);
    double tenDiff = 0.0;
    for (std::size_t k = 0; k < low.tension.size(); ++k)
        tenDiff = std::max(tenDiff, std::abs(full.tension[k] - low.tension[k]) /
                           std::max(std::abs(low.tension[k]), 1.0));
    check("规范自由度：索力仍唯一（相对）", tenDiff, 0.0, 1e-5);

    const DofMap3D dm = buildDofMap3D(m);
    std::vector<double> diffs;
    for (int n = 1; n <= 5; ++n) {
        const int d = dm.dof(n, 2);              // z（0-based）
        diffs.push_back(low.camber[d] - full.camber[d]);
    }
    const double mean = std::accumulate(diffs.begin(), diffs.end(), 0.0) /
                        static_cast<double>(diffs.size());
    double dev = 0.0, curDiff = 0.0;
    for (int n = 1; n <= 5; ++n) {
        const int d = dm.dof(n, 2);
        dev = std::max(dev, std::abs(diffs[n - 1] - mean));
        curDiff = std::max(curDiff,
            std::abs((low.camber[d] + low.u[d]) - (full.camber[d] + full.u[d])));
    }
    check("规范自由度：camber 差是均匀平移（偏差/平移量）",
          dev / std::max(std::abs(mean), 1e-30), 0.0, 1e-3);
    check("规范自由度：平移量显著（确实存在零空间）",
          std::abs(mean) > 1e-3 ? 1.0 : 0.0, 1.0, 0.0);
    check("规范自由度：当前成桥构形一致 (m)", curDiff, 0.0, 1e-8);
    std::printf("  两路径 camber 差均值=%.6f m（刚体平移），最大偏差=%.3e m；"
                "成桥构形差=%.3e m\n", mean, dev, curDiff);
}

}  // namespace

int main() {
    setupConsoleUTF8();
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("========== P3-3b 完整形状敏感度 K_shape 验收 oracle ==========\n");
    try {
        const double L = 30.0, H = 40.0, P = 1.0e6;
        Model3D m = makeModel(L, H, 3.45e10, 6.0, 4.0, 1.95e11, 5.0e-3, P);
        // 给三根索下料长度（刚性支承连续梁法目标索力），使体系可解。
        setTargetCableL0(m, P);
        const DofMap3D dm = buildDofMap3D(m);
        const std::vector<char> isT = translationDofs(m, dm);

        NonlinearOptions3D nl;
        nl.loadSteps = 2; nl.maxIter = 40;
        nl.tolForce = 1e-11; nl.tolDisp = 1e-12;
        NLTangent3D tan;
        const NonlinearResult3D r0 = solveModel3DNonlinear(m, nl, &tan);
        check("基准：非线性求解收敛", r0.converged ? 1.0 : 0.0, 1.0, 0.0);

        ShapeSensitivityResult3D ss;
        testResidualAndShapeTerm(m, dm, isT, tan, r0, ss, P);
        testDirectionalAgainstFD(m, dm, isT, r0, ss, nl);
        testThreePaths(m);
        testGaugeFreedom(L, H, P);
    } catch (const std::exception& e) {
        std::printf("\n[ERROR] 异常: %s\n", e.what());
        ++g_fail;
    }

    std::printf("\n====================================\n");
    std::printf("通过 %d 项，失败 %d 项\n", g_pass, g_fail);
    if (g_fail != 0) {
        std::printf("结果：FAIL —— 完整形状敏感度未通过验收\n");
        return 1;
    }
    std::printf("结果：ALL PASS —— 完整形状敏感度通过验收\n");
    return 0;
}
