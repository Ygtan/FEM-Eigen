#pragma once
// ============================================================================
// 成桥目标索力优化（带约束），以及结果文件的 [OPT] 段格式化。
//
// 与 `cableTensionsRigidSupport` / `cableTensionsMinBendingEnergy` 的区别：
//   两者是「按一个固定准则直接算一组索力」；本文件是「给定可行域，求一组
//   更优的索力」。RIGID 只解一次线性静力取支反力，MINENERGY 解无约束二次
//   规划（只有一个下限主动集）。这里补上三件它们都没有的东西：
//     · 目标可选 `min max|M主梁|`（切比雪夫，IRLS 逼近）或最小弯曲能量；
//     · 塔顶水平力平衡（逐塔一条 ΣT·cos 约束，用罚函数并入二次型）；
//     · 索力上下界（Tmax 是新增能力，MINENERGY 只有下限）。
//
// 力学全部复用内核已有积木：裸结构线性刚度 K、单位索力对 Ld、
// `Factorization3D::solveMulti` 给出 U=K⁻¹Ld、恒载位移 u0=K⁻¹Fdl、
// `beam3DLocalEndForces` + `beamEndForcesAt` 给出单位索力下的梁端弯矩。
// **不重跑非线性求解，也不在 Python 侧另算一套索力分配**。
//
// 有意不抽公共函数去改 MINENERGY：默认路径（OPTIMIZE=0）必须逐位不变，
// 复制约 80 行装配逻辑比冒着改动已验证数值路径的风险更划算。
// ============================================================================

#include "fem_3d_nl.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// 内核所有类型都在 `fem` 命名空间里（fem_3d.hpp:88 / fem_3d_nl.hpp:88）。
// 本头文件同样放进去，调用方靠 `using namespace fem;` 取用，与其它
// *_output_3d.hpp / fem_3d_nl.hpp 的口径一致。
namespace fem {

struct CableOptOptions3D {
    enum class Objective { MinBendingEnergy, MinMaxMoment };
    Objective objective = Objective::MinMaxMoment;
    bool towerBalance = true;       // 逐塔 ΣT·cos = 0
    // 平衡是**软约束**，权重必须相对弯矩项配平（见 addBalance）：1.0 表示
    // 与弯矩项同量级。实测权重过大时优化器只顾 ΣT·cos=0、把主梁弯矩推高，
    // 成桥态反而变差；0.1 能压住塔顶不平衡又不压倒弯矩目标。
    double balanceScale = 0.1;      // 平衡项相对弯矩项的权重（无量纲）
    double tMin = 0.0;              // 索力绝对下限 N；0 = 只用比例下限
    double tMax = 0.0;              // 索力绝对上限 N；0 = 用程序默认上界
    double minTensionFraction = 0.02;
    int maxIter = 400;              // MINENERGY 单轮投影梯度上限
    double tol = 1e-7;              // 相对步长 / 目标 / 弯矩峰值收敛判据
    int irlsSteps = 4;              // MinMaxMoment 的 IRLS 轮数（p = 2,4,8,16）
    int innerIter = 120;            // 每轮 IRLS 的投影梯度上限
    double roundTol = 1e-5;         // IRLS 相邻轮解的相对变化收敛判据
    double regularization = 1e-9;   // 相对正则，保证 Q 严格凸
    bool verbose = false;
};

struct CableOptStep3D {
    int irls = 0;                   // IRLS 轮次（1..）；MINENERGY 时为 0
    int iter = 0;                   // 该轮内投影梯度迭代
    double objective = 0.0;         // 二次目标值（含正则与平衡罚）
    double momentPeak = 0.0;        // max|My主梁|（N·m）
    double balance = 0.0;           // ||H·T||（N，塔顶水平合力）
    double step = 0.0;              // 相对步长
};

struct CableOptResult3D {
    bool converged = false;
    std::string status = "ok";      // ok | not_converged | <error code>
    std::string error;              // 非空 = 本次没有可用优化结果
    std::string objective = "MINMAXM";
    std::string start = "RIGID";
    std::vector<double> T0, T, tLo, tHi;
    int nc = 0, nm = 0, nb = 0;
    int iterations = 0, irlsSteps = 0;
    std::string convergedBy = "";   // projected_gradient | irls_round_stable
    // MPEAK0 / MPEAK 是**裸结构线性化影响矩阵**口径的目标量，只用于优化；
    // 成桥态实际弯矩由 cable_solver 从本次 [GIRDEREND] 回填到
    // momentPeakFinal。二者不是同一个量，段里分别标 KIND，报告必须并排
    // 写出，不能拿线性化目标冒充成桥态效应。未优化基线由调用方另跑一次
    // OPTIMIZE=0 取 [GIRDEREND] max|My| 对拍，内核不自己再造一份。
    double momentPeak0 = 0.0, momentPeak = 0.0;
    double momentPeakFinal = 0.0;    // 本次（优化后）成桥态 [GIRDEREND] max|My|
    double balance0 = 0.0, balance = 0.0;
    double tMaxUsed = 0.0;
    std::string tMaxSource = "DEFAULT";
    std::vector<CableOptStep3D> steps;
};

namespace cable_opt_detail {

inline std::vector<double> toStd(const Eigen::VectorXd& v) {
    return std::vector<double>(v.data(), v.data() + v.size());
}

}  // namespace cable_opt_detail

inline CableOptResult3D optimizeCableTensions3D(
    const Model3D& m, const std::vector<CableAnchor3D>& cables,
    const CableOptOptions3D& opt = {},
    const std::vector<double>& givenT = {}) {
    using cable_opt_detail::toStd;
    CableOptResult3D out;
    const int nc = static_cast<int>(cables.size());
    out.nc = nc;
    out.objective = (opt.objective == CableOptOptions3D::Objective::MinMaxMoment)
                        ? "MINMAXM" : "MINENERGY";
    out.irlsSteps = (out.objective == "MINMAXM") ? std::max(1, opt.irlsSteps) : 0;
    if (nc == 0) {
        out.status = out.error = "no_cables";
        return out;
    }

    // ---- ① 裸结构：去掉待定索，孤立节点补约束（与 MINENERGY 同口径） ----
    Model3D mc = m;
    {
        std::vector<Element3D> keep;
        keep.reserve(mc.elements.size());
        for (const auto& e : mc.elements) {
            bool isTgt = false;
            for (const auto& c : cables)
                if (c.elem == e.id) { isTgt = true; break; }
            if (!isTgt) keep.push_back(e);
        }
        mc.elements = keep;
    }
    {
        std::map<int, bool> connected;
        for (const auto& e : mc.elements) { connected[e.n1] = true; connected[e.n2] = true; }
        const DofMap3D dm0 = buildDofMap3D(mc);
        for (const auto& [id, p] : mc.nodes) {
            (void)p;
            if (connected[id]) continue;
            for (int k = 0; k < dm0.count(id); ++k) {
                bool hasSpring = false;
                for (const auto& c : mc.constraints)
                    if (c.node == id && c.dof == k + 1 && c.k > 0.0) {
                        hasSpring = true; break;
                    }
                if (!hasSpring) mc.constraints.push_back({id, k + 1, 0.0, 0.0});
            }
        }
    }
    // assembleStiffness3D 的签名要非 const DofMap3D&（fem_3d.hpp:787）。
    DofMap3D dm = buildDofMap3D(mc);
    const int nd = dm.nDofs;
    SparseMat3D K = assembleStiffness3D(mc, dm);
    Eigen::VectorXd F0 = Eigen::VectorXd::Zero(nd);
    for (const auto& [node, dof, val] : mc.loads) {
        const int d = dm.safeDof(node, dof - 1);
        if (d >= 0) F0(d) += val;
    }
    for (const auto& c : mc.constraints) {
        if (c.k <= 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0) F0(d) += c.k * c.value;
    }
    std::vector<int> cdofs;
    std::vector<double> cvals;
    for (const auto& c : mc.constraints) {
        if (c.k > 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        auto it = std::find(cdofs.begin(), cdofs.end(), d);
        if (it == cdofs.end()) { cdofs.push_back(d); cvals.push_back(c.value); }
        else cvals[std::distance(cdofs.begin(), it)] = c.value;
    }
    Eigen::VectorXd Fdl = F0;
    applyDirichlet(K, Fdl, cdofs, cvals);
    const auto isC = [&](int d) {
        return std::find(cdofs.begin(), cdofs.end(), d) != cdofs.end();
    };

    // ---- ② 单位索力对 Ld、索方向、塔侧节点 ----
    Eigen::MatrixXd Ld = Eigen::MatrixXd::Zero(nd, nc);
    std::vector<Eigen::Vector3d> ev(nc);
    std::vector<int> bNode(nc, 0);
    for (int k = 0; k < nc; ++k) {
        const Element3D* e = nullptr;
        for (const auto& q : m.elements)
            if (q.id == cables[k].elem) { e = &q; break; }
        if (!e) { out.status = out.error = "cable_element_missing"; return out; }
        const int A = cables[k].anchorNode;
        const int B = (e->n1 == A) ? e->n2 : e->n1;
        if (A != e->n1 && A != e->n2) {
            out.status = out.error = "anchor_not_on_cable";
            return out;
        }
        const Vec3 d3 = m.nodes.at(B) - m.nodes.at(A);
        const double L = d3.norm();
        if (L < 1e-12) { out.status = out.error = "degenerate_cable"; return out; }
        ev[k] = Eigen::Vector3d(d3.x(), d3.y(), d3.z()) / L;
        bNode[k] = B;
        const int sA = dm.dofStart[dm.index.at(A)];
        const int sB = dm.dofStart[dm.index.at(B)];
        for (int a = 0; a < 3; ++a) {
            const int dA = sA + a, dB = sB + a;
            if (!isC(dA)) Ld(dA, k) = +ev[k][a];
            if (!isC(dB)) Ld(dB, k) = -ev[k][a];
        }
    }

    Factorization3D fac(K);
    // 去索后若为机构，柔度矩阵无意义（MINENERGY 用同一判据，这里也照做）。
    {
        double eworst = 0.0;
        for (int probe = 0; probe < 2; ++probe) {
            Eigen::VectorXd r = Eigen::VectorXd::Zero(nd);
            for (int d = 0; d < nd; ++d) {
                if (isC(d)) continue;
                r(d) = (probe == 0) ? std::sin(1.7 * d + 0.3)
                                    : std::cos(0.31 * d + 1.1);
            }
            if (r.norm() < 1e-12) continue;
            const Eigen::VectorXd z = fac.solve(K * r);
            eworst = std::max(eworst, (z - r).norm() / r.norm());
        }
        if (!(eworst < 1e-4)) {
            out.status = out.error = "bare_structure_is_mechanism";
            return out;
        }
    }
    const Eigen::VectorXd u0 = fac.solve(Fdl);
    const Eigen::MatrixXd U = fac.solveMulti(Ld);     // nd×nc：单位索力工况位移
    Eigen::MatrixXd Fm = Ld.transpose() * U;
    Fm = 0.5 * (Fm + Fm.transpose());
    const Eigen::VectorXd g = Ld.transpose() * u0;    // 恒载 × 单位索力交叉项

    // ---- ③ 起点：给定索力优先，否则刚性支承连续梁法 --------------------
    // givenT 来自 *CABLE 段用户显式给的目标索力。全给齐且都为正时用它，
    // 优化在"用户方案附近"找可行改进；否则退回 RIGID 解析起点
    // （已有 11P/8、5P/16 解析 oracle）。start 如实写进 [OPT] 段，
    // 不能一边用 RIGID 一边标 GIVEN。
    Eigen::VectorXd T0 = Eigen::VectorXd::Zero(nc);
    const bool hasGiven = (givenT.size() == static_cast<std::size_t>(nc)) &&
        std::all_of(givenT.begin(), givenT.end(),
                    [](double v) { return v > 0.0; });
    if (hasGiven) {
        for (int k = 0; k < nc; ++k) T0(k) = givenT[k];
        out.start = "GIVEN";
    } else {
        try {
            const std::vector<double> Tr = cableTensionsRigidSupport(m, cables);
            for (int k = 0; k < nc; ++k) T0(k) = std::max(Tr[k], 0.0);
        } catch (const std::exception& ex) {
            (void)ex;
            out.status = out.error = "start_rigid_failed";
            return out;
        }
        out.start = "RIGID";
    }
    const double t0max = T0.maxCoeff();
    if (!(t0max > 0.0)) {
        out.status = out.error = "start_all_zero";
        return out;
    }

    // ---- ④ 上下界 ----
    const double floorAbs = std::max(opt.tMin, opt.minTensionFraction * t0max);
    const Eigen::VectorXd lo = Eigen::VectorXd::Constant(nc, std::max(0.0, floorAbs));
    // 不给上界时，实测优化器会把索力集中到少数索上（48 索算例从
    // 0.68~1.40 MN 被推到 0.14~6.67 MN），工程上不可用。程序默认上界取
    // 1.5 × RIGID 起点最大索力，属**程序口径不是规范限值**；需要别的用
    // TMAX= 显式给，段里会写明来源。
    out.tMaxSource = (opt.tMax > 0.0) ? "USER" : "DEFAULT_1.5_RIGID";
    out.tMaxUsed = (opt.tMax > 0.0) ? opt.tMax : 1.5 * t0max;
    Eigen::VectorXd hi = Eigen::VectorXd::Constant(nc, out.tMaxUsed);
    for (int k = 0; k < nc; ++k)
        hi(k) = std::max(hi(k), lo(k) * (1.0 + 1e-9) + 1e-6);
    Eigen::VectorXd T = T0.cwiseMax(lo).cwiseMin(hi);
    out.T0 = toStd(T0);
    out.tLo = toStd(lo);
    out.tHi = toStd(hi);

    // ---- ⑤ 主梁弯矩影响矩阵 A 与恒载基准 M0 ----
    // 只取「水平梁」：塔柱/斜撑不是主梁，混进来会让 min max|M| 变成压塔弯矩。
    struct BeamRow { const Element3D* e; int end; };
    std::vector<BeamRow> rows;
    for (const auto& e : mc.elements) {
        if (e.type != "BEAM") continue;
        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
        const double lh = std::hypot(d.x(), d.y());
        if (std::abs(d.z()) > 0.2 * std::max(lh, 1e-9)) continue;
        rows.push_back({&e, 1});
        rows.push_back({&e, 2});
    }
    const int nm = static_cast<int>(rows.size());
    out.nm = nm;
    if (nm == 0 && out.objective == "MINMAXM") {
        out.status = out.error = "no_girder_moments";
        return out;
    }
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nm, nc);
    Eigen::VectorXd M0 = Eigen::VectorXd::Zero(nm);
    const std::vector<double> u0s = toStd(u0);
    std::vector<std::vector<double>> ucols(nc);
    for (int k = 0; k < nc; ++k) {
        const Eigen::VectorXd uk = U.col(k);
        ucols[k] = toStd(uk);
    }
    for (int i = 0; i < nm; ++i) {
        const Element3D& e = *rows[i].e;
        M0(i) = beamEndForcesAt(
            beam3DLocalEndForces(mc, dm, e, u0s), 12, rows[i].end).My;
        for (int k = 0; k < nc; ++k)
            A(i, k) = beamEndForcesAt(
                beam3DLocalEndForces(mc, dm, e, ucols[k]), 12, rows[i].end).My;
    }

    // ---- ⑥ 塔顶水平平衡 H：按塔侧节点 x 分组，逐塔一行 ----
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(0, nc);
    std::vector<int> group(nc, 0);
    if (opt.towerBalance) {
        std::vector<double> towerX;
        for (int k = 0; k < nc; ++k) {
            const double xb = m.nodes.at(bNode[k]).x();
            int gi = -1;
            for (int j = 0; j < static_cast<int>(towerX.size()); ++j)
                if (std::abs(xb - towerX[j]) < 0.5) { gi = j; break; }
            if (gi < 0) { gi = static_cast<int>(towerX.size()); towerX.push_back(xb); }
            group[k] = gi;
        }
        const int nb = static_cast<int>(towerX.size());
        H = Eigen::MatrixXd::Zero(nb, nc);
        for (int k = 0; k < nc; ++k)
            // 塔侧节点受 −T·e（锚点 A 是 +T·e），水平分量 −T·e_x。
            H(group[k], k) = -ev[k].x();
        // 全零行（该塔所有索都竖直）不构成约束，去掉。用临时矩阵搬，
        // 避免 `H.row(i) = H.row(j)` 的别名问题。
        Eigen::MatrixXd Hkeep = Eigen::MatrixXd::Zero(nb, nc);
        int keep = 0;
        for (int j = 0; j < nb; ++j)
            if (H.row(j).cwiseAbs().maxCoeff() > 1e-12)
                Hkeep.row(keep++) = H.row(j);
        H = Hkeep.topRows(keep);
        out.nb = keep;
    }

    auto momentPeak = [&](const Eigen::VectorXd& Tv) -> double {
        if (nm == 0) return 0.0;
        return (M0 + A * Tv).cwiseAbs().maxCoeff();
    };
    auto balanceNorm = [&](const Eigen::VectorXd& Tv) -> double {
        if (out.nb == 0) return 0.0;
        return (H * Tv).norm();
    };
    out.momentPeak0 = momentPeak(T0);
    out.balance0 = balanceNorm(T0);

    // ---- ⑦ 无量纲化：优化变量 x = T / Tscale ----------------------------
    // 弯矩项 ~(N·m)²、平衡项 ~N²、变量本身 ~1e6 N，直接拼进同一个二次型
    // 会差十几个数量级；旧写法还把平衡罚乘上弯矩项的 qscale，罚项彻底
    // 淹没目标，实测出现"收敛了但主梁弯矩反而变大"。全部归一到 O(1)：
    //   Ã = (Tscale/Mscale)·A,  m̃0 = M0/Mscale
    //   H̃ = (Tscale/Bscale)·H
    //   f(x) = ‖m̃0 + Ãx‖²_W/(2·nm) + balanceScale·‖H̃x‖²/(2·nb)
    // 箱约束、步长、收敛判据都在 x 空间；输出时 T = Tscale·x 还原成 N。
    // 目标值本身是无量纲的，[OPT] 段以 OBJECTIVE_KIND=NORMALIZED 标注。
    const double Tscale = std::max(t0max, 1.0);
    const double Mscale = std::max(out.momentPeak0, 1.0);
    const Eigen::MatrixXd Ahat = A * (Tscale / Mscale);
    const Eigen::VectorXd m0hat = M0 / Mscale;
    // H 是无量纲 cos 分量，直接用；平衡残差 BAL 仍由真实 T = Tscale·X 算。
    const Eigen::VectorXd lox = lo / Tscale, hix = hi / Tscale;
    Eigen::VectorXd X = T / Tscale;
    const double wDenom = static_cast<double>(std::max(1, nm));
    const double hDenom = static_cast<double>(std::max(1, out.nb));
    const bool hasBalance = opt.towerBalance && out.nb > 0;
    const Eigen::MatrixXd HtH = H.transpose() * H;
    const double hscale = std::max(HtH.cwiseAbs().maxCoeff(), 1e-30);

    // 平衡罚必须**相对弯矩项配平**：两者量纲不同，直接加会让平衡项
    // 比弯矩项大几百倍，优化器只顾 ΣT·cos=0 而把弯矩推高。
    // qscale 用**固定的基准弯矩项**尺度（不含 IRLS 权重），否则 IRLS
    // 每轮 w 变化会让平衡项相对权重逐轮跳变，切比雪夫目标被推坏。
    const double qscaleFixed = std::max(
        (Ahat.transpose() * Ahat / wDenom).cwiseAbs().maxCoeff(), 1e-30);
    auto addBalance = [&](Eigen::MatrixXd& Q) {
        if (!hasBalance) return;
        Q += (opt.balanceScale * qscaleFixed / (hscale * hDenom)) * HtH;
    };
    auto makeMomentQP = [&](const Eigen::VectorXd& w) {
        const Eigen::MatrixXd WA = w.asDiagonal() * Ahat;
        Eigen::MatrixXd Q = (Ahat.transpose() * WA) / wDenom;
        Eigen::VectorXd q = (Ahat.transpose() * (w.asDiagonal() * m0hat)) / wDenom;
        Q += Eigen::MatrixXd::Identity(nc, nc) * opt.regularization;
        addBalance(Q);
        return std::make_pair(Q, q);
    };
    // MINENERGY 用裸结构弯曲应变能柔度 Fm、交叉项 g，同样按 Tscale 归一：
    //   Qe = Fm·Tscale²/Escale,  qe = g·Tscale/Escale
    const double Escale = std::max(Fm.cwiseAbs().maxCoeff() * Tscale * Tscale, 1.0);
    auto makeEnergyQP = [&]() {
        Eigen::MatrixXd Q = Fm * (Tscale * Tscale / Escale);
        Eigen::VectorXd q = g * (Tscale / Escale);
        Q += Eigen::MatrixXd::Identity(nc, nc) * opt.regularization;
        addBalance(Q);
        return std::make_pair(Q, q);
    };
    auto objectiveValue = [](const Eigen::VectorXd& Xv,
                             const Eigen::MatrixXd& Q,
                             const Eigen::VectorXd& q) -> double {
        return 0.5 * Xv.dot(Q * Xv) + q.dot(Xv);
    };

    // 投影梯度 + Jacobi 对角预条件 + Armijo 回溯，迭代变量是无量纲 X。
    // 试过无约束 Newton 方向：箱约束下全步长会一步把大半索力推到上界、
    // 弯矩反而暴涨，故不用。无量纲化后 Q 条件数回到可控范围，步长不再
    // 被 1e18 量级的罚项压死。这是**修量纲**而不是放宽收敛判据。
    auto solveQP = [&](const Eigen::MatrixXd& Q, const Eigen::VectorXd& q,
                       Eigen::VectorXd& Xv, int irls, int iterCap) -> bool {
        const Eigen::VectorXd diag = Q.diagonal().cwiseMax(1e-30);
        bool converged = false;
        double prevMpeak = momentPeak(Tscale * Xv);
        for (int it = 1; it <= iterCap; ++it) {
            const Eigen::VectorXd grad = Q * Xv + q;
            const Eigen::VectorXd d = -grad.cwiseQuotient(diag);  // 预条件下降方向
            const double denom = d.dot(Q * d);
            // 通用精确线搜索：d 不再等于 -grad，分子必须是 -(d·grad)。
            double alpha = (denom > 0.0)
                               ? (-d.dot(grad) / denom)
                               : 1.0 / std::max(1.0, Q.cwiseAbs().maxCoeff());
            const double J0 = objectiveValue(Xv, Q, q);
            Eigen::VectorXd Xt;
            bool ok = false;
            for (int ls = 0; ls < 40; ++ls) {
                Xt = (Xv + alpha * d).cwiseMax(lox).cwiseMin(hix);
                // 标准 Armijo：实际下降必须达到一阶预测下降的 1e-4。
                // 旧判据 `J0 - 1e-14·|J0|` 只要求"任意微降"，投影后大幅
                // 撞到箱边界的坏步也被放行，第一步就把弯矩峰值推高一倍。
                const double predicted = grad.dot(Xt - Xv);   // < 0
                if (predicted < 0.0 &&
                    objectiveValue(Xt, Q, q) <= J0 + 1e-4 * predicted) {
                    ok = true;
                    break;
                }
                alpha *= 0.5;
            }
            if (!ok) break;                       // 已到投影驻点
            const double step = (Xt - Xv).norm() / std::max(1.0, Xv.norm());
            Xv = Xt;
            ++out.iterations;
            const double Jt = objectiveValue(Xv, Q, q);
            const double mpeak = momentPeak(Tscale * Xv);
            const double mrel = std::abs(mpeak - prevMpeak)
                                / std::max(1.0, std::abs(prevMpeak));
            const double jrel = std::abs(Jt - J0) / std::max(1.0, std::abs(J0));
            prevMpeak = mpeak;
            const bool done = step < opt.tol
                              || (it > 5 && (mrel < opt.tol || jrel < opt.tol));
            if (it == 1 || it % 20 == 0 || done) {
                CableOptStep3D s;
                s.irls = irls;
                s.iter = it;
                s.objective = Jt;
                s.momentPeak = mpeak;
                s.balance = balanceNorm(Tscale * Xv);
                s.step = step;
                out.steps.push_back(s);
            }
            if (done) { converged = true; break; }
        }
        return converged;
    };

    bool converged = false;
    if (out.objective == "MINMAXM") {
        // IRLS：p = 2 → 4 → 8 → 16，权重 ∝ |M|^(p−2) 逐轮压平弯矩峰值。
        // 轮间以 roundTol 判断解是否稳定（IRLS 的标准停止准则），稳定即
        // 提前收束，不必把 4 轮 × innerIter 全跑满。
        for (int r = 0; r < out.irlsSteps; ++r) {
            const Eigen::VectorXd Xprev = X;
            const double p = 2.0 * std::pow(2.0, r);
            const Eigen::VectorXd Mv = M0 + A * (Tscale * X);
            const double scale = std::max(Mv.cwiseAbs().maxCoeff(), 1e-9);
            Eigen::VectorXd w(nm);
            for (int i = 0; i < nm; ++i)
                w(i) = std::max(std::pow(std::abs(Mv(i)) / scale, p - 2.0), 1e-6);
            w /= std::max(w.mean(), 1e-30);
            const auto Qq = makeMomentQP(w);
            const bool inner = solveQP(Qq.first, Qq.second, X, r + 1,
                                       opt.innerIter);
            const double round = (X - Xprev).norm()
                                 / std::max(1.0, Xprev.norm());
            // 只看**最后一轮**的收敛状态：某一轮达标、后续轮又漂移，
            // 不能拿旧轮的结果冒充整体收敛。
            converged = inner;
            out.convergedBy = inner ? "projected_gradient" : "";
            if (round < opt.roundTol) {
                // 轮间解稳定是 IRLS 的标准停止准则：再升 p 已无实质变化。
                converged = true;
                out.convergedBy = "irls_round_stable";
                break;
            }
        }
    } else {
        const auto Qq = makeEnergyQP();
        converged = solveQP(Qq.first, Qq.second, X, 0, opt.maxIter);
        if (converged) out.convergedBy = "projected_gradient";
    }
    T = Tscale * X;                              // 还原成 N
    out.converged = converged;
    out.status = converged ? "ok" : "not_converged";
    out.T = toStd(T);
    out.momentPeak = momentPeak(T);
    out.balance = balanceNorm(T);
    // 守卫：线性化目标都不降反升时，绝不把 T* 写回 *CABLE 冒充优化结果。
    // 未收敛但目标确实改善仍可采纳，STATUS 如实保留 not_converged。
    if (!(out.momentPeak <= out.momentPeak0 * (1.0 + 1e-6)))
        out.status = "objective_not_improved";
    if (opt.verbose)
        std::printf("      [opt] objective=%s nc=%d nm=%d nb=%d niter=%d "
                    "Mpeak %.6g -> %.6g  bal %.6g -> %.6g  status=%s\n",
                    out.objective.c_str(), nc, nm, out.nb, out.iterations,
                    out.momentPeak0, out.momentPeak, out.balance0, out.balance,
                    out.status.c_str());
    return out;
}

// ---------------------------------------------------------------------------
// [OPT] 段：只格式化内核已经算出的量，不做任何力学判断。
// ---------------------------------------------------------------------------
inline std::string cableOptStatus3D(const CableOptResult3D& r) {
    if (!r.error.empty()) return r.error;
    return r.status.empty() ? "ok" : r.status;
}

inline std::string formatCableOptSection3D(const CableOptResult3D& r,
                                           const char* header = "[OPT]",
                                           bool applied = false) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << header << "\n";
    ss << "KIND=constrained_cable_force STATUS=" << cableOptStatus3D(r) << "\n";
    if (!r.error.empty()) {
        ss << "ERROR=" << r.error << "\n";
        return ss.str();
    }
    // I 行的 objective 是**无量纲归一化目标**（x=T/Tscale），不是 N·m，
    // 只能看相对下降，不能与 MPEAK*/BAL* 混算或相减。
    ss << "OBJECTIVE=" << r.objective << " START=" << r.start
       << " OBJECTIVE_KIND=NORMALIZED\n";
    ss << "NCABLES=" << r.nc << " NMOMENTS=" << r.nm
       << " NBALANCE=" << r.nb << " IRLS=" << r.irlsSteps << "\n";
    ss << "NITER=" << r.iterations << " CONVERGED=" << (r.converged ? 1 : 0)
       << " CONVBY=" << (r.convergedBy.empty() ? "none" : r.convergedBy)
       // APPLIED=0 表示这组 T* 没被写回 *CABLE（目标未改善或求解失败），
       // 段里的 T 行只是过程量，不能拿去当设计索力。
       << " APPLIED=" << (applied ? 1 : 0) << "\n";
    ss << "TMIN=" << (r.tLo.empty() ? 0.0 : r.tLo[0])
       << " TMAX=" << r.tMaxUsed << " TMAXSRC=" << r.tMaxSource << "\n";
    // 两组弯矩峰值**分属不同口径**，段里各写一行 KIND：
    //   LINEAR_BARE  = 索去掉后的裸结构线性化影响矩阵，仅优化目标用；
    //   GIRDEREND_FINAL = 本次成桥态 [GIRDEREND] 实际 max|My|。
    // 实测可能反向变化，绝不能拿前者冒充后者，也不能互相相减当"改善量"。
    ss << "MPEAK_KIND=LINEAR_BARE\n";
    ss << "MPEAK0=" << r.momentPeak0 << " MPEAK=" << r.momentPeak
       << " BAL0=" << r.balance0 << " BAL=" << r.balance << "\n";
    // 只有 APPLIED=1 才回填了成桥态 [GIRDEREND] 峰值；未采纳时写
    // NOT_APPLICABLE 而不是 0，避免下游把 0 当成"弯矩被优化没了"。
    ss << "MPEAK_FINAL_KIND="
       << (applied ? "GIRDEREND_FINAL" : "NOT_APPLICABLE") << "\n";
    ss << "MPEAK_FINAL=" << (applied ? r.momentPeakFinal : 0.0) << "\n";
    for (const auto& s : r.steps)
        ss << "I " << s.irls << " " << s.iter << " " << s.objective << " "
           << s.momentPeak << " " << s.balance << " " << s.step << "\n";
    for (int k = 0; k < r.nc; ++k)
        ss << "T " << (k + 1) << " " << r.T0[k] << " " << r.T[k] << " "
           << r.tLo[k] << " " << r.tHi[k] << "\n";
    return ss.str();
}

inline std::string cableOptErrorSection3D(const char* header,
                                          const std::string& code) {
    std::ostringstream ss;
    ss << header << "\n"
       << "KIND=constrained_cable_force STATUS=" << code << "\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

}  // namespace fem
