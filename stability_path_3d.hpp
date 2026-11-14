// stability_path_3d.hpp — 三维共转结构的球面弧长法平衡路径跟踪
//
// 独立 opt-in 能力，不改 fem_3d_nl.hpp::solveModel3DNonlinear（仍是荷载
// 控制 Newton）、不改 NonlinearOptions3D 默认值、不接入 cable_solver.cpp
// 主链路。因此 B 类规格里「主链路尚未接入非线性稳定全过程」的如实边界
// 仍然成立；本文件只服务于独立的科研级对标探针 snap_probe.cpp。
//
// 力学口径：完全复用 fem_3d_nl.hpp 的 CR 单元内力与解析一致切线
// （crMakeElemData / crElementInternal），支持 TRUSS / CABLE / BEAM /
// BEAM7。唯一新增的是把荷载因子 lambda 当作未知量的增广系统：
//
//     R(u,lam) = lam * Fref - Fint(u) = 0                    (nd 个方程)
//     du^T du + psi^2 * dlam^2 = ds^2                        (Crisfield 球面约束)
//
// 每个校正迭代解两次线性系统（同一次 LU 分解、两个右端项）：
//     K * duR = R ,   K * duF = Fref ,   du = duR + dlam * duF
// 代入弧长约束得 dlam 的二次方程 a*dlam^2 + b*dlam + c = 0，按「增量与
// 上一收敛步增量夹角余弦最大」选根，避免极限点处回退。预测步用当前
// 切线方向给 dlam 初值，符号沿用上一步的切线方向。
//
// 越过极限点后 lambda 自然下降、可穿过 0 变为负值（反向荷载），这正是
// 荷载控制 Newton 做不到的下降段 / snap-through 全过程。
#pragma once

#include "fem_3d_nl.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {

struct ArcLengthOptions3D {
    double ds0 = 0.02;          // 初始弧长半径
    double dsMin = 1.0e-7;      // 弧长下限（低于此认为不可继续）
    double dsMax = 0.2;         // 弧长上限
    double psi = 1.0;           // 弧长约束里 lambda 增量的权重
    double tolForce = 1.0e-7;   // 残差相对 |Fref| 的收敛容差
    int maxSteps = 400;         // 最大平衡点数
    int maxIter = 30;           // 每个平衡点最大校正迭代
    int ctrlNode = 0;           // 控制节点（终止判据用）
    int ctrlDof = 2;            // 控制自由度，1~6
    double ctrlTarget = 1.0;    // |u_ctrl| 达到此值即终止
    double lambdaLimit = 1.0e9; // |lambda| 超过此值即终止
    // 轴向单元应变度量：0 = 工程应变 N=EA(L-L0)/L0（fem_3d_nl.hpp 口径），
    // 1 = 对数应变 N=EA*ln(L/L0)（ANSYS LINK180 KEYOPT(2)=1 口径）。
    // 两者都是共转、都用解析一致切线；仅本构不同，用于把 ANSYS 的应变
    // 度量差异与弧长算法差异严格分离。默认 0，不影响任何既有行为。
    int strainMeasure = 0;
};

struct ArcLengthPoint3D {
    int step = 0;
    double lambda = 0.0;
    double ctrlDisp = 0.0;      // u(ctrlDof)，带符号
    double residual = 0.0;      // 收敛时的自由自由度残差范数
    double ds = 0.0;            // 本步实际弧长半径
    double pivotRatio = 1.0;    // 切线奇异指标：min|diag(U)| / max|diag(U)|
    int iterations = 0;
    bool converged = false;
};

struct ArcLengthResult3D {
    bool converged = false;
    int steps = 0;
    double lambda = 0.0;
    std::vector<double> u;                    // 终点位移向量
    std::vector<ArcLengthPoint3D> points;    // 逐平衡点路径
    std::string error;                       // 失败原因（空 = 正常）
};

inline ArcLengthResult3D solveArcLengthPath3D(const Model3D& m,
                                              const ArcLengthOptions3D& opt = {}) {
    ArcLengthResult3D out;
    const DofMap3D dm = buildDofMap3D(m);
    const int nd = dm.nDofs;
    if (nd <= 0) {
        out.error = "empty_dof_map";
        return out;
    }

    // ---- 参考构形单元数据（与 solveModel3DNonlinear 完全同源）----
    std::vector<CRElemData> eds;
    eds.reserve(m.elements.size());
    for (const auto& e : m.elements) {
        if (!m.nodes.count(e.n1) || !m.nodes.count(e.n2))
            throw std::runtime_error("单元 " + std::to_string(e.id) + " 引用了不存在的节点");
        eds.push_back(crMakeElemData(m, e, dm));
    }

    // ---- 约束 / 弹性支承 ----
    std::vector<int> cdofs;
    std::vector<double> cvals;
    struct Spring { int dof; double k; double val; };
    std::vector<Spring> springs;
    for (const auto& c : m.constraints) {
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        if (c.k > 0.0) { springs.push_back({d, c.k, c.value}); continue; }
        auto it = std::find(cdofs.begin(), cdofs.end(), d);
        if (it == cdofs.end()) { cdofs.push_back(d); cvals.push_back(c.value); }
        else cvals[std::distance(cdofs.begin(), it)] = c.value;
    }

    // ---- 参考外荷载向量 Fref ----
    Eigen::VectorXd Fref = Eigen::VectorXd::Zero(nd);
    for (const auto& ld : m.loads) {
        const int d = dm.safeDof(std::get<0>(ld), std::get<1>(ld) - 1);
        if (d >= 0) Fref(d) += std::get<2>(ld);
    }
    for (const auto& s : springs) Fref(s.dof) += s.k * s.val;
    for (const auto& e : m.elements) {
        if (e.type != "CABLE" || e.w <= 0.0) continue;
        const double Lc = (e.L0 > 0.0) ? e.L0
                                       : (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm();
        const double half = -0.5 * e.w * Lc;
        const int d1 = dm.safeDof(e.n1, 2), d2 = dm.safeDof(e.n2, 2);
        if (d1 >= 0) Fref(d1) += half;
        if (d2 >= 0) Fref(d2) += half;
    }
    const double Fnorm = Fref.norm();
    if (Fnorm <= 0.0) {
        out.error = "zero_reference_load";
        return out;
    }

    // 约束自由度：增量恒为 0，故 K 行/列清零置 1、右端项清零。
    auto impose = [&](Eigen::MatrixXd& K, Eigen::VectorXd& v) {
        for (int d : cdofs) {
            K.row(d).setZero();
            K.col(d).setZero();
            K(d, d) = 1.0;
            v(d) = 0.0;
        }
    };
    Eigen::VectorXd FrefD = Fref;
    for (int d : cdofs) FrefD(d) = 0.0;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(nd);
    for (std::size_t i = 0; i < cdofs.size(); ++i) u(cdofs[i]) = cvals[i];

    // ---- 组装 Fint 与一致切线 K（稠密，探针规模足够）----
    auto assemble = [&](const Eigen::VectorXd& uu,
                        Eigen::VectorXd& Fint, Eigen::MatrixXd& K) {
        Fint = Eigen::VectorXd::Zero(nd);
        K = Eigen::MatrixXd::Zero(nd, nd);
        Eigen::Matrix<double, 14, 1> ue, fe;
        Eigen::Matrix<double, 14, 14> Ke;
        for (const auto& ed : eds) {
            if (opt.strainMeasure == 1 && ed.n == 6 && ed.type == "TRUSS") {
                // 对数应变轴向杆：N = EA*ln(L/L0)，dN/dL = EA/L。
                // 切线与 fem_3d_nl.hpp 的工程应变杆同构，仅 k1 不同。
                const Vec3 u1(uu(ed.gdofs[0]), uu(ed.gdofs[1]), uu(ed.gdofs[2]));
                const Vec3 u2(uu(ed.gdofs[3]), uu(ed.gdofs[4]), uu(ed.gdofs[5]));
                const Vec3 lv = (ed.p2 + u2) - (ed.p1 + u1);
                const double Ln = std::max(lv.norm(), 1.0e-12);
                const Vec3 c = lv / Ln;
                const double EA = ed.Eeff * ed.A;
                const double N = EA * std::log(Ln / ed.L0) - EA * ed.istrnEff;
                const Eigen::Matrix3d Gm =
                    (EA / Ln) * (c * c.transpose())
                    + (N / Ln) * (Eigen::Matrix3d::Identity() - c * c.transpose());
                Eigen::Matrix<double, 6, 6> K6;
                K6.block<3, 3>(0, 0) = Gm;
                K6.block<3, 3>(0, 3) = -Gm;
                K6.block<3, 3>(3, 0) = -Gm;
                K6.block<3, 3>(3, 3) = Gm;
                const Vec3 f1 = -N * c, f2 = N * c;
                for (int k = 0; k < 3; ++k) {
                    Fint(ed.gdofs[k]) += f1[k];
                    Fint(ed.gdofs[3 + k]) += f2[k];
                }
                for (int r = 0; r < 6; ++r)
                    for (int cc = 0; cc < 6; ++cc)
                        K(ed.gdofs[r], ed.gdofs[cc]) += K6(r, cc);
                continue;
            }
            ue.setZero();
            for (int k = 0; k < ed.n; ++k) ue[k] = uu(ed.gdofs[k]);
            crElementInternal(ed, ue, fe, true, &Ke);
            for (int k = 0; k < ed.n; ++k) Fint(ed.gdofs[k]) += fe[k];
            for (int r = 0; r < ed.n; ++r)
                for (int c = 0; c < ed.n; ++c)
                    K(ed.gdofs[r], ed.gdofs[c]) += Ke(r, c);
        }
        for (const auto& s : springs) {
            Fint(s.dof) += s.k * uu(s.dof);
            K(s.dof, s.dof) += s.k;
        }
    };

    const int ctrl = dm.safeDof(opt.ctrlNode, opt.ctrlDof - 1);
    if (ctrl < 0) {
        out.error = "invalid_control_dof";
        return out;
    }

    Eigen::PartialPivLU<Eigen::MatrixXd> lu;
    Eigen::VectorXd Fint, R, Rd, duR, duF;
    Eigen::MatrixXd K, rhs, X;

    // 组装 + Dirichlet + 分解（右端项单独施加）
    auto buildSystem = [&](const Eigen::VectorXd& uu, double& piv) {
        assemble(uu, Fint, K);
        Eigen::VectorXd zero = Eigen::VectorXd::Zero(nd);
        impose(K, zero);
        lu.compute(K);
        const Eigen::VectorXd dg = lu.matrixLU().diagonal().cwiseAbs();
        const double mx = dg.maxCoeff();
        piv = (mx > 0.0) ? dg.minCoeff() / mx : 0.0;
    };

    double ds = opt.ds0;
    double lambda = 0.0;
    Eigen::VectorXd prevDu = Eigen::VectorXd::Zero(nd);
    bool havePrev = false;

    ArcLengthPoint3D p0;
    p0.step = 0;
    p0.lambda = 0.0;
    p0.ctrlDisp = u(ctrl);
    p0.ds = ds;
    p0.converged = true;
    out.points.push_back(p0);

    for (int step = 1; step <= opt.maxSteps; ++step) {
        double piv = 1.0;
        buildSystem(u, piv);

        // ---- 预测步：沿当前切线方向给 dlam ----
        duF = lu.solve(FrefD);
        const double den = std::sqrt(duF.squaredNorm() + opt.psi * opt.psi);
        double sign = 1.0;
        if (havePrev) {
            const double proj = prevDu.dot(duF);
            sign = (proj >= 0.0) ? 1.0 : -1.0;
        }
        const double dlam0 = sign * ds / den;
        Eigen::VectorXd duAcc = dlam0 * duF;
        double dlamAcc = dlam0;

        Eigen::VectorXd uT = u + duAcc;
        double lamT = lambda + dlam0;
        for (int d : cdofs) uT(d) = u(d);   // 约束自由度不参与增量

        // ---- 校正步：增广系统的 Newton 迭代 ----
        bool conv = false;
        int iter = 0;
        double resNorm = 0.0;
        for (; iter < opt.maxIter; ++iter) {
            assemble(uT, Fint, K);
            R = lamT * Fref - Fint;
            Rd = R;
            impose(K, Rd);
            lu.compute(K);
            const Eigen::VectorXd dg = lu.matrixLU().diagonal().cwiseAbs();
            const double mx = dg.maxCoeff();
            piv = (mx > 0.0) ? dg.minCoeff() / mx : 0.0;

            rhs = Eigen::MatrixXd(nd, 2);
            rhs.col(0) = Rd;
            rhs.col(1) = FrefD;
            X = lu.solve(rhs);
            duR = X.col(0);
            duF = X.col(1);

            const Eigen::VectorXd base = duAcc + duR;
            const double a = duF.squaredNorm() + opt.psi * opt.psi;
            const double b = 2.0 * base.dot(duF)
                           + 2.0 * opt.psi * opt.psi * dlamAcc;
            const double c = base.squaredNorm()
                           + opt.psi * opt.psi * dlamAcc * dlamAcc - ds * ds;
            double disc = b * b - 4.0 * a * c;
            if (disc < 0.0) disc = 0.0;    // 数值上略负时夹紧到相切
            const double sq = std::sqrt(disc);
            const double r1 = (-b + sq) / (2.0 * a);
            const double r2 = (-b - sq) / (2.0 * a);
            const Eigen::VectorXd cand1 = base + r1 * duF;
            const Eigen::VectorXd cand2 = base + r2 * duF;
            const Eigen::VectorXd& ref = havePrev ? prevDu : duAcc;
            const double s1 = cand1.dot(ref), s2 = cand2.dot(ref);
            const double dd = (s1 >= s2) ? r1 : r2;
            const Eigen::VectorXd duNew = (s1 >= s2) ? cand1 : cand2;
            const Eigen::VectorXd du = duNew - duAcc;

            uT += du;
            lamT += dd;
            duAcc = duNew;
            dlamAcc += dd;
            for (int d : cdofs) uT(d) = u(d);

            resNorm = Rd.norm();
            if (resNorm <= opt.tolForce * Fnorm) { conv = true; ++iter; break; }
        }

        if (!conv) {
            ds *= 0.5;
            if (ds < opt.dsMin) {
                out.error = "arc_length_step_collapsed";
                break;
            }
            --step;   // 重试同一平衡点编号
            continue;
        }

        u = uT;
        lambda = lamT;
        prevDu = duAcc;
        havePrev = true;

        ArcLengthPoint3D p;
        p.step = step;
        p.lambda = lambda;
        p.ctrlDisp = u(ctrl);
        p.residual = resNorm;
        p.ds = ds;
        p.pivotRatio = piv;
        p.iterations = iter;
        p.converged = true;
        out.points.push_back(p);
        out.steps = step;

        if (std::abs(u(ctrl)) >= opt.ctrlTarget ||
            std::abs(lambda) >= opt.lambdaLimit) {
            out.converged = true;
            break;
        }
        // 轻微自适应：成功步略微放大弧长，不成功步已折半
        ds = std::min(opt.dsMax, std::max(opt.ds0 * 0.25, ds * 1.05));
    }

    if (out.error.empty() && !out.converged && out.steps >= opt.maxSteps)
        out.error = "max_steps_reached";
    out.lambda = lambda;
    out.u.assign(u.data(), u.data() + nd);
    return out;
}

}  // namespace fem
