// dynamic_newmark_3d.hpp — 三维共转结构的 Newmark-β 隐式动力时程（opt-in）
//
// 第二十一次力学⑤：动力 snap-through。本文件是**独立新模块**，与
// stability_path_3d.hpp（静力弧长）、fem_3d_nl.hpp（静力 CR Newton）完全
// 解耦：不改任何既有求解器 / 选项默认值 / CLI。默认不调用即零污染。
//
// 力学口径：
//   M·a + C·v + Fint(u) = Fext(t),   v = u̇, a = ü
//   Newmark 平均加速度（β=1/4, γ=1/2 默认；β,γ 可 opt-in 覆盖）：
//     a_{n+1} = a0·(u_{n+1}−u_n) − a2·v_n − a3·a_n
//     v_{n+1} = a1·(u_{n+1}−u_n) + a4·v_n + a5·a_n
//     a0=1/(βdt²), a1=γ/(βdt), a2=1/(βdt), a3=1/(2β)−1,
//     a4=1−γ/β, a5=dt·(1−γ/(2β))
//   每步用 Newton 解非线性有效方程：
//     (a0·M + a1·C + Kt(u))·Δu = Fext(t_{n+1}) − Fint(u)
//                                  − M·a_{n+1} − C·v_{n+1}
//   质量矩阵复用 fem_3d.hpp::assembleMassMatrix3D（一致质量），
//   lumped 模式按行和集中（row-sum lumping，保持总质量）。
//   阻尼为 Rayleigh 形式 C = cMass·M + cStiff·Kt(u_n)（每步起点冻结）。
//
// 荷载时程 loadType：
//   0 = 阶跃（t>0 起 Fext=Fref，动力 snap-through 经典激励）
//   1 = 线性爬升 Fext=Fref·t/tEnd
//   2 = 半正弦 Fext=Fref·sin(π·t/tEnd)
//
// 能量：动能 T=0.5·vᵀMv；应变能 U 用梯形增量积分（与单元类型无关）；
// 外功 W 同样梯形积分。无阻尼时 E=T+U−W 近似守恒（Newmark 对线性系统
// 严格能量守恒，非线性下为二阶精度，供趋势校验而非严格守恒锚点）。
#pragma once

#include "fem_3d.hpp"
#include "fem_3d_nl.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {

struct NewmarkOptions3D {
    double dt = 0.0;            // 时间步长（>0 必需）
    double tEnd = 0.0;          // 终止时刻（>0 必需）
    double beta = 0.25;         // Newmark β
    double gamma = 0.5;         // Newmark γ
    int massMode = 0;           // 0=一致质量（默认） 1=行和集中质量
    double cMass = 0.0;         // 质量比例阻尼
    double cStiff = 0.0;        // 刚度比例阻尼（每步起点冻结 Kt）
    int loadType = 0;           // 0=阶跃 1=爬升 2=半正弦
    double tolForce = 1.0e-7;   // 残差相对 |Fref| 容差
    int maxIter = 30;           // 每步最大 Newton 迭代
    int ctrlNode = 2;           // 控制节点（诊断 / snap 判据）
    int ctrlDof = 2;            // 控制自由度 1~6
    double snapThreshold = 0.0; // >0 时首次越过该位移即记录 snap 时刻
    int strainMeasure = 0;      // 轴向单元应变度量（与弧长口径一致）
};

struct NewmarkPoint3D {
    int step = 0;
    double t = 0.0;
    double lambda = 0.0;        // Fext(t)/Fref
    double ctrlDisp = 0.0;      // 带符号控制位移
    double ctrlVel = 0.0;
    double ctrlAcc = 0.0;
    double residual = 0.0;
    int iterations = 0;
    double kinetic = 0.0;
    double strain = 0.0;
    double externalWork = 0.0;
    bool converged = false;
};

struct NewmarkResult3D {
    bool converged = false;
    int steps = 0;
    double snapTime = -1.0;     // 首次越过 snapThreshold 的时刻（<0 = 未越过）
    double snapDisp = 0.0;
    double snapLambda = 0.0;
    double massTotal = 0.0;     // 一致质量矩阵迹 / 3（诊断）
    std::vector<double> u;      // 终点位移
    std::vector<double> v;      // 终点速度
    std::vector<double> a;      // 终点加速度
    std::vector<NewmarkPoint3D> points;
    std::string error;          // 空 = 正常
};

// 一致质量 → 稠密；massMode=1 时按行和集中（保持每方向总质量）。
inline Eigen::MatrixXd denseMassMatrix3D(const Model3D& m, const DofMap3D& dm,
                                         int massMode) {
    const SparseMat3D Ms = assembleMassMatrix3D(m, dm);
    Eigen::MatrixXd Md = Eigen::MatrixXd::Zero(dm.nDofs, dm.nDofs);
    for (int k = 0; k < Ms.outerSize(); ++k)
        for (SparseMat3D::InnerIterator it(Ms, k); it; ++it)
            Md(it.row(), it.col()) += it.value();
    if (massMode == 0) return Md;
    Eigen::MatrixXd Ml = Eigen::MatrixXd::Zero(dm.nDofs, dm.nDofs);
    for (int i = 0; i < dm.nDofs; ++i) {
        const double s = Md.row(i).sum();
        Ml(i, i) = (std::isfinite(s) && s > 0.0) ? s : Md(i, i);
    }
    return Ml;
}

inline NewmarkResult3D solveNewmarkDynamic3D(const Model3D& m,
                                             const NewmarkOptions3D& opt) {
    NewmarkResult3D out;
    if (!(opt.dt > 0.0)) { out.error = "invalid_dt"; return out; }
    if (!(opt.tEnd > 0.0)) { out.error = "invalid_tend"; return out; }
    if (!(opt.beta > 0.0)) { out.error = "invalid_beta"; return out; }
    if (!(opt.gamma > 0.0)) { out.error = "invalid_gamma"; return out; }

    const DofMap3D dm = buildDofMap3D(m);
    const int nd = dm.nDofs;
    if (nd <= 0) { out.error = "empty_dof_map"; return out; }

    // ---- 单元数据（与静力 CR 求解同源）----
    std::vector<CRElemData> eds;
    eds.reserve(m.elements.size());
    for (const auto& e : m.elements) {
        if (!m.nodes.count(e.n1) || !m.nodes.count(e.n2))
            throw std::runtime_error("动力单元 " + std::to_string(e.id) + " 引用不存在节点");
        eds.push_back(crMakeElemData(m, e, dm));
    }

    // ---- 约束 / 弹性支承 ----
    std::vector<int> cdofs;
    std::vector<double> cvals;
    struct DS { int dof; double k, val; };
    std::vector<DS> springs;
    for (const auto& c : m.constraints) {
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        if (c.k > 0.0) { springs.push_back({d, c.k, c.value}); continue; }
        auto it = std::find(cdofs.begin(), cdofs.end(), d);
        if (it == cdofs.end()) { cdofs.push_back(d); cvals.push_back(c.value); }
        else cvals[std::distance(cdofs.begin(), it)] = c.value;
    }

    // ---- 参考荷载 Fref（与静力弧长同口径）----
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
    if (!(Fnorm > 0.0)) { out.error = "zero_reference_load"; return out; }

    // ---- 质量矩阵 ----
    const Eigen::MatrixXd M = denseMassMatrix3D(m, dm, opt.massMode);
    double trace = 0.0;
    for (int i = 0; i < nd; ++i) trace += M(i, i);
    out.massTotal = trace / 3.0;
    if (!(trace > 0.0)) { out.error = "zero_mass_matrix"; return out; }

    // ---- Fint / Kt 装配（复用 CR 单元的一致切线）----
    auto assemble = [&](const Eigen::VectorXd& uu, Eigen::VectorXd& Fint,
                        Eigen::MatrixXd& Kt) {
        Fint = Eigen::VectorXd::Zero(nd);
        Kt = Eigen::MatrixXd::Zero(nd, nd);
        Eigen::Matrix<double, 14, 1> ue, fe;
        Eigen::Matrix<double, 14, 14> Ke;
        for (const auto& ed : eds) {
            if (opt.strainMeasure == 1 && ed.n == 6 && ed.type == "TRUSS") {
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
                K6.block<3, 3>(0, 0) = Gm;  K6.block<3, 3>(0, 3) = -Gm;
                K6.block<3, 3>(3, 0) = -Gm; K6.block<3, 3>(3, 3) = Gm;
                const Vec3 f1 = -N * c, f2 = N * c;
                for (int k = 0; k < 3; ++k) {
                    Fint(ed.gdofs[k]) += f1[k];
                    Fint(ed.gdofs[3 + k]) += f2[k];
                }
                for (int r = 0; r < 6; ++r)
                    for (int cc = 0; cc < 6; ++cc)
                        Kt(ed.gdofs[r], ed.gdofs[cc]) += K6(r, cc);
                continue;
            }
            ue.setZero();
            for (int k = 0; k < ed.n; ++k) ue[k] = uu(ed.gdofs[k]);
            crElementInternal(ed, ue, fe, true, &Ke);
            for (int k = 0; k < ed.n; ++k) Fint(ed.gdofs[k]) += fe[k];
            for (int r = 0; r < ed.n; ++r)
                for (int c = 0; c < ed.n; ++c)
                    Kt(ed.gdofs[r], ed.gdofs[c]) += Ke(r, c);
        }
        for (const auto& s : springs) {
            Fint(s.dof) += s.k * uu(s.dof);
            Kt(s.dof, s.dof) += s.k;
        }
    };

    auto impose = [&](Eigen::MatrixXd& A, Eigen::VectorXd& b) {
        for (int d : cdofs) {
            A.row(d).setZero(); A.col(d).setZero(); A(d, d) = 1.0;
            b(d) = 0.0;
        }
    };

    const int ctrl = dm.safeDof(opt.ctrlNode, opt.ctrlDof - 1);
    if (ctrl < 0) { out.error = "invalid_control_dof"; return out; }

    // ---- 初态：u=0, v=0, a0 = M⁻¹(Fext(0+)−Fint(0)) ----
    Eigen::VectorXd u = Eigen::VectorXd::Zero(nd);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(nd);
    Eigen::VectorXd a = Eigen::VectorXd::Zero(nd);
    for (std::size_t i = 0; i < cdofs.size(); ++i) u(cdofs[i]) = cvals[i];
    Eigen::VectorXd Fint, R;
    Eigen::MatrixXd Kt;
    assemble(u, Fint, Kt);
    auto lambdaAt = [&](double t) -> double {
        const double T = opt.tEnd;
        if (opt.loadType == 1) return std::min(1.0, std::max(0.0, t / T));
        if (opt.loadType == 2) {
            constexpr double pi = 3.14159265358979323846;
            return (t >= 0.0 && t <= T) ? std::sin(pi * t / T) : 0.0;
        }
        return (t > 0.0) ? 1.0 : 0.0;   // 阶跃
    };
    {
        const double lam0 = lambdaAt(opt.dt);   // 阶跃在第一步已施加
        R = lam0 * Fref - Fint;
        for (int d : cdofs) R(d) = 0.0;
        Eigen::MatrixXd A = M;
        Eigen::VectorXd b = R;
        impose(A, b);
        a = A.ldlt().solve(b);
        for (int d : cdofs) a(d) = 0.0;
    }

    const double a0 = 1.0 / (opt.beta * opt.dt * opt.dt);
    const double a1 = opt.gamma / (opt.beta * opt.dt);
    const double a2 = 1.0 / (opt.beta * opt.dt);
    const double a3 = 1.0 / (2.0 * opt.beta) - 1.0;
    const double a4 = 1.0 - opt.gamma / opt.beta;
    const double a5 = opt.dt * (1.0 - opt.gamma / (2.0 * opt.beta));

    double strainE = 0.0, workE = 0.0;
    Eigen::VectorXd FintOld = Fint;
    double lamOld = 0.0;

    NewmarkPoint3D p0;
    p0.step = 0; p0.t = 0.0; p0.lambda = 0.0;
    p0.ctrlDisp = u(ctrl); p0.ctrlVel = 0.0; p0.ctrlAcc = a(ctrl);
    p0.kinetic = 0.5 * v.dot(M * v);
    p0.strain = 0.0; p0.externalWork = 0.0;
    p0.converged = true;
    out.points.push_back(p0);

    const int nSteps = static_cast<int>(std::ceil(opt.tEnd / opt.dt));
    Eigen::MatrixXd Keff(nd, nd), Cmat(nd, nd), A(nd, nd);
    Eigen::VectorXd b(nd), aNew(nd), vNew(nd), Fext(nd);
    bool allConv = true;

    for (int step = 1; step <= nSteps; ++step) {
        const double t = step * opt.dt;
        const double lam = lambdaAt(t);
        Fext = lam * Fref;

        // Rayleigh 阻尼：Kt 在步起点冻结（保证 Keff 在迭代中对称定常）
        Cmat = opt.cMass * M + opt.cStiff * Kt;
        A.noalias() = a0 * M + a1 * Cmat;

        // Newmark 经典预测步：u_pred = u + dt·v + dt²(0.5−β)·a
        // （此时 a_pred=0、v_pred=v+dt(1−γ)a，与 Newmark 递推自洽）
        Eigen::VectorXd uTrial = u + opt.dt * v
                               + (opt.dt * opt.dt * (0.5 - opt.beta)) * a;
        for (int d : cdofs) uTrial(d) = u(d);
        Eigen::VectorXd FiTrial = FintOld;
        Eigen::MatrixXd KtTrial = Kt;
        bool conv = false;
        int iter = 0;
        double resNorm = 0.0;
        while (iter < opt.maxIter) {
            aNew = a0 * (uTrial - u) - a2 * v - a3 * a;
            vNew = a1 * (uTrial - u) + a4 * v + a5 * a;
            assemble(uTrial, FiTrial, KtTrial);
            R = Fext - FiTrial - M * aNew - Cmat * vNew;
            Keff = A + KtTrial;
            b = R;
            impose(Keff, b);   // 约束 DOF 清零：残差只统计自由自由度（排除支反力）
            resNorm = b.norm();   // 当前迭代点自由残差（标准 Newton 收敛判据）
            const Eigen::VectorXd du = Keff.ldlt().solve(b);
            uTrial += du;
            for (int d : cdofs) uTrial(d) = u(d);
            ++iter;
            if (resNorm <= opt.tolForce * Fnorm) { conv = true; break; }
        }
        if (!conv) {
            out.error = "newton_not_converged";
            allConv = false;
            break;
        }

        aNew = a0 * (uTrial - u) - a2 * v - a3 * a;
        vNew = a1 * (uTrial - u) + a4 * v + a5 * a;
        Eigen::VectorXd FiNew;
        Eigen::MatrixXd KtNew;
        assemble(uTrial, FiNew, KtNew);
        // 能量：梯形增量
        strainE += 0.5 * (FiNew + FintOld).dot(uTrial - u);
        workE += 0.5 * (lam + lamOld) * Fref.dot(uTrial - u);

        u = uTrial; v = vNew; a = aNew;
        Fint = FiNew; Kt = KtNew; FintOld = FiNew; lamOld = lam;

        NewmarkPoint3D p;
        p.step = step; p.t = t; p.lambda = lam;
        p.ctrlDisp = u(ctrl); p.ctrlVel = v(ctrl); p.ctrlAcc = a(ctrl);
        p.residual = resNorm; p.iterations = iter;
        p.kinetic = 0.5 * v.dot(M * v);
        p.strain = strainE; p.externalWork = workE;
        p.converged = true;
        out.points.push_back(p);

        if (opt.snapThreshold > 0.0 && out.snapTime < 0.0 &&
            std::abs(u(ctrl)) >= opt.snapThreshold) {
            out.snapTime = t;
            out.snapDisp = u(ctrl);
            out.snapLambda = lam;
        }
    }

    out.steps = static_cast<int>(out.points.size()) - 1;
    out.converged = allConv && out.steps == nSteps;
    out.u.assign(u.data(), u.data() + nd);
    out.v.assign(v.data(), v.data() + nd);
    out.a.assign(a.data(), a.data() + nd);
    return out;
}

}  // namespace fem
