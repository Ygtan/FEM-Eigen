// stability_path_3d.hpp — 三维共转结构的球面弧长法平衡路径跟踪
//
// 独立 opt-in 能力，不改 fem_3d_nl.hpp::solveModel3DNonlinear（仍是荷载
// 控制 Newton）、不改 NonlinearOptions3D 默认值。第十七次起以
// `STABILITY=1` 接入 cable_solver.exe 的**参考构形**独立模式；第十九次 ①
// 起另以 `STABASIS=COMPLETION` 接入**成桥态起点**模式（见下方
// ArcLengthState3D）。两种基准都是 opt-in，默认全关、默认行为逐位不变。
//
// 力学口径：完全复用 fem_3d_nl.hpp 的 CR 单元内力与解析一致切线
// （crMakeElemData / crElementInternal），支持 TRUSS / CABLE / BEAM /
// BEAM7。唯一新增的是把荷载因子 lambda 当作未知量的增广系统：
//
//     参考构形基准：R(u,lam) = lam * Fref - Fint(u) = 0,   u(0) = 0
//     成桥态基准：  R(u,lam) = F0 + lam * Finc - Fint(u) = 0, u(0) = u0
//     du^T du + psi^2 * dlam^2 = ds^2         (Crisfield 球面约束，恒对 du)
//
// 成桥态基准里 F0 是**不随 λ 缩放**的恒定荷载（成桥态恒载 + 弹性支承
// 预置 + 索自重），Finc 是 `*STABILITY_LOAD` 给出的增量荷载；两者必须
// 分列，不能合并成「一个 λ 放大全部荷载」。λ 也因此只放大增量荷载，
// 与参考构形基准的 λ（总参考荷载放大倍数）**不是同一个量**。
//
// 路径约束方法（ArcLengthOptions3D::pathMethod，opt-in）：
//   0 = Crisfield 球面弧长（默认，逐位不变）：
//       du^T du + psi^2 dlam^2 = ds^2，代入后 dλ 是二次方程，
//       按「与上一收敛步增量夹角余弦最大」选根。
//   1 = Wempner–Riks 正交轨迹（本项新增，仅 STAGE 起点入口放开）：
//       每轮校正增量 duR + dλ·duF 与上一收敛步增量正交，
//       即 duR·prevDu + dλ·(duF·prevDu + psi²·prevDlam) = 0，
//       代入后 dλ 是**线性方程**，无根可选、无夹角启发式。
//       两者跟踪同一平衡路径 R(u,λ)=0，只是参数化不同，可作
//       「同问题、双方法」交叉验证。
//   第一步（无 prev）一律走球面约束启动，保证球面默认路径逐位不变。
//   注意：本项目的 RIKS 是正交型，与 ANSYS ARCLEN / ABAQUS *STATIC,RIKS
//   的球面型不同源，不可当作同一算法对拍。
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

// ==================== 第二十一次力学②：线性屈曲模态缺陷场 ====================
//
// 缺陷源用线性屈曲模态 (K + λ·Kg)·φ = 0 的最低模态，而不是随便给的几何扰动：
// 模态本身来自参考荷载下的线性轴力几何刚度 Kg，因此"哪个方向最容易失稳"
// 由结构自己决定，符合缺陷敏感性分析的通行做法（Koiter 初始缺陷假设）。
//
// 本函数**不改动** fem_3d.hpp::computeBucklingResult3D（它只返回 factors，
// 被 example_3d_grillage 等多处调用，逐位不能漂移）。这里用同一套装配函数
// 重新求一次稠密特征问题，只为拿模态向量。规模是探针级（桁架/浅拱），
// 稠密求解足够；Spectra 路径不在本 opt-in 能力范围内。
struct BucklingModeShape3D {
    std::vector<double> phi;   // 全局位移向量（buildDofMap3D 序）
    double lambda = 0.0;       // 线性屈曲荷载因子 λ = −1/μ
    std::string error;         // 空 = 正常
};

inline BucklingModeShape3D computeBucklingModeShape3D(const Model3D& m,
                                                       int modeIndex = 0) {
    BucklingModeShape3D out;
    if (modeIndex < 0) {
        out.error = "invalid_mode_index";
        return out;
    }
    DofMap3D dm;
    const SparseMat3D K = assembleStiffness3D(m, dm);
    const int nd = dm.nDofs;
    const std::vector<double> U = solveModel3D(m);   // 参考荷载下线性位移 → 轴力
    const SparseMat3D Kg = assembleGeometricStiffness3D(m, dm, U);
    const std::vector<int> free = freeDofs3D(m, dm, K);
    const int nf = static_cast<int>(free.size());
    if (nf < 2) {
        out.error = "insufficient_free_dofs";
        return out;
    }
    if (Kg.nonZeros() == 0) {
        out.error = "no_geometric_stiffness";
        return out;
    }

    Eigen::MatrixXd Kd(nf, nf), Kgd(nf, nf);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            Kd(a, b) = K.coeff(free[a], free[b]);
            Kgd(a, b) = Kg.coeff(free[a], free[b]);
        }
    Eigen::LLT<Eigen::MatrixXd> llt(Kd);
    if (llt.info() != Eigen::Success) {
        out.error = "cholesky_failed";
        return out;
    }
    const Eigen::MatrixXd L = llt.matrixL();
    const Eigen::MatrixXd Linv = L.inverse();
    // A = L^{-1}·Kg·L^{-T}；A·v = μ·v，φ = L^{-T}·v 满足 Kg·φ = μ·K·φ，
    // 即 (K + λ·Kg)·φ = 0，λ = −1/μ（μ<0 才是受压屈曲模态）。
    const Eigen::MatrixXd A = Linv * Kgd * Linv.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(A);
    if (ses.info() != Eigen::Success) {
        out.error = "eigensolver_failed";
        return out;
    }
    std::vector<int> buck;   // ses 特征值升序：buck[0] 最负 μ → 最低 λ_cr
    for (int i = 0; i < ses.eigenvalues().size(); ++i) {
        const double mu = ses.eigenvalues()(i);
        if (mu < -1e-14 && std::isfinite(mu)) buck.push_back(i);
    }
    if (modeIndex >= static_cast<int>(buck.size())) {
        out.error = "mode_index_out_of_range";
        return out;
    }
    const int sel = buck[modeIndex];
    out.lambda = -1.0 / ses.eigenvalues()(sel);

    Eigen::VectorXd phiFree = Linv.transpose() * ses.eigenvectors().col(sel);
    // 符号与归一化固定：|φ| 最大分量恒为 −1（特征向量符号由求解器决定，
    // 不固定会导致同一模态正负号漂移，缺陷方向不可复现）。对 von Mises
    // 对称模态即"顶点向下"，amp>0 使顶点下移，物理方向明确。
    int iMax = 0;
    for (int i = 1; i < nf; ++i)
        if (std::abs(phiFree(i)) > std::abs(phiFree(iMax))) iMax = i;
    const double mx = std::abs(phiFree(iMax));
    if (!(mx > 0.0) || !std::isfinite(mx)) {
        out.error = "degenerate_mode";
        return out;
    }
    if (phiFree(iMax) > 0.0) phiFree = -phiFree;
    phiFree /= mx;

    out.phi.assign(nd, 0.0);
    for (int i = 0; i < nf; ++i) out.phi[free[i]] = phiFree(i);
    return out;
}

// 把屈曲模态按幅值施加到**当前参考构形**的节点坐标，形成初始缺陷场。
// 只改平动 DOF；屈曲模态定义在自由自由度上，约束节点分量本来就是 0，
// 因此不会破坏边界条件。L0（无应力长度）保持不变——这是"初始几何偏差"
// 而非"重新设计结构"，缺陷构形下会有初始内力，正是敏感性来源。
// updateL0=true（本 opt-in 采用）：缺陷是**无应力**初始几何偏差——把
// TRUSS/CABLE 的无应力长度同步改成缺陷后弦长，缺陷构形本身在 λ=0 就是
// 平衡态。这对应"制造误差 / 初始弯曲"型缺陷（构件成形后无残余应力），
// 也是缺陷敏感性分析的标准简化：稳定路径从平衡态起步，不会因 L0 不变
// 引入大初始内力把切线刚度推到奇异（实测 L0 不变时 amp≥0.005 即失稳）。
// updateL0=false：保留残余应力型缺陷（L0 不变），路径起点带初始内力，
// 仅供将来扩展，不在本项验收范围内。
inline bool applyImperfectionToNodes(Model3D& m, const DofMap3D& dm,
                                     const std::vector<double>& phi, double amp,
                                     bool updateL0) {
    if (static_cast<int>(phi.size()) != dm.nDofs) return false;
    for (auto& kv : m.nodes) {
        const int nid = kv.first;
        Vec3& p = kv.second;
        for (int k = 0; k < 3; ++k) {
            const int d = dm.safeDof(nid, k);
            if (d >= 0 && d < static_cast<int>(phi.size())) p(k) += amp * phi[d];
        }
    }
    if (updateL0) {
        for (auto& e : m.elements) {
            if ((e.type != "TRUSS" && e.type != "CABLE") ||
                !m.nodes.count(e.n1) || !m.nodes.count(e.n2)) continue;
            e.L0 = (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm();
        }
    }
    return true;
}

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
    // 路径约束：0 = Crisfield 球面弧长（默认，所有既有输出逐位不变），
    // 1 = Wempner–Riks 正交轨迹（opt-in，入口层仅 STAGE 起点放开）。
    // 两者复用同一残差装配 / 一致切线 / LU 分解，只在 dλ 的约束方程不同。
    int pathMethod = 0;
    // 第二十一次力学①：位移控制 / 线搜索（opt-in，默认 0 逐位不变）。
    // controlMode=0（默认）：沿用上面的弧长 / Riks 约束，所有既有输出不变。
    // controlMode=1（位移控制）：每步**控制自由度的位移增量固定为 ds**，
    //   λ 为未知量，由线性约束 du_ctrl = ds 唯一确定（与 Riks 同构，
    //   只是约束从「与上步增量正交」换成「控制 DOF 位移增量给定」）：
    //     dλ = -duR(ctrl) / duF(ctrl)
    //   这是越过极值点 / snap-through 的经典位移控制，天然不依赖 λ 单调。
    // lineSearch=1：校正增量里对残差修正 duR 做 Armijo 回溯（缩放 alpha 后
    //   重算 dλ 保持控制位移约束），抑制跳跃失稳时的残差过冲；正常步
    //   alpha=1 立即接受，结果与不开线搜索一致。仅 controlMode=1 时生效。
    int controlMode = 0;
    double ctrlDispStep = 0.0;   // 位移控制每步控制 DOF 的目标位移增量（>0 生效，否则回退 ds）
    int lineSearch = 0;
    double lsArmijo = 1.0e-4;    // Armijo 充分下降系数 c1
    double lsAlphaMin = 1.0e-3;  // 回溯步长下限
    int lsMaxCut = 12;           // 最多回溯次数
    // 第二十一次力学④：分岔点检测与分支切换（opt-in，默认 0/0.05/1/0 逐位不变）。
    // bifurcate=0：整段检测/切换逻辑跳过，所有既有输出逐字节不变。
    // bifurcate=1：每个收敛点先看 pivotRatio（已逐点输出但此前未使用），
    //   低于 bifPivot 时再用当前一致切线自由子块的最小特征值确认临界点；
    //   确认后沿 computeBucklingModeShape3D(m, bifMode) 的模态施加 bifAmp
    //   扰动到收敛位移 u，后续弧长步从扰动构形继续 → 跟踪分岔后分支。
    // 扰动只改位移，不改 L0 / 荷载 / 求解器算法；多阶分岔只切换一次。
    int bifurcate = 0;
    double bifPivot = 0.05;      // pivotRatio 触发阈值
    int bifMode = 1;             // 扰动模态序号（0=最低阶，1=次低阶/反对称）
    double bifAmp = 0.0;         // 扰动幅值（长度；0 = 检测但不切换，负例守卫用）
    double bifEigTol = 1.0e-3;   // |μ_min|/max|diag(Kff)| 确认阈值
};

// ---- 第十九次 ①：成桥态起点的初态 ------------------------------------
// 三个向量都按 buildDofMap3D(m) 的自由度序排列；**全空** = 参考构形基准，
// 与第十七 / 十八次逐位一致（这是第 19/20/21 项 CTest 不漂移的保证）。
// 只要给了一个，三个都必须给齐且长度等于 nd，否则返回
// invalid_initial_state，绝不静默截断或补 0。
struct ArcLengthState3D {
    std::vector<double> u0;    // 成桥态收敛位移（路径起点）
    std::vector<double> F0;    // 恒定荷载，不随 λ 缩放
    std::vector<double> Finc;  // 增量参考荷载，λ 只放大这一份
};

struct ArcLengthPoint3D {
    int step = 0;
    double lambda = 0.0;
    double ctrlDisp = 0.0;      // u(ctrlDof)，带符号
    double residual = 0.0;      // 收敛时的自由自由度残差范数
    double ds = 0.0;            // 本步实际弧长半径
    double pivotRatio = 1.0;    // 切线奇异指标：min|diag(U)| / max|diag(U)|
    int iterations = 0;
    int lineCuts = 0;           // 第二十一次力学①：本步线搜索回溯次数（默认恒 0）
    double alpha = 1.0;         // 第二十一次力学①：本步接受的步长（默认恒 1.0）
    // 第二十一次力学④：分岔诊断，默认全 0（既有 benchmark 不读 → 不漂移）。
    int bifDetected = 0;        // 本步满足分岔候选（pivot 低 + 一致切线特征值过零）
    int bifSwitched = 0;        // 本步实际施加了模态扰动
    double bifEigRatio = 0.0;   // |μ_min| / max|diag(Kff)|（确认前的临界度）
    double bifModalAmp = 0.0;   // 当前位移在扰动模态上的投影幅值（分支进入深度）
    bool converged = false;
};

struct ArcLengthResult3D {
    bool converged = false;
    int steps = 0;
    double lambda = 0.0;
    std::vector<double> u;                    // 终点位移向量
    std::vector<ArcLengthPoint3D> points;    // 逐平衡点路径
    std::string error;                       // 失败原因（空 = 正常）
    double res0 = 0.0;   // 第十九次 ①：起点闭合残差 |F0 − Fint(u0)|（自由 dof）
    // 第二十一次力学④：分岔切换汇总，默认全 0（既有 benchmark 不读 → 不漂移）。
    int bifDetected = 0;        // 是否检测到分岔候选
    int bifSwitched = 0;        // 是否实际切换分支
    int bifStep = 0;            // 切换发生的平衡点步号
    double bifLambda = 0.0;     // 切换时刻的荷载因子
    double bifCtrlDisp = 0.0;   // 切换时刻的控制位移
    double bifMinEigRatio = 0.0;// 切换时刻的临界度 |μ_min|/max|diag(Kff)|
    double bifPivot = 0.0;      // 切换时刻的 pivotRatio（轻量入口记录）
};

inline ArcLengthResult3D solveArcLengthPath3D(const Model3D& m,
                                              const ArcLengthOptions3D& opt = {},
                                              const ArcLengthState3D& st = {}) {
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

    // ---- 参考外荷载向量 Fref -------------------------------------------
    // 默认（st 全空）：与 fem_3d_nl.hpp::solveModel3DNonlinear 的 Fext
    // 逐项同源 —— 节点荷载 + 弹性支承预置 + CABLE 自重两端各半。
    // 成桥态基准：F0 走 st.F0、Finc 走 st.Finc，二者由调用方按同一口径
    // 装配；F0 不随 λ 缩放，λ 只放大 Finc。
    const bool hasState = !st.u0.empty() || !st.F0.empty() || !st.Finc.empty();
    if (hasState &&
        (static_cast<int>(st.u0.size()) != nd ||
         static_cast<int>(st.F0.size()) != nd ||
         static_cast<int>(st.Finc.size()) != nd)) {
        out.error = "invalid_initial_state";
        return out;
    }
    Eigen::VectorXd Fref = Eigen::VectorXd::Zero(nd);
    if (!st.Finc.empty()) {
        for (int i = 0; i < nd; ++i) Fref(i) = st.Finc[static_cast<std::size_t>(i)];
    } else {
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
    }
    Eigen::VectorXd F0 = Eigen::VectorXd::Zero(nd);
    if (!st.F0.empty())
        for (int i = 0; i < nd; ++i) F0(i) = st.F0[static_cast<std::size_t>(i)];
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
    if (!st.u0.empty())
        for (int i = 0; i < nd; ++i) u(i) = st.u0[static_cast<std::size_t>(i)];
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

    // 段内 CTRL 一律是**相对起点的控制位移增量**：参考构形基准起点为 0，
    // 与既有输出逐位一致；成桥态基准起点为 u0(ctrl)，避免把成桥态已有
    // 位移误当增量，也避免终止判据在第一步就误触发。
    const double ctrlDisp0 = u(ctrl);
    // 起点闭合残差 |F0 − Fint(u0)|：成桥态 u0 应是 F0 的平衡解，这里如实
    // 记录差值而不是补 0，供 STATE RES0= 与下游自洽对拍。
    if (hasState) {
        Eigen::VectorXd Fint0, R0;
        Eigen::MatrixXd K0tmp;
        assemble(u, Fint0, K0tmp);
        R0 = F0 - Fint0;
        for (int d : cdofs) R0(d) = 0.0;
        out.res0 = R0.norm();
    }

    // 第二十一次力学④：当前一致切线自由子块最小特征值（分岔确认用）。
    // 只在 bifurcate=1 且 pivotRatio 触阈时调用；默认零开销、零污染。
    auto minTangentEig = [&](const Eigen::VectorXd& uu, double& ratio) -> bool {
        Eigen::VectorXd Fi;
        Eigen::MatrixXd Kt;
        assemble(uu, Fi, Kt);   // 未施加 Dirichlet 的原始一致切线
        std::vector<int> freeList;
        freeList.reserve(nd);
        for (int i = 0; i < nd; ++i) {
            bool fixed = false;
            for (int d : cdofs) if (d == i) { fixed = true; break; }
            if (!fixed) freeList.push_back(i);
        }
        const int nf = static_cast<int>(freeList.size());
        if (nf < 1) return false;
        Eigen::MatrixXd Kff(nf, nf);
        double maxDiag = 0.0;
        for (int a = 0; a < nf; ++a) {
            maxDiag = std::max(maxDiag, std::abs(Kt(freeList[a], freeList[a])));
            for (int b = 0; b < nf; ++b)
                Kff(a, b) = Kt(freeList[a], freeList[b]);
        }
        if (!(maxDiag > 0.0)) return false;
        Kff = 0.5 * (Kff + Kff.transpose());   // 消去装配舍入的非对称
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Kff);
        if (ses.info() != Eigen::Success) return false;
        ratio = std::abs(ses.eigenvalues()(0)) / maxDiag;
        return std::isfinite(ratio);
    };
    bool bifTriggered = false;   // 多阶分岔只检测/切换一次
    // 扰动模态在 opt-in 时只求一次，供切换与逐点投影共用（φ 已按最大分量归一化）。
    Eigen::VectorXd bifPhi;
    bool bifPhiOk = false;
    if (opt.bifurcate != 0) {
        const fem::BucklingModeShape3D bk0 =
            computeBucklingModeShape3D(m, opt.bifMode);
        if (bk0.error.empty() && static_cast<int>(bk0.phi.size()) == nd) {
            bifPhi.resize(nd);
            for (int i = 0; i < nd; ++i) bifPhi(i) = bk0.phi[i];
            bifPhiOk = true;
        }
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
    double prevDlam = 0.0;   // RIKS 正交约束需要上一步的 λ 增量
    bool havePrev = false;
    // 位移控制方向在第一步按荷载响应符号锁定，越过极限点后保持不变，
    // 否则控制位移会在极值点两侧来回翻转（切线变号）。
    double ctrlSignLocked = 0.0;

    ArcLengthPoint3D p0;
    p0.step = 0;
    p0.lambda = 0.0;
    p0.ctrlDisp = u(ctrl) - ctrlDisp0;
    p0.residual = out.res0;   // 成桥态基准 = 起点闭合残差；参考构形 = 0
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
        double dlam0 = sign * ds / den;
        // 第二十一次力学①：位移控制预测步让控制 DOF 沿荷载响应方向精确推进
        // ctrlDispStep（有符号目标）。CTRLSTEP 给的是幅值，符号取 duF(ctrl) 的
        // 方向，保证向下荷载对应向下位移；未给时幅值回退 ds。
        double ctrlTargetStep = 0.0;   // 本步控制位移的**有符号**目标增量
        if (opt.controlMode == 1) {
            const double duFctrl = duF(ctrl);
            if (std::abs(duFctrl) <= 1.0e-14) {
                out.error = "singular_control_dof";
                break;
            }
            const double mag = (opt.ctrlDispStep > 0.0) ? opt.ctrlDispStep : ds;
            if (ctrlSignLocked == 0.0)
                ctrlSignLocked = (duFctrl >= 0.0 ? 1.0 : -1.0);
            ctrlTargetStep = ctrlSignLocked * mag;
            dlam0 = ctrlTargetStep / duFctrl;
        }
        Eigen::VectorXd duAcc = dlam0 * duF;
        double dlamAcc = dlam0;

        Eigen::VectorXd uT = u + duAcc;
        double lamT = lambda + dlam0;
        for (int d : cdofs) uT(d) = u(d);   // 约束自由度不参与增量

        // ---- 校正步：增广系统的 Newton 迭代 ----
        bool conv = false;
        int iter = 0;
        double resNorm = 0.0;
        int lineCutsHere = 0;
        double alphaHere = 1.0;
        for (; iter < opt.maxIter; ++iter) {
            assemble(uT, Fint, K);
            R = F0 + lamT * Fref - Fint;
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
            double dd = 0.0;
            Eigen::VectorXd duNew = base;
            if (opt.controlMode == 1) {
                // 第二十一次力学①：位移控制。每步把控制 DOF 相对步起点 u 的
                // 位移增量锁定为 ctrlDispStep，λ 由线性约束唯一确定：
                //   base(ctrl) + dλ·duF(ctrl) = ctrlDispStep
                //   dλ = (ctrlDispStep - base(ctrl)) / duF(ctrl)
                // 预测步已令 duAcc(ctrl)=ctrlDispStep，故迭代中 duAcc(ctrl)
                // 恒等于该值；整步增量 duFull = duR + ddFull·duF 满足约束。
                // 线搜索只缩放「整步校正」alpha·duFull / alpha·ddFull，约束
                // 仍成立（因 duAcc(ctrl)=ctrlDispStep），故不改变平衡路径。
                const double duFctrl = duF(ctrl);
                if (std::abs(duFctrl) <= 1.0e-14) {
                    out.error = "singular_control_dof";
                    break;
                }
                const double ctrlTarget = ctrlTargetStep;
                const double ddFull = (ctrlTarget - base(ctrl)) / duFctrl;
                const Eigen::VectorXd duFull = duR + ddFull * duF;
                if (opt.lineSearch) {
                    // Armijo：phi(a)=||R(uT+a·duFull, lamT+a·ddFull)||，
                    // 接受 phi(a) <= (1 - c1·a)·phi(0)。整步通常直接满足，
                    // alpha=1、0 次回溯，与不开线搜索逐位一致。
                    auto phi = [&](double a) {
                        Eigen::VectorXd uu = uT + a * duFull;
                        for (int d : cdofs) uu(d) = u(d);
                        Eigen::VectorXd Fi, Rr;
                        Eigen::MatrixXd Kt;
                        assemble(uu, Fi, Kt);
                        Rr = F0 + (lamT + a * ddFull) * Fref - Fi;
                        for (int d : cdofs) Rr(d) = 0.0;
                        return Rr.norm();
                    };
                    const double phi0 = Rd.norm();
                    double alpha = 1.0;
                    lineCutsHere = 0;
                    while (lineCutsHere < opt.lsMaxCut &&
                           phi(alpha) > (1.0 - opt.lsArmijo * alpha) * phi0 &&
                           alpha > opt.lsAlphaMin) {
                        alpha *= 0.5;
                        ++lineCutsHere;
                    }
                    alphaHere = alpha;
                    dd = alpha * ddFull;
                    duNew = duAcc + alpha * duFull;
                } else {
                    dd = ddFull;
                    duNew = duAcc + duFull;
                }
                for (int d : cdofs) duNew(d) = 0.0;
            } else if (opt.pathMethod == 1 && havePrev) {
                // Wempner–Riks 正交轨迹：**本轮校正增量** δu = duR + dλ·duF
                // 与上一收敛步增量正交：
                //   duR·prevDu + dλ·(duF·prevDu + psi²·prevDlam) = 0
                // dλ 由线性方程唯一确定，无需求根、无夹角选根。
                // 注意约束必须落在**校正增量**上，不能落在「本步总增量」
                // duAcc 上：预测步沿当前切线走，恰是正交面法向，若约束
                // 总增量会把预测点全部投影掉，路径原地不动（实测卡死）。
                const double w2 = opt.psi * opt.psi;
                const double den = duF.dot(prevDu) + w2 * prevDlam;
                dd = (std::abs(den) > 1.0e-30)
                   ? (-duR.dot(prevDu) / den) : 0.0;
                duNew = base + dd * duF;
            } else {
                // 默认（含 RIKS 第一步启动）：Crisfield 球面弧长，逐位不变。
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
                dd = (s1 >= s2) ? r1 : r2;
                duNew = (s1 >= s2) ? cand1 : cand2;
            }
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
        prevDlam = dlamAcc;
        havePrev = true;

        // ---- 第二十一次力学④：分岔点检测与分支切换（opt-in）----
        // 用已逐点输出但此前未使用的 pivotRatio 做轻量入口，再用当前一致
        // 切线自由子块的最小特征值确认（一致切线奇异 ⇔ 临界点）。确认后
        // 沿指定线性屈曲模态施加位移扰动，后续弧长步从扰动构形继续，
        // 即跟踪分岔后分支。bifurcate=0 时整段跳过，逐位不变。
        int bifDetHere = 0, bifSwHere = 0;
        double bifEigHere = 0.0;
        if (opt.bifurcate != 0 && !bifTriggered && piv < opt.bifPivot) {
            double eigRatio = 0.0;
            if (minTangentEig(u, eigRatio) && eigRatio < opt.bifEigTol) {
                bifDetHere = 1;
                bifEigHere = eigRatio;
                bifTriggered = true;   // 首次候选即锁定，避免多阶重复触发
                out.bifDetected = 1;
                out.bifStep = step;
                out.bifLambda = lambda;
                out.bifCtrlDisp = u(ctrl) - ctrlDisp0;
                out.bifMinEigRatio = eigRatio;
                out.bifPivot = piv;
                if (bifPhiOk && std::abs(opt.bifAmp) > 0.0) {
                    for (int i = 0; i < nd; ++i)
                        u(i) += opt.bifAmp * bifPhi(i);
                    for (std::size_t i = 0; i < cdofs.size(); ++i)
                        u(cdofs[i]) = cvals[i];   // 约束 DOF 不参与扰动
                    bifSwHere = 1;
                    out.bifSwitched = 1;
                    ds = std::min(ds, 0.5 * opt.ds0); // 扰动后缩小步长助 Newton 收敛
                }
            }
        }

        ArcLengthPoint3D p;
        p.step = step;
        p.lambda = lambda;
        p.ctrlDisp = u(ctrl) - ctrlDisp0;
        p.residual = resNorm;
        p.ds = ds;
        p.pivotRatio = piv;
        p.iterations = iter;
        p.lineCuts = lineCutsHere;   // 第二十一次力学①：默认恒 0
        p.alpha = alphaHere;         // 第二十一次力学①：默认恒 1.0
        p.bifDetected = bifDetHere;  // 第二十一次力学④：默认恒 0
        p.bifSwitched = bifSwHere;   // 第二十一次力学④：默认恒 0
        p.bifEigRatio = bifEigHere;  // 第二十一次力学④：默认恒 0
        if (bifPhiOk) {
            double num = 0.0, den = 0.0;
            for (int i = 0; i < nd; ++i) {
                num += u(i) * bifPhi(i);
                den += bifPhi(i) * bifPhi(i);
            }
            p.bifModalAmp = (den > 0.0) ? num / den : 0.0;
        }
        p.converged = true;
        out.points.push_back(p);
        out.steps = step;

        if (std::abs(u(ctrl) - ctrlDisp0) >= opt.ctrlTarget ||
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
