// ============================================================
// fem_3d_nl.hpp —— 三维几何非线性内核（共转列式 Co-Rotational, CR）
//
// 定位：大位移 / 大转动 / 小应变（桥梁成桥态、索承体系、施工过程）。
// 与线性内核（fem_3d.hpp）是并列关系：线性路径不动，非线性路径新开。
//
// ----------------------------------------------------------------
// 一、为什么必须用共转列式
//   桥梁成桥态的非线性来源有三：
//     (1) 缆索垂度（索自身）—— 由索单元处理（P3-2）；
//     (2) 梁柱效应 P-Δ（轴力改变弯曲刚度）—— 由几何刚度项处理；
//     (3) 大转动（悬索桥加劲梁、施工过程中的大位移）—— 由 CR 的刚体分离处理。
//   直接用 T.L./U.L. 列式，转动项与应变项耦合，小转动下退化不干净；
//   CR 的核心思想是**把刚体运动从位移场里精确剔除**，只在随体坐标系里
//   用线性单元刚度算出内力，因此：
//     - 刚体运动（任意大小）下内力恒为零（test_3d_nl 里有专门回归）；
//     - 大转动下精度只受"局部变形是否小"影响，而不受转动大小影响；
//     - 退化到线性时与 fem_3d.hpp 的线性梁完全一致。
//
// 二、列式（Crisfield / Rankin 的单元无关 CR 方法）
//   状态变量就是位移向量 u（含三个转角分量，采用**空间增量约定**：
//   转角向量 ψ 满足 节点三元组 T = expMapRot(ψ)·E0，E0 为单元参考三元组）。
//   这样状态更新无路径依赖，三元组每次迭代由 u 重算，Newton 不必存储
//   转动历史——对"施工阶段 + 成桥态"流程尤其重要。
//
//   1) 当前弦：      l = (p2+u2) - (p1+u1)， Ln = |l|， e1 = l/Ln
//   2) 刚体三元组：  e2 = 正交化( ½(T1.col(1) + T2.col(1)) 在 ⊥e1 平面上的投影 )，
//                    e3 = e1 × e2，  Rr = [e1, e2, e3]
//      （"平分向量"法：两节点三元组各贡献一半，保证刚体运动下 Rr 精确跟随）
//   3) 局部变形自由度（12）：
//        d_l = [ 0, 0, 0, log(RrᵀT1), Ln - L0, 0, 0, log(RrᵀT2) ]
//      即节点 1 平动恒为 0（随体坐标原点），节点 2 只有轴向伸长，
//      横向位移恒为 0 —— 这是 CR 对刚体平动的剔除。
//   4) 局部内力：    q_l = K_l · d_l（K_l 取参考长度 L0 的线性梁刚度）
//   5) 全局内力：    f = Bᵀ q_l， B = ∂d_l/∂d_g（12×12，解析式见下）
//   6) 一致切线：    K_t = Bᵀ K_l B + K_σ， K_σ = ∂(Bᵀ q_l)/∂u
//
// 三、B 矩阵的解析式（本文件的关键推导）
//   记 S(a)v ≡ a × v。对节点三元组施加空间小增量 δ₁、δ₂，节点平动
//   δu₁、δu₂，可证刚体三元组的转动增量为
//        δ_r = e1 × δe1 + ½((δ₁+δ₂)·e1) e1， δe1 = (I - e1e1ᵀ)(δu₂-δu₁)/Ln
//   （推导：A = ½(T₁+T₂) 的转动增量是 ½(δ₁+δ₂)；把 A 的第一列对齐到弦
//    所需的修正转动恰好抵消其中的 ⊥ 分量，只留下弦转动与轴向自旋。）
//   于是 δθ_l,i = Rrᵀ(δ_i - δ_r)，得到 B 的四个 3×3 块：
//        平动块  A   = Rrᵀ S(e1)/Ln
//        转动块  K11 = Rrᵀ(I - ½ e1e1ᵀ)， K12 = -½ Rrᵀ e1e1ᵀ
//        B(θ_l1) = [ +A , K11 , -A , K12 ]
//        B(θ_l2) = [ +A , K12 , -A , K11 ]
//        B(ΔL)   = [ -e1ᵀ , 0 , +e1ᵀ , 0 ]
//   自检（test_3d_nl 中均有对应用例）：
//     · 刚体转动 δ₁=δ₂=δ（⊥e1）且节点随刚体平动 → δθ_l = 0；
//     · 纯扭转 δ₁=δ₂=αe1 → δθ_l = 0；
//     · 仅节点 2 三元组转动、无平动 → θ_l = [0, δ]（应变能 = 2EIδ²/L 线性解）。
//
// 四、关于当前切线的实现方式
//   K_σ 目前用**单元级中心差分**求得（每单元 2n 次内力重算，n = 12/14）。
//   理由：CR 的 K_σ 解析式含大量有限转动相关项，写错不易察觉；而差分切线
//   与内力严格共轭，Newton 收敛阶与解析切线相同（二次）。
//   代价是单元形成时间约为解析式的 6~10 倍——但本内核的瓶颈在稀疏分解，
//   单元形成不是主导项。
//   TODO(P3-1b)：补解析 K_σ，并用本文件的差分切线作黄金基线回归。
// ============================================================
#pragma once

#include "fem_3d.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {

constexpr double CR_PI = 3.14159265358979323846;

// ==================== 有限旋转工具 ====================

// 反对称矩阵：S(a)·v ≡ a × v
inline Eigen::Matrix3d skew3(const Vec3& v) {
    Eigen::Matrix3d S;
    S << 0.0, -v.z(), v.y(),
         v.z(), 0.0, -v.x(),
         -v.y(), v.x(), 0.0;
    return S;
}

// vee 算子：反对称矩阵 → 轴向量（vee(S(a)) = a）
inline Vec3 vee3(const Eigen::Matrix3d& S) {
    return 0.5 * Vec3(S(2, 1) - S(1, 2),
                      S(0, 2) - S(2, 0),
                      S(1, 0) - S(0, 1));
}

// 旋转向量 → 旋转矩阵（Rodrigues；小角度走级数避免除零）
inline Eigen::Matrix3d expMapRot(const Vec3& psi) {
    const double th = psi.norm();
    const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    if (th < 1e-8) {
        const Eigen::Matrix3d S = skew3(psi);
        return I3 + S + 0.5 * S * S;      // 保留二阶项，误差 O(th³)
    }
    const Vec3 a = psi / th;
    const double c = std::cos(th), s = std::sin(th);
    return c * I3 + (1.0 - c) * (a * a.transpose()) + s * skew3(a);
}

// 旋转矩阵 → 旋转向量（含 θ→0 与 θ→π 两端的稳定处理）
//
// 关键：必须用 atan2(sinθ, cosθ) 求转角，**不能用 acos(cosθ)**。
// acos 在 θ→0 时条件数发散：c = 1 - θ²/2 在双精度下的 1e-16 舍入，
// 会被 √(2ε) 放大成 1e-8 量级的转角误差。CR 单元的局部转角正是小量，
// 这个放大效应曾让单元级差分切线的相对误差达到 1e-5（见 test_3d_nl）。
// atan2 在两端都良态：θ→0 时 ≈ 输入本身，θ→π/2 时也无放大。
inline Vec3 logMapRot(const Eigen::Matrix3d& R) {
    const Vec3 w = vee3(R - R.transpose());        // = 2 sinθ · a
    const double sn2 = 0.5 * w.norm();             // = |sinθ|
    const double c = 0.5 * (R.trace() - 1.0);      // = cosθ
    if (sn2 > 1e-12) {
        const double th = std::atan2(sn2, c);      // ∈ (0, π)，两端良态
        return (th / (2.0 * sn2)) * w;             // = θ · a
    }
    if (c > 0.0) return 0.5 * w;                   // θ ≈ 0：sinθ ≈ θ
    // θ ≈ π：sinθ ≈ 0，轴由 a aᵀ = (R + I)/2 恢复（取模最大的那一列）
    const Eigen::Matrix3d B = R + Eigen::Matrix3d::Identity();
    int col = 0;
    double best = -1.0;
    for (int j = 0; j < 3; ++j) {
        const double nrm = B.col(j).norm();
        if (nrm > best) { best = nrm; col = j; }
    }
    Vec3 a = B.col(col);
    if (a.norm() < 1e-12) return Vec3::Zero();
    a.normalize();
    // 注：θ 恰为 π 时轴的符号不可辨识。CR 单元局部转角是小量，正常不进入此分支。
    return CR_PI * a;
}

// 单元参考三元组：与线性内核 beam3DTransformation 内部的 R 完全一致
// （列 = 局部轴在全局坐标下的方向），保证退化到线性时两条路径逐位一致。
inline Eigen::Matrix3d beamRefTriad(const Vec3& p1, const Vec3& p2) {
    Vec3 ex = p2 - p1;
    const double L = std::max(ex.norm(), 1e-12);
    ex /= L;
    Vec3 ez(0.0, 0.0, 1.0);
    if (std::abs(ex.dot(ez)) > 0.999) ez = Vec3(0.0, 1.0, 0.0);
    Vec3 ey = ez.cross(ex); ey.normalize();
    ez = ex.cross(ey); ez.normalize();
    Eigen::Matrix3d R;
    R.col(0) = ex; R.col(1) = ey; R.col(2) = ez;
    return R;
}

// ==================== CR 梁单元状态 ====================

struct CRBeamState {
    Vec3 e1 = Vec3::Zero();                                 // 当前弦方向
    double Ln = 0.0;                                        // 当前弦长
    Eigen::Matrix3d Rr = Eigen::Matrix3d::Identity();       // 刚体三元组
    Eigen::Matrix<double, 12, 12> B = Eigen::Matrix<double, 12, 12>::Zero();
    Eigen::Matrix<double, 12, 1> d = Eigen::Matrix<double, 12, 1>::Zero();
    Eigen::Matrix<double, 12, 1> q = Eigen::Matrix<double, 12, 1>::Zero();
    Eigen::Matrix<double, 12, 1> f = Eigen::Matrix<double, 12, 1>::Zero();
};

// 由「当前节点位置 + 当前节点三元组 + 参考长度 + 局部刚度」算出 CR 状态。
inline CRBeamState crBeamState(const Vec3& p1c, const Vec3& p2c,
                               const Eigen::Matrix3d& T1, const Eigen::Matrix3d& T2,
                               double L0,
                               const Eigen::Matrix<double, 12, 12>& Kl) {
    CRBeamState st;

    // ---- 当前弦 ----
    const Vec3 l = p2c - p1c;
    st.Ln = std::max(l.norm(), 1e-12);
    st.e1 = l / st.Ln;

    // ---- 刚体三元组（平分向量法）----
    const Vec3 qv = 0.5 * (T1.col(1) + T2.col(1));
    Vec3 e2 = qv - st.e1 * qv.dot(st.e1);
    if (e2.norm() < 1e-8) {                       // 退化：平分向量与弦共线
        Vec3 gz(0.0, 0.0, 1.0);
        if (std::abs(st.e1.dot(gz)) > 0.999) gz = Vec3(0.0, 1.0, 0.0);
        e2 = gz - st.e1 * gz.dot(st.e1);
    }
    e2.normalize();
    const Vec3 e3 = st.e1.cross(e2);
    st.Rr.col(0) = st.e1; st.Rr.col(1) = e2; st.Rr.col(2) = e3;

    // ---- 局部转角（节点三元组相对刚体系）----
    const Vec3 t1 = logMapRot(st.Rr.transpose() * T1);
    const Vec3 t2 = logMapRot(st.Rr.transpose() * T2);

    // ---- 局部变形自由度（12）；横向平动恒为 0，节点 2 只有轴向伸长 ----
    st.d.setZero();
    st.d[3] = t1.x(); st.d[4] = t1.y(); st.d[5] = t1.z();
    st.d[6] = st.Ln - L0;
    st.d[9] = t2.x(); st.d[10] = t2.y(); st.d[11] = t2.z();

    // ---- 局部内力 ----
    st.q = Kl * st.d;

    // ---- B 矩阵（∂d_l/∂d_g，列序 [u1 θ1 u2 θ2]）----
    const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d Rt = st.Rr.transpose();
    const Eigen::Matrix3d A = (Rt * skew3(st.e1)) / st.Ln;      // 平动 → 局部转角
    const Eigen::Matrix3d H = 0.5 * (st.e1 * st.e1.transpose());
    const Eigen::Matrix3d K11 = Rt * (I3 - H);
    const Eigen::Matrix3d K12 = -(Rt * H);

    st.B.setZero();
    st.B.block<3, 3>(3, 0) = A;
    st.B.block<3, 3>(3, 3) = K11;
    st.B.block<3, 3>(3, 6) = -A;
    st.B.block<3, 3>(3, 9) = K12;
    st.B.block<1, 3>(6, 0) = -st.e1.transpose();
    st.B.block<1, 3>(6, 6) = st.e1.transpose();
    st.B.block<3, 3>(9, 0) = A;
    st.B.block<3, 3>(9, 3) = K12;
    st.B.block<3, 3>(9, 6) = -A;
    st.B.block<3, 3>(9, 9) = K11;

    // ---- 全局内力（虚功共轭：f = Bᵀ q）----
    st.f = st.B.transpose() * st.q;
    return st;
}

// ==================== 单元级内力 / 一致切线 ====================

// 单元位移自由度数：BEAM=12，BEAM7=14，其余（TRUSS/CABLE/STRUT）=6
inline int crElemDofs(const std::string& type) {
    if (type == "BEAM") return 12;
    if (type == "BEAM7") return 14;
    return 6;
}

// 单元内部数据（由参考构形预先算好，迭代中不变）
struct CRElemData {
    int id = 0;
    std::string type;
    int n = 6;
    std::vector<int> gdofs;
    Vec3 p1 = Vec3::Zero(), p2 = Vec3::Zero();
    Eigen::Matrix3d E0 = Eigen::Matrix3d::Identity();
    double Eeff = 0.0, A = 0.0, Iy = 0.0, Iz = 0.0, J = 0.0, G = 0.0, Iw = 0.0;
    double L0 = 0.0;              // 无应力长度
    double istrnEff = 0.0;        // 初应变（含收缩）
    Eigen::Matrix<double, 12, 12> Kl12 = Eigen::Matrix<double, 12, 12>::Zero();
    Eigen::Matrix<double, 14, 14> Kl14 = Eigen::Matrix<double, 14, 14>::Zero();
};

// 组装单元内部数据（参考构形 + 局部刚度 + 徐变/收缩折减）
inline CRElemData crMakeElemData(const Model3D& m, const Element3D& e, const DofMap3D& dm) {
    CRElemData d;
    d.id = e.id;
    d.type = e.type;
    d.n = crElemDofs(e.type);
    d.gdofs = elementGdofs3D(dm, e);
    d.p1 = m.nodes.at(e.n1);
    d.p2 = m.nodes.at(e.n2);
    d.E0 = beamRefTriad(d.p1, d.p2);
    d.Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
    d.A = e.A; d.Iy = e.Iy; d.Iz = e.Iz; d.J = e.J; d.G = e.G; d.Iw = e.Iw;
    const double Lg = std::max((d.p2 - d.p1).norm(), 1e-12);
    d.L0 = (e.L0 > 0.0) ? e.L0 : Lg;
    d.istrnEff = e.istrn - m.shrinkageEcs;
    if (d.n == 12 || d.n == 14) {
        d.Kl12 = beam3DLocalStiffness(d.Eeff, e.A, e.Iy, e.Iz, e.J, e.G, Lg);
    }
    if (d.n == 14) {
        d.Kl14 = beam7LocalStiffness(d.Eeff, e.A, e.Iy, e.Iz, e.J, e.G, e.Iw, Lg);
    }
    return d;
}

// 单元内力（可选：一致切线刚度、局部内力、轴力）。
// 统一用 14 维容器承载（12/6 自由度取前 n 个分量）。
inline void crElementInternal(const CRElemData& ed,
                              const Eigen::Matrix<double, 14, 1>& ue,
                              Eigen::Matrix<double, 14, 1>& fe,
                              bool wantTangent,
                              Eigen::Matrix<double, 14, 14>* Ke = nullptr,
                              Eigen::Matrix<double, 14, 1>* ql = nullptr,
                              double* axial = nullptr) {
    // 由单元位移向量算全局内力（含 12/14 两种梁与轴向单元）
    auto eval = [&](const Eigen::Matrix<double, 14, 1>& v) -> Eigen::Matrix<double, 14, 1> {
        Eigen::Matrix<double, 14, 1> f = Eigen::Matrix<double, 14, 1>::Zero();

        if (ed.n == 6) {
            // ---------------- 轴向单元（TRUSS / CABLE / STRUT）----------------
            const Vec3 u1 = v.segment<3>(0), u2 = v.segment<3>(3);
            const Vec3 p1c = ed.p1 + u1, p2c = ed.p2 + u2;
            const Vec3 l = p2c - p1c;
            const double Ln = std::max(l.norm(), 1e-12);
            const Vec3 c = l / Ln;
            const double kA = ed.Eeff * ed.A / ed.L0;
            // 机械应变 = 工程应变 - 初应变
            double N = kA * (Ln - ed.L0) - ed.Eeff * ed.A * ed.istrnEff;
            const bool tensionOnly = (ed.type == "CABLE");
            const double Nuse = (tensionOnly && N < 0.0) ? 0.0 : N;   // 索：松弛时不传力
            f.segment<3>(0) = -Nuse * c;
            f.segment<3>(3) = Nuse * c;
            if (ql) { (*ql).setZero(); (*ql)[0] = -Nuse; (*ql)[1] = Nuse; }
            if (axial) *axial = Nuse;
            return f;
        }

        // ---------------- 梁单元 ----------------
        const bool b7 = (ed.n == 14);
        const int off = b7 ? 7 : 6;
        const Vec3 u1 = v.segment<3>(0);
        const Vec3 th1 = v.segment<3>(3);
        const Vec3 u2 = v.segment<3>(off);
        const Vec3 th2 = v.segment<3>(off + 3);
        const Eigen::Matrix3d T1 = expMapRot(th1) * ed.E0;
        const Eigen::Matrix3d T2 = expMapRot(th2) * ed.E0;
        CRBeamState st = crBeamState(ed.p1 + u1, ed.p2 + u2, T1, T2, ed.L0, ed.Kl12);

        if (!b7) {
            Eigen::Matrix<double, 12, 1> d0 = Eigen::Matrix<double, 12, 1>::Zero();
            d0[6] = ed.istrnEff * ed.L0;
            const Eigen::Matrix<double, 12, 1> q = st.q - ed.Kl12 * d0;
            f.head<12>() = st.B.transpose() * q;
            if (ql) { (*ql).setZero(); ql->head<12>() = q; }
            if (axial) *axial = -q[0];          // q[0] = -N（受拉为正）
            return f;
        }

        // BEAM7：翘曲 φ 是绕局部 x 轴的标量，随体坐标下不变号
        Eigen::Matrix<double, 14, 1> d14 = Eigen::Matrix<double, 14, 1>::Zero();
        for (int i = 0; i < 6; ++i) { d14[i] = st.d[i]; d14[7 + i] = st.d[6 + i]; }
        d14[6] = v[6];
        d14[13] = v[13];
        Eigen::Matrix<double, 14, 1> d0 = Eigen::Matrix<double, 14, 1>::Zero();
        d0[7] = ed.istrnEff * ed.L0;
        const Eigen::Matrix<double, 14, 1> q = ed.Kl14 * (d14 - d0);

        Eigen::Matrix<double, 14, 14> B14 = Eigen::Matrix<double, 14, 14>::Zero();
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                B14(i, j) = st.B(i, j);
                B14(i, 7 + j) = st.B(i, 6 + j);
                B14(7 + i, j) = st.B(6 + i, j);
                B14(7 + i, 7 + j) = st.B(6 + i, 6 + j);
            }
        }
        B14(6, 6) = 1.0;
        B14(13, 13) = 1.0;
        const Eigen::Matrix<double, 14, 1> f14 = B14.transpose() * q;
        f.head<14>() = f14;
        if (ql) *ql = q;
        if (axial) *axial = -q[0];
        return f;
    };

    fe = eval(ue);

    // ---------------- 一致切线：单元级中心差分 ----------------
    // 扰动直接加在自由度上（转角按空间增量约定经 expMapRot 施加），
    // 因此差分切线与内力严格共轭，Newton 保持二次收敛。
    if (wantTangent && Ke) {
        const double h = 1e-5;   // 差分步长：1e-6 时转角自由度上噪声过大，1e-5 兼顾截断与舍入
        Eigen::Matrix<double, 14, 14> Kt = Eigen::Matrix<double, 14, 14>::Zero();
        for (int j = 0; j < ed.n; ++j) {
            Eigen::Matrix<double, 14, 1> vp = ue, vm = ue;
            vp[j] += h;
            vm[j] -= h;
            const Eigen::Matrix<double, 14, 1> fp = eval(vp);
            const Eigen::Matrix<double, 14, 1> fm = eval(vm);
            Kt.block(0, j, ed.n, 1) = (fp - fm).head(ed.n) / (2.0 * h);
        }
        // 差分带来的轻微非对称 → 对称化（保证稀疏 LDLT 行为一致）
        Kt.topLeftCorner(ed.n, ed.n) =
            0.5 * (Kt.topLeftCorner(ed.n, ed.n) + Kt.topLeftCorner(ed.n, ed.n).transpose());
        *Ke = Kt;
    }
}

// ==================== 非线性求解器 ====================

struct NonlinearOptions3D {
    int loadSteps = 10;          // 荷载步数（1 = 一步加载）
    int maxIter = 30;            // 每荷载步最大 Newton 迭代
    double tolForce = 1e-6;      // 残差相对外荷载范数
    double tolDisp = 1e-9;       // 位移增量相对范数
    bool verbose = false;
};

struct NonlinearResult3D {
    bool converged = false;
    int totalIterations = 0;
    int steps = 0;
    std::vector<double> u;                              // 自由度位移（含转角向量）
    std::vector<int> elemIds;                           // 单元 id（与下列数组同序）
    std::vector<std::vector<double>> localForces;       // 单元局部端力
    std::vector<double> axial;                          // 轴力（拉为正）
    std::vector<double> Ln;                             // 变形后弦长
};

// 三维几何非线性静力求解（Newton + 荷载步，共转列式）。
// 只使用节点荷载（3D 内核无单元荷载面荷载需预先转为节点力）。
inline NonlinearResult3D solveModel3DNonlinear(const Model3D& m,
                                               const NonlinearOptions3D& opt = {}) {
    NonlinearResult3D out;
    const DofMap3D dm = buildDofMap3D(m);
    const int nd = dm.nDofs;

    // ---- 单元数据（参考构形）----
    std::vector<CRElemData> eds;
    eds.reserve(m.elements.size());
    for (const auto& e : m.elements) {
        if (!m.nodes.count(e.n1) || !m.nodes.count(e.n2))
            throw std::runtime_error("单元 " + std::to_string(e.id) + " 引用了不存在的节点");
        eds.push_back(crMakeElemData(m, e, dm));
    }

    // ---- 约束：硬约束（含支座沉降）与弹性支承分类 ----
    std::vector<char> isC(nd, 0);
    std::vector<int> cdofs;
    std::vector<double> cvals;
    struct Spring { int dof; double k; double val; };
    std::vector<Spring> springs;
    for (const auto& c : m.constraints) {
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;                    // 纯轴向节点没有该转动自由度
        if (c.k > 0.0) { springs.push_back({d, c.k, c.value}); continue; }
        isC[d] = 1;
        auto it = std::find(cdofs.begin(), cdofs.end(), d);
        if (it == cdofs.end()) { cdofs.push_back(d); cvals.push_back(c.value); }
        else cvals[std::distance(cdofs.begin(), it)] = c.value;
    }

    // ---- 外荷载（节点荷载；初应变/收缩已计入单元内力，不再进外荷载）----
    Eigen::VectorXd Fext = Eigen::VectorXd::Zero(nd);
    for (const auto& [node, dof, val] : m.loads) {
        const int d = dm.safeDof(node, dof - 1);
        if (d >= 0) Fext(d) += val;
    }
    for (const auto& s : springs) Fext(s.dof) += s.k * s.val;

    // ---- 状态：位移（约束自由度直接置为指定值）----
    Eigen::VectorXd u = Eigen::VectorXd::Zero(nd);
    for (std::size_t i = 0; i < cdofs.size(); ++i) u(cdofs[i]) = cvals[i];

    // 单元位移提取缓冲
    auto gather = [&](const CRElemData& ed, Eigen::Matrix<double, 14, 1>& ue) {
        ue.setZero();
        for (int k = 0; k < ed.n; ++k) ue[k] = u(ed.gdofs[k]);
    };

    // 组装 F_int 与 K_t（K_t 可跳过）
    auto assemble = [&](Eigen::VectorXd& Fint, SparseMat3D& K, bool wantK) {
        Fint = Eigen::VectorXd::Zero(nd);
        std::vector<Eigen::Triplet<double>> trips;
        Eigen::Matrix<double, 14, 1> ue;
        Eigen::Matrix<double, 14, 1> fe;
        Eigen::Matrix<double, 14, 14> Ke;
        for (const auto& ed : eds) {
            gather(ed, ue);
            crElementInternal(ed, ue, fe, wantK, wantK ? &Ke : nullptr);
            for (int k = 0; k < ed.n; ++k) Fint(ed.gdofs[k]) += fe[k];
            if (wantK) {
                const double scale = std::max(1.0, Ke.topLeftCorner(ed.n, ed.n).cwiseAbs().maxCoeff());
                for (int r = 0; r < ed.n; ++r)
                    for (int c = 0; c < ed.n; ++c) {
                        const double v = Ke(r, c);
                        if (std::abs(v) > 1e-12 * scale)
                            trips.emplace_back(ed.gdofs[r], ed.gdofs[c], v);
                    }
            }
        }
        // 弹性支承：内力侧为 k·u（等效外荷载 k·value 已计入 Fext）
        for (const auto& s : springs) Fint(s.dof) += s.k * u(s.dof);
        if (wantK) {
            // 注意：K 在这里才被构造，弹簧的对角项必须放在构造之后写入，
            // 否则会写到默认构造的 0×0 矩阵上。
            K = SparseMat3D(nd, nd);
            K.setFromTriplets(trips.begin(), trips.end());
            for (const auto& s : springs) K.coeffRef(s.dof, s.dof) += s.k;
            K.makeCompressed();
        }
    };

    // ---- 初始切线：顺带做零刚度自由度诊断（沿用线性内核的严格检查）----
    {
        Eigen::VectorXd Fint0;
        SparseMat3D K0;
        assemble(Fint0, K0, true);
        for (int d = 0; d < nd; ++d) {
            if (isC[d] || K0.coeff(d, d) != 0.0) continue;
            int owner = -1;
            for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
                if (d >= dm.dofStart[i] && d < dm.dofStart[i] + dm.dofCount[i]) {
                    owner = dm.nodeIds[i];
                    break;
                }
            }
            throw std::runtime_error(
                "非线性三维求解：零刚度自由度——节点 " + std::to_string(owner) +
                " 的第 " + std::to_string(d - dm.dofStart[dm.index.at(owner)] + 1) +
                " 个自由度未约束且未被任何单元连接");
        }
    }

    // ---- Newton + 荷载步 ----
    const int nSteps = std::max(1, opt.loadSteps);
    const double Fnorm = Fext.norm();
    bool ok = true;
    for (int istep = 1; istep <= nSteps && ok; ++istep) {
        const double lam = static_cast<double>(istep) / nSteps;
        const Eigen::VectorXd Ftarget = lam * Fext;
        for (int iter = 0; iter < opt.maxIter; ++iter) {
            Eigen::VectorXd Fint;
            SparseMat3D K;
            assemble(Fint, K, true);
            Eigen::VectorXd R = Ftarget - Fint;
            // 约束自由度：增量恒为 0（u 已含指定位移）
            std::vector<double> zeros(cdofs.size(), 0.0);
            applyDirichlet(K, R, cdofs, zeros);
            const Eigen::VectorXd du = solveSparse3D(K, R);

            // 收敛判据（只看自由自由度）
            double rn = 0.0, fn = 0.0, dn = 0.0;
            for (int d = 0; d < nd; ++d) {
                if (isC[d]) continue;
                rn += R(d) * R(d);
                fn += Ftarget(d) * Ftarget(d);
                dn += du(d) * du(d);
            }
            u += du;
            for (std::size_t i = 0; i < cdofs.size(); ++i) u(cdofs[i]) = cvals[i];
            ++out.totalIterations;
            if (opt.verbose) {
                std::printf("  步 %d/%d 迭代 %d: |R|=%.3e  |du|=%.3e\n",
                            istep, nSteps, iter + 1,
                            std::sqrt(rn), std::sqrt(dn));
            }
            const double cref = std::sqrt(fn);
            if (std::sqrt(rn) <= opt.tolForce * std::max(cref, 1e-12) ||
                std::sqrt(dn) <= opt.tolDisp * std::max(1.0, u.norm())) {
                break;
            }
            if (iter == opt.maxIter - 1) ok = false;   // 未收敛
        }
        out.steps = istep;
    }
    out.converged = ok;
    if (!ok) throw std::runtime_error("非线性求解未收敛：达到最大迭代次数");

    // ---- 结果输出 ----
    out.u = toStdVec(u);
    Eigen::Matrix<double, 14, 1> ue, fe, q;
    for (const auto& ed : eds) {
        gather(ed, ue);
        q.setZero();
        double ax = 0.0;
        crElementInternal(ed, ue, fe, false, nullptr, &q, &ax);
        out.elemIds.push_back(ed.id);
        out.axial.push_back(ax);
        // 变形后弦长
        const int off = (ed.n == 14) ? 7 : ((ed.n == 12) ? 6 : 3);
        const Vec3 p1c = ed.p1 + ue.segment<3>(0);
        const Vec3 p2c = ed.p2 + ue.segment<3>(off);
        out.Ln.push_back((p2c - p1c).norm());
        std::vector<double> qv;
        const int nq = (ed.n == 14) ? 14 : ((ed.n == 12) ? 12 : 2);
        for (int i = 0; i < nq; ++i) qv.push_back(q[i]);
        out.localForces.push_back(qv);
    }
    return out;
}

}  // namespace fem
