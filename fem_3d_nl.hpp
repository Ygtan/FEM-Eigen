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
// 四、关于切线的实现方式（已从差分改为解析）
//   早期版本用**单元级中心差分**求一致切线（每单元 2n 次内力重算），理由是
//   CR 的几何刚度解析式容易写错。实测证明这个取舍不成立：单元级微基准
//   （bench/nl_bench）显示差分把单元成本放大 26.7 倍、占单元耗时 96%，
//   而稀疏分解只占小头——于是"耗时随单元数线性增长"而不是随自由度超线性。
//   现在切线走解析式（crBeamTangent12 / 索的 G 块），差分路径保留为
//   tangentByFD = true，作**黄金基线**回归（test_3d_nl 逐元素对照）。
//   实测：单元切线 20.2 µs → 3.1 µs（切线本身 8.5x），单次非线性求解、
//   索长迭代、成桥态闭环整体快 3.3~3.4x；解析式与差分一致到 4e-11。
//
//   推导中踩过的三个坑（都写进了对应函数注释）：
//     ① f = Bᵀq 的全导 = Bᵀ·Kl·∂d/∂u + (∂Bᵀ/∂u)q，δq 项只属于前者；
//        重复计入会让切线正好放大近 2 倍。
//     ② 节点转角扰动是**分量**扰动：∂T/∂ψ_j = S(T·E0ᵀ·J_r(ψ)e_j)·T，
//        需要的是 J_r 本身（不是 J_r⁻¹，也不是忽略）。忽略得 1.6e-3 误差。
//     ③ (I-H) 与 Rr **不可交换**：(I-H)Rr q ≠ Rr(I-H)q，旋转块会错 13%。
// ============================================================
#pragma once

#include "fem_3d.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
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

// 对数映射的**左** Jacobian 逆：δR = S(v)R  ⇒  δt = J_l⁻¹(t)·v
//
//      J_l⁻¹(t) = I - ½S(t) + (1/θ² - (1+cosθ)/(2θ·sinθ))·S(t)²
//
// 符号必须是 **-½S(t)**（左扰动）：本内核的节点转角按**空间增量**约定施加
// （T = expMap(ψ)·E0），扰动在左侧，对应 J_l。数值验证见 test_3d_nl 的
// "∂d/∂u 与差分一致"用例——这里若取 +½S(t) 会立刻被它抓出来。
// θ→0 用级数（1/12 + θ²/720），θ→π 取极限值 1/θ²（CR 局部转角不会到 π）。
inline Eigen::Matrix3d logMapLeftJacInv(const Vec3& t) {
    const double th2 = t.squaredNorm();
    const double th = std::sqrt(th2);
    const Eigen::Matrix3d St = skew3(t);
    const Eigen::Matrix3d St2 = St * St;
    const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    // 局部转角是**变形**转角（共转列式已剔除刚体转动），恒为小量，
    // 因此级数阈值取 0.3 rad：c2 = 1/12 + θ²/720 在 θ=0.3 时的下一项
    // θ⁴/30240 ≈ 2.7e-7，相对 c2 只 3e-6，且省掉每单元 2 次 sin/cos。
    if (th < 0.3) return I3 - 0.5 * St + (1.0 / 12.0 + th2 / 720.0) * St2;
    const double sn = std::sin(th);
    const double c2 = (std::abs(sn) < 1e-9)
                          ? 1.0 / th2
                          : 1.0 / th2 - (1.0 + std::cos(th)) / (2.0 * th * sn);
    return I3 - 0.5 * St + c2 * St2;
}

// 对数映射的**右** Jacobian：Exp(ψ + δ) = Exp(ψ)·Exp(J_r(ψ)·δ)
//
//      J_r(t) = I - ((1-cosθ)/θ²)·S(t) + ((θ-sinθ)/θ³)·S(t)²
//
// 用途：节点转角的自由度是旋转向量的**分量**，差分与 Newton 都是直接扰动
// 分量，故 ∂T/∂ψ_j = S(w)·T，其中 w = T·E0ᵀ·J_r(ψ)·e_j —— 左乘向量既不是
// 简单的 e_j，也不是 T·J_r⁻¹e_j（两个都试过：忽略得 1.6e-3 误差，用逆得
// 3.1e-3，用 J_r 本身才降到 1e-9 以内）。见 test_3d_nl 的对照用例。
inline Eigen::Matrix3d logMapRightJac(const Vec3& t) {
    const double th2 = t.squaredNorm();
    const double th = std::sqrt(th2);
    const Eigen::Matrix3d St = skew3(t);
    const Eigen::Matrix3d St2 = St * St;
    const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    if (th < 0.3) {
        const double c1 = 0.5 - th2 / 24.0;
        const double c2 = 1.0 / 6.0 - th2 / 120.0;
        return I3 - c1 * St + c2 * St2;
    }
    const double c1 = (1.0 - std::cos(th)) / th2;
    const double c2 = (th - std::sin(th)) / (th2 * th);
    return I3 - c1 * St + c2 * St2;
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

// ==================== 解析一致切线（CR 梁）====================
//
// 单元级微基准（bench/nl_bench）实测：差分切线要把单元成本放大 26.7 倍，
// 占单元总耗时的 96%（每单元 2n 次内力评估，每次含 4 个超越函数）。
// 这里给出解析式，把切线成本压到「1 次内力评估 + O(n³) 矩阵运算」。
//
// 把 B 的三个块代入 f = Bᵀq，按 [p1 ψ1 p2 ψ2] 分块展开，得显式结构
// （记 w = q_t1 + q_t2，H = ½e1e1ᵀ，q_ax = q[6]；Rr 只以**左乘**出现，
//  故 f 对 Rr 是线性的，对 Rr 的全部依赖可以收进一个刚体自旋 ω）：
//     f_p1 = -S(e1)·Rr·w/Ln - e1·q_ax
//     f_ψ1 = Rr·[(I-H)q_t1 - H q_t2]
//     f_p2 = -f_p1
//     f_ψ2 = Rr·[-H q_t1 + (I-H)q_t2]
// 于是 δf = Bᵀ·Kl·Bd（材料项）+ ∂f/∂(e1,ω,δq,δLn)（几何项），
// 其中 Bd = ∂d/∂u 是局部变形的精确导数（含对数映射的 Jacobian）。
//
// ω 的三个分量：δe1 = ω×e1 给出垂直于 e1 的两个；第三个由
//     ω·e1 = δe2·e3,   e2 = P·qv/|P·qv|,  P = I - e1e1ᵀ
// 补上（δe2 的 e1 分量对真旋转是冗余的，不必用）。
//
// 退化情形（平分向量与弦共线，nQ→0）与 crBeamState 的兜底分支同源，
// 此时 ω 无定义，返回 false 由调用方退回差分。

// 刚体三元组的自旋 ω（δRr = S(ω)Rr）
inline Vec3 crRigidSpin(const CRBeamState& st, const Vec3& qv, double nQ,
                        const Vec3& de1, const Vec3& dqv) {
    const Vec3 e3 = st.Rr.col(2);
    const double e1qv = st.e1.dot(qv);
    const double wa = (dqv.dot(e3) - de1.dot(e3) * e1qv) / std::max(nQ, 1e-12);
    return st.e1.cross(de1) + wa * st.e1;
}

// 单元切向的运动学辅助量（避免重复算 qv / P / Jacobian）
struct CRTangentKin {
    Vec3 e1 = Vec3::Zero();
    Vec3 qv = Vec3::Zero();
    double nQ = 0.0;
    double Ln = 0.0;
    Eigen::Matrix3d P = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Rt = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d J1 = Eigen::Matrix3d::Identity();    // 左 Jacobian 逆 @ t1
    Eigen::Matrix3d J2 = Eigen::Matrix3d::Identity();    // 左 Jacobian 逆 @ t2
    Eigen::Matrix3d Jp1 = Eigen::Matrix3d::Identity();   // 右 Jacobian @ ψ1
    Eigen::Matrix3d Jp2 = Eigen::Matrix3d::Identity();   // 右 Jacobian @ ψ2
    // 分量扰动 ψ_j += h 对应的**等效空间转角向量** w = T·E0ᵀ·J_r(ψ)e_j
    // （T = exp(ψ)E0 ⇒ ∂T/∂ψ_j = S(w)T）。参考态下 w = e_j，退回常规形式。
    Eigen::Matrix3d W1 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d W2 = Eigen::Matrix3d::Identity();
};

inline bool crTangentKin(const CRBeamState& st, const Eigen::Matrix3d& T1,
                         const Eigen::Matrix3d& T2, const Eigen::Matrix3d& E0,
                         const Vec3& psi1, const Vec3& psi2, CRTangentKin& k) {
    k.e1 = st.e1;
    k.Ln = st.Ln;
    k.qv = 0.5 * (T1.col(1) + T2.col(1));
    k.nQ = (k.qv - k.e1 * k.e1.dot(k.qv)).norm();
    if (k.nQ < 1e-8) return false;              // 与 crBeamState 退化分支一致
    k.P = Eigen::Matrix3d::Identity() - k.e1 * k.e1.transpose();
    k.Rt = st.Rr.transpose();
    k.J1 = logMapLeftJacInv(st.d.segment<3>(3));
    k.J2 = logMapLeftJacInv(st.d.segment<3>(9));
    k.Jp1 = logMapRightJac(psi1);
    k.Jp2 = logMapRightJac(psi2);
    k.W1 = T1 * (E0.transpose() * k.Jp1);
    k.W2 = T2 * (E0.transpose() * k.Jp2);
    return true;
}

// 第 j 个全局自由度方向上的局部运动学增量（de1, ω, δLn 与 δd）。
// 转角方向上的 δψ_j 先换算成"等效空间转角向量" T·J_r⁻¹(ψ)e_j 再往下传。
inline void crTangentDir(const CRBeamState& st, const CRTangentKin& k,
                         const Eigen::Matrix3d& T1, const Eigen::Matrix3d& T2,
                         int j, Vec3& de1, Vec3& om, double& dLn,
                         Eigen::Matrix<double, 12, 1>& dd) {
    de1.setZero(); om.setZero(); dLn = 0.0;
    dd.setZero();
    Vec3 dp1 = Vec3::Zero(), dp2 = Vec3::Zero(), dqv = Vec3::Zero();
    if (j < 3)            { de1 = -k.P.col(j) / k.Ln;          dLn = -k.e1[j]; }
    else if (j < 6)       { dp1 = k.W1.col(j - 3);
                            dqv = -0.5 * T1.col(1).cross(dp1); }
    else if (j < 9)       { de1 = k.P.col(j - 6) / k.Ln;       dLn = k.e1[j - 6]; }
    else                  { dp2 = k.W2.col(j - 9);
                            dqv = -0.5 * T2.col(1).cross(dp2); }

    om = crRigidSpin(st, k.qv, k.nQ, de1, dqv);
    dd.segment<3>(3) = k.J1 * (k.Rt * (dp1 - om));
    dd[6] = dLn;
    dd.segment<3>(9) = k.J2 * (k.Rt * (dp2 - om));
}

// ∂d/∂u 的解析式（12×12；行 = 局部自由度 d，列 = [u1 θ1 u2 θ2]）
inline bool crBeamLocalDerivative(const CRBeamState& st, const Eigen::Matrix3d& T1,
                                  const Eigen::Matrix3d& T2, const Eigen::Matrix3d& E0,
                                  const Vec3& psi1, const Vec3& psi2,
                                  Eigen::Matrix<double, 12, 12>& Bd) {
    CRTangentKin k;
    if (!crTangentKin(st, T1, T2, E0, psi1, psi2, k)) return false;
    Bd.setZero();
    Bd.block<1, 3>(6, 0) = -k.e1.transpose();
    Bd.block<1, 3>(6, 6) = k.e1.transpose();
    for (int j = 0; j < 12; ++j) {
        Vec3 de1, om;
        double dLn;
        Eigen::Matrix<double, 12, 1> dd;
        crTangentDir(st, k, T1, T2, j, de1, om, dLn, dd);
        Bd.col(j) = dd;
    }
    return true;
}

// 解析一致切线 K = Bᵀ·Kl·Bd + K_g（12×12）。退化时返回 false。
// Kg/Bd 非空时一并输出（BEAM7 的 14 维组装要用）。
inline bool crBeamTangent12(const CRBeamState& st, const Eigen::Matrix3d& T1,
                            const Eigen::Matrix3d& T2, const Eigen::Matrix3d& E0,
                            const Vec3& psi1, const Vec3& psi2,
                            const Eigen::Matrix<double, 12, 12>& Kl,
                            const Eigen::Matrix<double, 12, 1>& q,
                            Eigen::Matrix<double, 12, 12>& K,
                            Eigen::Matrix<double, 12, 12>* KgOut = nullptr,
                            Eigen::Matrix<double, 12, 12>* BdOut = nullptr) {
    CRTangentKin k;
    if (!crTangentKin(st, T1, T2, E0, psi1, psi2, k)) return false;

    // 几何项用到的常量（一次算好，逐方向复用）
    const Vec3 qa = q.segment<3>(3), qb = q.segment<3>(9);
    const Vec3 w = qa + qb;
    const double qax = q[6];
    const Eigen::Matrix3d H = 0.5 * (k.e1 * k.e1.transpose());
    const Vec3 Rw = st.Rr * w;
    // f_ψ1 = Rr·q_t1 - H·Rr·w，f_ψ2 = Rr·q_t2 - H·Rr·w
    // 注意 (I-H) 与 Rr **不可交换**：(I-H)Rr·q_t1 ≠ Rr·(I-H)·q_t1。
    // （第一版按 Rr(I-H)q_t1 写，旋转块误差 13%，被"固定 q 只扰动 B"的
    //   数值校验抓出来——平动块因为 Rr 直接乘向量而恰好正确。）
    const Vec3 y1 = st.Rr * qa;
    const Vec3 y2 = st.Rr * qb;

    K.setZero();
    Eigen::Matrix<double, 12, 12> Kg = Eigen::Matrix<double, 12, 12>::Zero();
    Eigen::Matrix<double, 12, 12> Bd = Eigen::Matrix<double, 12, 12>::Zero();
    Bd.block<1, 3>(6, 0) = -k.e1.transpose();
    Bd.block<1, 3>(6, 6) = k.e1.transpose();

    for (int j = 0; j < 12; ++j) {
        Vec3 de1, om;
        double dLn;
        Eigen::Matrix<double, 12, 1> dd;
        crTangentDir(st, k, T1, T2, j, de1, om, dLn, dd);
        Bd.col(j) = dd;

        // ---- 几何项：显式结构对 (e1, ω, δLn) 求导 ----
        // 注意：**不能**把 ∂f/∂q·δq 的项放进来——那部分正是材料项 Bᵀ·Kl·Bd
        // 已经覆盖的内容（f = Bᵀq 对 u 的全导 = Bᵀ Kl Bd + (∂Bᵀ/∂u)q），
        // 重复计入会让切线正好放大近 2 倍（第一版就踩了这个坑）。
        // 这里只保留 B 显式依赖 e1 与 Rr 的项。
        // δH 作用在 Rr·w 上（不是 w 上）
        const Vec3 dHRw = 0.5 * ((de1.dot(Rw)) * k.e1 + (k.e1.dot(Rw)) * de1);
        const Vec3 hTerm = H * (om.cross(Rw));

        // f_p1 及其负 = f_p2
        const Vec3 dfp1 = Rw.cross(de1) / k.Ln
                        + k.e1.cross(Rw.cross(om)) / k.Ln
                        + (k.e1.cross(Rw)) * (dLn / (k.Ln * k.Ln))
                        - de1 * qax;
        Kg.block<3, 1>(0, j) += dfp1;
        Kg.block<3, 1>(6, j) -= dfp1;

        // f_ψ1 / f_ψ2：δ(Rr q) = -（Rr q)×ω，δH·Rw，以及 -H·δ(Rr w)
        Kg.block<3, 1>(3, j) = -(y1.cross(om)) - dHRw - hTerm;
        Kg.block<3, 1>(9, j) = -(y2.cross(om)) - dHRw - hTerm;
    }
    // ---- 材料项一次性成块：Bᵀ·Kl·Bd ----
    // 逐列做两个 12×1 的矩阵-向量乘（每列 288 次乘加，SIMD 利用率极低）会
    // 吃掉大半时间；改成两次 12×12 整块乘（Eigen 走分块 + 向量化）快得多。
    const Eigen::Matrix<double, 12, 12> Km = st.B.transpose() * (Kl * Bd);
    K = Km + Kg;
    // 对称化：f = Bᵀq 中 B ≠ ∂d/∂u（B 是一阶近似），故严格 Jacobian 微不对称；
    // 与旧差分路径保持一致（那边也做了对称化，LDLT 只取三角）。
    K = 0.5 * (K + K.transpose());
    if (KgOut) *KgOut = Kg;
    if (BdOut) *BdOut = Bd;
    return true;
}

// ==================== 缆索垂度（sag）====================
//
// 力-伸长关系与它的反解（cableTensionFromChord）现在放在 fem_3d.hpp 里，
// 因为**线性施工阶段调索也要用它**：施工阶段把"目标索力"翻译成"无应力
// 索长"，成桥态找形再把这个 L0 交给非线性求解——两边必须是同一个式子，
// 否则阶段索力与成桥态索力没法互相对账。本文件通过 #include "fem_3d.hpp"
// 直接使用，不再重复定义。

// Ernst 等效弹性模量（= 上述力-伸长关系的切线 dT/dLn·L0/A，与内力严格共轭）
inline double cableErnstModulus(double E, double A, double w, double Lh, double T,
                                double L0 = 0.0) {
    if (w <= 0.0 || Lh <= 1e-12 || T <= 0.0) return E;
    const double Lr = (L0 > 0.0) ? L0 : Lh;      // 未给 L0 时退化为标准 Ernst 式
    return E / (1.0 + w * w * Lh * Lh * Lh * E * A / (12.0 * T * T * T * Lr));
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
    double w = 0.0;               // 缆索单位长度重量 N/m（>0 启用垂度效应）
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
    d.w = e.w;
    // 正的收缩应变与正的初应变相反：初应变伸长，收缩缩短。
    // 阶段循环把模型收缩写进 e.ecs，且 mc.shrinkageEcs 保持 0，避免两个字段各算一遍。
    d.istrnEff = e.istrn - m.shrinkageEcs - e.ecs;
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
                              double* axial = nullptr,
                              bool tangentByFD = false) {
    bool tangentOk = false;      // 解析切线是否可用（退化时为 false，退回差分）
    // 由单元位移向量算全局内力（含 12/14 两种梁与轴向单元）；
    // Kout 非空时同时给出**解析**一致切线
    auto eval = [&](const Eigen::Matrix<double, 14, 1>& v,
                    Eigen::Matrix<double, 14, 14>* Kout) -> Eigen::Matrix<double, 14, 1> {
        Eigen::Matrix<double, 14, 1> f = Eigen::Matrix<double, 14, 1>::Zero();

        if (ed.n == 6) {
            // ---------------- 轴向单元（TRUSS / CABLE / STRUT）----------------
            const Vec3 u1 = v.segment<3>(0), u2 = v.segment<3>(3);
            const Vec3 p1c = ed.p1 + u1, p2c = ed.p2 + u2;
            const Vec3 l = p2c - p1c;
            const double Ln = std::max(l.norm(), 1e-12);
            const Vec3 c = l / Ln;
            const double Lh = std::sqrt(l.x() * l.x() + l.y() * l.y());
            const double EA = ed.Eeff * ed.A;
            const bool tensionOnly = (ed.type == "CABLE");
            const bool sag = tensionOnly && ed.w > 0.0;
            double N = 0.0;
            if (sag) {
                // 缆索垂度：弦长 → 索力走浅垂度弹性悬链线的精确反解
                // （N = N(Ln, Lh)，下面给出它的两个偏导，二者共同构成切线）
                N = cableTensionFromChord(Ln, ed.L0, EA, ed.w, Lh);
            } else {
                N = EA / ed.L0 * (Ln - ed.L0);
            }
            N -= EA * ed.istrnEff;                   // 初应变（可作为索力调整手段）
            const bool slack = (tensionOnly && N < 0.0);
            if (slack) N = 0.0;                      // 索：松弛时不传力
            f.segment<3>(0) = -N * c;
            f.segment<3>(3) = N * c;
            if (ql) { (*ql).setZero(); (*ql)[0] = -N; (*ql)[1] = N; }
            if (axial) *axial = N;

            // ---- 解析切线 ----
            // f = [-N·c; N·c]，N = N(Ln, Lh)，Ln、Lh 都只通过 l = p2c-p1c 依赖位移：
            //     a = k1·c + k2·g,  k1 = ∂N/∂Ln, k2 = ∂N/∂Lh, g = (lx,ly,0)/Lh
            //     K = [[G, -G], [-G, G]],   G = c·aᵀ - (N/Ln)·(I - c cᵀ)
            // 松弛索不传力 ⇒ 切线为零（与差分的活动集一致）。
            if (Kout) {
                Kout->setZero();
                if (!slack) {
                    double k1 = EA / ed.L0, k2 = 0.0;
                    if (sag && Lh > 1e-12) {
                        // 隐式求导：g(T,Ln,Lh) = L0 + L0T/EA - cSag/T² - Ln = 0
                        //   ∂T/∂Ln = 1/(L0/EA + 2cSag/T³)  ← 就是 Ernst 切线 Eeq·A/L0
                        //   ∂T/∂Lh = k1·w²Lh²/(8T²)
                        const double Tstar = N + EA * ed.istrnEff;   // 未扣初应变的索力
                        const double cSag = ed.w * ed.w * Lh * Lh * Lh / 24.0;
                        if (Tstar > 1e-6) {
                            k1 = 1.0 / (ed.L0 / EA + 2.0 * cSag / (Tstar * Tstar * Tstar));
                            k2 = k1 * (ed.w * ed.w * Lh * Lh / (8.0 * Tstar * Tstar));
                        }
                    }
                    const Vec3 g = (k2 != 0.0 && Lh > 1e-12)
                                       ? Vec3(l.x() / Lh, l.y() / Lh, 0.0)
                                       : Vec3::Zero();
                    const Vec3 a = k1 * c + k2 * g;
                    const Eigen::Matrix3d Pm =
                        Eigen::Matrix3d::Identity() - c * c.transpose();
                    const Eigen::Matrix3d G = c * a.transpose() + (N / Ln) * Pm;
                    Kout->block<3, 3>(0, 0) = G;
                    Kout->block<3, 3>(0, 3) = -G;
                    Kout->block<3, 3>(3, 0) = -G;
                    Kout->block<3, 3>(3, 3) = G;
                }
                tangentOk = true;
            }
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
            if (Kout) {
                Eigen::Matrix<double, 12, 12> K12;
                if (crBeamTangent12(st, T1, T2, ed.E0, th1, th2, ed.Kl12, q, K12)) {
                    Kout->setZero();
                    Kout->topLeftCorner(12, 12) = K12;
                    tangentOk = true;
                }
            }
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
        if (Kout) {
            // 12 维的几何刚度与 ∂d/∂u（B14 对位移的依赖与 B 完全同源），
            // 材料项必须用 14 维的 Kl14（含翘曲—扭转耦合，不能靠嵌入代替）。
            Eigen::Matrix<double, 12, 1> q12;
            for (int i = 0; i < 6; ++i) { q12[i] = q[i]; q12[6 + i] = q[7 + i]; }
            Eigen::Matrix<double, 12, 12> Ktmp, Kg12, Bd12;
            if (crBeamTangent12(st, T1, T2, ed.E0, th1, th2, ed.Kl12, q12, Ktmp, &Kg12, &Bd12)) {
                const auto idx = [](int i) { return (i < 6) ? i : i + 1; };
                Eigen::Matrix<double, 14, 14> Bd14 = Eigen::Matrix<double, 14, 14>::Zero();
                for (int i = 0; i < 12; ++i)
                    for (int j = 0; j < 12; ++j) Bd14(idx(i), idx(j)) = Bd12(i, j);
                Bd14(6, 6) = 1.0;
                Bd14(13, 13) = 1.0;
                Eigen::Matrix<double, 14, 14> K14 = B14.transpose() * (ed.Kl14 * Bd14);
                for (int i = 0; i < 12; ++i)
                    for (int j = 0; j < 12; ++j) K14(idx(i), idx(j)) += Kg12(i, j);
                K14 = 0.5 * (K14 + K14.transpose());
                Kout->setZero();
                *Kout = K14;
                tangentOk = true;
            }
        }
        return f;
    };

    // ---------------- 一致切线 ----------------
    // 默认走解析式（见上文 crBeamTangent12 / 轴向单元的 G 块）：成本约等于
    // 1 次内力评估 + O(n³) 矩阵运算，而差分要 2n 次内力评估（实测 26.7 倍）。
    // tangentByFD = true 保留原差分路径，作为解析式的**黄金基线**（回归测试
    // 逐元素对照）与退化情形的兜底。
    const bool wantAnalytic = wantTangent && Ke && !tangentByFD;
    Eigen::Matrix<double, 14, 14> Kt = Eigen::Matrix<double, 14, 14>::Zero();
    fe = eval(ue, wantAnalytic ? &Kt : nullptr);

    if (wantTangent && Ke) {
        if (wantAnalytic && tangentOk) {
            *Ke = Kt;
        } else {
            // 单元级中心差分：扰动直接加在自由度上（转角按空间增量约定经
            // expMapRot 施加），因此与内力严格共轭。
            const double h = 1e-5;   // 1e-6 时转角自由度上噪声过大，1e-5 兼顾截断与舍入
            for (int j = 0; j < ed.n; ++j) {
                Eigen::Matrix<double, 14, 1> vp = ue, vm = ue;
                vp[j] += h;
                vm[j] -= h;
                const Eigen::Matrix<double, 14, 1> fp = eval(vp, nullptr);
                const Eigen::Matrix<double, 14, 1> fm = eval(vm, nullptr);
                Kt.block(0, j, ed.n, 1) = (fp - fm).head(ed.n) / (2.0 * h);
            }
            // 差分带来的轻微非对称 → 对称化（保证稀疏 LDLT 行为一致）
            Kt.topLeftCorner(ed.n, ed.n) =
                0.5 * (Kt.topLeftCorner(ed.n, ed.n) + Kt.topLeftCorner(ed.n, ed.n).transpose());
            *Ke = Kt;
        }
    }
}

// ==================== 非线性求解器 ====================

struct NonlinearOptions3D {
    int loadSteps = 10;          // 荷载步数（1 = 一步加载）
    int maxIter = 30;            // 每荷载步最大 Newton 迭代
    double tolForce = 1e-6;      // 残差相对外荷载范数
    double tolDisp = 1e-9;       // 位移增量相对范数
    bool verbose = false;
    // Newton 迭代内复用符号分解（只做数值分解）。默认开；关掉可用于对照计时。
    bool reuseFactorization = true;
    // 线性求解后端。默认 SimplicialLDLT；装了 SuiteSparse 时可显式选 CHOLMOD
    // （大规模稀疏上通常快数倍——这是把 P1-7 里"从未接到求解路径"的 CHOLMOD
    // 真正用起来的地方）。不可用时 Factorization3D 会自动回退。
    SolverType solver = SolverType::SimplicialLDLT;
};

// ---- 内层 Newton 过程留痕（只记录，不参与任何判据 / 步长 / 收敛计算）----
// 与外层索长迭代 / 几何闭合的 IterationTrace3D 是两层不同的迭代：外层改模型
// 参数（L0 或制造线形 ΔX），内层在给定模型上做共转列式 Newton + 荷载步求平衡。
// 这一层以前只打印在 verbose 控制台里，结果文件里完全没有；暴露它是为了能
// 从计算书直接看出「本次平衡是否真收敛、残差按什么量级降下来」。
struct NewtonStep3D {
    int loadStep = 0;        // 荷载步 1..loadSteps
    int iter = 0;            // 步内 Newton 迭代 1..maxIter
    double lambda = 0.0;     // 荷载因子 λ = loadStep / loadSteps
    double residual = 0.0;   // |R|（自由自由度，力）
    double du = 0.0;         // |du|（自由自由度，位移）
    double ref = 0.0;        // |Ftarget|，残差判据的参考范数
    bool refactor = false;   // 本迭代是否重建了符号分解（稀疏模式变化 / 首次）
    bool converged = false;  // 本迭代满足收敛判据，随后 break
};

struct NewtonTrace3D {
    int loadSteps = 0;
    int maxIter = 0;
    double tolForce = 0.0, tolDisp = 0.0;
    int totalIterations = 0;
    int refactorizations = 0;
    bool converged = false;
    double residualFirst = 0.0, residualLast = 0.0;
    std::vector<NewtonStep3D> steps;
};

struct NonlinearResult3D {
    bool converged = false;
    int totalIterations = 0;
    int steps = 0;
    NewtonTrace3D trace;                                // 内层 Newton 逐迭代留痕
    std::vector<double> u;                              // 自由度位移（含转角向量）
    std::vector<int> elemIds;                           // 单元 id（与下列数组同序）
    std::vector<std::vector<double>> localForces;       // 单元局部端力
    std::vector<double> axial;                          // 轴力（拉为正）
    std::vector<double> Ln;                             // 变形后弦长
};

// 收敛构形的切线（供成桥态几何迭代、影响矩阵复用分解等使用）。
// K 已施加 Dirichlet（约束行/列清零、对角置 1），可直接分解回代。
struct NLTangent3D {
    SparseMat3D K;
    std::vector<int> cdofs;
};

// 三维几何非线性静力求解（Newton + 荷载步，共转列式）。
// 只使用节点荷载（3D 内核无单元荷载面荷载需预先转为节点力）。
// tangent 非空时，额外返回收敛构形的一致切线（用于敏感度/影响矩阵）。
inline NonlinearResult3D solveModel3DNonlinear(const Model3D& m,
                                               const NonlinearOptions3D& opt = {},
                                               NLTangent3D* tangent = nullptr) {
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
    // 缆索自重：两端各承担 w·L0/2（向下）。
    // 物理上索把自重传给锚固点，弦向只传索力，故这里按外荷载施加。
    for (const auto& e : m.elements) {
        if (e.type != "CABLE" || e.w <= 0.0) continue;
        const double Lc = (e.L0 > 0.0) ? e.L0
                                       : (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm();
        const double half = -0.5 * e.w * Lc;
        const int d1 = dm.safeDof(e.n1, 2);      // dof 3 = z（竖向）
        const int d2 = dm.safeDof(e.n2, 2);
        if (d1 >= 0) Fext(d1) += half;
        if (d2 >= 0) Fext(d2) += half;
    }

    // Newton 迭代用的可复用分解（整个求解只 compute 一次，其余 refactorize）
    std::unique_ptr<Factorization3D> fac;

    // ---- 状态：位移（约束自由度直接置为指定值）----
    Eigen::VectorXd u = Eigen::VectorXd::Zero(nd);
    for (std::size_t i = 0; i < cdofs.size(); ++i) u(cdofs[i]) = cvals[i];

    // 单元位移提取缓冲
    auto gather = [&](const CRElemData& ed, Eigen::Matrix<double, 14, 1>& ue) {
        ue.setZero();
        for (int k = 0; k < ed.n; ++k) ue[k] = u(ed.gdofs[k]);
    };

    // 组装 F_int 与 K_t（K_t 可跳过）
    //
    // 稀疏模式"单调增长"：位置一旦进入 pattern 就不再移除，只有新增时才
    // 触发重新分解。这是 Newton 迭代里复用符号分解的前提——Eigen 的
    // factorize() 要求矩阵稀疏结构逐次不变，若模式抖动就会读到错误的符号
    // 分解（表现为 Newton 突然不收敛）。
    //
    // 关键坑：pattern 不能在 u=0 时一次性定死。几何刚度 Kσ 正比于内力，
    // 在参考态恒为零，它的非零位置此刻全部落在阈值以下；一旦结构受力，
    // Kσ 出现而 pattern 里没有 → 切线缺一大块 → 收敛失败（实测即如此）。
    // 单调增长让 Kσ 的结构在第一次受力时自动补进来，之后就稳定了。
    std::vector<std::vector<std::pair<int, int>>> pat;
    std::vector<std::vector<char>> inPat;
    bool patChanged = false;

    auto assemble = [&](Eigen::VectorXd& Fint, SparseMat3D& K, bool wantK) {
        Fint = Eigen::VectorXd::Zero(nd);
        std::vector<Eigen::Triplet<double>> trips;
        Eigen::Matrix<double, 14, 1> ue;
        Eigen::Matrix<double, 14, 1> fe;
        Eigen::Matrix<double, 14, 14> Ke;
        if (wantK && pat.empty()) {
            pat.resize(eds.size());
            inPat.assign(eds.size(), std::vector<char>(196, 0));
        }
        for (std::size_t ei = 0; ei < eds.size(); ++ei) {
            const auto& ed = eds[ei];
            gather(ed, ue);
            crElementInternal(ed, ue, fe, wantK, wantK ? &Ke : nullptr);
            for (int k = 0; k < ed.n; ++k) Fint(ed.gdofs[k]) += fe[k];
            if (!wantK) continue;
            const double scale = std::max(1.0, Ke.topLeftCorner(ed.n, ed.n).cwiseAbs().maxCoeff());
            for (int r = 0; r < ed.n; ++r)
                for (int c = 0; c < ed.n; ++c) {
                    const double v = Ke(r, c);
                    const int p = r * 14 + c;
                    if (inPat[ei][p]) {
                        trips.emplace_back(ed.gdofs[r], ed.gdofs[c], v);
                    } else if (std::abs(v) > 1e-12 * scale) {
                        inPat[ei][p] = 1;
                        pat[ei].push_back({r, c});
                        trips.emplace_back(ed.gdofs[r], ed.gdofs[c], v);
                        patChanged = true;
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
    out.trace.loadSteps = nSteps;
    out.trace.maxIter = opt.maxIter;
    out.trace.tolForce = opt.tolForce;
    out.trace.tolDisp = opt.tolDisp;
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
            // 复用符号分解：矩阵结构不变，只有数值变 → 只做数值分解，
            // 省掉重复的符号分解与列排序（大模型上通常是 factorize 的大头）。
            const bool didRefactor = !fac || patChanged || !opt.reuseFactorization;
            if (didRefactor) {
                fac.reset(new Factorization3D(K, opt.solver));
                patChanged = false;
                ++out.trace.refactorizations;
            } else {
                fac->refactorize(K);
            }
            const Eigen::VectorXd du = fac->solve(R);

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
            NewtonStep3D ns;
            ns.loadStep = istep;
            ns.iter = iter + 1;
            ns.lambda = lam;
            ns.residual = std::sqrt(rn);
            ns.du = std::sqrt(dn);
            ns.ref = std::sqrt(fn);
            ns.refactor = didRefactor;
            ns.converged = (ns.residual <= opt.tolForce * std::max(ns.ref, 1e-12) ||
                            ns.du <= opt.tolDisp * std::max(1.0, u.norm()));
            out.trace.steps.push_back(ns);   // 只留痕；判据数值与改动前逐位一致
            if (opt.verbose) {
                std::printf("  步 %d/%d 迭代 %d: |R|=%.3e  |du|=%.3e\n",
                            istep, nSteps, iter + 1, ns.residual, ns.du);
            }
            if (ns.converged) break;
            if (iter == opt.maxIter - 1) ok = false;   // 未收敛
        }
        out.steps = istep;
    }
    out.trace.converged = ok;
    out.trace.totalIterations = out.totalIterations;
    if (!out.trace.steps.empty()) {
        out.trace.residualFirst = out.trace.steps.front().residual;
        out.trace.residualLast = out.trace.steps.back().residual;
    }
    out.converged = ok;
    if (!ok) throw std::runtime_error("非线性求解未收敛：达到最大迭代次数");

    // ---- 收敛构形的一致切线（可选导出）----
    // 成桥态几何迭代用它做低秩 Newton；影响矩阵 J = M·D 由
    // cableSensitivity3D 复用**这一次分解** + nc 次回代闭式给出
    // （P3-2b 已实现，不再是待办），不用为每根索扰动重跑一次完整非线性
    // 求解。J 冻结复用直到残差停滞才重建，重建次数见 [LENGTH] 的
    // SENS JREBUILD。
    if (tangent) {
        Eigen::VectorXd Fint;
        SparseMat3D Kt;
        assemble(Fint, Kt, true);
        Eigen::VectorXd dummy = Eigen::VectorXd::Zero(nd);
        std::vector<double> zeros(cdofs.size(), 0.0);
        applyDirichlet(Kt, dummy, cdofs, zeros);
        tangent->K = Kt;
        tangent->cdofs = cdofs;
    }

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

// ==================== 索体系的低秩敏感度算子 ====================
//
// 索承结构里，凡是"动索参数 → 看索力怎么变"的问题，敏感度都是**低秩**的：
//
//     B   (nd×nc)：列 k = 索 k 的轴向方向向量 [-e_k; +e_k]
//     Kc  (nc)   ：dT_k/dLn_k（Ernst 切线，含垂度软化）
//     D   (nc)   ：dT_k/dL0_k（固定弦长，纯标量关系）
//     W   = K⁻¹B （一次分解，nc 次回代）
//     M   = I - Kc·(BᵀW)      （nc×nc）
//
// 有了它，下面两件事都不再需要"每根索扰动一次、重跑一次完整非线性求解"：
//
//   · 索长迭代的影响矩阵   J = ∂T_i/∂L0_j = M·D
//       推导：扰动 L0_j → 索 j 内力变化 D_j·ΔL0_j → 不平衡力 -D_j·b_j·ΔL0_j
//             → du = K⁻¹(-D_j b_j ΔL0_j) → 索 i 索力变化 -Kc_i b_iᵀ K⁻¹ b_j D_j ΔL0_j
//             加上 i=j 的直接项 D_i·ΔL0_j，整理即得 J = (I - Kc BᵀK⁻¹B)·D = M·D
//
//   · 成桥态几何迭代       (I + G)⁻¹，G = -K⁻¹B Kc Bᵀ（Woodbury，见 P3-3）
//
// 两者共用同一次矩阵分解：nc 根索的算价从 O(nc) 次非线性求解降到
// 1 次分解 + nc 次回代。对百根以上拉索的实桥，这是数量级的差别。

struct CableSensitivity3D {
    Eigen::MatrixXd B;      // nd×nc
    Eigen::VectorXd Kc;     // nc：dT/dLn
    Eigen::VectorXd D;      // nc：dT/dL0
    Eigen::MatrixXd W;      // nd×nc：K⁻¹B
    Eigen::MatrixXd M;      // nc×nc：I - Kc·(BᵀW)
};

inline CableSensitivity3D cableSensitivity3D(const Model3D& m, const DofMap3D& dm,
                                             const std::vector<int>& cableIds,
                                             const std::vector<double>& u,
                                             const std::vector<double>& tension,
                                             const std::vector<char>& isT,
                                             const NLTangent3D& tan) {
    const int nd = dm.nDofs;
    const int nc = static_cast<int>(cableIds.size());
    CableSensitivity3D s;
    s.B = Eigen::MatrixXd::Zero(nd, nc);
    s.Kc = Eigen::VectorXd::Zero(nc);
    s.D = Eigen::VectorXd::Zero(nc);
    for (int k = 0; k < nc; ++k) {
        const Element3D* e = nullptr;
        for (const auto& q : m.elements) if (q.id == cableIds[k]) { e = &q; break; }
        if (!e) throw std::runtime_error("索敏感度：找不到单元 " + std::to_string(cableIds[k]));
        const int s1 = dm.dofStart[dm.index.at(e->n1)];
        const int s2 = dm.dofStart[dm.index.at(e->n2)];
        Vec3 p1 = m.nodes.at(e->n1), p2 = m.nodes.at(e->n2);
        for (int a = 0; a < 3; ++a) { p1[a] += u[s1 + a]; p2[a] += u[s2 + a]; }
        const Vec3 d3 = p2 - p1;
        const double Ln = d3.norm();
        if (Ln < 1e-12) throw std::runtime_error("索敏感度：索 "
                                                 + std::to_string(cableIds[k]) + " 弦长退化");
        const Vec3 ev = d3 / Ln;
        const double Lh = std::sqrt(d3.x() * d3.x() + d3.y() * d3.y());
        const double L0 = (e->L0 > 0.0) ? e->L0 : Ln;
        for (int a = 0; a < 3; ++a) {
            if (isT[s1 + a]) s.B(s1 + a, k) = -ev[a];
            if (isT[s2 + a]) s.B(s2 + a, k) = +ev[a];
        }
        const double EA = e->E * e->A;
        // 松弛索（T ≤ 0）不提供刚度，也不参与敏感度——否则会给 Newton 一个
        // "这根索还能拉"的错误方向。
        s.Kc(k) = (tension[k] > 0.0)
                      ? cableErnstModulus(e->E, e->A, e->w, Lh, tension[k], L0) * e->A / L0
                      : 0.0;
        // D = ∂T/∂L0（固定弦长）：只是标量关系，中心差分代价可忽略且不会推错
        const double h = 1e-7 * L0;
        s.D(k) = (cableTensionFromChord(Ln, L0 + h, EA, e->w, Lh)
                  - cableTensionFromChord(Ln, L0 - h, EA, e->w, Lh)) / (2.0 * h);
    }
    Factorization3D fac(tan.K);
    s.W = fac.solveMulti(s.B);
    for (int d = 0; d < nd; ++d) if (!isT[d]) s.W.row(d).setZero();
    s.M = Eigen::MatrixXd::Identity(nc, nc) - s.Kc.asDiagonal() * (s.B.transpose() * s.W);
    if (getenv("FEM_NL_DIAG")) {
        // λ_k = Kc_k·b_kᵀK⁻¹b_k ∈ (0,1)：→1 表示"索刚度远大于结构其余刚度"，
        // (I+G) 的对应特征值 1-λ → 0，Newton 步会被放大 1/(1-λ) 倍。
        double lmin = 1e30, lmax = -1e30;
        for (int k = 0; k < nc; ++k) {
            const double lam = s.Kc(k) * (s.B.transpose() * s.W)(k, k);
            lmin = std::min(lmin, lam); lmax = std::max(lmax, lam);
        }
        const Eigen::VectorXcd ev = s.M.eigenvalues();
        double emin = 1e30;
        for (int k = 0; k < ev.size(); ++k) emin = std::min(emin, std::abs(ev(k)));
        std::printf("  [diag] λ∈[%.4g, %.4g]，|M 最小特征值|=%.4g\n", lmin, lmax, emin);
    }
    return s;
}

// ==================== 成桥态：无应力索长迭代 ====================
//
// 成桥态的目标：在设计线形下恒载由索力平衡，各索达到目标索力 T*。
// 索是"做短一点"才张得紧的构件，因此控制量不是索力而是**无应力索长 L0**
// ——它同时也是施工下料的直接依据。
//
// 未知量与方程：每根索一个未知量 L0、一个目标 T*，故迭代对 T(L0) 求根：
//   ① 初值取"目标状态估计"：
//        L0 = (L* + s(T*)) / (1 + T*/(EA))，  s(T) = w²Lh³/(24T²) 为垂度缩短
//      （L* = 设计弦长；即先假设"结构上正好落在设计几何时索力为 T*"）
//   ② 完整非线性求解（含恒载 + 索自重），取出该索的索力 Tn；
//   ③ 影响矩阵 Newton 步：ΔL0 = -J⁻¹·(Tn - T*)，  J_ij = ∂T_i/∂L0_j
//      J = M·D 由 cableSensitivity3D 闭式给出：复用收敛构形的一次切线
//      分解 + nc 次回代，而不是 nc 次扰动重跑非线性求解（详见实现处注释）。
//      J 冻结复用，残差停滞才重建。
//      索与索之间通过结构位移耦合，若只用对角斜率，收敛率≈耦合系数，
//      实测会退化到数十轮；完整 J 才恢复二次收敛。
//
// 为什么必须对 T(L0) 求根，而不能"跟着弦长走"（L0 ← L0 + (Ln - L*)）：
//   弦长与索力在柔性体系里是**反号**的。以竖索吊悬臂梁为例：索力偏大 →
//   梁端被拉高 → 弦变短 → 上式继续缩短 L0 → 索力更大 → 发散。
//   而对 T(L0)，串联弹簧论证给出 dT/dL0 = -EA/(L0 + EA/k_s) < 0 恒成立
//   （k_s 为结构沿索轴方向的局部刚度），故 Newton 方向恒正确、必收敛。
//
// 收敛后：
//   · 索力 = 目标成桥索力 T*；
//   · 位移 u 是**物理上正确的成桥变形**：本模型以设计线形为参考构形（梁在
//     该构形下无内力），故恒载下梁必然变形，u ≠ 0。-u 是预抬高的**一阶估计**。
//     要真正让结构"恒载下正好落到设计线形"，用 P3-3 的 solveCompletionState3D：
//     它先在**设计弦长**上闭式反算 L0（本函数是在变形后的弦长上迭代，二者
//     差别正是成桥变形量），再解几何闭合方程 ΔX + u(X0 + ΔX) = 0 得制造线形。
//     直接拿 -u 当预抬高会过冲：移动节点改变索弦长 → 改变索力 → 又改变位移，
//     反馈环增益 = Kc/(kb + Kc)，索-梁刚度同量级时 ≈ 0.5，不动点迭代数十轮。
//     （成桥态梁内应有的弯矩由"无应力曲率"承担，属无应力状态法的后续扩展，
//       不影响索长迭代本身。）
//   · L0 可直接用于下料与正装分析。

// ==================== 迭代过程留痕（只记录，不改数值）====================
// 段格式见 cable_length_output_3d.hpp。这两个结构由索长迭代 / 制造线形
// 几何闭合在循环里顺手填，绝不参与判据、步长或收敛计算——"暴露过程量"
// 不允许改变既有数值口径。
struct IterationStep3D {
    int iter = 0;                  // 从 1 起
    double err = 0.0;              // 控制误差：索长迭代=max 相对索力偏差；几何闭合=max|ΔX+u| (m)
    double geomErr = 0.0;          // max|成桥弦长 − 设计弦长| (m)
    double stepMax = 0.0;          // 本轮最大修正量 |ΔL0| / |ΔX| (m)
    double stepRelMax = 0.0;       // stepMax / 对应参考长度
    double alpha = 1.0;            // 线搜索实际步长因子（索长迭代恒 1，靠 stepCap 限幅）
    int lineSearch = 0;            // 本轮线搜索评估次数
    bool rebuilt = false;          // 本轮是否重建影响矩阵 / 索敏感度
    bool capped = false;           // 修正是否被步长上限截断
    bool fallback = false;         // 是否退化为不动点步（Newton 方向失效）
    bool converged = false;        // 本轮命中收敛判据
};

struct IterationTrace3D {
    std::string kind;              // "cable_length" | "erection_geometry"
    std::string errKind;           // "tension_rel" | "closure_m"
    int maxIter = 0;
    double tol = 0.0;
    bool converged = false;
    std::vector<IterationStep3D> steps;
    std::vector<int> elems;
    std::vector<double> targetTension;   // N（仅索长迭代）
    Eigen::MatrixXd J;                   // 索长迭代最后一次影响矩阵 (nc×nc)
    bool jAvailable = false;
    int jRebuilds = 0;
    double sensDiagMin = 0.0, sensDiagMax = 0.0;  // 影响/敏感度矩阵对角范围
};

struct CableTarget3D {
    int elem = 0;
    double tension = 0.0;        // 目标成桥索力（N，拉为正）
};

struct CableLengthOptions3D {
    int maxIter = 30;
    double tolTension = 1e-6;    // 索力相对容差（唯一收敛判据）
    double stepCap = 0.1;        // 单轮 L0 步长上限（相对 L0），防野步长
    NonlinearOptions3D nl;       // 每轮内部使用的非线性求解选项
};

struct CableLengthResult3D {
    bool converged = false;
    int iterations = 0;
    std::vector<int> elems;
    std::vector<double> L0;          // 收敛后的无应力索长（施工下料用）
    std::vector<double> chord;       // 最终弦长
    std::vector<double> tension;     // 最终索力
    std::vector<double> Ltarget;     // 设计弦长（仅用于报告成桥变形量级）
    double maxGeomErr = 0.0;         // max|弦长 - 设计弦长|：成桥变形量（报告，不作判据）
    double maxTensionErr = 0.0;
    std::vector<double> u;           // 成桥态位移；其反号 = 需要的预抬高量
    IterationTrace3D trace;          // 逐轮过程留痕（只记录，不参与计算）
    NewtonTrace3D lastNewton;        // 收敛那一轮的内层 Newton 过程（透传给 CLI）
};

inline CableLengthResult3D iterateCableUnstrainedLengths(
        Model3D m,
        const std::vector<CableTarget3D>& targets,
        const CableLengthOptions3D& opt = {}) {
    CableLengthResult3D out;
    if (targets.empty()) return out;

    // ---- 索单元指针、设计弦长、水平投影、EA ----
    const std::size_t nc = targets.size();
    std::vector<Element3D*> ep(nc, nullptr);
    std::vector<double> Ltar(nc, 0.0), Lh(nc, 0.0), EA(nc, 0.0);
    for (std::size_t k = 0; k < nc; ++k) {
        Element3D* e = nullptr;
        for (auto& q : m.elements) if (q.id == targets[k].elem) { e = &q; break; }
        if (!e) throw std::runtime_error("索长迭代：找不到单元 " + std::to_string(targets[k].elem));
        if (targets[k].tension <= 0.0)
            throw std::runtime_error("索长迭代：目标索力必须为正（单元 "
                                     + std::to_string(targets[k].elem) + "）");
        ep[k] = e;
        const Vec3 p1 = m.nodes.at(e->n1), p2 = m.nodes.at(e->n2);
        const Vec3 d = p2 - p1;
        Ltar[k] = d.norm();
        Lh[k] = std::sqrt(d.x() * d.x() + d.y() * d.y());
        EA[k] = e->E * e->A;
        out.elems.push_back(targets[k].elem);
        out.trace.targetTension.push_back(targets[k].tension);
        // 初始 L0：目标状态估计（含垂度缩短 s(T*) = w²Lh³/(24T*²)）
        double sag = 0.0;
        if (e->w > 0.0 && Lh[k] > 1e-12)
            sag = e->w * e->w * Lh[k] * Lh[k] * Lh[k]
                  / (24.0 * targets[k].tension * targets[k].tension);
        e->L0 = (Ltar[k] + sag) / (1.0 + targets[k].tension / EA[k]);
    }
    out.Ltarget = Ltar;
    out.trace.kind = "cable_length";
    out.trace.errKind = "tension_rel";
    out.trace.maxIter = opt.maxIter;
    out.trace.tol = opt.tolTension;
    out.trace.elems = out.elems;

    // ---- 影响矩阵 Newton：对 T(L0) 求根 ----
    // 索与索之间通过结构位移耦合。若只用对角斜率，收敛率≈耦合系数（本仓库
    // 双索算例实测 0.7，40 轮才收敛到 1e-6），因此必须建立**完整影响矩阵**
    // J_ij = ∂T_i/∂L0_j。
    //
    // J 不再靠"每根索扰动一次 + 重跑一次完整非线性求解"去凑，而是闭式给出：
    //     J = M·D，  M = I - Kc·(BᵀK⁻¹B)，  D = diag(dT/dL0)
    // 推导见 cableSensitivity3D。代价 = 1 次矩阵分解 + nc 次回代，而不是
    // nc 次非线性求解（每含若干 Newton 迭代 + 若干次分解）。
    // 之后冻结复用 J——它随 L0 的变化很小，Newton 保持二次收敛；
    // 若某轮收敛停滞（残差降不到一半）就重建一次。
    const DofMap3D dm = buildDofMap3D(m);
    const int nd = dm.nDofs;
    std::vector<char> isC(nd, 0), isT(nd, 0);
    for (const auto& c : m.constraints) {
        if (c.k > 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0) isC[d] = 1;
    }
    for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
        const int cnt = std::min(dm.dofCount[i], 3);
        for (int a = 0; a < cnt; ++a) isT[dm.dofStart[i] + a] = 1;
    }
    for (int d = 0; d < nd; ++d) if (isC[d]) isT[d] = 0;

    std::vector<int> ids(nc, 0);
    for (std::size_t k = 0; k < nc; ++k) ids[k] = targets[k].elem;

    const int ncv = static_cast<int>(nc);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(ncv, ncv);
    Eigen::VectorXd res = Eigen::VectorXd::Zero(ncv);
    std::vector<double> chord(nc, 0.0), tens(nc, 0.0);
    bool needJ = true;
    double terrPrev = 0.0;
    for (int it = 0; it < opt.maxIter; ++it) {
        NLTangent3D tan;
        const NonlinearResult3D r =
            solveModel3DNonlinear(m, opt.nl, needJ ? &tan : nullptr);
        chord.assign(nc, 0.0);
        tens.assign(nc, 0.0);
        for (std::size_t k = 0; k < nc; ++k) {
            std::size_t idx = 0;
            while (idx < r.elemIds.size() && r.elemIds[idx] != targets[k].elem) ++idx;
            if (idx >= r.elemIds.size())
                throw std::runtime_error("索长迭代：结果缺少单元 "
                                         + std::to_string(targets[k].elem));
            chord[k] = r.Ln[idx];
            tens[k] = r.axial[idx];
        }
        out.iterations = it + 1;
        out.u = r.u;
        out.lastNewton = r.trace;

        IterationStep3D step;       // 逐轮留痕，不参与任何判据/步长计算
        step.iter = it + 1;
        double gerr = 0.0, terr = 0.0;
        for (std::size_t k = 0; k < nc; ++k) {
            gerr = std::max(gerr, std::abs(chord[k] - Ltar[k]));
            terr = std::max(terr, std::abs(tens[k] - targets[k].tension) / targets[k].tension);
        }
        out.chord = chord;
        out.tension = tens;
        out.maxGeomErr = gerr;
        out.maxTensionErr = terr;
        step.err = terr;
        step.geomErr = gerr;
        if (terr <= opt.tolTension) {
            step.converged = true;
            out.trace.steps.push_back(step);
            out.converged = true;
            break;
        }

        // 影响矩阵：J = M·D（一次分解 + nc 次回代）
        step.rebuilt = needJ;
        if (needJ) {
            const CableSensitivity3D s =
                cableSensitivity3D(m, dm, ids, r.u, tens, isT, tan);
            J = s.M * s.D.asDiagonal();
            out.trace.J = J;
            out.trace.jAvailable = true;
            ++out.trace.jRebuilds;
        }

        // Newton 步 + 步长保护（J 奇异时退化为对角步）
        for (std::size_t k = 0; k < nc; ++k) res(k) = tens[k] - targets[k].tension;
        const Eigen::VectorXd dL = J.fullPivLu().solve(-res);
        needJ = (it > 0 && terrPrev > 0.0 && terr > 0.5 * terrPrev);   // 停滞则重建
        terrPrev = terr;
        for (std::size_t k = 0; k < nc; ++k) {
            const int kk = static_cast<int>(k);
            double s = dL(kk);
            if (!std::isfinite(s))
                s = (J(kk, kk) < -1e-12) ? -res(kk) / J(kk, kk) : 0.0;
            const double cap = opt.stepCap * ep[k]->L0;
            if (s > cap) { s = cap; step.capped = true; }
            if (s < -cap) { s = -cap; step.capped = true; }
            step.stepMax = std::max(step.stepMax, std::abs(s));
            step.stepRelMax = std::max(step.stepRelMax,
                                       std::abs(s) / std::max(ep[k]->L0, 1e-12));
            ep[k]->L0 = std::max(ep[k]->L0 + s, 1e-6);
        }
        out.trace.steps.push_back(step);
    }

    out.trace.converged = out.converged;
    if (out.trace.jAvailable && out.trace.J.rows() > 0) {
        double dmin = out.trace.J(0, 0), dmax = dmin;
        for (int k = 1; k < out.trace.J.rows(); ++k) {
            dmin = std::min(dmin, out.trace.J(k, k));
            dmax = std::max(dmax, out.trace.J(k, k));
        }
        out.trace.sensDiagMin = dmin;
        out.trace.sensDiagMax = dmax;
    }

    // ---- 输出最终 L0 ----
    out.L0.clear();
    for (const auto& t : targets) {
        for (const auto& e : m.elements)
            if (e.id == t.elem) { out.L0.push_back(e.L0); break; }
    }
    return out;
}

// ==================== 成桥态闭环（P3-3）====================
//
// 成桥态要交付两件东西，也正是施工最需要的两件：
//   · 下料索长 L0 —— 工厂按此长度制索（无应力长度，可直接正装分析）；
//   · 制造/安装线形 Xm —— 现场按此线形拼梁挂索，与设计线形之差即预抬高。
//
// 关键结论：**这两个量是解耦的**，不必联立迭代。
//
//   L0 —— 成桥位置被定义为设计线形 X0，故成桥弦长必为设计弦长 L*，索力由
//         力-伸长关系单值决定，于是 L0 可**闭式反算**，一次都不用解结构：
//              L0 = (L* + s(T*)) / (1 + T*/(EA))，  s(T*) = w²Lh³/(24T*²)
//
//   Xm —— 梁/塔在制造线形下无应力（参考构形 = 制造线形），成桥时恒载 + 索力
//         使其变形到 X0，即几何闭合方程
//              F(ΔX) = ΔX + u(X0 + ΔX) = 0
//
// 这个拆分把"索力 ↔ 几何"的强反馈消掉了。常规做法先在一个**会变形**的构形上
// 调索力，再把位移反号当预抬高，二者互相打架：闭环增益 = Kc/(kb + Kc)，
// 索-梁刚度同量级时 ≈ 0.5，不动点迭代数十轮，单次代入必过冲（P3-2 实测）。
// 本方案 L0 一步到位，几何迭代只处理剩下的弱耦合，Newton 2~3 轮收敛。
//
// 几何闭合方程的 Jacobian 是良态的，但**条件数随索刚度增大而恶化**：
//      ∂u/∂ΔX = -K⁻¹K_索   ⇒   I + ∂u/∂ΔX = K⁻¹K_梁塔
// （K = K_梁塔 + K_索）。索越刚、梁塔相对越柔，K⁻¹K_梁塔 越小 → Newton 步被
// 放大 1/|λ_min| 倍（斜拉桥算例实测 λ_min ≈ 0.03，步长放大 30 倍）。因此
// 几何迭代**必须**带回溯线搜索：取整步会直接反弹发散。
//
// 这个"只计索"的 Jacobian 究竟够不够？实测（46 索斜拉桥、在成桥构形上逐列
// 对拍：扰动一个节点坐标、重解、与预测的 dF/dΔX 比较）：
//     主跨跨中梁竖向 1.1%、索锚点梁竖向 1.9%、塔顶竖向 0.4%
// 即已经准到 1~2%。早前注释里"G 的误差约 20%"是在**畸变索力**下测的，属于
// 另一个问题（索力确定法），不是 Jacobian 的锅。
// 残下这 1~2% 来自梁的**参考轴**项：移动节点同时改变参考长度 L0，轴力随之
// 变化，量级 N/L（与梁弯曲刚度 12EI/L³ 相比约 1%，恰好对上实测）。
// 它使收敛从二次退化为线性（每轮 ×0.12，46 索算例 10 轮收敛到 1e-9）。
// 补齐这一项要精确的 ∂F_int/∂X，而且一旦 K_X 变成满秩，(I+G) 不再是低秩形式、
// Woodbury 直接失效，必须改成每次分解 (K - K_X)——在 145 ms 的成桥态里只能
// 省下 4 轮，不划算。因此这里保留低秩 Newton + 线搜索 + 不动点兜底。
// （试过的弯路：把 K_X 猜成轴力几何刚度 (N/L)[[I,-I],[-I,I]]，结果比不加更差
//   70 倍——说明真修正项的符号/结构都不平凡，必须严格推导，不能拍脑袋。）
//
// 最后一条工程判断：预抬高量应当**小**（厘米～分米级）。若迭代给出的 ΔX 达到
// 米级，说明目标索力本身不合理（梁被强迫变形到设计线形会产生巨大内力），
// 应当回到索力的确定方法，而不是靠几何迭代硬凑——刚性支承连续梁法在"支座与
// 索锚点紧邻"时会给出畸变的反力（实测 min/max 相差 17 倍），这是它的已知局限。
//
// 几何迭代的敏感度是**低秩**的：共转梁单元对节点坐标的平移/转动严格不变
// （刚体不变性，test_3d_nl 已钉死到 1e-16），所以 ΔX 只通过索产生不平衡力：
//
//     δLn_k = b_kᵀΔX，  b_k = [-e_k; +e_k]
//     Δu    = -K⁻¹·Σ_k Kc_k·b_k·δLn_k
//     G     = ∂u/∂ΔX = -K⁻¹·B·Kc·Bᵀ      （秩 ≤ nc，nc = 索数）
//
// 于是 (I + G) 的 Newton 步可以用 Woodbury 精确求解，只需 nc 次**线性回代**
// （复用同一次分解）+ 一个 nc×nc 求逆，而不是 nd 维矩阵运算：
//
//     W = K⁻¹B，  M = I - Kc·(BᵀW)
//     δ = -F - W·M⁻¹·(Kc·(BᵀF))
//
// 对百根以上拉索的实桥，这正是把成桥态从"几十次非线性求解"降到"几次"的关键。
//
// 目标索力 T* 若未给定，用**最小弯曲能量法**自动确定（见 ②b）：令全桥弯曲应变
// 能最小，等价于"每根索两端锚点沿索向的相对位移为零"的变形协调方程。相比刚性
// 支承连续梁法（把每个锚点的竖向位移各自压成 0），它在锚点紧邻支座时不会给出
// 畸变反力，索力沿跨径平滑分布。

// ① 由设计弦长与目标索力闭式反算无应力（下料）索长。（定义见 fem_3d.hpp：
//    施工阶段调索也要用同一个式子，故放在公共头文件里。）

// ② 刚性支承连续梁法：索 → 锚点处刚性竖向支承，反力沿索向分解即得目标索力。
struct CableAnchor3D {
    int elem = 0;
    int anchorNode = 0;        // 梁侧锚点（在该节点加刚性竖向支承）
};

inline std::vector<double> cableTensionsRigidSupport(const Model3D& m,
                                                     const std::vector<CableAnchor3D>& cables,
                                                     int verticalDof = 3) {
    const std::size_t nc = cables.size();
    std::vector<double> T(nc, 0.0);
    if (nc == 0) return T;

    // 去掉索，改为锚点处的刚性竖向支承
    Model3D mc = m;
    std::vector<Element3D> keep;
    keep.reserve(mc.elements.size());
    for (const auto& e : mc.elements) {
        bool isTgt = false;
        for (const auto& c : cables) if (c.elem == e.id) { isTgt = true; break; }
        if (!isTgt) keep.push_back(e);
    }
    mc.elements = keep;
    for (const auto& c : cables)
        mc.constraints.push_back({c.anchorNode, verticalDof, 0.0, 0.0});
    constrainDisconnected3D(mc);      // 塔顶锚点等孤立节点：补齐约束，避免奇异

    DofMap3D dm = buildDofMap3D(mc);
    const int nd = dm.nDofs;
    const SparseMat3D K0 = assembleStiffness3D(mc, dm);     // 未施加 Dirichlet
    Eigen::VectorXd F0 = Eigen::VectorXd::Zero(nd);
    for (const auto& [node, dof, val] : mc.loads) {
        const int d = dm.safeDof(node, dof - 1);
        if (d >= 0) F0(d) += val;
    }
    std::vector<int> cdofs;
    std::vector<double> cvals;
    for (const auto& c : mc.constraints) {
        if (c.k > 0.0) continue;                    // 弹性支承不参与
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        auto it = std::find(cdofs.begin(), cdofs.end(), d);
        if (it == cdofs.end()) { cdofs.push_back(d); cvals.push_back(c.value); }
        else cvals[std::distance(cdofs.begin(), it)] = c.value;
    }
    SparseMat3D K = K0;
    Eigen::VectorXd F = F0;
    applyDirichlet(K, F, cdofs, cvals);
    const Eigen::VectorXd u = solveSparse3D(K, F);

    // 支承反力：K0·u = F0 + R  ⇒  R = K0·u - F0（沿坐标正向为正）
    const Eigen::VectorXd R = K0 * u - F0;

    for (std::size_t k = 0; k < nc; ++k) {
        const int d = dm.safeDof(cables[k].anchorNode, verticalDof - 1);
        if (d < 0) throw std::runtime_error("刚性支承连续梁法：锚点 "
                                            + std::to_string(cables[k].anchorNode) + " 无竖向自由度");
        const Element3D* e = nullptr;
        for (const auto& q : m.elements) if (q.id == cables[k].elem) { e = &q; break; }
        if (!e) throw std::runtime_error("刚性支承连续梁法：找不到索单元 "
                                         + std::to_string(cables[k].elem));
        const int other = (e->n1 == cables[k].anchorNode) ? e->n2 : e->n1;
        // 索对锚点的力沿"锚点 → 另一端"方向，其竖向分量 = T·ea.z
        const Vec3 ea = (m.nodes.at(other) - m.nodes.at(cables[k].anchorNode)).normalized();
        if (std::abs(ea.z()) < 1e-6)
            throw std::runtime_error("刚性支承连续梁法：索 "
                                     + std::to_string(cables[k].elem) + " 近水平，竖向分量退化");
        T[k] = R(d) / ea.z();
    }
    return T;
}

// ②b 最小弯曲能量法（斜拉桥成桥索力的经典优化方法，本内核默认）
//
// 目标函数：恒载 + 索力作用下全桥（主梁 + 桥塔）的弯曲应变能
//      U(T) = ∫ M²/(2EI) ds ,     T = 各索索力（作为赘余力）
//
// 关键实现技巧：**不必真的去积分弯矩图**。由虚功原理，U 对 T 的一、二阶导数
// 都能用裸结构（去掉全部索）的**柔度**表示：
//      f_i = 索 i 的单位力对（锚点 +e_i、另端 -e_i 的自平衡力系）
//      u_j = K⁻¹ f_j ,   u_0 = 裸结构在恒载下的位移
//      F_ij = f_i·u_j ,  g_i = f_i·u_0
//      ⇒  U(T) = U_0 + ½TᵀFT + Tᵀg ,   ∂U/∂T = 0  ⇒  F·T = -g
// 算价 = 1 次分解 + (nc+1) 次回代（与索长迭代共用同一套低秩思路），
// 对 100 根索的实桥也是毫秒级——而不是"nc 次完整求解 + 逐个积分弯矩"。
//
// 这个方程的力学含义很清楚：f_i·u = -(u_B - u_A)·e_i = -δ_i，故 F·T = -g 就是
//      δ_i(T) = 0   —— 每根索两端锚点**沿索向**的相对位移为零（索视为不可伸长）
// 即力法的变形协调方程。它和"刚性竖向支承法"的区别正在这里：
//   · 刚性支承法把锚点的**竖向位移**各自压成 0：锚点紧邻支座时该点竖向刚度
//     极大，反力按局部刚度比例分配，实测 min/max 相差 17 倍（畸变）；
//   · 最小弯曲能量法约束的是**沿索向的相对位移**（允许锚点随塔顶一起位移），
//     力由变形协调而非局部刚度决定，索力沿跨径平滑分布。
//
// 计入索自身伸长时（rigidCables = false），协调方程变为
//      δ_i(T) = T_i·L_i/(EA_i)   ⇒   F 的对角加 L_i/(EA_i)
// 等价于"把索换成刚度 EA/L 的弹性支承"，可用线性模型的弹性支承反力独立校核。
//
// 已知局限（工程上公认）：跨中短索（紧邻塔柱）可能算出零/负索力——那里梁本身
// 刚度大、不需要索。本实现用**主动集**把这类索固定在下限（默认 0）后重解其余
// 索，而不是让它们参与方程制造噪声。实桥若要求所有索都受力，应改用影响矩阵法
// 或直接指定索力。
struct CableForceOptions3D {
    bool rigidCables = true;     // true = 经典最小弯曲能量法（索不可伸长）
    double minTension = 0.0;     // 索力下限（绝对值）：低于此值的索钉在该值
    double minTensionFraction = 0.02;  // 保力下限 = 该比例 × 最大索力（短索不为 0）
    int maxActiveSet = 40;       // 主动集迭代上限
    bool checkStability = true;  // 检查"去掉索后结构是否仍稳定"（见下）
    bool verbose = false;
};

inline std::vector<double> cableTensionsMinBendingEnergy(const Model3D& m,
                                                         const std::vector<CableAnchor3D>& cables,
                                                         const CableForceOptions3D& opt = {}) {
    const std::size_t nc = cables.size();
    if (nc == 0) return {};

    // ---- 裸结构：去掉全部待定索；只由索连接的孤立节点补齐约束 ----
    Model3D mc = m;
    {
        std::vector<Element3D> keep;
        keep.reserve(mc.elements.size());
        for (const auto& e : mc.elements) {
            bool isTgt = false;
            for (const auto& c : cables) if (c.elem == e.id) { isTgt = true; break; }
            if (!isTgt) keep.push_back(e);
        }
        mc.elements = keep;
    }
    // 只由索连接的孤立节点（例如去索后的浮式塔顶）：补齐约束避免矩阵奇异。
    // **不能直接调 constrainDisconnected3D** —— 它会把有弹簧支承的节点也硬约束
    // 掉，于是"塔顶柔性"变成"塔顶刚性"，弹簧刚度彻底失效，索力退化成刚性支承
    // 法的结果（实测：本该 947 kN 的索被算成 2400 kN = P，恰好等于把塔当无限刚）。
    // 这里只补"没有弹簧支承"的那些自由度。
    {
        std::map<int, bool> connected;
        for (const auto& e : mc.elements) { connected[e.n1] = true; connected[e.n2] = true; }
        const DofMap3D dm0 = buildDofMap3D(mc);
        for (const auto& [id, x] : mc.nodes) {
            (void)x;
            if (connected[id]) continue;
            for (int k = 0; k < dm0.count(id); ++k) {
                bool hasSpring = false;
                for (const auto& c : mc.constraints)
                    if (c.node == id && c.dof == k + 1 && c.k > 0.0) { hasSpring = true; break; }
                if (!hasSpring) mc.constraints.push_back({id, k + 1, 0.0, 0.0});
            }
        }
    }

    DofMap3D dm = buildDofMap3D(mc);
    const int nd = dm.nDofs;
    SparseMat3D K = assembleStiffness3D(mc, dm);
    Eigen::VectorXd F0 = Eigen::VectorXd::Zero(nd);
    for (const auto& [node, dof, val] : mc.loads) {
        const int d = dm.safeDof(node, dof - 1);
        if (d >= 0) F0(d) += val;
    }
    // 弹性支承（k>0）：刚度已由 assembleStiffness3D 计入（见其"弹性支承留在
    // 刚度里"一段），这里**只补 k·value 的等效荷载项**。重复再加一次 k 会让
    // 弹簧刚度翻倍——实测「悬臂梁 + 塔顶弹簧」算例的索力被算成 1358 kN 而
    // 闭式解是 947 kN（误差 43%）。F0 的 k·value 项 assembleStiffness3D 不管，
    // 必须在这里补上（弹性支承带非零强制位移时才有）。
    for (const auto& c : mc.constraints) {
        if (c.k <= 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        F0(d) += c.k * c.value;
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

    // ---- 单位力对矩阵 Ld（nd×nc）：列 k = 锚点 +e、另端 -e ----
    const auto isC = [&](int d) {
        return std::find(cdofs.begin(), cdofs.end(), d) != cdofs.end();
    };
    Eigen::MatrixXd Ld = Eigen::MatrixXd::Zero(nd, static_cast<int>(nc));
    std::vector<double> chord(nc, 0.0), EA(nc, 0.0);
    for (std::size_t k = 0; k < nc; ++k) {
        const Element3D* e = nullptr;
        for (const auto& q : m.elements) if (q.id == cables[k].elem) { e = &q; break; }
        if (!e) throw std::runtime_error("最小弯曲能量法：找不到索单元 "
                                         + std::to_string(cables[k].elem));
        const int A = cables[k].anchorNode;
        const int B = (e->n1 == A) ? e->n2 : e->n1;
        if (A != e->n1 && A != e->n2)
            throw std::runtime_error("最小弯曲能量法：锚点 " + std::to_string(A)
                                     + " 不是索单元 " + std::to_string(e->id) + " 的端点");
        const Vec3 d3 = m.nodes.at(B) - m.nodes.at(A);
        const double L = d3.norm();
        if (L < 1e-12) throw std::runtime_error("最小弯曲能量法：索 "
                                                + std::to_string(e->id) + " 长度退化");
        chord[k] = L;
        EA[k] = e->E * e->A;
        const Vec3 ev = d3 / L;
        const int sA = dm.dofStart[dm.index.at(A)];
        const int sB = dm.dofStart[dm.index.at(B)];
        // 单位**拉力**对：锚点被拉向另一端（+e），另一端被拉向锚点（-e）。
        // 符号不能反：写成"压力对"时整个解会反号，索力被解成负值后被主动集
        // 钉在下限 0 上（第一版就是这样，T 全是 0）。
        for (int a = 0; a < 3; ++a) {
            const int dA = sA + a, dB = sB + a;
            if (!isC(dA)) Ld(dA, static_cast<int>(k)) = +ev[a];
            if (!isC(dB)) Ld(dB, static_cast<int>(k)) = -ev[a];
        }
    }

    // ---- 柔度矩阵：一次分解 + (nc+1) 次回代 ----
    Factorization3D fac(K);

    // 裸结构稳定性自检：最小弯曲能量法的前提是"去掉索后结构仍稳定"。
    // 这在索承体系里并非自动成立——若主梁的竖向支承全由索提供（纯索承、
    // 无支座），去索后梁在竖向是**机构**，柔度矩阵毫无意义，解出来的索力是
    // 舍入噪声（实测：三索算例的端索被解成 0）。用"随机向量 → K⁻¹ → 回代"
    // 检验分解是否为真逆：机构会沿零能模式放大到天文数字。
    if (opt.checkStability && nd > 0) {
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
        if (!(eworst < 1e-4))
            throw std::runtime_error(
                "最小弯曲能量法：去掉索后结构接近机构（回代相对误差 "
                + std::to_string(eworst) + "），柔度矩阵不可用。"
                "该体系的主梁竖向支承由索提供（纯索承），请改用刚性支承连续梁法"
                "（CompletionStateOptions3D::CableForceMethod::RigidSupport）或直接指定索力。");
    }

    const Eigen::VectorXd u0 = fac.solve(Fdl);
    const Eigen::MatrixXd U = fac.solveMulti(Ld);
    Eigen::MatrixXd Fm = Ld.transpose() * U;
    Fm = 0.5 * (Fm + Fm.transpose());                  // 消掉舍入造成的非对称
    Eigen::VectorXd g = Ld.transpose() * u0;
    if (!opt.rigidCables)
        for (std::size_t k = 0; k < nc; ++k) Fm(static_cast<int>(k), static_cast<int>(k))
                                                += chord[k] / EA[k];

    // 锚点自由度全被约束时该索无法定力（柔度恒为 0）——直接报错，而不是给个伪解
    {
        double dmax = 0.0;
        for (std::size_t k = 0; k < nc; ++k) dmax = std::max(dmax, std::abs(Fm(static_cast<int>(k), static_cast<int>(k))));
        for (std::size_t k = 0; k < nc; ++k) {
            if (std::abs(Fm(static_cast<int>(k), static_cast<int>(k))) <= 1e-14 * std::max(dmax, 1e-30))
                throw std::runtime_error("最小弯曲能量法：索 " + std::to_string(cables[k].elem)
                                         + " 的锚点自由度被完全约束，无法由变形协调定力");
        }
    }

    // ---- 主动集：解 F·T = -g；低于下限的索钉住后重解其余 ----
    //
    // 下限的来历：跨中短索（紧邻塔柱）在最小弯曲能量下可能算出 0 或**负**索力
    // ——那里梁本身刚度大、确实"不需要"这根索。但实桥不许索松弛（风振、疲劳、
    // 索力测点、换索），工程上统一按"短索保力"处理：给一个下限，把这根索从
    // 方程中移出、力固定在下限。下限取**最大索力的固定比例**（默认 2%），
    // 这样与桥的规模、单位无关。
    auto runActiveSet = [&](double floor, std::vector<double>& Tv) {
        Tv.assign(nc, 0.0);
        std::vector<char> fx(nc, 0);
        int npin = 0;
        for (int pass = 0; pass < opt.maxActiveSet; ++pass) {
            std::vector<int> freeIdx;
            for (std::size_t k = 0; k < nc; ++k) if (!fx[k]) freeIdx.push_back(static_cast<int>(k));
            if (freeIdx.empty()) break;
            const int nf = static_cast<int>(freeIdx.size());
            Eigen::MatrixXd Fr(nf, nf);
            Eigen::VectorXd gr(nf);
            for (int i = 0; i < nf; ++i) {
                const int ki = freeIdx[i];
                double bi = -g(ki);
                for (std::size_t k = 0; k < nc; ++k)
                    if (fx[k]) bi -= Fm(ki, static_cast<int>(k)) * Tv[k];
                gr(i) = bi;
                for (int j = 0; j < nf; ++j) Fr(i, j) = Fm(ki, freeIdx[j]);
            }
            const Eigen::VectorXd Ts = Fr.ldlt().solve(gr);
            for (int i = 0; i < nf; ++i) Tv[freeIdx[i]] = Ts(i);
            int added = 0;
            for (int i = 0; i < nf; ++i) {
                const int ki = freeIdx[i];
                if (!(Tv[ki] > floor)) { Tv[ki] = floor; fx[ki] = 1; ++added; }
            }
            npin += added;
            if (added == 0) break;
        }
        return npin;
    };

    std::vector<double> T;
    int pinned = runActiveSet(opt.minTension, T);
    if (opt.minTensionFraction > 0.0) {
        double Tmax0 = 0.0;
        for (double t : T) Tmax0 = std::max(Tmax0, t);
        const double floor = std::max(opt.minTension, opt.minTensionFraction * Tmax0);
        if (floor > opt.minTension) {
            // 只需在"确有索低于该下限"时才重跑
            bool need = false;
            for (double t : T) if (t < floor) { need = true; break; }
            if (need) pinned = runActiveSet(floor, T);
        }
    }
    if (opt.verbose) {
        double tmin = 1e30, tmax = -1e30, sum = 0.0;
        for (std::size_t k = 0; k < nc; ++k) {
            tmin = std::min(tmin, T[k]);
            tmax = std::max(tmax, T[k]);
            sum += T[k];
        }
        std::printf("  [最小弯曲能量] T ∈ [%.4g, %.4g] kN，ΣT = %.4g kN，"
                    "保力下限钉住 %d 根（下限 %.4g kN）\n",
                    tmin / 1e3, tmax / 1e3, sum / 1e3, pinned, tmin / 1e3);
    }
    return T;
}

// ③ 制造线形（预抬高）迭代：解几何闭合 ΔX + u(X0 + ΔX) = 0。
struct ErectionGeometryOptions3D {
    int maxIter = 20;
    double tolClosure = 1e-9;    // 几何闭合容差（m）
    double stepCap = 2.0;        // 单轮节点坐标变化上限（m），防野步长
    int maxLineSearch = 8;       // 回溯线搜索次数
    double lineSearchFactor = 0.4;   // 每次回溯的衰减系数
    // 每轮刷新索敏感度（W=K⁻¹B 与 Kc）。关掉能省一次分解，但 W 被冻在误差最大
    // 的那一轮上，收敛会明显变慢（实测 8 轮 → 10 轮，且每轮收敛率更差）。
    bool refreshSensitivity = true;
    bool verbose = false;
    NonlinearOptions3D nl;
};

struct ErectionGeometryResult3D {
    bool converged = false;
    int iterations = 0;
    std::vector<int> cableIds;
    std::map<int, Vec3> erected;    // 制造/安装线形节点坐标
    std::vector<double> camber;     // 预抬高 = erected - design（按自由度向量）
    std::vector<double> u;          // 成桥位移（相对制造线形）
    std::vector<double> tension;
    std::vector<double> chord;
    double maxClosureErr = 0.0;
    IterationTrace3D trace;          // 逐轮过程留痕（只记录，不参与计算）
    NewtonTrace3D lastNewton;        // 最后一次 evalGeom 的内层 Newton 过程
};

inline ErectionGeometryResult3D iterateErectionGeometry(Model3D m,
                                                        const std::vector<int>& cableIds,
                                                        const ErectionGeometryOptions3D& opt = {}) {
    ErectionGeometryResult3D out;
    if (cableIds.empty()) throw std::runtime_error("制造线形迭代：未指定索单元");
    out.cableIds = cableIds;
    out.trace.kind = "erection_geometry";
    out.trace.errKind = "closure_m";
    out.trace.maxIter = opt.maxIter;
    out.trace.tol = opt.tolClosure;
    out.trace.elems = cableIds;

    const DofMap3D dm = buildDofMap3D(m);
    const int nd = dm.nDofs;
    const std::map<int, Vec3> X0 = m.nodes;             // 设计线形

    // 约束自由度（不可动）与平动自由度（几何闭合只针对位置，不含转角）：
    // 制造线形是节点坐标，梁端转角由变形自行决定，不参与闭合方程。
    // **弹性支承（k>0）的自由度必须排除**：弹簧的平衡本身决定 u = T/k ≠ 0，
    // 这个沉降是物理的（基础/支座沉降），不能靠预抬高"抬掉"——若把它写进闭合
    // 方程，迭代会不断要求 u→0，而弹簧力只认位移，最终把索力压到零（实测）。
    std::vector<char> isC(nd, 0), isT(nd, 0), isSpring(nd, 0);
    for (const auto& c : m.constraints) {
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        if (c.k > 0.0) { isSpring[d] = 1; continue; }
        isC[d] = 1;
    }
    for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
        const int cnt = std::min(dm.dofCount[i], 3);
        for (int a = 0; a < cnt; ++a) isT[dm.dofStart[i] + a] = 1;
    }
    for (int d = 0; d < nd; ++d) if (isC[d] || isSpring[d]) isT[d] = 0;

    const std::size_t nc = cableIds.size();
    std::vector<const Element3D*> ep(nc, nullptr);
    for (std::size_t k = 0; k < nc; ++k) {
        for (const auto& q : m.elements) if (q.id == cableIds[k]) { ep[k] = &q; break; }
        if (!ep[k]) throw std::runtime_error("制造线形迭代：找不到索单元 "
                                             + std::to_string(cableIds[k]));
    }
    // 设计弦长（只用于过程留痕里的"成桥弦长偏差"，不参与闭合方程）
    std::vector<double> Ltar(nc, 0.0);
    for (std::size_t k = 0; k < nc; ++k)
        Ltar[k] = (X0.at(ep[k]->n2) - X0.at(ep[k]->n1)).norm();

    // 一次"应用几何 → 非线性求解 → 取闭合残差"
    struct GeomEval {
        double err = 0.0;
        std::vector<double> u, tension, chord;
        Eigen::VectorXd Fv;
        NewtonTrace3D newton;       // 本次评估的内层 Newton 过程（只留痕）
    };
    auto evalGeom = [&](const std::vector<double>& dx, NLTangent3D* tan) {
        // 把 ΔX 落到节点坐标上（参考构形随之改变 = 制造线形无应力）
        for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
            const int s = dm.dofStart[i];
            Vec3 p = X0.at(dm.nodeIds[i]);
            const int cnt = std::min(dm.dofCount[i], 3);
            for (int a = 0; a < cnt; ++a) p[a] += dx[s + a];
            m.nodes[dm.nodeIds[i]] = p;
        }
        const NonlinearResult3D r = solveModel3DNonlinear(m, opt.nl, tan);
        GeomEval g;
        g.newton = r.trace;         // 每次评估自带内层过程，由主循环挑选收敛那轮
        g.u = r.u;
        g.chord.assign(nc, 0.0);
        g.tension.assign(nc, 0.0);
        for (std::size_t k = 0; k < nc; ++k) {
            std::size_t idx = 0;
            while (idx < r.elemIds.size() && r.elemIds[idx] != cableIds[k]) ++idx;
            if (idx >= r.elemIds.size())
                throw std::runtime_error("制造线形迭代：结果缺少索单元 "
                                         + std::to_string(cableIds[k]));
            g.chord[k] = r.Ln[idx];
            g.tension[k] = r.axial[idx];
        }
        g.Fv = Eigen::VectorXd::Zero(nd);
        double err = 0.0;
        for (int d = 0; d < nd; ++d) {
            if (!isT[d]) continue;
            g.Fv(d) = dx[d] + r.u[d];
            err = std::max(err, std::abs(g.Fv(d)));
        }
        g.err = err;
        return g;
    };

    std::vector<double> dX(nd, 0.0);                    // 制造线形相对设计线形的偏移
    CableSensitivity3D sens;                            // 冻结复用（G 随 ΔX 变化很小）
    bool haveS = false;
    double errPrev = 0.0;
    GeomEval cur;

    for (int it = 0; it < opt.maxIter; ++it) {
        const std::vector<double> dXprev = dX;   // 只用于过程留痕的步长统计
        NLTangent3D tan;
        const bool needS = !haveS || opt.refreshSensitivity ||
                       (it > 0 && errPrev > 0.0 && cur.err > 0.5 * errPrev);
        cur = evalGeom(dX, needS ? &tan : nullptr);
        out.lastNewton = cur.newton;   // 收敛轮（或末轮）的内层 Newton 过程
        out.iterations = it + 1;
        out.u = cur.u;
        out.tension = cur.tension;
        out.chord = cur.chord;
        out.maxClosureErr = cur.err;

        IterationStep3D step;       // 逐轮留痕，不参与任何判据/步长计算
        step.iter = it + 1;
        step.err = cur.err;
        for (std::size_t k = 0; k < nc; ++k)
            step.geomErr = std::max(step.geomErr, std::abs(cur.chord[k] - Ltar[k]));
        step.rebuilt = needS;
        if (opt.verbose)
            std::printf("  几何迭代 %d：闭合误差 = %.3e m\n", it + 1, cur.err);
        if (cur.err <= opt.tolClosure) {
            step.converged = true;
            out.trace.steps.push_back(step);
            out.converged = true;
            break;
        }
        if (needS) {
            sens = cableSensitivity3D(m, dm, cableIds, cur.u, cur.tension, isT, tan);
            haveS = true;
            ++out.trace.jRebuilds;
            if (sens.M.rows() > 0) {
                double dmin = sens.M(0, 0), dmax = dmin;
                for (int q = 1; q < sens.M.rows(); ++q) {
                    dmin = std::min(dmin, sens.M(q, q));
                    dmax = std::max(dmax, sens.M(q, q));
                }
                out.trace.sensDiagMin = dmin;
                out.trace.sensDiagMax = dmax;
            }
        }
        errPrev = cur.err;

        // ---- 低秩 Newton（Woodbury）：共用索敏感度算子 ----
        const Eigen::VectorXd y =
            sens.M.fullPivLu().solve(sens.Kc.asDiagonal() * (sens.B.transpose() * cur.Fv));
        Eigen::VectorXd delta = -cur.Fv - sens.W * y;
        for (int d = 0; d < nd; ++d) if (!isT[d]) delta(d) = 0.0;

        // ---- 回溯线搜索 ----
        // Newton 步在远离解时会严重超调：(I+G) 的最小特征值实测可达 0.03，
        // 步长被放大 30 倍，直接取整步会让误差反弹震荡（实桥算例实测）。
        // Armijo 条件按**实际**闭合误差衰减步长，是这里唯一可靠的判据。
        double mx = 0.0;
        for (std::size_t i = 0; i < dm.nodeIds.size(); ++i) {
            const int s = dm.dofStart[i];
            mx = std::max(mx, std::sqrt(delta(s) * delta(s) + delta(s + 1) * delta(s + 1)
                                        + delta(s + 2) * delta(s + 2)));
        }
        step.capped = (mx > opt.stepCap);
        double alpha = (mx > opt.stepCap) ? opt.stepCap / mx : 1.0;
        bool accepted = false;
        GeomEval bestG = cur;
        std::vector<double> bestDx = dX;
        double bestErr = cur.err;
        for (int ls = 0; ls < opt.maxLineSearch; ++ls) {
            ++step.lineSearch;
            std::vector<double> trial(nd, 0.0);
            for (int d = 0; d < nd; ++d) trial[d] = dX[d] + alpha * delta(d);
            const GeomEval g = evalGeom(trial, nullptr);
            if (g.err < bestErr) { bestErr = g.err; bestDx = trial; bestG = g; }
            if (opt.verbose)
                std::printf("      线搜索 α=%.4g → %.3e m\n", alpha, g.err);
            if (g.err <= opt.tolClosure || g.err <= (1.0 - 1.0e-4 * alpha) * cur.err) {
                dX = trial;
                cur = g;
                accepted = true;
                step.alpha = alpha;
                break;
            }
            alpha *= opt.lineSearchFactor;
        }
        if (!accepted) {
            step.fallback = true;
            step.alpha = 1.0;       // 不动点步是全量替换，不是带 α 的增量
            // Newton 方向失效时的兜底：退化为不动点步 ΔX ← -u。
            //
            // 上面的 G 只计入了索（ΔX 只通过索的弦长产生不平衡力），而梁的
            // 参考形状同样参与：移动锚点会改变梁的参考轴方向，在"空间增量
            // 约定"下节点三元组跟着转、局部转角因此改变一阶量——实测扰动无
            // 索节点也能引起位移变化，且量级与索项相当（斜拉桥算例里 G 的
            // 误差约 20%）。完整计入需构造 K_shape = ∂F_int/∂X|_{u fixed}
            // 并令 G_full = -K⁻¹K_shape；K_shape 是满秩 nd×nd 算子，
            // 不能再走 Woodbury 低秩公式。实现规格与验收 oracle 见
            // USER_MANUAL.md §9.8 / CHANGELOG 2026-11-09（TODO(P3-3b)）。
            // 2026-11-09 用真实斜拉桥 COMPLETION 算例复核：6 轮闭合误差
            // 0.158 → 1.4e-9 m、逐轮收敛比约 0.025，每轮都是 Newton 步
            // 被线搜索接受、零次 FALLBACK，故暂不改默认数值路径。
            //
            // 不动点步 ΔX ← -u 的收敛率 = |G| 的谱半径 < 1（结构稳定时），
            // 不需要任何敏感度，是这里唯一的收敛保证。
            std::vector<double> fp(nd, 0.0);
            for (int d = 0; d < nd; ++d) if (isT[d]) fp[d] = -cur.u[d];
            const GeomEval g = evalGeom(fp, nullptr);
            if (opt.verbose)
                std::printf("      Newton 步失效 → 不动点步 → %.3e m\n", g.err);
            if (g.err < bestErr) { dX = fp; cur = g; }
            else { dX = bestDx; cur = bestG; }   // 取试过的最优点
        }
        double smax = 0.0;
        for (int d = 0; d < nd; ++d)
            smax = std::max(smax, std::abs(dX[d] - dXprev[d]));
        step.stepMax = smax;
        // 几何闭合没有单一参考长度，stepRelMax 留给 0，Python 侧不展示该列。
        out.trace.steps.push_back(step);
    }
    out.trace.converged = out.converged;

    out.camber = dX;
    out.erected = m.nodes;
    return out;
}

// ④ 组合：一次给出下料索长 + 制造线形（施工要的两件套）。
struct CompletionCable3D {
    int elem = 0;
    int anchorNode = 0;       // 梁侧锚点（自动定索力时用）
    double tension = 0.0;     // 目标成桥索力；<=0 表示用刚性支承连续梁法自动确定
};

// ---- 目标索力：给定则用；否则按所选方法自动确定 ----
struct CompletionStateOptions3D {
    int maxIter = 20;
    double tolClosure = 1e-9;    // 几何闭合容差（m）
    double stepCap = 2.0;
    bool verbose = false;
    // 自动定索力的方法：默认最小弯曲能量法（见 ②b）。
    // RigidSupport 保留作对照——它在锚点紧邻支座时会给出畸变索力，只适合
    // 无塔、竖索、锚点远离支座的简单体系。
    enum class CableForceMethod { MinBendingEnergy, RigidSupport };
    CableForceMethod cableForce = CableForceMethod::MinBendingEnergy;
    CableForceOptions3D minEnergy;      // 最小弯曲能量法的选项
    NonlinearOptions3D nl;
};

struct CompletionStateResult3D {
    bool converged = false;
    int iterations = 0;
    std::vector<int> elems;
    std::vector<double> L0;          // 下料索长（无应力长度）
    std::vector<double> tension;     // 成桥索力
    std::vector<double> chord;       // 成桥弦长（应回到设计弦长）
    std::vector<double> Ltarget;     // 设计弦长
    std::map<int, Vec3> erected;     // 制造/安装线形
    std::vector<double> camber;      // 预抬高 = erected - design（按自由度向量）
    std::vector<double> u;           // 成桥位移（相对制造线形）
    double maxClosureErr = 0.0;
    double maxTensionErr = 0.0;
    IterationTrace3D trace;          // 透传制造线形几何闭合的逐轮过程
    NewtonTrace3D lastNewton;        // 透传收敛轮的内层 Newton 过程
};

inline CompletionStateResult3D solveCompletionState3D(Model3D m,
                                                      const std::vector<CompletionCable3D>& cables,
                                                      const CompletionStateOptions3D& opt = {}) {
    CompletionStateResult3D out;
    if (cables.empty()) throw std::runtime_error("成桥态求解：未指定索");
    const std::size_t nc = cables.size();
    out.elems.resize(nc);
    for (std::size_t k = 0; k < nc; ++k) out.elems[k] = cables[k].elem;

    // ---- 目标索力：给定则用；否则按所选方法自动定 ----
    std::vector<double> Tst(nc, 0.0);
    std::vector<int> needAuto;
    for (std::size_t k = 0; k < nc; ++k) {
        if (cables[k].tension > 0.0) Tst[k] = cables[k].tension;
        else needAuto.push_back(static_cast<int>(k));
    }
    if (!needAuto.empty()) {
        std::vector<CableAnchor3D> ca;
        for (int k : needAuto) ca.push_back({cables[k].elem, cables[k].anchorNode});
        std::vector<double> Ta;
        if (opt.cableForce == CompletionStateOptions3D::CableForceMethod::RigidSupport) {
            Ta = cableTensionsRigidSupport(m, ca);
        } else {
            CableForceOptions3D fo = opt.minEnergy;
            fo.verbose = opt.minEnergy.verbose || opt.verbose;
            Ta = cableTensionsMinBendingEnergy(m, ca, fo);
        }
        for (std::size_t j = 0; j < needAuto.size(); ++j) {
            Tst[needAuto[j]] = Ta[j];
            // 索力非正 ⇒ 下料索长反算无解（索做多长都张不紧）。必须报错而不是
            // 塞一个数进去：这意味着该索在此恒载下不参与受力（短索紧邻塔柱时
            // 常见），应当指定索力或用影响矩阵法。
            if (Ta[j] <= 0.0)
                throw std::runtime_error(
                    "成桥态求解：自动定出的索力非正（索 " + std::to_string(cables[needAuto[j]].elem)
                    + " = " + std::to_string(Ta[j] / 1e3) + " kN）。"
                    "请指定该索的目标索力，或改用影响矩阵法/最小弯曲能量法加索力下限。");
        }
    }

    // ---- 下料索长：由设计弦长 + 目标索力闭式反算（无需任何结构迭代）----
    std::vector<double> Ltar(nc, 0.0);
    for (std::size_t k = 0; k < nc; ++k) {
        Element3D* e = nullptr;
        for (auto& q : m.elements) if (q.id == cables[k].elem) { e = &q; break; }
        if (!e) throw std::runtime_error("成桥态求解：找不到索单元 "
                                         + std::to_string(cables[k].elem));
        const Vec3 d3 = m.nodes.at(e->n2) - m.nodes.at(e->n1);
        Ltar[k] = d3.norm();
        const double Lh = std::sqrt(d3.x() * d3.x() + d3.y() * d3.y());
        e->L0 = cableUnstrainedLengthForTension(Ltar[k], Tst[k], e->E * e->A, e->w, Lh);
    }
    out.Ltarget = Ltar;

    // ---- 制造线形：几何闭合迭代 ----
    std::vector<int> ids(nc);
    for (std::size_t k = 0; k < nc; ++k) ids[k] = cables[k].elem;
    ErectionGeometryOptions3D eo;
    eo.maxIter = opt.maxIter;
    eo.tolClosure = opt.tolClosure;
    eo.stepCap = opt.stepCap;
    eo.verbose = opt.verbose;
    eo.nl = opt.nl;
    const ErectionGeometryResult3D eg = iterateErectionGeometry(m, ids, eo);

    out.trace = eg.trace;              // 制造线形几何闭合的逐轮过程
    out.lastNewton = eg.lastNewton;    // 收敛轮的内层 Newton 过程
    out.trace.targetTension = Tst;     // 给定或自动定出的目标成桥索力
    out.iterations = eg.iterations;
    out.converged = eg.converged;
    out.maxClosureErr = eg.maxClosureErr;
    out.tension = eg.tension;
    out.chord = eg.chord;
    out.camber = eg.camber;
    out.erected = eg.erected;
    out.u = eg.u;
    out.L0.assign(nc, 0.0);
    for (std::size_t k = 0; k < nc; ++k) {
        for (const auto& e : m.elements)
            if (e.id == cables[k].elem) { out.L0[k] = e.L0; break; }
        out.maxTensionErr = std::max(out.maxTensionErr,
                                     std::abs(eg.tension[k] - Tst[k]) / Tst[k]);
    }
    return out;
}

// 预应力等效荷载，公式与线性阶段 solveModel3DStages 相同：
// 直线束 = 两端轴力 Pe·c 与偏心弯矩 Pe·e·(ẑ×c)；
// 曲线束再加两端同向竖向力 4Pe·s/L 和均布 q=−8Pe·s/L² 的节点化。
// 非线性阶段每次都重算全量平衡，这些荷载必须在张拉之后一直保留。
inline std::vector<std::tuple<int, int, double>> prestressEquivalentLoads(
        const Model3D& m, const std::vector<Element3D>& active,
        const Model3D::Prestress3D& ps, const std::string& ptag) {
    auto elemIn = [&](int id) -> const Element3D* {
        for (const auto& e : active) if (e.id == id) return &e;
        return nullptr;
    };
    if (ps.elems.empty())
        throw std::runtime_error(ptag + " 没有 ELEM 列表");
    if (ps.P <= 0.0)
        throw std::runtime_error(ptag + " 的张拉力必须为正");
    if (ps.loss <= 0.0 || ps.loss > 1.0)
        throw std::runtime_error(ptag + " 的 LOSS 必须在 (0,1]");
    std::map<int, int> deg;
    for (int eid : ps.elems) {
        const Element3D* e = elemIn(eid);
        if (!e)
            throw std::runtime_error(ptag + " 引用了本阶段未激活的单元 "
                                     + std::to_string(eid));
        if (e->type != "BEAM" && e->type != "BEAM7")
            throw std::runtime_error(ptag + " 只能施加到梁单元");
        deg[e->n1]++; deg[e->n2]++;
    }
    int head = 0, tail = 0;
    for (const auto& [nd, c] : deg)
        if (c == 1) { if (head == 0) head = nd; else tail = nd; }
    if (head == 0 || tail == 0 || deg.size() != ps.elems.size() + 1)
        throw std::runtime_error(ptag + " 的单元不构成一条首尾相接的开链");
    const Vec3 dv = m.nodes.at(tail) - m.nodes.at(head);
    const double L = dv.norm();
    if (L <= 1e-12)
        throw std::runtime_error(ptag + " 的锚固端节点重合");
    const double ax = dv.x() / L, ay = dv.y() / L, az = dv.z() / L;
    const double Pe = ps.P * ps.loss;
    std::vector<std::tuple<int, int, double>> loads;
    const double mx = Pe * ps.e * (-ay);
    const double my = Pe * ps.e * ax;
    loads.push_back({head, 1, Pe * ax});
    loads.push_back({head, 2, Pe * ay});
    loads.push_back({head, 3, Pe * az});
    loads.push_back({head, 4, mx});
    loads.push_back({head, 5, my});
    loads.push_back({tail, 1, -Pe * ax});
    loads.push_back({tail, 2, -Pe * ay});
    loads.push_back({tail, 3, -Pe * az});
    loads.push_back({tail, 4, -mx});
    loads.push_back({tail, 5, -my});
    if (ps.s == 0.0) return loads;

    std::map<int, std::vector<int>> adj;
    for (int eid : ps.elems) {
        const Element3D* e = elemIn(eid);
        adj[e->n1].push_back(e->n2);
        adj[e->n2].push_back(e->n1);
    }
    std::vector<int> chain{head};
    int cur = head, prev = -1;
    while (cur != tail) {
        int nxt = -1;
        for (int nb : adj[cur]) if (nb != prev) { nxt = nb; break; }
        if (nxt < 0) throw std::runtime_error(ptag + " 的单元链不连通");
        chain.push_back(nxt);
        prev = cur;
        cur = nxt;
    }
    auto dist = [&](int a, int b) {
        return (m.nodes.at(b) - m.nodes.at(a)).norm();
    };
    for (int nd : chain) {
        if (nd == head || nd == tail) continue;
        const Vec3 w = m.nodes.at(nd) - m.nodes.at(head);
        const double cx = w.y() * dv.z() - w.z() * dv.y();
        const double cy = w.z() * dv.x() - w.x() * dv.z();
        const double cz = w.x() * dv.y() - w.y() * dv.x();
        const double perp = std::sqrt(cx * cx + cy * cy + cz * cz) / L;
        if (perp > 1e-6 * L)
            throw std::runtime_error(ptag + " 是曲线筋，单元必须共线");
    }
    const double tz = 4.0 * Pe * ps.s / L;
    const double qz = -8.0 * Pe * ps.s / (L * L);
    loads.push_back({head, 3, tz});
    loads.push_back({tail, 3, tz});
    for (std::size_t k = 0; k < chain.size(); ++k) {
        double half = 0.0;
        if (k == 0) half = dist(chain[0], chain[1]) / 2.0;
        else if (k + 1 == chain.size()) half = dist(chain[k - 1], chain[k]) / 2.0;
        else half = (dist(chain[k - 1], chain[k]) + dist(chain[k], chain[k + 1])) / 2.0;
        loads.push_back({chain[k], 3, qz * half});
    }
    return loads;
}

// 非线性施工阶段。每个阶段是当前已激活结构上的**全量平衡**，不是线性内核的
// 荷载增量。位移不能再按阶段相加。TENSION= 用无应力索长迭代命中目标力；
// 没有 TENSION= 的阶段只做一次共转求解，已有索的 L0 保持不变。
// 预应力等效荷载从张拉阶段起一直保留在外荷载里。
// 与线性内核 compute3DBeamEnvelope 同一套截面内力约定：N 拉为正，My 下侧受拉为正。
struct BeamEnd3D {
    double N = 0, Vy = 0, Vz = 0, T = 0, My = 0, Mz = 0;
};

struct NLStageResult3D {
    int id = 0;
    bool converged = false;
    NewtonTrace3D newton;   // 本阶段全量平衡的内层 Newton 过程（每阶段一段）
    DofMap3D dm;
    std::vector<double> u;
    std::map<int, double> cableForce;
    std::map<int, double> cableL0;
    // 键是单元号，值是 i 端、j 端。只含本阶段激活的梁。
    std::map<int, std::pair<BeamEnd3D, BeamEnd3D>> beamEnds;
    // 仅最后一个含索阶段：切线刚度上的活载梁端包络。无索模型不填，避免和线弹性包络重复。
    std::vector<LiveBeamEnvelope3D> liveBeam;
    std::vector<StageLiveEnvelope3D> liveDisp;
    LaneImpact3D liveImpact;
    bool liveTangent = false;
};

inline void fillTangentGirderLive(const Model3D& mc, const NLTangent3D& tangent,
                                  NLStageResult3D& row);

inline std::vector<NLStageResult3D> solveModel3DStagesNonlinear(
        const Model3D& m, NonlinearOptions3D nl = {}) {
    std::vector<char> active(m.elements.size(), 0);
    std::map<int, double> L0state;
    std::map<std::pair<int, int>, double> loadSum;
    // 束号 → 等效节点荷载。后一阶段不写这束时荷载仍在，全量平衡才不会把预加力丢掉。
    std::map<int, std::vector<std::tuple<int, int, double>>> psKeep;
    // 徐变只记扣除等效荷载后的应力合力。下一次增量乘的是这个力，不是含徐变位移的内力。
    // 节点力必须是单元实际交到节点上的力。用线性坐标变换去乘局部力，会和共转列式差一截，
    // 弯矩会逐次变大。等效荷载一旦加上就留在后面的全量平衡里。
    std::map<int, std::array<double, 12>> creepStress;      // 局部，用于回报内力
    std::map<int, std::array<double, 12>> creepApplied;     // 已作为外荷载的局部合力
    std::map<int, std::array<double, 12>> creepForceG;      // 局部应力对应的全局节点力
    std::map<int, std::array<double, 12>> creepLoadG;       // 已施加的全局徐变节点力
    double prevCreepTime = 0.0;
    // 单元第一次激活的龄期。规范收缩只计这之后新发展的应变，不把浇筑前的收缩放进平衡。
    std::map<int, double> shrinkBirth;
    int psSeq = 0;
    std::vector<NLStageResult3D> out;

    for (const auto& st : m.stages) {
        for (int id : st.activate) {
            bool found = false;
            for (std::size_t i = 0; i < m.elements.size(); ++i)
                if (m.elements[i].id == id) { active[i] = 1; found = true; break; }
            if (!found)
                throw std::runtime_error("非线性阶段激活了不存在单元 " + std::to_string(id));
        }
        for (int id : st.deactivate) {
            bool found = false;
            for (std::size_t i = 0; i < m.elements.size(); ++i)
                if (m.elements[i].id == id) { active[i] = 0; found = true; break; }
            if (!found)
                throw std::runtime_error("非线性阶段钝化了不存在单元 " + std::to_string(id));
        }
        for (const auto& ld : st.loads)
            loadSum[{std::get<0>(ld), std::get<1>(ld)}] += std::get<2>(ld);
        for (std::size_t i = 0; i < m.elements.size(); ++i) {
            if (!active[i]) continue;
            const int id = m.elements[i].id;
            if (!shrinkBirth.count(id))
                shrinkBirth[id] = st.time > 0.0 ? st.time : m.creepT0;
        }

        Model3D mc;
        mc.nodes = m.nodes;
        mc.constraints = m.constraints;
        // 一次徐变：刚度用 E/(1+χφ)。等效预加力不随 E 变，静定梁位移按 1+χφ 放大。
        mc.creepPhi = m.creepPhi;
        mc.creepChi = m.creepChi;
        for (const auto& [key, val] : loadSum)
            mc.loads.emplace_back(key.first, key.second, val);
        for (std::size_t i = 0; i < m.elements.size(); ++i) {
            if (!active[i]) continue;
            Element3D e = m.elements[i];
            // 收缩是初应变，每个全量平衡都在，但不是再加一份荷载。
            // 正值缩短。后面的徐变阶段仍用同一个 e.ecs，不会再缩短一个 ε·L。
            if (m.shrinkageEcs != 0.0) {
                e.ecs += m.shrinkageEcs;
            } else if (e.type == "BEAM" && st.time > 0.0
                       && m.creepRh > 0.0 && m.creepFcm > 0.0) {
                // 手写收缩应变优先。没写时，混凝土梁按附录 C.1 随龄期缩短。索不收缩。
                // 模型上的理论厚度优先；否则用该截面自己的 h=2A/u。
                const double h = m.creepH > 0.0 ? m.creepH : e.hTheory;
                if (h <= 0.0) {
                    mc.elements.push_back(e);
                    continue;
                }
                const double birth = shrinkBirth.count(e.id) ? shrinkBirth[e.id] : m.creepT0;
                const double ts = m.shrinkageTs > 0.0 ? m.shrinkageTs : m.creepT0;
                const double eps = jtg3362ShrinkageStrain(
                    st.time, birth, ts, m.creepRh, h, m.creepFcm);
                if (eps != 0.0) e.ecs += eps;
            }
            if (e.type == "CABLE") {
                if (!L0state.count(e.id))
                    L0state[e.id] = (e.L0 > 0.0)
                                        ? e.L0
                                        : (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm();
                e.L0 = L0state[e.id];
            }
            mc.elements.push_back(e);
        }
        for (const auto& [eid, val] : st.unstrained) {
            if (val <= 0.0)
                throw std::runtime_error("非线性阶段 L0= 必须为正");
            bool found = false;
            for (auto& e : mc.elements)
                if (e.id == eid && e.type == "CABLE") {
                    e.L0 = val;
                    L0state[eid] = val;
                    found = true;
                    break;
                }
            if (!found)
                throw std::runtime_error("非线性阶段 L0= 指向未激活的索 " + std::to_string(eid));
        }
        const std::string stag = "*STAGE " + std::to_string(st.id);
        for (const auto& ps : st.prestresses) {
            // 束号按阶段内顺序从 1 编，跨阶段会撞号。全量平衡用独立序号，两阶段各张一束才都留得住。
            psKeep[++psSeq] = prestressEquivalentLoads(
                m, mc.elements, ps, stag + " 的 PRESTRESS= 束 " + std::to_string(ps.id));
        }
        for (const auto& [id, lds] : psKeep) {
            (void)id;
            for (const auto& ld : lds) mc.loads.push_back(ld);
        }
        // 有龄期时 Δφ=φ(t)−φ(上一阶段)。第一阶段还没有已积累的应力，增量为 0。
        // 手写 DPHI 只在没写 TIME 时使用，避免同一阶段乘两次。
        const bool codeCreep = m.creepRh > 0.0 && m.creepFcm > 0.0
            && (m.creepH > 0.0 || std::any_of(mc.elements.begin(), mc.elements.end(),
                    [](const Element3D& el) { return el.hTheory > 0.0; }));
        const bool ageCreep = m.creepPhi0 > 0.0 || codeCreep;
        double dphi = st.creepDphi;
        if (st.time > 0.0 && ageCreep && !codeCreep) {
            const double phiNow = creepCoefficientAt3D(m, st.time);
            const double phiPrev = prevCreepTime > 0.0
                                       ? creepCoefficientAt3D(m, prevCreepTime) : 0.0;
            dphi = prevCreepTime > 0.0 ? (phiNow - phiPrev) : 0.0;
        }
        auto dphiAt = [&](double hElem) {
            if (!(st.time > 0.0 && codeCreep)) return dphi;
            const double phiNow = creepCoefficientAt3D(m, st.time, hElem);
            const double phiPrev = prevCreepTime > 0.0
                                       ? creepCoefficientAt3D(m, prevCreepTime, hElem) : 0.0;
            return prevCreepTime > 0.0 ? (phiNow - phiPrev) : 0.0;
        };
        for (const auto& e : mc.elements) {
            if (e.type != "BEAM") continue;
            const double dphiE = codeCreep ? dphiAt(e.hTheory) : dphi;
            if (dphiE == 0.0) continue;
            const auto itG = creepForceG.find(e.id);
            const auto itL = creepStress.find(e.id);
            if (itG == creepForceG.end() || itL == creepStress.end()) continue;
            auto& g = creepLoadG[e.id];
            auto& loc = creepApplied[e.id];
            for (int i = 0; i < 12; ++i) {
                g[static_cast<std::size_t>(i)] +=
                    dphiE * itG->second[static_cast<std::size_t>(i)];
                loc[static_cast<std::size_t>(i)] +=
                    dphiE * itL->second[static_cast<std::size_t>(i)];
            }
        }
        for (const auto& e : mc.elements) {
            if (e.type != "BEAM") continue;
            const auto it = creepLoadG.find(e.id);
            if (it == creepLoadG.end()) continue;
            bool any = false;
            for (double v : it->second)
                if (v != 0.0) any = true;
            if (!any) continue;
            for (int k = 0; k < 6; ++k) {
                mc.loads.push_back({e.n1, k + 1, it->second[static_cast<std::size_t>(k)]});
                mc.loads.push_back({e.n2, k + 1, it->second[static_cast<std::size_t>(6 + k)]});
            }
        }
        constrainDisconnected3D(mc);

        std::vector<CableTarget3D> targets;
        for (const auto& [eid, val] : st.tensions) {
            if (val <= 0.0)
                throw std::runtime_error("非线性阶段 TENSION= 必须为正");
            bool found = false;
            for (const auto& e : mc.elements)
                if (e.id == eid && e.type == "CABLE") { found = true; break; }
            if (!found)
                throw std::runtime_error("非线性阶段 TENSION= 指向未激活的索 " + std::to_string(eid));
            targets.push_back(CableTarget3D{eid, val});
        }

        bool tuned = true;
        if (!targets.empty()) {
            CableLengthOptions3D co;
            co.nl = nl;
            co.maxIter = 12;
            co.tolTension = 1e-4;
            const CableLengthResult3D cl = iterateCableUnstrainedLengths(mc, targets, co);
            tuned = cl.converged;
            for (std::size_t k = 0; k < cl.elems.size(); ++k) {
                L0state[cl.elems[k]] = cl.L0[k];
                for (auto& e : mc.elements)
                    if (e.id == cl.elems[k]) e.L0 = cl.L0[k];
            }
        }

        const bool lastStage = (&st == &m.stages.back());
        bool hasCable = false;
        for (const auto& e : mc.elements)
            if (e.type == "CABLE") hasCable = true;
        NLTangent3D tangent;
        const NonlinearResult3D r = solveModel3DNonlinear(
            mc, nl, (lastStage && hasCable && !m.lanes.empty()) ? &tangent : nullptr);
        NLStageResult3D row;
        row.id = st.id;
        row.newton = r.trace;
        row.converged = tuned && r.converged;
        row.dm = buildDofMap3D(mc);
        row.u = r.u;
        row.cableL0 = L0state;
        std::map<int, std::string> kind;
        for (const auto& e : mc.elements) kind[e.id] = e.type;
        for (std::size_t k = 0; k < r.elemIds.size(); ++k) {
            const int id = r.elemIds[k];
            const std::string& tp = kind[id];
            if (tp == "CABLE") {
                row.cableForce[id] = r.axial[k];
                continue;
            }
            if (tp != "BEAM" && tp != "BEAM7") continue;
            const auto& q = r.localForces[k];
            if (tp == "BEAM7" && q.size() >= 14) {
                BeamEnd3D i{-q[0], q[1], q[2], -q[3], -q[4], -q[5]};
                BeamEnd3D j{q[7], -q[8], -q[9], q[10], q[11], q[12]};
                row.beamEnds[id] = {i, j};
            } else if (q.size() >= 12) {
                std::array<double, 12> qs{};
                for (int i = 0; i < 12; ++i) qs[static_cast<std::size_t>(i)] = q[i];
                const auto applied = creepApplied.find(id);
                if (applied != creepApplied.end()) {
                    for (int i = 0; i < 12; ++i)
                        qs[static_cast<std::size_t>(i)] -= applied->second[static_cast<std::size_t>(i)];
                }
                creepStress[id] = qs;
                const Element3D* el = nullptr;
                for (const auto& e : mc.elements)
                    if (e.id == id) { el = &e; break; }
                if (el) {
                    const CRElemData ed = crMakeElemData(mc, *el, row.dm);
                    Eigen::Matrix<double, 14, 1> ue = Eigen::Matrix<double, 14, 1>::Zero();
                    Eigen::Matrix<double, 14, 1> fe = Eigen::Matrix<double, 14, 1>::Zero();
                    for (int i = 0; i < ed.n; ++i) {
                        const int d = ed.gdofs[static_cast<std::size_t>(i)];
                        if (d >= 0 && d < static_cast<int>(row.u.size())) ue[i] = row.u[d];
                    }
                    crElementInternal(ed, ue, fe, false);
                    std::array<double, 12> fg{};
                    const auto held = creepLoadG.find(id);
                    for (int i = 0; i < 12; ++i) {
                        fg[static_cast<std::size_t>(i)] = fe[i];
                        if (held != creepLoadG.end())
                            fg[static_cast<std::size_t>(i)] -= held->second[static_cast<std::size_t>(i)];
                    }
                    creepForceG[id] = fg;
                }
                BeamEnd3D i{-qs[0], qs[1], qs[2], -qs[3], -qs[4], -qs[5]};
                BeamEnd3D j{qs[6], -qs[7], -qs[8], qs[9], qs[10], qs[11]};
                row.beamEnds[id] = {i, j};
            }
        }
        if (st.time > 0.0) prevCreepTime = st.time;
        if (lastStage && hasCable && !m.lanes.empty() && r.converged && tangent.K.rows() > 0) {
            mc.lanes = m.lanes;
            mc.influences = m.influences;
            fillTangentGirderLive(mc, tangent, row);
        }
        out.push_back(std::move(row));
    }
    return out;
}

inline void fillTangentGirderLive(const Model3D& mc, const NLTangent3D& tangent,
                                  NLStageResult3D& row) {
    std::vector<const Element3D*> tgt;
    for (const auto& e : mc.elements) {
        if (e.type != "BEAM" && e.type != "BEAM7") continue;
        const Vec3 d = mc.nodes.at(e.n2) - mc.nodes.at(e.n1);
        if (std::abs(d.x()) < std::max(std::abs(d.y()), std::abs(d.z()))) continue;
        tgt.push_back(&e);
    }
    std::vector<char> isC(tangent.K.rows(), 0);
    for (int d : tangent.cdofs)
        if (d >= 0 && d < static_cast<int>(isC.size())) isC[d] = 1;
    const Factorization3D fac(tangent.K);
    const DofMap3D dm = buildDofMap3D(mc);
    bool wantAuto = false;
    for (const auto& lane : mc.lanes)
        if (lane.impact == 1.0) wantAuto = true;
    // 基频用这套切线，不用线弹性 EA/L。没有梁时仍要算挠度和冲击系数。
    row.liveImpact = wantAuto ? resolveLaneImpact3D(mc, &tangent.K) : LaneImpact3D{};
    row.liveTangent = true;
    if (!tgt.empty()) {
    for (const auto& lane : mc.lanes) {
        if (lane.nodes.empty()) continue;
        Model3D::Lane3D laneUse = lane;
        if (lane.impact == 1.0 && row.liveImpact.automatic)
            laneUse.impact = row.liveImpact.factor;
        std::vector<std::vector<double>> eta(tgt.size() * 4,
                                             std::vector<double>(lane.nodes.size(), 0.0));
        scanTangentLane3D(fac, dm, isC, lane.nodes, lane.dof,
            [&](std::size_t k, const std::vector<double>& U) {
                if (k >= lane.nodes.size()) return;
                for (std::size_t t = 0; t < tgt.size(); ++t) {
                    const Element3D& e = *tgt[t];
                    const bool is7 = (e.type == "BEAM7");
                    const int dpe = is7 ? 14 : 12;
                    const Eigen::VectorXd fl =
                        is7 ? beam7LocalEndForces(mc, dm, e, U)
                            : beam3DLocalEndForces(mc, dm, e, U);
                    for (int end = 1; end <= 2; ++end) {
                        const BeamEndForces bf = beamEndForcesAt(fl, dpe, end);
                        eta[t * 4 + (end - 1) * 2 + 0][k] = bf.My;
                        eta[t * 4 + (end - 1) * 2 + 1][k] = bf.Vz;
                    }
                }
            });
        for (std::size_t t = 0; t < tgt.size(); ++t) {
            for (int end = 1; end <= 2; ++end) {
                const auto eMy = envelope3DFromInfluence(
                    mc, laneUse, eta[t * 4 + (end - 1) * 2 + 0]);
                const auto eVz = envelope3DFromInfluence(
                    mc, laneUse, eta[t * 4 + (end - 1) * 2 + 1]);
                LiveBeamEnvelope3D slot;
                slot.elem = tgt[t]->id;
                slot.end = end;
                slot.MyMin = eMy.second;
                slot.MyMax = eMy.first;
                slot.VzMin = eVz.second;
                slot.VzMax = eVz.first;
                row.liveBeam.push_back(slot);
            }
        }
    }
    }
    if (mc.influences.empty()) return;
    const auto& spec = mc.influences.front();
    const int dResp = dm.safeDof(spec.node, spec.dof - 1);
    if (dResp < 0) return;
    double maxAll = -1e300, minAll = 1e300;
    for (const auto& lane : mc.lanes) {
        if (lane.nodes.empty()) continue;
        Model3D::Lane3D laneUse = lane;
        laneUse.impact = 1.0;   // 位移包络不计冲击（表 4.3.2 挠度限值，冲击只放大内力）
        std::vector<double> eta(lane.nodes.size(), 0.0);
        scanTangentLane3D(fac, dm, isC, lane.nodes, lane.dof,
            [&](std::size_t k, const std::vector<double>& U) {
                if (k >= eta.size()) return;
                if (dResp >= 0 && dResp < static_cast<int>(U.size()))
                    eta[k] = U[dResp];
            });
        const auto env = envelope3DFromInfluence(mc, laneUse, eta);
        maxAll = std::max(maxAll, env.first);
        minAll = std::min(minAll, env.second);
    }
    if (maxAll > -1e299)
        row.liveDisp.push_back({row.id, maxAll, minAll});
}

}  // namespace fem
