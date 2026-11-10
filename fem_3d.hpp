#pragma once

// 复用 fem.hpp 中的稀疏求解器与 Dirichlet 处理（applyDirichlet / solve）。
#include "fem.hpp"
#include "code_loads.hpp"
#include "prestress_loss.hpp"
#include "section.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace fem {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using SparseMat3D = Eigen::SparseMatrix<double>;

struct Model3D;
inline std::vector<double> solveModel3D(const Model3D& m);

// 三维欧拉-伯努利梁局部刚度矩阵，12 自由度：
// [u1,v1,w1,θx1,θy1,θz1,u2,v2,w2,θx2,θy2,θz2]
inline Eigen::Matrix<double, 12, 12> beam3DLocalStiffness(
    double E, double A, double Iy, double Iz, double J, double G, double L) {
    if (L <= 0.0) throw std::runtime_error("三维梁单元长度为零");
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    const double EAL = E * A / L;
    const double kz = E * Iz / (L * L * L);
    const double ky = E * Iy / (L * L * L);
    const double GJL = G * J / L;

    K(0, 0) = EAL;  K(0, 6) = -EAL;
    K(6, 0) = -EAL; K(6, 6) = EAL;

    K(1, 1) = 12 * kz; K(1, 5) = 6 * kz * L; K(1, 7) = -12 * kz; K(1, 11) = 6 * kz * L;
    K(5, 1) = 6 * kz * L; K(5, 5) = 4 * kz * L * L; K(5, 7) = -6 * kz * L; K(5, 11) = 2 * kz * L * L;
    K(7, 1) = -12 * kz; K(7, 5) = -6 * kz * L; K(7, 7) = 12 * kz; K(7, 11) = -6 * kz * L;
    K(11, 1) = 6 * kz * L; K(11, 5) = 2 * kz * L * L; K(11, 7) = -6 * kz * L; K(11, 11) = 4 * kz * L * L;

    K(2, 2) = 12 * ky; K(2, 4) = -6 * ky * L; K(2, 8) = -12 * ky; K(2, 10) = -6 * ky * L;
    K(4, 2) = -6 * ky * L; K(4, 4) = 4 * ky * L * L; K(4, 8) = 6 * ky * L; K(4, 10) = 2 * ky * L * L;
    K(8, 2) = -12 * ky; K(8, 4) = 6 * ky * L; K(8, 8) = 12 * ky; K(8, 10) = 6 * ky * L;
    K(10, 2) = -6 * ky * L; K(10, 4) = 2 * ky * L * L; K(10, 8) = 6 * ky * L; K(10, 10) = 4 * ky * L * L;

    K(3, 3) = GJL; K(3, 9) = -GJL;
    K(9, 3) = -GJL; K(9, 9) = GJL;
    return K;
}

inline Eigen::Matrix<double, 12, 12> beam3DTransformation(const Vec3& p1, const Vec3& p2) {
    Vec3 ex = p2 - p1;
    const double L = std::max(ex.norm(), 1e-12);
    ex /= L;
    Vec3 ez(0.0, 0.0, 1.0);
    if (std::abs(ex.dot(ez)) > 0.999) ez = Vec3(0.0, 1.0, 0.0);
    Vec3 ey = ez.cross(ex); ey.normalize();
    ez = ex.cross(ey); ez.normalize();
    Eigen::Matrix3d R;
    R.col(0)=ex; R.col(1)=ey; R.col(2)=ez;
    Eigen::Matrix<double, 12, 12> T = Eigen::Matrix<double, 12, 12>::Zero();
    // 方向余弦矩阵 λ = R^T 将全局位移映射到局部：u_local = λ·u_global。
    // （修正：此前误用 R，把变换矩阵转置了，导致非轴线对齐梁刚度算错。）
    for (int b=0; b<4; ++b) T.block<3,3>(3*b,3*b)=R.transpose();
    return T;
}

// 7 自由度（含翘曲）薄壁梁局部刚度，14 自由度：
// [u1,v1,w1,θx1,θy1,θz1,φ1, u2,v2,w2,θx2,θy2,θz2,φ2]
// φ = dθx/dx 为翘曲率。扭转项 = 圣维南 GJ·∫(φ)² + 翘曲 E·Iw·∫(φ')²。
inline Eigen::Matrix<double, 14, 14> beam7LocalStiffness(
    double E, double A, double Iy, double Iz, double J, double G, double Iw, double L) {
    if (L <= 0.0) throw std::runtime_error("7 自由度梁单元长度为零");
    Eigen::Matrix<double, 14, 14> K = Eigen::Matrix<double, 14, 14>::Zero();
    const double EAL = E * A / L;
    const double kz = E * Iz / (L * L * L);
    const double ky = E * Iy / (L * L * L);

    K(0, 0) = EAL;  K(0, 7) = -EAL;
    K(7, 0) = -EAL; K(7, 7) = EAL;

    // 绕 z 轴弯曲（v,θz）：1,5,8,12
    K(1, 1) = 12*kz; K(1, 5) = 6*kz*L; K(1, 8) = -12*kz; K(1, 12) = 6*kz*L;
    K(5, 1) = 6*kz*L; K(5, 5) = 4*kz*L*L; K(5, 8) = -6*kz*L; K(5, 12) = 2*kz*L*L;
    K(8, 1) = -12*kz; K(8, 5) = -6*kz*L; K(8, 8) = 12*kz; K(8, 12) = -6*kz*L;
    K(12, 1) = 6*kz*L; K(12, 5) = 2*kz*L*L; K(12, 8) = -6*kz*L; K(12, 12) = 4*kz*L*L;

    // 绕 y 轴弯曲（w,θy）：2,4,9,11（耦合符号相反）
    K(2, 2) = 12*ky; K(2, 4) = -6*ky*L; K(2, 9) = -12*ky; K(2, 11) = -6*ky*L;
    K(4, 2) = -6*ky*L; K(4, 4) = 4*ky*L*L; K(4, 9) = 6*ky*L; K(4, 11) = 2*ky*L*L;
    K(9, 2) = -12*ky; K(9, 4) = 6*ky*L; K(9, 9) = 12*ky; K(9, 11) = 6*ky*L;
    K(11, 2) = -6*ky*L; K(11, 4) = 2*ky*L*L; K(11, 9) = 6*ky*L; K(11, 11) = 4*ky*L*L;

    // 扭转（θx,φ）：3,6,10,13。翘曲项按梁弯曲刚度，圣维南项按 ∫N'^T N'。
    const double kw = E * Iw / (L * L * L);
    const double ks = G * J / L;
    K(3, 3)   += 12*kw;      K(3, 6)   += 6*kw*L;    K(3, 10)  += -12*kw;     K(3, 13)  += 6*kw*L;
    K(6, 3)   += 6*kw*L;     K(6, 6)   += 4*kw*L*L;  K(6, 10)  += -6*kw*L;    K(6, 13)  += 2*kw*L*L;
    K(10, 3)  += -12*kw;     K(10, 6)  += -6*kw*L;   K(10, 10) += 12*kw;      K(10, 13) += -6*kw*L;
    K(13, 3)  += 6*kw*L;     K(13, 6)  += 2*kw*L*L;  K(13, 10) += -6*kw*L;    K(13, 13) += 4*kw*L*L;

    K(3, 3)   += 6.0/5.0*ks;   K(3, 6)   += L/10.0*ks;      K(3, 10)  += -6.0/5.0*ks;  K(3, 13)  += L/10.0*ks;
    K(6, 3)   += L/10.0*ks;    K(6, 6)   += 2.0*L*L/15.0*ks; K(6, 10)  += -L/10.0*ks;   K(6, 13)  += -L*L/30.0*ks;
    K(10, 3)  += -6.0/5.0*ks;  K(10, 6)  += -L/10.0*ks;     K(10, 10) += 6.0/5.0*ks;   K(10, 13) += -L/10.0*ks;
    K(13, 3)  += L/10.0*ks;    K(13, 6)  += -L*L/30.0*ks;   K(13, 10) += -L/10.0*ks;   K(13, 13) += 2.0*L*L/15.0*ks;
    return K;
}

// 7 自由度梁变换（全局->局部）。每节点：平动/转动用 λ=R^T，翘曲率 φ 不参与旋转。
inline Eigen::Matrix<double, 14, 14> beam7Transformation(const Vec3& p1, const Vec3& p2) {
    Vec3 ex = p2 - p1;
    const double L = std::max(ex.norm(), 1e-12);
    ex /= L;
    Vec3 ez(0.0, 0.0, 1.0);
    if (std::abs(ex.dot(ez)) > 0.999) ez = Vec3(0.0, 1.0, 0.0);
    Vec3 ey = ez.cross(ex); ey.normalize();
    ez = ex.cross(ey); ez.normalize();
    Eigen::Matrix3d R;
    R.col(0)=ex; R.col(1)=ey; R.col(2)=ez;
    const Eigen::Matrix3d lam = R.transpose();
    Eigen::Matrix<double, 14, 14> T = Eigen::Matrix<double, 14, 14>::Zero();
    for (int n = 0; n < 2; ++n) {
        const int b = 7 * n;
        T.block<3, 3>(b, b) = lam;         // u,v,w
        T.block<3, 3>(b + 3, b + 3) = lam; // θx,θy,θz
        T(b + 6, b + 6) = 1.0;             // φ
    }
    return T;
}

inline Eigen::Matrix<double, 6, 6> truss3DStiffness(double E, double A, const Vec3& p1, const Vec3& p2) {
    const Vec3 d = p2 - p1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("三维桁架单元长度为零");
    const Vec3 c = d / L;
    const Mat3 K = (E * A / L) * (c * c.transpose());
    Eigen::Matrix<double, 6, 6> Ke = Eigen::Matrix<double, 6, 6>::Zero();
    Ke.block<3, 3>(0, 0) = K;
    Ke.block<3, 3>(0, 3) = -K;
    Ke.block<3, 3>(3, 0) = -K;
    Ke.block<3, 3>(3, 3) = K;
    return Ke;
}

// ==================== 缆索垂度：索力 ↔ 弦长 ====================
//
// 这一段放在线性内核里，而不是非线性内核里，是为了让**两条路径共用同一
// 个式子**：施工阶段调索（线性、增量法）要把"目标索力"翻译成"无应力索长"，
// 成桥态找形（几何非线性）要用同一个关系反算 L0。两边各写一份的话，迟早
// 会漂移成两个不同的物理模型——而阶段索力与成桥态索力是要互相校核的。
//
// 缆索不是直杆：自重让它垂成悬链线，弦的伸长因此"看起来"更软。浅垂度
// （f/L < 1/8，斜拉桥拉索与悬索桥主缆均满足）下的力-伸长关系：
//
//      Ln = L0·(1 + T/(EA)) - w²·Lh³/(24·T²)
//
// 其中 Ln = 弦长，L0 = 无应力（弧）长度，Lh = 弦的水平投影，w = 单位长重度。
// 第二项即垂度缩短量 s(T) = w²Lh³/(24T²)（抛物线悬链线结果）：索绷得越紧，
// 垂度越小、弦越长；索放松则垂度增大、弦缩短。
//
// 对该关系求导即得**切线刚度**（注意分母里的 L0 不能省）：
//
//      dT/dLn = 1 / (L0/(EA) + w²Lh³/(12T³))
//             = (EA/L0) / (1 + w²·Lh³·EA/(12·T³·L0)) = Eeq·A/L0
//      Eeq    = E / (1 + w²·Lh³·E·A/(12·T³·L0))   ← 这正是 Ernst 等效弹模
//
// 标准 Ernst 式写作 E/(1 + (wL)²EA/(12T³))，那里默认 L ≈ L0 ≈ Lh；
// 本内核保留分母中的 L0，使"内力—切线"严格共轭（二者只差 O(弹性应变)）。
//
// 注意：上式是**固定水平投影 Lh** 下的偏导，正是 Ernst 式的定义（Ernst 推导
// 时把水平跨径 L 当常量）。单元实际切线由差分对节点位置求得，会额外带上
// ∂/∂Lh 项（端点沿弦向移动时 Lh 同时变），量级仅 O(垂度²)；它与内力严格
// 共轭，才是 Newton 二次收敛的保证。Ernst 模量在这里的定位是"工程上可读
// 的等效弹模"，而不是切线的替代品。

// 由弦长反解索力（浅垂度弹性悬链线）。T>0 唯一解（g 单调）。
inline double cableTensionFromChord(double Ln, double L0, double EA, double w, double Lh) {
    if (w <= 0.0 || Lh <= 1e-12) {
        return EA * (Ln - L0) / L0;          // 竖直索或未给重度：退化为直杆
    }
    const double cSag = w * w * Lh * Lh * Lh / 24.0;    // s(T) = cSag / T²
    // 初值：张紧时用直杆估计；松弛时用抛物线垂度下的水平力估计
    double T = EA * (Ln - L0) / L0;
    if (T < 1e-6) {
        const double slack = std::max(L0 - Ln, 1e-9);
        T = std::max(w * Lh * Lh / (8.0 * std::sqrt(3.0 * Lh * slack / 8.0)), 1e-3);
    }
    // g(T) = L0 + L0·T/EA - cSag/T² - Ln；g'(T) = L0/EA + 2cSag/T³ > 0
    for (int it = 0; it < 60; ++it) {
        const double g = L0 + L0 * T / EA - cSag / (T * T) - Ln;
        const double gp = L0 / EA + 2.0 * cSag / (T * T * T);
        const double step = g / gp;
        T -= step;
        if (T < 1e-6) T = 1e-6;             // 防止 Newton 步长跨越奇点
        if (std::abs(step) <= 1e-12 * std::max(1.0, T)) break;
    }
    return T;
}

// ① 由设计弦长与目标索力闭式反算无应力（下料）索长。
//    即 Ln = L0(1 + T/EA) - s(T) 对 L0 求解，s(T) = w²Lh³/(24T²)。
inline double cableUnstrainedLengthForTension(double chordTarget, double T,
                                              double EA, double w, double Lh) {
    if (T <= 0.0) throw std::runtime_error("反算下料索长：目标索力必须为正");
    double sag = 0.0;
    if (w > 0.0 && Lh > 1e-12) sag = w * w * Lh * Lh * Lh / (24.0 * T * T);
    return (chordTarget + sag) / (1.0 + T / EA);
}

// 变截面梁单元：截面属性沿梁轴在控制点间线性插值，取 [s1,s2] 中点处截面
// 作为该单元的等效常截面（与 2D 变截面“单元中点截面”近似一致）。
inline BeamElement beamElementFromSection(const VariableSection& vs,
                                          double s1, double s2,
                                          const Vec3& p1, const Vec3& p2,
                                          double E, double G) {
    const SectionProperties sec = elementMidSection(vs, s1, s2);
    BeamElement b;
    b.E = E; b.G = G;
    b.A = sec.A; b.Iy = sec.Iy; b.Iz = sec.Iz; b.J = sec.J;
    b.x0 = p1; b.x1 = p2;
    return b;
}

struct Node3D {
    int id = 0;
    Vec3 x = Vec3::Zero();
};

struct Element3D {
    int id = 0;
    std::string type = "TRUSS";
    double E = 0.0, A = 0.0;
    double Iy = 0.0, Iz = 0.0, J = 0.0, G = 0.0;
    double Iw = 0.0;                 // 扇形惯性矩/翘曲常数（BEAM7 约束扭转用）
    double dens = 0.0, addmass = 0.0;   // 密度 kg/m^3、附加线质量 kg/m（模态/自重）
    double istrn = 0.0;
    // 无应力长度（m）：CABLE/索长迭代用；<=0 表示取几何长度。
    // 索的"无应力索长"是成桥态找形的核心量——给定索力与几何即可反算 L0，
    // 给定 L0 才能在非线性求解里得到正确的成桥态（见索长迭代）。
    double L0 = 0.0;
    // 缆索单位长度重量（N/m）：>0 时启用垂度效应（浅垂度弹性悬链线 +
    // Ernst 等效模量）。斜拉桥拉索、悬索桥主缆的垂度非线性由此体现。
    double w = 0.0;
    double Wtop = 0.0, Wbot = 0.0;
    double Wy = 0.0;   // 绕 y 轴截面模量（水平弯曲 My 用）
    double Wt = 0.0;   // 扭转截面模量 J/r_max（扭转剪应力 τ=T/Wt 用）
    int n1 = 0, n2 = 0;
    // 本阶段才计入的收缩应变。增量法里全局 shrinkageEcs 不能每个阶段都加，
    // 只在单元第一次参与受力时写入这里。
    double ecs = 0.0;
    // 理论厚度 mm。0 表示没有。模型上的 creepH>0 时以模型值为准。
    double hTheory = 0.0;
};

// 三维约束：节点 + 自由度(1~6) + 指定位移值。
// value 默认 0；非零即支座沉降 / 强制位移。
// （旧实现只有 (node, dof) 两个字段，且求解时硬编码 F(d)=0，
//  三维模型完全无法计算支座沉降工况。）
struct Constraint3D {
    int node = 0;
    int dof = 1;
    double value = 0.0;
    double k = 0.0;   // 弹性支承刚度（>0 为弹簧，0 为刚性约束）
};

struct Model3D {
    std::map<int, Vec3> nodes;
    std::vector<Element3D> elements;
    std::vector<Constraint3D> constraints;
    std::vector<std::tuple<int, int, double>> loads; // node, dof, value

    // ---- 一条预应力筋（等效荷载法，直线 + 二次抛物线曲线）----
    //
    // 与 2D 的 `*PRESTRESS` 不同：这里**不生成力筋单元**，只把它翻译成等效
    // 节点荷载（见 solveModel3DStages 里的 applyPrestress）：
    //   1. 轴向压力 Pe，沿梁轴（两端一对，压短梁）；
    //   2. 锚固端偏心弯矩 Pe·e，绕横轴（竖向偏心 e 产生的附加弯矩）；
    //   3. 曲线筋（s≠0）的跨内反向均布力：两端锚固切线竖向力 4Pe·s/L（两端
    //      同向，指"束凹侧"）＋ 竖向均布 q = −8Pe·s/L² 节点化（首尾各半跨、
    //      中间节点全跨）。这两项合起来静力自平衡（竖向合力=0），是曲线筋
    //      "上拱"的唯一来源。s=0 时这两项均为 0，退化为直线筋。
    //
    // 约定：
    //   P   = 张拉控制力（N，正）；
    //   e   = 端部偏心（m，**形心上方为正**，即 +z 方向）；直线筋时即常偏心；
    //   s   = 矢高 = 跨中偏心 − 端部偏心（m，向上为正）。s<0 → 跨中筋更低 →
    //         产生向上的均布力 → **上拱**；s=0 → 直线筋。曲线筋的抛物线定义
    //         在首尾锚固端连线上（矢高在弦中点），故要求链共线；
    //   loss = 有效/张拉比（0<loss≤1）。两种给了法互斥：
    //          LOSS=<比值> 是「总损失打包成一个数」的简化写法；
    //          LOSSREF=<编号> 引用 *PRESTRESS_LOSS 表，按 JTG 3362-2018 §6.2
    //          分项算 σl1~σl6 后由 σpe/σcon 得到有效比。解析完成后统一写回
    //          loss，所以下面的等效荷载与结果段都只认 loss 这一个量。
    // 直线筋（s=0）**不产生上拱**——它只给轴压 + 端弯矩，跨内无反弯力。
    struct Prestress3D {
        int id = 0;
        std::vector<int> elems;      // 连续梁单元（按节点顺序，首尾即锚固端）
        double P = 0.0;
        double e = 0.0;
        double s = 0.0;              // 矢高（m，向上为正；0 = 直线筋）
        double loss = 1.0;
        bool lossGiven = false;      // 是否显式写过 LOSS=（用于与 LOSSREF 互斥）
        int lossRef = 0;             // >0：引用 m.prestressLosses[lossRef]
        PrestressLossResult lossDetail;   // 分项明细（LOSSREF 路径才有）
        bool hasLossDetail = false;
    };

    struct Stage3D {
        int id = 0;
        std::vector<int> activate;
        std::vector<int> deactivate;
        std::vector<std::tuple<int, int, double>> loads;
        // ---- 本阶段的索张拉（施工阶段"调索"的唯一入口）----
        //
        // tensions:   单元号 → 目标索力 N（正为张拉）。求解时迭代无应力索长
        //             L0，使本阶段末索力命中目标。
        // unstrained: 单元号 → 无应力索长 m。直接给定下料长度，不做迭代——
        //             L0 本身就是最终答案，再迭代只会把它改掉。
        // 两者对同一单元**互斥**：同时给出直接报错，因为"迭代命中目标"与
        // "索长已给定、力是结果"是相反的语义，不能取其一了事。
        std::vector<std::pair<int, double>> tensions;
        std::vector<std::pair<int, double>> unstrained;
        // 本阶段张拉的预应力束。预应力的效应是**状态量**（张拉一次、永久存在），
        // 但阶段是增量法：等效荷载只在本阶段施加一次，产生的变形与内力通过
        // 累计位移跨阶段延续——这与索力"状态量、不累加"是同一套逻辑。
        std::vector<Prestress3D> prestresses;
        // 相对上一阶段已积累内力的徐变增量 Δφ。只作用在已有应力上，
        // 本阶段新荷载仍用弹性模量。静定结构的附加位移约为 Δφ 倍既有位移。
        double creepDphi = 0.0;
        // 阶段末龄期（天）。>0 且模型给了 PHI0 时，Δφ 改由 φ(t)−φ(t_prev) 计算，
        // 不再使用 creepDphi。本阶段新加的荷载要到下一龄期才开始徐变。
        double time = 0.0;
    };
    std::vector<Stage3D> stages;
    // *PRESTRESS_LOSS 表：编号 → JTG 3362-2018 §6.2 分项损失参数。
    // 与 2D 共用 prestress_loss.hpp 的结构与公式；*STAGE 的 PRESTRESS= 用
    // LOSSREF=<编号> 引用，同一张表可被多束复用。
    std::map<int, PrestressLossInput> prestressLosses;
    double creepPhi = 0.0;
    double creepChi = 0.8;
    // 徐变曲线 φ = φ0·[(t−t0)/(βH+t−t0)]^power。φ0=0 表示不按龄期。
    double creepPhi0 = 0.0;
    double creepT0 = 7.0;
    double creepBetaH = 700.0;
    double creepPower = 0.3;
    // 三者都大于 0 时，名义徐变系数和 βH 改按 JTG 3362 附录 C.2.1，不再用上面的 PHI0。
    // RH 为百分数，H 为理论厚度 mm，FCM 为圆柱体平均强度 MPa。
    double creepRh = 0.0;
    double creepH = 0.0;
    double creepFcm = 0.0;
    // 干燥开始龄期（天）。0 表示与加载龄期 creepT0 相同。
    // 只影响规范收缩：参与受力之前已经发生的收缩不进入平衡。
    double shrinkageTs = 0.0;
    double shrinkageEcs = 0.0;

    // 与二维 creepCoefficientAt 的简化曲线相同。t≤t0 时为 0。
    // 多龄期用相邻阶段的差值，不把整段 φ 一次乘到全部应力上。

    struct Influence3D {
        int node = 0;
        int dof = 1;
    };
    std::vector<Influence3D> influences;

    struct ElementInfluence3D {
        int elem = 0;
    };
    std::vector<ElementInfluence3D> elementInfluences;

    struct Lane3D {
        std::vector<int> nodes;
        double q = 0.0;
        double p = 0.0;
        int dof = 2;
        double reduction = 1.0;
        double impact = 1.0;
        bool autoCode = false;   // CODE=JTG_D60：按跨径自动生成车道荷载
        // 横向布载（横向分布系数法，JTG D60 多车道）：单梁模型下横向影响面
        // 退化为「横向分布系数」，即车道数 × 多车道横向折减 × 偏载增大。
        // 默认值 (1, 1.0) 等价于旧的单车道近似，向后兼容。
        int laneCount = 1;         // 横向设计车道数
        double eccentricity = 1.0; // 偏载增大系数（横向偏载对纵向内力的放大）
    };
    std::vector<Lane3D> lanes;

    struct LaneSurface3D {
        std::vector<int> xNodes;
        std::vector<int> yNodes;
        double q = 0.0;
        double p = 0.0;
        int dof = 2;
        double reduction = 1.0;
        int laneCount = 1;
        double impact = 1.0;
        bool autoCode = false;
    };
    LaneSurface3D laneSurface;
};

// 与二维 creepCoefficientAt 的简化曲线相同。t≤t0 时 φ=0。
// 多龄期取相邻阶段的差值，不把从加载日起的整段 φ 一次乘到后来才加上的应力上。
// hMm 是构件自己的理论厚度。模型写了 creepH 时全桥用那一个数。
inline double creepCoefficientAt3D(const Model3D& m, double timeDays, double hMm = 0.0) {
    const double h = m.creepH > 0.0 ? m.creepH : hMm;
    double phi0 = m.creepPhi0;
    double betaH = m.creepBetaH;
    double power = m.creepPower;
    if (m.creepRh > 0.0 && h > 0.0 && m.creepFcm > 0.0) {
        const Jtg3362Creep c = jtg3362NominalCreep(
            m.creepRh, h, m.creepFcm, m.creepT0);
        phi0 = c.phi0;
        betaH = c.betaH;
        power = 0.3;
    }
    if (phi0 <= 0.0) return 0.0;
    const double dt = std::max(0.0, timeDays - m.creepT0);
    const double den = betaH + dt;
    if (den <= 0.0) return 0.0;
    return phi0 * std::pow(dt / den, power);
}

// ==================== 自由度映射 ====================
//
// 按节点分配自由度数：
//   - 与 BEAM 单元相连的节点：6 个自由度 [u,v,w,θx,θy,θz]
//   - 仅由轴向单元（TRUSS/CABLE/STRUT）连接的节点：3 个平动自由度
//
// 旧实现用"只要存在一个 BEAM 就给所有节点 6 个自由度"的一刀切方式，
// 会让纯桁架节点产生恒为零的转动自由度行/列，导致总体矩阵奇异。
using DofMap3D = DofMap;   // 统一自由度映射见 fem.hpp

inline DofMap3D buildDofMap3D(const Model3D& m) {
    DofMap3D dm;
    std::map<int, int> need;
    for (const auto& [id, x] : m.nodes) need[id] = 3;
    for (const auto& e : m.elements) {
        if (!m.nodes.count(e.n1) || !m.nodes.count(e.n2)) continue;
        if (e.type == "BEAM") {
            need[e.n1] = std::max(need[e.n1], 6);
            need[e.n2] = std::max(need[e.n2], 6);
        } else if (e.type == "BEAM7") {
            need[e.n1] = 7;   // 6 + 翘曲自由度 φ
            need[e.n2] = 7;
        }
    }
    int start = 0;
    for (const auto& [id, x] : m.nodes) {
        dm.index[id] = static_cast<int>(dm.nodeIds.size());
        dm.nodeIds.push_back(id);
        dm.dofStart.push_back(start);
        dm.dofCount.push_back(need[id]);
        start += need[id];
    }
    dm.nDofs = start;
    return dm;
}

// 单元全局自由度编号：梁 12 个、BEAM7 14 个、轴向单元 6 个（各取前 3 个平动自由度）。
inline std::vector<int> elementGdofs3D(const DofMap3D& dm, const Element3D& e) {
    std::vector<int> g;
    if (e.type == "BEAM") {
        g.resize(12);
        for (int k = 0; k < 6; ++k) {
            g[k]     = dm.dof(e.n1, k);
            g[6 + k] = dm.dof(e.n2, k);
        }
    } else if (e.type == "BEAM7") {
        g.resize(14);
        for (int k = 0; k < 7; ++k) {
            g[k]      = dm.dof(e.n1, k);
            g[7 + k]  = dm.dof(e.n2, k);
        }
    } else {
        g.resize(6);
        for (int k = 0; k < 3; ++k) {
            g[k]     = dm.dof(e.n1, k);
            g[3 + k] = dm.dof(e.n2, k);
        }
    }
    return g;
}

inline std::vector<double> toStdVec(const Eigen::VectorXd& v) {
    return std::vector<double>(v.data(), v.data() + v.size());
}

// ==================== 装配与求解 ====================

// 三维线性系统的装配结果。影响线/影响面扫描可复用同一份 K 的分解。
struct Assembled3D {
    DofMap3D dm;
    SparseMat3D K0;                 // 未施加约束（可用于反力 K0*U - F0）
    Eigen::VectorXd F0;
    SparseMat3D K;                  // 已施加 Dirichlet，可直接分解
    Eigen::VectorXd F;              // 已施加 Dirichlet
    std::vector<int> cdofs;
    std::vector<double> cvals;
    std::vector<char> isConstrained;
};

// 稀疏直接法求解。首选 SimplicialLDLT（对称正定），失败时回退 SparseLU。
// 旧实现把稀疏矩阵转成 Eigen::MatrixXd 后用 FullPivLU 稠密求解，
// 内存 O(n²)、时间 O(n³)：6000 自由度需 288 MB，是三维分析的规模天花板。
inline Eigen::VectorXd solveSparse3D(const SparseMat3D& K, const Eigen::VectorXd& F) {
    const double fNorm = std::max(1.0, F.norm());

    Eigen::VectorXd U;
    Eigen::SimplicialLDLT<SparseMat3D> ldlt;
    ldlt.compute(K);
    if (ldlt.info() == Eigen::Success) {
        U = ldlt.solve(F);
        if (ldlt.info() != Eigen::Success || !U.allFinite()) U.resize(0);
    }
    if (U.size() == 0) {
        Eigen::SparseLU<SparseMat3D, Eigen::COLAMDOrdering<int>> lu;
        lu.compute(K);
        if (lu.info() != Eigen::Success)
            throw std::runtime_error("三维求解失败：矩阵奇异，请检查约束是否充分");
        U = lu.solve(F);
        if (lu.info() != Eigen::Success || !U.allFinite())
            throw std::runtime_error("三维求解失败：矩阵奇异，请检查约束是否充分");
    }

    // 残差分两级：明显错误直接抛异常；中等残差只警告——
    // 细长梁单元链的条件数按 ~N⁴ 增长，残差 1e-5 量级属于模型本身病态，
    // 并不是求解错误，不应因此中断计算。
    const double rel = (K * U - F).norm() / fNorm;
    if (rel > 1e-3)
        throw std::runtime_error("三维求解：相对残差过大 (" + std::to_string(rel) +
                                 ")，矩阵可能奇异或约束不足");
    if (rel > 1e-6)
        std::cerr << "[警告] 三维求解相对残差 " << rel
                  << "：模型可能病态（细长梁单元过多或刚度比悬殊）\n";
    return U;
}

// 可复用的稀疏分解：一次分解，多次回代。
//
// 默认 SimplicialLDLT，失败自动回退 SparseLU（保持旧行为）。
// 显式传 SolverType::CHOLMOD 时走 SuiteSparse——fem.hpp 里那套
// analyzePattern/factorize 分离的 CHOLMOD 封装此前只在 example_cholmod.cpp
// 里演示过一次，从来没有接到真正的求解路径上（P1-7）。CHOLMOD 不可用
// （未装 SuiteSparse）或矩阵不适用时，这里会静默回退到 LDLT/ LU。
class Factorization3D {
public:
    explicit Factorization3D(const SparseMat3D& K,
                             SolverType type = SolverType::SimplicialLDLT) {
        compute(K, type);
    }

    void compute(const SparseMat3D& K,
                 SolverType type = SolverType::SimplicialLDLT) {
        n_ = static_cast<int>(K.rows());
        type_ = type;
        useCholmod_ = false;
        useLu_ = false;
        ok_ = false;

        if (type == SolverType::CHOLMOD) {
            try {
                cholmod_.analyze(K);
                cholmod_.factorize(K);
                useCholmod_ = true;
                ok_ = true;
                return;
            } catch (const std::exception&) {
                useCholmod_ = false;   // 未启用或非正定 → 回退
            }
        }

        ldlt_.compute(K);
        if (ldlt_.info() == Eigen::Success) {
            ok_ = true;
            return;
        }
        useLu_ = true;
        lu_.compute(K);
        ok_ = (lu_.info() == Eigen::Success);
    }

    // 结构不变、只有数值变化时只做数值分解，省掉重复的符号分解与排序。
    void refactorize(const SparseMat3D& K) {
        if (!ok_) { compute(K, type_); return; }
        if (useCholmod_) { cholmod_.factorize(K); return; }
        if (useLu_) {
            lu_.factorize(K);
            ok_ = (lu_.info() == Eigen::Success);
            return;
        }
        ldlt_.factorize(K);
        if (ldlt_.info() != Eigen::Success) {
            useLu_ = true;
            lu_.compute(K);
            ok_ = (lu_.info() == Eigen::Success);
        }
    }

    SolverType type() const { return type_; }
    bool usingCholmod() const { return useCholmod_; }

    Eigen::VectorXd solve(const Eigen::VectorXd& b) const {
        if (!ok_) throw std::runtime_error("三维求解：矩阵分解失败");
        // 三种分解的 solve() 返回不同的表达式类型，不能放在同一个 ?: 里。
        Eigen::VectorXd x;
        if (useCholmod_) x = cholmod_.solve(b);
        else if (useLu_) x = lu_.solve(b);
        else             x = ldlt_.solve(b);
        if (!x.allFinite())
            throw std::runtime_error("三维求解失败：结果含 NaN/Inf，矩阵可能奇异");
        return x;
    }

    // 多右端项：B 的每一列一个右端项，返回同尺寸解矩阵。
    // CHOLMOD 原生支持多右端项；Eigen 后端逐列回代。
    Eigen::MatrixXd solveMulti(const Eigen::MatrixXd& B) const {
        if (useCholmod_) return cholmod_.solveMulti(B);
        Eigen::MatrixXd X(n_, B.cols());
        for (int c = 0; c < static_cast<int>(B.cols()); ++c)
            X.col(c) = solve(Eigen::VectorXd(B.col(c)));
        return X;
    }

private:
    Eigen::SimplicialLDLT<SparseMat3D> ldlt_;
    Eigen::SparseLU<SparseMat3D, Eigen::COLAMDOrdering<int>> lu_;
    ReusableSolver cholmod_{SolverType::CHOLMOD};
    SolverType type_ = SolverType::SimplicialLDLT;
    bool useCholmod_ = false;
    bool useLu_ = false;
    bool ok_ = false;
    int n_ = 0;
};

inline Assembled3D assembleModel3D(const Model3D& m) {
    Assembled3D a;
    a.dm = buildDofMap3D(m);
    const int nd = a.dm.nDofs;

    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        if (!m.nodes.count(e.n1) || !m.nodes.count(e.n2))
            throw std::runtime_error("单元 " + std::to_string(e.id) + " 引用了不存在的节点");
        const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
        const std::vector<int> g = elementGdofs3D(a.dm, e);

        if (e.type == "BEAM") {
            const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
            const double L = (p2 - p1).norm();
            const auto Kl = beam3DLocalStiffness(Eeff, e.A, e.Iy, e.Iz, e.J, e.G, L);
            const auto T = beam3DTransformation(p1, p2);
            const auto Ke = T.transpose() * Kl * T;
            // 坐标变换会引入 1e-17 量级的数值噪声，按相对阈值过滤以控制 nnz。
            const double scale = std::max(1.0, Ke.cwiseAbs().maxCoeff());
            for (int r = 0; r < 12; ++r)
                for (int c = 0; c < 12; ++c)
                    if (std::abs(Ke(r, c)) > 1e-12 * scale)
                        trips.emplace_back(g[r], g[c], Ke(r, c));
        } else if (e.type == "BEAM7") {
            const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
            const double L = (p2 - p1).norm();
            const auto Kl = beam7LocalStiffness(Eeff, e.A, e.Iy, e.Iz, e.J, e.G, e.Iw, L);
            const auto T = beam7Transformation(p1, p2);
            const auto Ke = T.transpose() * Kl * T;
            const double scale = std::max(1.0, Ke.cwiseAbs().maxCoeff());
            for (int r = 0; r < 14; ++r)
                for (int c = 0; c < 14; ++c)
                    if (std::abs(Ke(r, c)) > 1e-12 * scale)
                        trips.emplace_back(g[r], g[c], Ke(r, c));
        } else {
            const auto Ke = truss3DStiffness(Eeff, e.A, m.nodes.at(e.n1), m.nodes.at(e.n2));
            const double scale = std::max(1.0, Ke.cwiseAbs().maxCoeff());
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c)
                    if (std::abs(Ke(r, c)) > 1e-12 * scale)
                        trips.emplace_back(g[r], g[c], Ke(r, c));
        }
    }

    SparseMat3D K(nd, nd);
    K.setFromTriplets(trips.begin(), trips.end());
    K.makeCompressed();
    a.K0 = K;

    Eigen::VectorXd F = Eigen::VectorXd::Zero(nd);
    for (const auto& [node, dof, val] : m.loads) {
        const int d = a.dm.safeDof(node, dof - 1);
        if (d >= 0) F(d) += val;
    }

    // 初应变 / 收缩等效节点荷载
    for (const auto& e : m.elements) {
        const double ecs = m.shrinkageEcs + e.ecs;
        if (e.istrn == 0.0 && ecs == 0.0) continue;
        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
        const double L = std::max(d.norm(), 1e-12);
        const Vec3 c = d / L;
        const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
        // 正的收缩应变与正的初应变一样，把两端拉向对方（缩短）。
        const double P = Eeff * e.A * (e.istrn + ecs);
        for (int k = 0; k < 3; ++k) {
            const int d1 = a.dm.safeDof(e.n1, k);
            const int d2 = a.dm.safeDof(e.n2, k);
            if (d1 >= 0) F(d1) += P * c(k);
            if (d2 >= 0) F(d2) -= P * c(k);
        }
    }
    a.F0 = F;

    // 弹性支承（k>0）：刚度对角加 k，等效荷载加 k·value；不作为硬约束。
    for (const auto& c : m.constraints) {
        if (c.k <= 0.0) continue;
        const int d = a.dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;
        K.coeffRef(d, d) += c.k;
        a.K0.coeffRef(d, d) += c.k;
        F(d) += c.k * c.value;
        a.F0(d) += c.k * c.value;
    }

    // 约束：支持非零指定位移（支座沉降）。
    a.isConstrained.assign(nd, 0);
    for (const auto& c : m.constraints) {
        const int d = a.dm.safeDof(c.node, c.dof - 1);
        if (d < 0) continue;   // 纯轴向节点没有该转动自由度，约束无意义
        if (c.k > 0.0) continue;   // 弹性支承不作硬约束
        a.isConstrained[d] = 1;
        const auto it = std::find(a.cdofs.begin(), a.cdofs.end(), d);
        if (it == a.cdofs.end()) {
            a.cdofs.push_back(d);
            a.cvals.push_back(c.value);
        } else {
            a.cvals[std::distance(a.cdofs.begin(), it)] = c.value;
        }
    }

    // 诊断：未约束且未被任何单元连接的自由度（零对角线）。
    for (int d = 0; d < nd; ++d) {
        if (a.isConstrained[d]) continue;
        if (K.coeff(d, d) == 0.0) {
            int owner = -1;
            for (std::size_t i = 0; i < a.dm.nodeIds.size(); ++i) {
                if (d >= a.dm.dofStart[i] && d < a.dm.dofStart[i] + a.dm.dofCount[i]) {
                    owner = a.dm.nodeIds[i];
                    break;
                }
            }
            throw std::runtime_error(
                "三维求解：零刚度自由度——节点 " + std::to_string(owner) +
                " 的第 " + std::to_string(d - a.dm.dofStart[a.dm.index.at(owner)] + 1) +
                " 个自由度未约束且未被任何单元连接");
        }
    }

    a.K = K;
    a.F = F;
    applyDirichlet(a.K, a.F, a.cdofs, a.cvals);
    return a;
}

inline std::vector<double> solveModel3D(const Model3D& m) {
    const Assembled3D a = assembleModel3D(m);
    return toStdVec(solveSparse3D(a.K, a.F));
}

// ==================== 3D 模态 / 屈曲 ====================

struct ModalResult3D {
    std::vector<double> frequencies;              // Hz
    std::vector<std::vector<double>> modes;       // 每个模态一个全局位移向量
};

// 3D 桁架一致质量矩阵（6 自由度 [u1,v1,w1,u2,v2,w2]），与 2D 约定一致。
inline Eigen::Matrix<double, 6, 6> truss3DConsistentMassMatrix(double rhoA, double L) {
    const double m = rhoA * L;
    const double c = m / 6.0;
    Eigen::Matrix<double, 6, 6> Me = Eigen::Matrix<double, 6, 6>::Zero();
    for (int d : {0, 1, 2}) {
        Me(d, d) = 2.0 * c; Me(d, d + 3) = c;
        Me(d + 3, d) = c; Me(d + 3, d + 3) = 2.0 * c;
    }
    return Me;
}

// 未施加 Dirichlet 的刚度矩阵（用于模态/屈曲），与 assembleModel3D 的 K0 一致，
// 但不做零刚度诊断——模态分析自行丢弃零刚度/无质量自由度。
inline SparseMat3D assembleStiffness3D(const Model3D& m, DofMap3D& dm) {
    dm = buildDofMap3D(m);
    const int nd = dm.nDofs;
    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
        const std::vector<int> g = elementGdofs3D(dm, e);
        if (e.type == "BEAM") {
            BeamElement b;
            b.E = Eeff; b.A = e.A; b.Iy = e.Iy; b.Iz = e.Iz; b.J = e.J; b.G = e.G;
            b.x0 = m.nodes.at(e.n1); b.x1 = m.nodes.at(e.n2);
            const DenseMat Ke = b.stiffness();
            const double scale = std::max(1.0, Ke.cwiseAbs().maxCoeff());
            for (int r = 0; r < 12; ++r)
                for (int c = 0; c < 12; ++c)
                    if (std::abs(Ke(r, c)) > 1e-12 * scale)
                        trips.emplace_back(g[r], g[c], Ke(r, c));
        } else {
            const DenseMat Ke = truss3DStiffness(Eeff, e.A, m.nodes.at(e.n1), m.nodes.at(e.n2));
            const double scale = std::max(1.0, Ke.cwiseAbs().maxCoeff());
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c)
                    if (std::abs(Ke(r, c)) > 1e-12 * scale)
                        trips.emplace_back(g[r], g[c], Ke(r, c));
        }
    }
    // 弹性支承留在刚度里。硬约束不进这张矩阵，由自由自由度筛选去掉。
    for (const auto& c : m.constraints) {
        if (c.k <= 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0) trips.emplace_back(d, d, c.k);
    }
    SparseMat3D K(nd, nd);
    K.setFromTriplets(trips.begin(), trips.end());
    K.makeCompressed();
    return K;
}

// 未施加 Dirichlet 的一致质量矩阵。线质量密度 = dens·A + addmass。
inline SparseMat3D assembleMassMatrix3D(const Model3D& m, const DofMap3D& dm) {
    const int nd = dm.nDofs;
    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        const double rhoA = e.dens * e.A + e.addmass;
        if (rhoA <= 0.0) continue;
        const std::vector<int> g = elementGdofs3D(dm, e);
        if (e.type == "BEAM") {
            BeamElement b;
            b.E = e.E; b.A = e.A; b.Iy = e.Iy; b.Iz = e.Iz; b.J = e.J; b.G = e.G;
            b.x0 = m.nodes.at(e.n1); b.x1 = m.nodes.at(e.n2);
            const DenseMat T = b.transformation();
            const DenseMat Me = T.transpose() * b.consistentMassMatrix(rhoA) * T;
            for (int r = 0; r < 12; ++r)
                for (int c = 0; c < 12; ++c)
                    if (Me(r, c) != 0.0) trips.emplace_back(g[r], g[c], Me(r, c));
        } else {
            const double L = (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm();
            const DenseMat Me = truss3DConsistentMassMatrix(rhoA, L);
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c)
                    if (Me(r, c) != 0.0) trips.emplace_back(g[r], g[c], Me(r, c));
        }
    }
    SparseMat3D M(nd, nd);
    M.setFromTriplets(trips.begin(), trips.end());
    M.makeCompressed();
    return M;
}

// 自由自由度：去掉约束 + 零刚度（+ 可选零质量）自由度。
inline std::vector<int> freeDofs3D(const Model3D& m, const DofMap3D& dm,
                                   const SparseMat3D& K, const SparseMat3D* M = nullptr) {
    const int nd = dm.nDofs;
    std::vector<char> isC(nd, 0);
    for (const auto& c : m.constraints) {
        if (c.k > 0.0) continue;   // 弹性支承不是固结，刚度已在 K 里
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0) isC[d] = 1;
    }
    std::vector<int> free;
    for (int d = 0; d < nd; ++d) {
        if (isC[d]) continue;
        if (K.coeff(d, d) == 0.0) continue;
        if (M && M->coeff(d, d) == 0.0) continue;
        free.push_back(d);
    }
    return free;
}

// Kt 非空时用这套刚度（成桥态切线，含 Ernst 与几何刚度），不再线弹性装配。
// 切线的索块可能因 ∂T/∂Lh 略不对称，稠密特征问题取对称部分。
inline ModalResult3D computeNaturalFrequencies3D(const Model3D& m, int nModes = 6,
                                                 const SparseMat3D* Kt = nullptr) {
    DofMap3D dm;
    SparseMat3D Kown;
    const SparseMat3D* Kp = Kt;
    if (!Kp) {
        Kown = assembleStiffness3D(m, dm);
        Kp = &Kown;
    } else {
        dm = buildDofMap3D(m);
    }
    const SparseMat3D& K = *Kp;
    const SparseMat3D M = assembleMassMatrix3D(m, dm);
    if (M.nonZeros() == 0)
        throw std::runtime_error("模态分析：全部单元质量为零，请在 *MP 给 DENS（或 *R 给 ADDMASS）");
    const std::vector<int> free = freeDofs3D(m, dm, K, &M);
    const int nf = static_cast<int>(free.size());
    ModalResult3D out;
    if (nf == 0) return out;
    const int n = std::min(nModes, nf);
    constexpr double pi = 3.14159265358979323846;

#ifdef FEM_HAS_SPECTRA
    const SparseModalResult sm = solveSparseModal(K, M, free, dm.nDofs, n);
    if (!sm.eigenvalues.empty()) {
        for (int i = 0; i < static_cast<int>(sm.eigenvalues.size()); ++i) {
            const double lam = sm.eigenvalues[i];
            const double f = (lam > 0.0) ? std::sqrt(lam) / (2.0 * pi) : 0.0;
            out.frequencies.push_back(f);
            out.modes.push_back(sm.modes[i]);
        }
        return out;
    }
    // 回退稠密（Spectra 未收敛或自由度太少）
#endif
    Eigen::MatrixXd Kd(nf, nf), Md(nf, nf);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            Kd(a, b) = 0.5 * (K.coeff(free[a], free[b]) + K.coeff(free[b], free[a]));
            Md(a, b) = M.coeff(free[a], free[b]);
        }
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Kd, Md);
    for (int i = 0; i < n; ++i) {
        const double lam = ges.eigenvalues()(i);
        const double f = (lam > 0.0) ? std::sqrt(lam) / (2.0 * pi) : 0.0;
        out.frequencies.push_back(f);
        std::vector<double> mode(dm.nDofs, 0.0);
        for (int a = 0; a < nf; ++a) mode[free[a]] = ges.eigenvectors()(a, i);
        out.modes.push_back(std::move(mode));
    }
    return out;
}

// 竖向弯曲基频：在频率升序的模态里，取梁节点竖向位移占平动能量过半的最低阶。
// 识别不出返回 -1。汽车冲击要用这一阶，不能拿纵向或扭转的最低频率。
inline double verticalBendingFrequency3D(const Model3D& m, const DofMap3D& dm,
                                         const ModalResult3D& modal) {
    std::vector<int> beams;
    for (const auto& e : m.elements) {
        if (e.type != "BEAM" && e.type != "BEAM7") continue;
        beams.push_back(e.n1);
        beams.push_back(e.n2);
    }
    if (beams.empty()) return -1.0;
    for (std::size_t i = 0; i < modal.frequencies.size(); ++i) {
        const double f = modal.frequencies[i];
        if (f <= 0.0) continue;
        const auto& mode = modal.modes[i];
        double w2 = 0.0, uv2 = 0.0;
        for (int nid : beams) {
            const int dw = dm.safeDof(nid, 2);
            const int du = dm.safeDof(nid, 0);
            const int dv = dm.safeDof(nid, 1);
            const double w = (dw >= 0 && dw < static_cast<int>(mode.size())) ? mode[dw] : 0.0;
            const double u = (du >= 0 && du < static_cast<int>(mode.size())) ? mode[du] : 0.0;
            const double v = (dv >= 0 && dv < static_cast<int>(mode.size())) ? mode[dv] : 0.0;
            w2 += w * w;
            uv2 += u * u + v * v;
        }
        if (w2 + uv2 <= 0.0) continue;
        if (w2 / (w2 + uv2) > 0.5) return f;
    }
    return -1.0;
}

// 由参考荷载下的线性位移组装几何刚度 Kg，(K + λ·Kg)·φ = 0 的输入之一。
//
// 两种单元各有各的 Kg，别只装 BEAM：
//   BEAM       —— 一致几何刚度（横向 v/w 与转角的耦合项 + 扭转项），
//                 局部 12×12 经 T 变换到全局；
//   TRUSS/CABLE—— Kg = (N/L)·[[P,-P],[-P,P]]，P = I - n⊗n 是垂直于杆轴
//                 的投影。轴向没有几何刚度（轴向刚度已在 K 里），
//                 只有横向分量随轴力增刚/软化。
// 早期版本只遍历 `e.type == "BEAM"`，于是纯桁架桥的 Kg 恒为零 —— 明明
// 内核有屈曲求解器，桁架却永远报"没有屈曲模态"。现在按实际单元类型补齐。
//
// 轴力 N（受拉为正）统一由参考位移反算。CABLE 在**线性**内核里本来就按
// TRUSS 处理（见 input_format_3d.md「阶段调索」），这里沿用同一物理模型；
// 几何非线性成桥态的切线刚度不走这条路径。BEAM7 的 14 自由度翘曲项没有
// 对应的几何刚度式，显式跳过而不是塞一个 12×12 进去错位。
inline SparseMat3D assembleGeometricStiffness3D(const Model3D& m,
                                                const DofMap3D& dm,
                                                const std::vector<double>& U) {
    auto uget = [&](int nid, int k) -> double {
        const int d = dm.safeDof(nid, k);
        return (d >= 0 && d < static_cast<int>(U.size())) ? U[d] : 0.0;
    };
    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
        const double L = d.norm();
        if (L <= 0.0) continue;
        const std::vector<int> g = elementGdofs3D(dm, e);
        if (e.type == "BEAM") {
            BeamElement b;
            b.E = e.E; b.A = e.A; b.Iy = e.Iy; b.Iz = e.Iz; b.J = e.J; b.G = e.G;
            b.x0 = m.nodes.at(e.n1); b.x1 = m.nodes.at(e.n2);
            const DenseMat T = b.transformation();
            Eigen::VectorXd ug(12);
            for (int k = 0; k < 6; ++k) {
                ug(k)     = uget(e.n1, k);
                ug(6 + k) = uget(e.n2, k);
            }
            const Eigen::VectorXd ul = T * ug;
            const double N = e.E * e.A / L * (ul(6) - ul(0));
            const DenseMat Kg = T.transpose() * b.geometricStiffness(N) * T;
            for (int r = 0; r < 12; ++r)
                for (int c = 0; c < 12; ++c)
                    if (Kg(r, c) != 0.0) trips.emplace_back(g[r], g[c], Kg(r, c));
        } else if (g.size() == 6) {
            const Vec3 n = d / L;
            const Vec3 u1(uget(e.n1, 0), uget(e.n1, 1), uget(e.n1, 2));
            const Vec3 u2(uget(e.n2, 0), uget(e.n2, 1), uget(e.n2, 2));
            const double N = e.E * e.A / L * n.dot(u2 - u1);
            const Mat3 P = Mat3::Identity() - n * n.transpose();
            const Mat3 K = (N / L) * P;
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) {
                    const double v = K(r, c);
                    if (v == 0.0) continue;
                    trips.emplace_back(g[r], g[c], v);
                    trips.emplace_back(g[r], g[c + 3], -v);
                    trips.emplace_back(g[r + 3], g[c], -v);
                    trips.emplace_back(g[r + 3], g[c + 3], v);
                }
        }
    }
    SparseMat3D Kg(dm.nDofs, dm.nDofs);
    Kg.setFromTriplets(trips.begin(), trips.end());
    Kg.makeCompressed();
    return Kg;
}

// 屈曲结果 + **诊断**。λ 为空时不能一律叫"求解失败"：三跨连续梁格在
// 竖向均布荷载下几乎没有轴力，Kg 里没有受压软化项，本来就不存在压屈模态。
// 把 Kg 的非零元数、是否存在受压对角、自由自由度数一并带出去，下游才能
// 区分「物理上无压杆」和「有压杆却没算出来」这两件完全不同的事。
struct BucklingResult3D {
    std::vector<double> factors;
    int kgNonZeros = 0;
    int nFree = 0;
    bool choleskyOk = true;
    bool compression = false;   // Kg 对角存在负值（有单元受压软化）
};

inline BucklingResult3D computeBucklingResult3D(const Model3D& m, int nModes = 3) {
    BucklingResult3D res;
    const std::vector<double> U = solveModel3D(m);   // 参考荷载下的线性位移
    DofMap3D dm;
    const SparseMat3D K = assembleStiffness3D(m, dm);
    const SparseMat3D Ksg = assembleGeometricStiffness3D(m, dm, U);
    res.kgNonZeros = static_cast<int>(Ksg.nonZeros());
    for (int k = 0; k < Ksg.outerSize(); ++k)
        for (SparseMat3D::InnerIterator it(Ksg, k); it; ++it)
            if (it.row() == it.col() && it.value() < -1e-12) res.compression = true;

    const std::vector<int> free = freeDofs3D(m, dm, K);
    res.nFree = static_cast<int>(free.size());
    const int nf = res.nFree;

    // Kg 恒为零（例如连续梁格只受竖向荷载、没有轴力）时 (K+λKg) 恒正定，
    // 不可能有屈曲模态。早退也避免稀疏求解无解后掉进 nf³ 的稠密回退。
    if (res.kgNonZeros == 0) return res;
    if (nf == 0) return res;

#ifdef FEM_HAS_SPECTRA
    const std::vector<double> sp = solveSparseBuckling(K, Ksg, free, dm.nDofs, nModes);
    if (!sp.empty()) {
        res.factors = sp;
        return res;
    }
    // 稀疏求解失败/未收敛时回退稠密
#endif

    Eigen::MatrixXd Kd(nf, nf), Kgd(nf, nf);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            Kd(a, b) = K.coeff(free[a], free[b]);
            Kgd(a, b) = Ksg.coeff(free[a], free[b]);
        }
    // (K + λ·Kg)·φ = 0 → A = L^{-1}·Kg·L^{-T}，λ = -1/μ（μ<0 为屈曲模态）。
    Eigen::LLT<Eigen::MatrixXd> llt(Kd);
    if (llt.info() != Eigen::Success) {
        res.choleskyOk = false;
        return res;
    }
    const Eigen::MatrixXd L = llt.matrixL();
    const Eigen::MatrixXd Linv = L.inverse();
    const Eigen::MatrixXd A = Linv * Kgd * Linv.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(A);
    const int n = std::min(nModes, nf);
    std::vector<double> mu;
    for (int i = 0; i < ses.eigenvalues().size() && static_cast<int>(mu.size()) < n; ++i) {
        const double mval = ses.eigenvalues()(i);
        if (mval < -1e-14 && std::isfinite(mval)) mu.push_back(mval);
    }
    for (double mval : mu) res.factors.push_back(-1.0 / mval);
    return res;
}

// 旧接口保持不变：example_3d_modal / modal_bench 只关心 λ 序列。
inline std::vector<double> computeBucklingFactors3D(const Model3D& m, int nModes = 3) {
    return computeBucklingResult3D(m, nModes).factors;
}

// ==================== 施工阶段 ====================

// 某阶段末的索状态。
//
// 与梁内力相反：**索力是状态量，不是增量**。梁内力在各阶段间累加
// （见 compute3DFinalStageBeamForces），索力不能累加——把各阶段的索力相加
// 得到的是"张拉动作的总和"，不是"桥上有多少索力"。所以这里记的是阶段末的
// 总量，接口上也就没有对应的"索力包络"函数。
struct CableStageState3D {
    int elem = 0;
    double L0 = 0.0;        // 本阶段末的无应力索长（跨阶段继承）
    double chord = 0.0;     // 本阶段末的弦长（含累积变形）
    double T = 0.0;         // 本阶段末的索力（含垂度，与非线性内核同口径）
    double target = 0.0;    // 本阶段的目标索力；0 = 本阶段未调这根索
    double dev = 0.0;       // |T - target| / target
    // 该索**最近一次**被明确给定的目标索力（跨阶段继承），以及当前偏离。
    // 施工阶段真正要回答的问题是"张好之后被后来的索带跑多少"，所以末阶段的
    // 每一根索都要能看到自己相对设计目标的漂移，而不只是本阶段调过的那几根。
    double Tdes = 0.0;
    double devDes = 0.0;
};

struct StageResult3D {
    int stage = 0;
    std::vector<double> U;
    // 本阶段的自由度映射。
    //
    // 阶段 U 是按"该阶段激活单元"的自由度编号来排的：某节点在阶段 1 只连桁架
    // （3 个自由度），阶段 2 接上梁后变 6 个，于是同一个节点在两个阶段里的
    // 全局自由度编号完全不同，总自由度数也不一样。
    // 包络函数若拿全模型的 DofMap3D 去索引阶段 U，轻则读错自由度（静默算错），
    // 重则下标越界（std::vector 断言失败直接 abort —— example_3d_beam_envelope
    // 原来就是这么崩的）。所以每个阶段必须自带自己的映射。
    DofMap3D dm;
    std::vector<CableStageState3D> cables;   // 空 = 本阶段没有索参与受力
    // 本阶段张拉的预应力束的有效预加力（束号 → 有效力 N）。用于前端校核
    // "有效预应力 vs 张拉控制应力"，以及分阶段看"哪批束在哪阶段张拉"。
    std::vector<std::pair<int, double>> prestressForces;
    // LOSSREF 路径的分项损失明细：(束号, σl1~σl6/σpe, σcon)。
    // 只给 LOSS= 或不写损失时为空——那时内核根本没做分项计算，输出段也
    // 不出现，不拿一堆 0 冒充「各项损失为 0」。
    std::vector<std::tuple<int, PrestressLossResult, double>> prestressLossDetails;
    // 本阶段施加的徐变等效杆端力（局部坐标，12 分量）。
    // 位移里已经含徐变变形，再按 K·Δu 取内力会把同一批应力算第二遍。
    // 报告用的截面内力要扣掉这一项：静定、均匀 Δφ 时应力增量应为 0。
    std::map<int, std::array<double, 12>> creepRelax;
};

// 给未参与任何激活单元的节点补齐约束（否则这些节点自由度零刚度）。
inline void constrainDisconnected3D(Model3D& mc) {
    const DofMap3D dm = buildDofMap3D(mc);
    std::map<int, bool> connected;
    for (const auto& e : mc.elements) { connected[e.n1] = true; connected[e.n2] = true; }
    for (const auto& [id, x] : mc.nodes) {
        if (connected[id]) continue;
        for (int k = 0; k < dm.count(id); ++k)
            mc.constraints.push_back({id, k + 1, 0.0});
    }
}

// 逐阶段求解（增量法）。
//
// 语义（这条决定了下面所有事）：**每个阶段解的是"本阶段荷载在当前构型上的
// 增量"**，不是"累计荷载"——`mc.loads = st.loads` 而不是累计。这是本文件
// 一贯的约定，compute3DFinalStageBeamForces 正是按"成桥态 = Σ 各阶段增量"
// 来累加的。谁把 mc.loads 改成累计，那边的求和就立刻双重计数。
//
// 索张拉在这里实现，原因是线性求解器把 CABLE 当桁架（assembleStiffness3D →
// truss3DStiffness），**根本不读 e.L0**：只改 L0 对结果毫无影响。所以"目标
// 索力"必须翻译成一对**初应变等效节点荷载**（把两端沿索向拉近 P），才能真正
// 进入方程。而 P 的增量才是本阶段该加的量：前面阶段已经张拉过的索，若每阶段
// 都把全量 P 再加一遍，就等于反复张拉。
inline Eigen::VectorXd beam3DLocalEndForces(const Model3D& m, const DofMap3D& dm,
                                            const Element3D& e,
                                            const std::vector<double>& U);

inline std::vector<StageResult3D> solveModel3DStages(const Model3D& m,
                                                     int tensionRounds = 6,
                                                     double tensionTol = 1e-4) {
    std::vector<char> active(m.elements.size(), 0);
    std::vector<StageResult3D> out;

    // 累积位移按**节点**存：各阶段的自由度编号不同（见 StageResult3D::dm），
    // 但"节点 n 沿 x 的位移"这个物理量跨阶段连续。索的弦长必须用累积变形后
    // 的位置——否则逐阶段调索会把前几阶段的伸长量当成本阶段新产生的，越调越大。
    std::map<int, std::array<double, 6>> uTotal;
    // 各索的状态，跨阶段继承：
    //   Tcur     当前（上阶段末）索力。张拉增量 ΔP = N目标 − Tcur 必须建立在
    //            这个状态上，不能每次从"参考几何 + L0"重算（理由见下面张拉段）。
    //   L0state  当前无应力索长。
    //   Tdes     最近一次明确给定的目标索力（写输出算漂移用）。
    std::map<int, double> Tcur, L0state, Tdes;
    // 已积累的截面应力合力（局部 12 分量），不是 K·u。
    // 用 std::array：Eigen 定长矩阵的默认构造不填零，map::operator[]
    // 会把未初始化内存加进下一阶段的等效荷载。
    std::map<int, std::array<double, 12>> creepForce;
    // 收缩等效荷载只在单元第一次激活时进增量方程，后面的阶段不再加 ε·L。
    std::set<int> shrinkageDone;

    for (const auto& st : m.stages) {
        for (int id : st.activate) {
            bool found = false;
            for (std::size_t i = 0; i < m.elements.size(); ++i)
                if (m.elements[i].id == id) { active[i] = 1; found = true; break; }
            if (!found) throw std::runtime_error("三维阶段激活了不存在单元");
        }
        for (int id : st.deactivate) {
            bool found = false;
            for (std::size_t i = 0; i < m.elements.size(); ++i)
                if (m.elements[i].id == id) { active[i] = 0; found = true; break; }
            if (!found) throw std::runtime_error("三维阶段钝化了不存在单元");
        }
        Model3D mc;
        mc.nodes = m.nodes;
        mc.constraints = m.constraints;
        mc.loads = st.loads;
        mc.creepPhi = m.creepPhi; mc.creepChi = m.creepChi; mc.shrinkageEcs = 0.0;
        for (std::size_t i = 0; i < m.elements.size(); ++i) {
            if (!active[i]) continue;
            Element3D e = m.elements[i];
            if (m.shrinkageEcs != 0.0 && !shrinkageDone.count(e.id)) {
                e.ecs = m.shrinkageEcs;
                shrinkageDone.insert(e.id);
            }
            mc.elements.push_back(e);
        }
        constrainDisconnected3D(mc);
        const DofMap3D dm = buildDofMap3D(mc);

        // ---- 索的几何工具 ----
        auto elemIn = [&mc](int id) -> const Element3D* {
            for (const auto& e : mc.elements) if (e.id == id) return &e;
            return nullptr;
        };
        auto horiz = [&m](const Element3D& e) {
            const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
            return std::sqrt(d.x() * d.x() + d.y() * d.y());
        };
        // 弦长 = 参考几何 + 已累积位移（+ 本阶段增量，若给了）。
        auto chordOf = [&m](const Element3D& e,
                            const std::map<int, std::array<double, 6>>& u0,
                            const std::vector<double>* dU,
                            const DofMap3D* dmap) {
            Vec3 a = m.nodes.at(e.n1), b = m.nodes.at(e.n2);
            const auto it1 = u0.find(e.n1), it2 = u0.find(e.n2);
            for (int k = 0; k < 3; ++k) {
                if (it1 != u0.end()) a[k] += it1->second[k];
                if (it2 != u0.end()) b[k] += it2->second[k];
                if (dU && dmap) {
                    const int d1 = dmap->safeDof(e.n1, k);
                    const int d2 = dmap->safeDof(e.n2, k);
                    if (d1 >= 0 && d1 < static_cast<int>(dU->size())) a[k] += (*dU)[d1];
                    if (d2 >= 0 && d2 < static_cast<int>(dU->size())) b[k] += (*dU)[d2];
                }
            }
            return (b - a).norm();
        };

        const std::string stag = "*STAGE " + std::to_string(st.id);
        // 本阶段各索要动的力：单元 → 本阶段末的**目标索力增量基准**。
        // tune 标记"需要迭代命中目标"的索（TENSION=）；L0= 是直接给索长，
        // 索力是结果而不是目标，故不进 tune。
        std::map<int, double> Ngoal;
        std::map<int, bool> tune;
        std::map<int, double> tgt;      // 本阶段显式给定的目标索力（写输出用）

        for (const auto& e : mc.elements) {
            if (e.type != "CABLE") continue;
            if (!Tcur.count(e.id)) Tcur[e.id] = 0.0;
            if (!L0state.count(e.id))
                L0state[e.id] = (e.L0 > 0.0
                                     ? e.L0
                                     : (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm());
        }

        // 显式下料长度：直接换 L0，索力按新索长在"当前弦长"上取值。
        for (const auto& [eid, val] : st.unstrained) {
            const Element3D* e = elemIn(eid);
            if (!e)
                throw std::runtime_error(stag + " 的 L0= 指向单元 " + std::to_string(eid)
                                         + "，该单元在本阶段未激活");
            if (e->type != "CABLE")
                throw std::runtime_error(stag + " 的 L0= 指向单元 " + std::to_string(eid)
                                         + "，它不是 CABLE");
            if (val <= 0.0)
                throw std::runtime_error(stag + " 的 L0= 必须为正（单元 "
                                         + std::to_string(eid) + "）");
            L0state[eid] = val;
            Ngoal[eid] = cableTensionFromChord(chordOf(*e, uTotal, nullptr, nullptr),
                                               val, e->E * e->A, e->w, horiz(*e));
        }

        // 目标索力：本阶段末必须正好是这个力，阶段内迭代。
        for (const auto& [eid, val] : st.tensions) {
            const Element3D* e = elemIn(eid);
            if (!e)
                throw std::runtime_error(stag + " 的 TENSION= 指向单元 " + std::to_string(eid)
                                         + "，该单元在本阶段未激活");
            if (e->type != "CABLE")
                throw std::runtime_error(stag + " 的 TENSION= 指向单元 " + std::to_string(eid)
                                         + "，它不是 CABLE（索力只对索单元有意义）");
            if (val <= 0.0)
                throw std::runtime_error(stag + " 的 TENSION= 必须为正（单元 "
                                         + std::to_string(eid)
                                         + "）；要放张请改用 L0= 给无应力索长");
            if (st.unstrained.end() != std::find_if(
                    st.unstrained.begin(), st.unstrained.end(),
                    [eid](const std::pair<int, double>& p) { return p.first == eid; }))
                throw std::runtime_error(stag + " 对单元 " + std::to_string(eid)
                                         + " 同时给了 TENSION= 与 L0=："
                                           "前者要迭代索长、后者直接给索长，语义冲突");
            Ngoal[eid] = val;
            tune[eid] = true;
            tgt[eid] = val;
            Tdes[eid] = val;      // 设计目标跨阶段继承，末阶段的漂移要能看出来
        }

        // ---- 预应力等效荷载：校验一次、预计算，循环内只做叠加 ----
        //
        // 为什么等效荷载而不是力筋单元：本内核没有 RIGID_PIN 刚臂，直线束的
        // 全部效应就是"两端一对轴压 + 两端一对偏心弯矩"；曲线束（s≠0）再叠加
        // 跨内反向均布力。轴向方向取**首尾节点连线**的单位向量。
        //
        // 曲线筋（抛物线，矢高 s 定义在首尾弦中点）的跨内等效荷载分两项，
        // 缺一不可、且必须同时施加才能静力自平衡（竖向合力 = 0）：
        //   1. 锚固端切线竖向力 4Pe·s/L：两端**同向**（s>0 即跨中筋更高时向上）；
        //   2. 竖向均布 q = −8Pe·s/L²：节点化到链上每个节点（首尾各半跨、
        //      中间节点全跨）。
        // 单单元特例下首尾两项恰好抵消；多单元时首尾节点净竖向力 ≠ 0，不能
        // 图省事只给中间节点——那会漏掉锚固竖向力、竖向合力不为 0，凭空多出
        // 支座反力。均布荷载节点化的精度依赖单元细分：单元内弯矩仍线性，跨内
        // 越细越接近抛物线的二次弯矩。
        //
        // 曲线筋要求链**共线**（抛物线定义在首尾弦上，折线束的线形无定义），
        // 不共线直接报错；直线筋（s=0）保持旧行为：按首尾连线近似轴向。
        // 弯矩向量 = Pe·e·(ẑ × c)，ẑ=(0,0,1) 为竖向。
        struct PsEquiv {
            int ni = 0, nj = 0;
            double fx = 0, fy = 0, fz = 0;   // Pe·c（i 端 +，j 端 −）
            double mx = 0, my = 0, mz = 0;   // 偏心弯矩（i 端 +，j 端 −）
            double tz = 0.0;                 // 锚固切线竖向力（i、j 端都 +tz）
            std::vector<std::pair<int, double>> distFz;  // 均布竖向荷载节点化
        };
        std::vector<PsEquiv> psLoads;
        for (const auto& ps : st.prestresses) {
            const std::string ptag = stag + " 的 PRESTRESS= 束 " + std::to_string(ps.id);
            if (ps.lossRef > 0 && !ps.hasLossDetail)
                throw std::runtime_error(ptag
                    + " 指定了 LOSSREF 却没解析分项损失：请用 readModel3D 读入 .dat，"
                      "或直接构造模型时先调用 resolvePrestressLosses3D(m)");
            if (ps.elems.empty())
                throw std::runtime_error(ptag + " 没有 ELEM 列表");
            if (ps.P <= 0.0)
                throw std::runtime_error(ptag + " 的张拉力必须为正");
            if (ps.loss <= 0.0 || ps.loss > 1.0)
                throw std::runtime_error(ptag + " 的 LOSS 必须在 (0,1]");
            // 链校验：都是 BEAM、都在本阶段激活、首尾相接成一条开链。
            std::map<int, int> deg;
            for (int eid : ps.elems) {
                const Element3D* e = elemIn(eid);
                if (!e)
                    throw std::runtime_error(ptag + " 引用了本阶段未激活的单元 "
                                             + std::to_string(eid));
                if (e->type != "BEAM" && e->type != "BEAM7")
                    throw std::runtime_error(ptag + " 只能施加到梁单元（单元 "
                                             + std::to_string(eid) + " 不是梁）");
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
            PsEquiv q;
            q.ni = head; q.nj = tail;
            q.fx = Pe * ax; q.fy = Pe * ay; q.fz = Pe * az;
            // ẑ × c = (0,0,1)×(ax,ay,az) = (−ay, ax, 0)；M = Pe·e·(ẑ × c)
            q.mx = Pe * ps.e * (-ay);
            q.my = Pe * ps.e * ( ax);
            q.mz = 0.0;
            if (ps.s != 0.0) {
                // ---- 曲线筋：先沿链排序节点，再做共线校验 + 均布节点化 ----
                std::map<int, std::vector<int>> adj;
                for (int eid : ps.elems) {
                    const Element3D* e = elemIn(eid);
                    adj[e->n1].push_back(e->n2);
                    adj[e->n2].push_back(e->n1);
                }
                std::vector<int> chain;
                chain.push_back(head);
                int cur = head, prev = -1;
                while (cur != tail) {
                    int nxt = -1;
                    for (int nb : adj[cur]) if (nb != prev) { nxt = nb; break; }
                    if (nxt < 0)
                        throw std::runtime_error(ptag + " 的单元链不连通");
                    chain.push_back(nxt);
                    prev = cur; cur = nxt;
                }
                auto dist = [&](int a, int b) {
                    return (m.nodes.at(b) - m.nodes.at(a)).norm();
                };
                // 共线校验：每个中间节点到首尾弦的垂距 < 1e-6·L。
                for (int nd : chain) {
                    if (nd == head || nd == tail) continue;
                    const Vec3 w = m.nodes.at(nd) - m.nodes.at(head);
                    const double cx = w.y() * dv.z() - w.z() * dv.y();
                    const double cy = w.z() * dv.x() - w.x() * dv.z();
                    const double cz = w.x() * dv.y() - w.y() * dv.x();
                    const double perp = std::sqrt(cx * cx + cy * cy + cz * cz) / L;
                    if (perp > 1e-6 * L)
                        throw std::runtime_error(
                            ptag + " 是曲线筋（S≠0），单元必须共线；节点 "
                            + std::to_string(nd) + " 偏离锚固连线 "
                            + std::to_string(perp) + " m");
                }
                // 锚固端切线竖向力（两端同向，s>0 时向上 +z）
                q.tz = 4.0 * Pe * ps.s / L;
                // 均布 q = −8Pe·s/L²，节点化：首尾各半跨、中间节点全跨。
                const double qz = -8.0 * Pe * ps.s / (L * L);
                for (std::size_t k = 0; k < chain.size(); ++k) {
                    double half = 0.0;
                    if (k == 0)
                        half = dist(chain[0], chain[1]) / 2.0;
                    else if (k + 1 == chain.size())
                        half = dist(chain[k - 1], chain[k]) / 2.0;
                    else
                        half = (dist(chain[k - 1], chain[k])
                                + dist(chain[k], chain[k + 1])) / 2.0;
                    q.distFz.push_back({chain[k], qz * half});
                }
            }
            psLoads.push_back(q);
        }

        // ---- 求解；有目标索力时在阶段内迭代 ----
        //
        // 张拉量必须按**索力状态**算，不能按"参考几何上的初张力"算：
        //     ΔP = N目标 − N当前(本阶段开始时)
        // 早期版本用的是 f(L_geom, L0新) − f(L_geom, L0旧)，在斜拉桥分批张拉
        // 里直接翻车——结构变形 0.5 m 时 f(L_geom, ·) 与真实索力能差近一倍，
        // 实测后期张拉的索把前期张好的索带到 +89%。
        // 迭代修的是"刚度响应"这一项：解出来 Δu 使弦长变化 ext，索力还会
        // 跟着变 EA/L·ext，所以要把这点残差从 ΔP 里扣回去。
        std::map<int, double> corr;
        std::vector<double> U;
        for (int round = 0; round <= tensionRounds; ++round) {
            Model3D mcs = mc;
            mcs.loads = st.loads;
            // 预应力等效荷载：每个 round 都要叠加（mcs 每次都重建，不叠的话
            // 第二次迭代起预应力就"消失"了）。端节点是梁端、必有 6 自由度，
            // 不会落到无效 dof 上；曲线筋的均布节点化会作用到链的中间节点，
            // 那些节点只加竖向力（dof 3），梁节点必有该自由度。
            for (const auto& q : psLoads) {
                mcs.loads.push_back({q.ni, 1,  q.fx});
                mcs.loads.push_back({q.ni, 2,  q.fy});
                mcs.loads.push_back({q.ni, 3,  q.fz});
                mcs.loads.push_back({q.ni, 4,  q.mx});
                mcs.loads.push_back({q.ni, 5,  q.my});
                mcs.loads.push_back({q.ni, 6,  q.mz});
                mcs.loads.push_back({q.nj, 1, -q.fx});
                mcs.loads.push_back({q.nj, 2, -q.fy});
                mcs.loads.push_back({q.nj, 3, -q.fz});
                mcs.loads.push_back({q.nj, 4, -q.mx});
                mcs.loads.push_back({q.nj, 5, -q.my});
                mcs.loads.push_back({q.nj, 6, -q.mz});
                if (q.tz != 0.0) {
                    mcs.loads.push_back({q.ni, 3, q.tz});
                    mcs.loads.push_back({q.nj, 3, q.tz});
                }
                for (const auto& [nd, f] : q.distFz)
                    mcs.loads.push_back({nd, 3, f});
            }
            if (st.creepDphi != 0.0) {
                for (const auto& e : mc.elements) {
                    if (e.type != "BEAM") continue;
                    const auto itc = creepForce.find(e.id);
                    if (itc == creepForce.end()) continue;
                    Eigen::Matrix<double, 12, 1> fl;
                    for (int i = 0; i < 12; ++i)
                        fl(i) = st.creepDphi * itc->second[static_cast<std::size_t>(i)];
                    const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
                    const auto Tg = beam3DTransformation(p1, p2);
                    // fl = Kl·(T·u) 时，T^T·fl = Ke·u，再求解得到 Δu = Δφ·u。
                    const Eigen::VectorXd fe = Tg.transpose() * fl;
                    for (int k = 0; k < 6; ++k) {
                        mcs.loads.push_back({e.n1, k + 1, fe(k)});
                        mcs.loads.push_back({e.n2, k + 1, fe(6 + k)});
                    }
                }
            }
            for (auto& e : mcs.elements) {
                if (e.type != "CABLE") continue;
                const auto it = Ngoal.find(e.id);
                if (it == Ngoal.end()) continue;       // 本阶段不动这根索
                const double dP = it->second + corr[e.id] - Tcur[e.id];
                if (std::abs(dP) <= 1e-9) continue;
                // 组装时用的是徐变折减后的模量：写成 istrn = ΔP/(Eeff·A) 才能
                // 保证等效节点力正好等于 ΔP，不受 creepPhi 影响。
                const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
                e.istrn += dP / (Eeff * e.A);
            }
            U = solveModel3D(mcs);
            if (tgt.empty()) break;              // 未调索：一次求解即结果
            bool needMore = false;
            // 注意：这里必须从 tgt（单元→目标索力）取目标值。曾经写成从
            // tune（单元→是否需要迭代，map<int,bool>）里结构化绑定，Ttar
            // 就成了 bool 的 true=1.0——迭代会把索力精确地拉到 1 N，
            // 而且编译器一声不响。
            for (const auto& [eid, Ttar] : tgt) {
                const Element3D* e = elemIn(eid);
                const double c0 = chordOf(*e, uTotal, nullptr, nullptr);
                const double c1 = chordOf(*e, uTotal, &U, &dm);
                const double Tnow = Tcur[eid] + (Ngoal[eid] + corr[eid] - Tcur[eid])
                                    + e->E * e->A / c0 * (c1 - c0);
                const double err = Tnow - Ttar;
                if (std::abs(err) > tensionTol * Ttar) {
                    corr[eid] -= err;
                    needMore = true;
                }
            }
            if (!needMore || round == tensionRounds) break;
        }

        // ---- 记录本阶段末各索的力（在提交位移之前算，弦长要用 c0/c1）----
        std::map<int, double> Tnew, chordEnd;
        for (const auto& e : mc.elements) {
            if (e.type != "CABLE") continue;
            const double c0 = chordOf(e, uTotal, nullptr, nullptr);
            const double c1 = chordOf(e, uTotal, &U, &dm);
            const double dP = Ngoal.count(e.id) ? (Ngoal[e.id] + corr[e.id] - Tcur[e.id]) : 0.0;
            chordEnd[e.id] = c1;
            Tnew[e.id] = Tcur[e.id] + dP + e.E * e.A / c0 * (c1 - c0);
            // 调过索的，把 L0 落成"与该状态一致"的下料长度：f(chord, L0) = T。
            // 这样下一阶段接手时，两种口径（跟踪力 / 弦长-索长关系）不会打架。
            const auto it = tune.find(e.id);
            if (it != tune.end())
                L0state[e.id] = cableUnstrainedLengthForTension(
                    c1, Tnew[e.id], e.E * e.A, e.w, horiz(e));
        }

        // ---- 提交状态（只在本阶段收敛后写，迭代过程中的中间值不落地）----
        for (const auto& [id, x] : m.nodes) {
            (void)x;
            const int cnt = dm.count(id);
            for (int k = 0; k < cnt; ++k) {
                const int d = dm.safeDof(id, k);
                if (d >= 0 && d < static_cast<int>(U.size())) uTotal[id][k] += U[d];
            }
        }
        for (const auto& [eid, v] : Tnew) Tcur[eid] = v;

        StageResult3D sr;
        sr.stage = st.id;
        sr.dm = dm;
        sr.U = U;
        for (const auto& e : mc.elements) {
            if (e.type != "CABLE") continue;
            CableStageState3D cs;
            cs.elem = e.id;
            cs.L0 = L0state[e.id];
            cs.chord = chordEnd[e.id];
            cs.T = Tnew[e.id];
            const auto it = tgt.find(e.id);
            if (it != tgt.end() && it->second > 0.0) {
                cs.target = it->second;
                cs.dev = std::abs(cs.T - cs.target) / cs.target;
            }
            const auto id2 = Tdes.find(e.id);
            if (id2 != Tdes.end() && id2->second > 0.0) {
                cs.Tdes = id2->second;
                cs.devDes = std::abs(cs.T - cs.Tdes) / cs.Tdes;
            }
            sr.cables.push_back(cs);
        }
        for (const auto& ps : st.prestresses) {
            sr.prestressForces.push_back({ps.id, ps.P * ps.loss});
            if (ps.hasLossDetail) {
                const auto it = m.prestressLosses.find(ps.lossRef);
                const double sc = (it != m.prestressLosses.end()) ? it->second.sigmaCon : 0.0;
                sr.prestressLossDetails.push_back({ps.id, ps.lossDetail, sc});
            }
        }
        for (const auto& e : mc.elements) {
            if (e.type != "BEAM") continue;
            const Eigen::VectorXd fel = beam3DLocalEndForces(m, dm, e, U);
            std::array<double, 12> prev{};
            const auto itc = creepForce.find(e.id);
            if (itc != creepForce.end()) prev = itc->second;
            // 收缩是初应变。自由缩短时 K·Δu 与等效荷载抵消，应力里不留这份力。
            const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
            const double Psh = (e.ecs == 0.0) ? 0.0 : Eeff * e.A * e.ecs;
            std::array<double, 12> stress = prev;
            std::array<double, 12> applied{};
            for (int i = 0; i < 12 && i < static_cast<int>(fel.size()); ++i) {
                const std::size_t k = static_cast<std::size_t>(i);
                applied[k] = st.creepDphi * prev[k];
                if (i == 0) applied[k] += Psh;
                if (i == 6) applied[k] -= Psh;
                // 应力合力增量 = K·Δu − 徐变等效荷载 − 收缩等效荷载。
                stress[k] = prev[k] + fel(i) - applied[k];
            }
            creepForce[e.id] = stress;
            if (st.creepDphi != 0.0 || Psh != 0.0) sr.creepRelax[e.id] = applied;
        }
        out.push_back(sr);
    }
    return out;
}

struct Envelope3D {
    int id = 0;
    double Nmin = 0.0, Nmax = 0.0;
};

inline std::vector<Envelope3D> compute3DStageEnvelope(const Model3D& m,
                                                      const std::vector<StageResult3D>& stages) {
    std::map<int, Envelope3D> env;
    for (const auto& sr : stages) {
        const DofMap3D& dm = sr.dm;   // 必须用本阶段的映射
        for (const auto& e : m.elements) {
            const int d1 = dm.safeDof(e.n1, 0), d2 = dm.safeDof(e.n2, 0);
            if (d1 < 0 || d2 < 0) continue;
            const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
            const double L = std::max(d.norm(), 1e-12);
            const Vec3 c = d / L;
            double dL = 0.0;
            for (int k = 0; k < 3; ++k)
                dL += (sr.U[d2 + k] - sr.U[d1 + k]) * c(k);
            const double N = e.E * e.A / L * dL;
            auto& r = env[e.id]; r.id = e.id;
            r.Nmin = std::min(r.Nmin, N); r.Nmax = std::max(r.Nmax, N);
        }
    }
    std::vector<Envelope3D> out;
    for (const auto& [id, r] : env) out.push_back(r);
    return out;
}

struct Camber3D {
    int node = 0, dof = 0;
    double camber = 0.0;
};

struct DisplacementEnvelope3D {
    int node = 0, dof = 0;
    double minVal = 0.0, maxVal = 0.0;
};

struct BeamEnvelope3D {
    int id = 0;
    double Nmin=0,Nmax=0,VyMin=0,VyMax=0,VzMin=0,VzMax=0;
    double Tmin=0,Tmax=0,MyMin=0,MyMax=0,MzMin=0,MzMax=0;
    double Bmin=0,Bmax=0;   // 双力矩（BEAM7 约束扭转）
    bool valid=false;
};

// 梁单元的局部杆端力：fl = Kl · (T · ug)
inline Eigen::VectorXd beam3DLocalEndForces(const Model3D& m, const DofMap3D& dm,
                                            const Element3D& e,
                                            const std::vector<double>& U) {
    Eigen::VectorXd ug(12);
    for (int k = 0; k < 6; ++k) {
        // 该阶段该节点可能只有 3 个平动自由度（梁尚未激活）→ 转动分量取 0。
        const int d1 = dm.safeDof(e.n1, k), d2 = dm.safeDof(e.n2, k);
        ug(k)     = (d1 >= 0 && d1 < static_cast<int>(U.size())) ? U[d1] : 0.0;
        ug(6 + k) = (d2 >= 0 && d2 < static_cast<int>(U.size())) ? U[d2] : 0.0;
    }
    const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
    const double L = (p2 - p1).norm();
    const auto T = beam3DTransformation(p1, p2);
    const auto Kl = beam3DLocalStiffness(e.E, e.A, e.Iy, e.Iz, e.J, e.G, L);
    return Kl * (T * ug);
}

// 7 自由度梁局部杆端力：fl = Kl · (T · ug)，14 分量含双力矩 B。
inline Eigen::VectorXd beam7LocalEndForces(const Model3D& m, const DofMap3D& dm,
                                           const Element3D& e,
                                           const std::vector<double>& U) {
    Eigen::VectorXd ug(14);
    for (int k = 0; k < 7; ++k) {
        const int d1 = dm.safeDof(e.n1, k), d2 = dm.safeDof(e.n2, k);
        ug(k)     = (d1 >= 0 && d1 < static_cast<int>(U.size())) ? U[d1] : 0.0;
        ug(7 + k) = (d2 >= 0 && d2 < static_cast<int>(U.size())) ? U[d2] : 0.0;
    }
    const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
    const double L = (p2 - p1).norm();
    const auto T = beam7Transformation(p1, p2);
    const auto Kl = beam7LocalStiffness(e.E, e.A, e.Iy, e.Iz, e.J, e.G, e.Iw, L);
    return Kl * (T * ug);
}

// 徐变阶段的 K·Δu 含徐变变形。截面应力要扣掉本阶段等效荷载。
inline Eigen::VectorXd withoutCreepRelax(
    const Eigen::VectorXd& fl,
    const std::map<int, std::array<double, 12>>& relax, int eid) {
    const auto it = relax.find(eid);
    if (it == relax.end()) return fl;
    Eigen::VectorXd out = fl;
    const int n = std::min(12, static_cast<int>(out.size()));
    for (int i = 0; i < n; ++i)
        out(i) -= it->second[static_cast<std::size_t>(i)];
    return out;
}

// 累积一个截面内力样本（首个样本初始化，避免 0 初值吞掉负值极值）。
static inline void accumBeamEnv(BeamEnvelope3D& r,
                                double N, double Vy, double Vz, double T,
                                double My, double Mz, double B = 0.0) {
    if (!r.valid) {
        r.Nmin = r.Nmax = N;
        r.VyMin = r.VyMax = Vy;
        r.VzMin = r.VzMax = Vz;
        r.Tmin = r.Tmax = T;
        r.MyMin = r.MyMax = My;
        r.MzMin = r.MzMax = Mz;
        r.Bmin = r.Bmax = B;
        r.valid = true;
        return;
    }
    r.Nmin = std::min(r.Nmin, N);   r.Nmax = std::max(r.Nmax, N);
    r.VyMin = std::min(r.VyMin, Vy); r.VyMax = std::max(r.VyMax, Vy);
    r.VzMin = std::min(r.VzMin, Vz); r.VzMax = std::max(r.VzMax, Vz);
    r.Tmin = std::min(r.Tmin, T);   r.Tmax = std::max(r.Tmax, T);
    r.MyMin = std::min(r.MyMin, My); r.MyMax = std::max(r.MyMax, My);
    r.MzMin = std::min(r.MzMin, Mz); r.MzMax = std::max(r.MzMax, Mz);
    r.Bmin = std::min(r.Bmin, B);   r.Bmax = std::max(r.Bmax, B);
}

inline std::vector<BeamEnvelope3D> compute3DBeamEnvelope(
    const Model3D& m, const std::vector<StageResult3D>& stages) {
    std::map<int, BeamEnvelope3D> env;
    for (const auto& sr : stages) {
        const DofMap3D& dm = sr.dm;   // 必须用本阶段的映射
        for (const auto& e : m.elements) {
            if (e.type != "BEAM" && e.type != "BEAM7") continue;
            // 本阶段该梁是否激活：两端节点都得有转动自由度（6 自由度）才算。
            if (dm.safeDof(e.n1, 5) < 0 || dm.safeDof(e.n2, 5) < 0) continue;
            auto& r = env[e.id]; r.id = e.id;
            if (e.type == "BEAM7") {
                const Eigen::VectorXd fl = withoutCreepRelax(
                    beam7LocalEndForces(m, dm, e, sr.U), sr.creepRelax, e.id);
                accumBeamEnv(r, -fl(0),  fl(1),  fl(2), -fl(3), -fl(4), -fl(5), -fl(6));   // i 端
                accumBeamEnv(r,  fl(7), -fl(8), -fl(9),  fl(10), fl(11), fl(12), fl(13)); // j 端
            } else {
                const Eigen::VectorXd fl = withoutCreepRelax(
                    beam3DLocalEndForces(m, dm, e, sr.U), sr.creepRelax, e.id);
                // 截面内力（成员端力，两端同一约定）：
                //   N(受拉正) = -fl(0) / +fl(6)
                //   Vy,Vz     = +fl(1) / -fl(7), +fl(2) / -fl(8)
                //   T(扭转)   = -fl(3) / +fl(9)
                //   My,Mz     = -fl(4)/+fl(10), -fl(5)/+fl(11)   (下侧受拉正)
                accumBeamEnv(r, -fl(0),  fl(1),  fl(2), -fl(3), -fl(4), -fl(5));  // i 端
                accumBeamEnv(r,  fl(6), -fl(7), -fl(8),  fl(9),  fl(10), fl(11)); // j 端
            }
        }
    }
    std::vector<BeamEnvelope3D> out;
    for (auto& [id, r] : env) out.push_back(r);
    return out;
}

// 成桥态（最后一个施工阶段）的梁端内力：运营期荷载组合的基准。
// 实现放在 beamEndForcesAt 之后（见前文 LiveBeamEnvelope3D 一节），
// 这里只留一条指向，避免把工具函数提前声明得七零八落。

// 成桥态位移：阶段包络里同样混着半成品构型，运营期挠度组合要用单一状态。
// 成桥态位移：同样是**各阶段累加**（理由见 compute3DFinalStageBeamForces）。
// 每个阶段按自己的自由度映射取值——同一节点在只连桁架的阶段有 3 个自由度、
// 接上梁后变 6 个，两个阶段的全局编号根本不同，拿最后阶段的映射去索引前一
// 阶段的 U 会读错自由度。
inline std::vector<DisplacementEnvelope3D> compute3DFinalStageDisplacement(
    const Model3D& m, const std::vector<StageResult3D>& stages) {
    std::vector<DisplacementEnvelope3D> out;
    if (stages.empty()) return out;
    const DofMap3D& dmLast = stages.back().dm;
    for (const auto& [id, x] : m.nodes) {
        const int nd = dmLast.count(id);
        for (int k = 0; k < nd; ++k) {
            double sum = 0.0;
            for (const auto& s : stages) {
                const int d = s.dm.safeDof(id, k);
                if (d < 0 || d >= static_cast<int>(s.U.size())) continue;
                sum += s.U[d];
            }
            DisplacementEnvelope3D e;
            e.node = id; e.dof = k + 1; e.minVal = sum; e.maxVal = sum;
            out.push_back(e);
        }
    }
    return out;
}

inline std::vector<DisplacementEnvelope3D> compute3DDisplacementEnvelope(
    const Model3D& m, const std::vector<StageResult3D>& stages) {    // 用"最后一个阶段"的映射决定输出哪些自由度；每个阶段再用它自己的映射取值。
    // 未激活阶段的自由度视为 0（该阶段此节点还没参与受力）。
    const DofMap3D dmLast = stages.empty() ? buildDofMap3D(m) : stages.back().dm;
    std::vector<DisplacementEnvelope3D> out;
    for (const auto& [id, x] : m.nodes) {
        const int nd = dmLast.count(id);
        for (int k = 0; k < nd; ++k) {
            double mn = 1e300, mx = -1e300;
            bool seen = false;
            for (const auto& s : stages) {
                const int d = s.dm.safeDof(id, k);
                if (d < 0 || d >= static_cast<int>(s.U.size())) continue;
                const double v = s.U[d];
                if (!seen) { mn = v; mx = v; seen = true; }
                else { mn = std::min(mn, v); mx = std::max(mx, v); }
            }
            if (!seen) { mn = 0.0; mx = 0.0; }
            out.push_back({id, k + 1, mn, mx});
        }
    }
    return out;
}

// ==================== 车道扫描（一次分解，多次回代）====================
//
// 旧实现对每个加载位置都重新装配并重新分解一次总体矩阵，而 K 完全相同、
// 只有右端项相差一个单位荷载。这里装配一次、分解一次，之后只做回代。
template <class F>
inline void scanLane3D(const Model3D& m, const Model3D::Lane3D& lane, F&& f,
                       SolverType type = SolverType::SimplicialLDLT) {
    Model3D bare = m;
    bare.loads.clear();                    // 影响线只反映单位荷载的增量
    const Assembled3D a = assembleModel3D(bare);
    const Factorization3D fac(a.K, type);
    const std::vector<double> Ubase = toStdVec(fac.solve(a.F));

    Eigen::VectorXd Fk(a.F.size());
    for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
        Fk = a.F;
        const int d = a.dm.safeDof(lane.nodes[k], lane.dof - 1);
        // 加载点若被约束，单位荷载直接由支座承担，对结构无影响。
        if (d >= 0 && !a.isConstrained[d]) Fk(d) += 1.0;
        f(k, toStdVec(fac.solve(Fk)), Ubase);
    }
}

inline std::vector<double> compute3DInfluence(const Model3D& m,
                                              const Model3D::Influence3D& spec,
                                              const Model3D::Lane3D& lane) {
    const DofMap3D dm = buildDofMap3D(m);
    const int dResp = dm.safeDof(spec.node, spec.dof - 1);
    if (dResp < 0)
        throw std::runtime_error("三维影响线：响应自由度不存在（节点 " +
                                std::to_string(spec.node) + " dof " +
                                std::to_string(spec.dof) + "）");

    std::vector<double> eta(lane.nodes.size(), 0.0);
    scanLane3D(m, lane,
               [&](std::size_t k, const std::vector<double>& U, const std::vector<double>& Ub) {
                   eta[k] = U[dResp] - Ub[dResp];
               });
    return eta;
}

inline std::vector<double> compute3DBeamAxialInfluence(const Model3D& m,
                                                       int elemId,
                                                       const Model3D::Lane3D& lane) {
    const Element3D* target = nullptr;
    for (const auto& e : m.elements) if (e.id == elemId) { target = &e; break; }
    if (!target) throw std::runtime_error("三维梁单元不存在");

    const DofMap3D dm = buildDofMap3D(m);
    const Vec3 d0 = m.nodes.at(target->n2) - m.nodes.at(target->n1);
    const double L = std::max(d0.norm(), 1e-12);
    const Vec3 c = d0 / L;

    auto axial = [&](const std::vector<double>& U) {
        double dl = 0.0;
        for (int k = 0; k < 3; ++k)
            dl += (U[dm.dof(target->n2, k)] - U[dm.dof(target->n1, k)]) * c(k);
        return target->E * target->A / L * dl;
    };

    std::vector<double> eta(lane.nodes.size(), 0.0);
    scanLane3D(m, lane,
               [&](std::size_t k, const std::vector<double>& U, const std::vector<double>& Ub) {
                   eta[k] = axial(U) - axial(Ub);
               });
    return eta;
}

// 梁端截面内力（12/14 自由度统一）：N 轴力、Vy/Vz 剪力、T 扭矩、
// My/Mz 弯矩、B 双力矩（仅 BEAM7）。
struct BeamEndForces {
    double N = 0.0, Vy = 0.0, Vz = 0.0, T = 0.0, My = 0.0, Mz = 0.0, B = 0.0;
};

// 从局部杆端力向量取出 i 端（end=1）或 j 端（end=2）的截面内力。
inline BeamEndForces beamEndForcesAt(const Eigen::VectorXd& fl, int dofsPerElem, int end) {
    BeamEndForces f;
    if (dofsPerElem == 14) {
        if (end == 1) {
            f.N = -fl(0); f.Vy = fl(1); f.Vz = fl(2); f.T = -fl(3);
            f.My = -fl(4); f.Mz = -fl(5); f.B = -fl(6);
        } else {
            f.N = fl(7); f.Vy = -fl(8); f.Vz = -fl(9); f.T = fl(10);
            f.My = fl(11); f.Mz = fl(12); f.B = fl(13);
        }
    } else {
        if (end == 1) {
            f.N = -fl(0); f.Vy = fl(1); f.Vz = fl(2); f.T = -fl(3);
            f.My = -fl(4); f.Mz = -fl(5);
        } else {
            f.N = fl(6); f.Vy = -fl(7); f.Vz = -fl(8); f.T = fl(9);
            f.My = fl(10); f.Mz = fl(11);
        }
    }
    return f;
}

// 梁单元内力/应力影响线的统一实现：extractor 从梁端截面内力取出关注分量。
inline std::vector<double> compute3DBeamResponseInfluence(
    const Model3D& m, int elemId, int end,
    const std::function<double(const BeamEndForces&, const Element3D&)>& extract,
    const Model3D::Lane3D& lane) {
    const Element3D* target = nullptr;
    for (const auto& e : m.elements) if (e.id == elemId) { target = &e; break; }
    if (!target || (target->type != "BEAM" && target->type != "BEAM7"))
        throw std::runtime_error("三维梁单元不存在: " + std::to_string(elemId));

    const bool is7 = (target->type == "BEAM7");
    const int dpe = is7 ? 14 : 12;
    const DofMap3D dm = buildDofMap3D(m);
    std::vector<double> eta(lane.nodes.size(), 0.0);
    scanLane3D(m, lane,
               [&](std::size_t k, const std::vector<double>& U, const std::vector<double>& Ub) {
                   const Eigen::VectorXd fl  = is7 ? beam7LocalEndForces(m, dm, *target, U)
                                                   : beam3DLocalEndForces(m, dm, *target, U);
                   const Eigen::VectorXd flb = is7 ? beam7LocalEndForces(m, dm, *target, Ub)
                                                   : beam3DLocalEndForces(m, dm, *target, Ub);
                   eta[k] = extract(beamEndForcesAt(fl, dpe, end), *target)
                          - extract(beamEndForcesAt(flb, dpe, end), *target);
               });
    return eta;
}

inline std::vector<double> compute3DBeamMomentInfluence(const Model3D& m, int elemId,
                                                        const std::string& comp, int end,
                                                        const Model3D::Lane3D& lane) {
    if (comp != "MY" && comp != "MZ")
        throw std::runtime_error("未知弯矩分量: " + comp);
    return compute3DBeamResponseInfluence(m, elemId, end,
        [comp](const BeamEndForces& f, const Element3D&) {
            return (comp == "MY") ? f.My : f.Mz;
        }, lane);
}

inline std::vector<double> compute3DBeamShearInfluence(const Model3D& m, int elemId,
                                                       const std::string& comp, int end,
                                                       const Model3D::Lane3D& lane) {
    if (comp != "VY" && comp != "VZ")
        throw std::runtime_error("未知剪力分量: " + comp);
    return compute3DBeamResponseInfluence(m, elemId, end,
        [comp](const BeamEndForces& f, const Element3D&) {
            return (comp == "VY") ? f.Vy : f.Vz;
        }, lane);
}

inline std::vector<double> compute3DBeamBimomentInfluence(const Model3D& m, int elemId, int end,
                                                          const Model3D::Lane3D& lane) {
    return compute3DBeamResponseInfluence(m, elemId, end,
        [](const BeamEndForces& f, const Element3D&) { return f.B; }, lane);
}

inline std::vector<double> compute3DBeamStressInfluence(const Model3D& m, int elemId,
                                                        const std::string& comp, int end,
                                                        const Model3D::Lane3D& lane) {
    return compute3DBeamResponseInfluence(m, elemId, end,
        [comp](const BeamEndForces& f, const Element3D& e) {
            if (comp == "TORSION") {
                if (e.Wt <= 0)
                    throw std::runtime_error("梁扭转剪应力需要 WT（扭转截面模量）");
                return f.T / e.Wt;   // τ = T / Wt
            }
            if (comp == "TOP" || comp == "BOT" || comp == "MAX" || comp == "MIN") {
                if (e.Wtop <= 0 || e.Wbot <= 0)
                    throw std::runtime_error("梁应力需要 WT/WB");
            }
            if (comp == "LEFT" || comp == "RIGHT" || comp == "MAX" || comp == "MIN") {
                if (e.Wy <= 0)
                    throw std::runtime_error("梁水平应力需要 WY");
            }
            if (comp == "TOP")   return f.N / e.A - f.Mz / e.Wtop;
            if (comp == "BOT")   return f.N / e.A + f.Mz / e.Wbot;
            if (comp == "LEFT")  return f.N / e.A - f.My / e.Wy;
            if (comp == "RIGHT") return f.N / e.A + f.My / e.Wy;
            // MAX/MIN：四角点组合应力
            const double s1 = f.N / e.A - f.Mz / e.Wtop - f.My / e.Wy;  // 左上
            const double s2 = f.N / e.A - f.Mz / e.Wtop + f.My / e.Wy;  // 右上
            const double s3 = f.N / e.A + f.Mz / e.Wbot - f.My / e.Wy;  // 左下
            const double s4 = f.N / e.A + f.Mz / e.Wbot + f.My / e.Wy;  // 右下
            if (comp == "MAX") return std::max(std::max(s1, s2), std::max(s3, s4));
            return std::min(std::min(s1, s2), std::min(s3, s4));   // MIN
        }, lane);
}

// ==================== 影响面 ====================

// 由纵向节点 xn 与横向节点 yn 确定网格交点处的节点号。
// 交点定义为坐标为 ( nodes[xn].x, nodes[yn].y ) 的节点。
//
// 旧实现的内层循环完全不使用 yn，"影响面"实际是同一条一维影响线被复制
// N 遍，横向维度从未生效，得到的包络值是错的。
inline int surfaceGridNode(const Model3D& m, int xn, int yn) {
    const auto itx = m.nodes.find(xn);
    const auto ity = m.nodes.find(yn);
    if (itx == m.nodes.end() || ity == m.nodes.end())
        throw std::runtime_error("LANE_SURFACE_3D 引用了不存在的节点");
    const double targetX = itx->second.x();
    const double targetY = ity->second.y();

    double span = 1.0;
    for (const auto& [id, p] : m.nodes)
        span = std::max(span, std::max(std::abs(p.x()), std::abs(p.y())));
    const double tol = 1e-6 * std::max(1.0, span);

    int best = -1;
    double bestD2 = 0.0;
    for (const auto& [id, p] : m.nodes) {
        const double dx = p.x() - targetX;
        const double dy = p.y() - targetY;
        const double d2 = dx * dx + dy * dy;
        if (best < 0 || d2 < bestD2) { bestD2 = d2; best = id; }
    }
    if (best < 0 || std::sqrt(bestD2) > tol)
        throw std::runtime_error(
            "LANE_SURFACE_3D 需要规则的节点网格：找不到坐标 (" +
            std::to_string(targetX) + ", " + std::to_string(targetY) + ") 处的节点");
    return best;
}

inline std::vector<std::tuple<int,int,double>> compute3DInfluenceSurface(
    const Model3D& m, const Model3D::Influence3D& spec,
    const Model3D::LaneSurface3D& surface) {
    const DofMap3D dm = buildDofMap3D(m);
    const int dResp = dm.safeDof(spec.node, spec.dof - 1);
    if (dResp < 0)
        throw std::runtime_error("三维影响面：响应自由度不存在");

    // 预计算每个 (纵向, 横向) 组合对应的实际加载节点。
    struct Cell { int xn, yn, node; };
    std::vector<Cell> cells;
    cells.reserve(surface.xNodes.size() * surface.yNodes.size());
    for (int yn : surface.yNodes)
        for (int xn : surface.xNodes)
            cells.push_back({xn, yn, surfaceGridNode(m, xn, yn)});

    Model3D bare = m;
    bare.loads.clear();
    const Assembled3D a = assembleModel3D(bare);
    const Factorization3D fac(a.K);
    const std::vector<double> Ubase = toStdVec(fac.solve(a.F));
    const double baseResp = Ubase[dResp];

    std::vector<std::tuple<int,int,double>> out;
    out.reserve(cells.size());
    Eigen::VectorXd Fk(a.F.size());
    for (const auto& cell : cells) {
        Fk = a.F;
        const int d = a.dm.safeDof(cell.node, surface.dof - 1);
        if (d >= 0 && !a.isConstrained[d]) Fk(d) += 1.0;
        const std::vector<double> U = toStdVec(fac.solve(Fk));
        out.push_back({cell.xn, cell.yn, U[dResp] - baseResp});
    }
    return out;
}

inline std::pair<double,double> envelope3DFromSurface(
    const Model3D::LaneSurface3D& surface,
    const std::vector<std::tuple<int,int,double>>& vals) {
    // 每个纵向位置取所有横向位置中的最不利值（影响面的标准用法）。
    std::map<int, double> perXMax, perXMin;
    std::map<int, bool> seen;
    for (auto [xn, yn, v] : vals) {
        if (!seen[xn]) {
            perXMax[xn] = v; perXMin[xn] = v; seen[xn] = true;
        } else {
            perXMax[xn] = std::max(perXMax[xn], v);
            perXMin[xn] = std::min(perXMin[xn], v);
        }
    }

    // 旧实现用 std::map 的默认值 0 做初值，负值影响线会被 max(0, v) 吞掉，
    // 负向车道荷载（如 Q=-10）算出的包络恒为 0。这里改用真实极值初始化。
    double pos = 0.0, neg = 0.0;
    double globalMax = -1e300, globalMin = 1e300;
    for (const auto& [xn, mx] : perXMax) {
        const double mn = perXMin[xn];
        const double vsel = (surface.q >= 0.0) ? mx : mn;   // 均布荷载取同向最不利
        if (vsel > 0.0) pos += vsel; else neg += vsel;
        globalMax = std::max(globalMax, mx);
        globalMin = std::min(globalMin, mn);
    }

    const double qMaxArea = (surface.q >= 0.0) ? pos : neg;
    const double qMinArea = (surface.q >= 0.0) ? neg : pos;
    const double pMaxOrd  = (surface.p >= 0.0) ? globalMax : globalMin;
    const double pMinOrd  = (surface.p >= 0.0) ? globalMin : globalMax;

    double factor = surface.reduction;
    // 仅在 reduction 为默认值 1.0 且 laneCount>1 时，按规范自动折减；
    // 用户显式指定了 REDUCTION= 则完全尊重用户输入，不再覆盖。
    if (factor == 1.0 && surface.laneCount > 1) {
        factor = jtgD60LaneReduction(surface.laneCount);
    }

    return {surface.impact * factor * (surface.q * qMaxArea + surface.p * pMaxOrd),
            surface.impact * factor * (surface.q * qMinArea + surface.p * pMinOrd)};
}

// 车道各加载位置的积分权重（梯形法则）。
//
// 为什么必须有这一步：线荷载 q 的效应是 ∫q·η dx，而影响线只在离散节点上
// 有值。直接 Ση 等于假设"每个节点代表单位长度"，量纲少了一个 m——
// q·Ση 得到的是 N 而不是 N·m，数值上偏小一个节点间距（本项目 2.5 m）。
// 实测：中跨跨中单元，活载弯矩 q 项 72 kN·m（错）对 181 kN·m（对）；
// 位移包络同样小 2.5 倍（−2.72 mm 应为 −6.8 mm）。
//
// 集中力 p 只在最不利的一个位置作用，取影响线最大竖标即可，不乘权重。
inline std::vector<double> laneWeights3D(const Model3D& m,
                                         const Model3D::Lane3D& lane) {
    const std::size_t n = lane.nodes.size();
    std::vector<double> w(n, 0.0);
    if (n == 0) return w;
    if (n == 1) { w[0] = 1.0; return w; }
    auto seg = [&](std::size_t i, std::size_t j) {
        const Vec3& a = m.nodes.at(lane.nodes[i]);
        const Vec3& b = m.nodes.at(lane.nodes[j]);
        const Vec3 d = b - a;
        return std::sqrt(d.x() * d.x() + d.y() * d.y() + d.z() * d.z());
    };
    // 端点只取梁内的半段。若再向外镜像半段，总长度变成 L 加两个半跨，
    // 简支跨中弯矩碰巧不受影响（支座处影响线为 0），悬臂固端弯矩会被放大。
    for (std::size_t i = 0; i < n; ++i) {
        const double left = (i > 0) ? 0.5 * seg(i - 1, i) : 0.0;
        const double right = (i + 1 < n) ? 0.5 * seg(i, i + 1) : 0.0;
        w[i] = left + right;
    }
    return w;
}

inline std::pair<double,double> envelope3DFromInfluence(const Model3D& m,
                                                        const Model3D::Lane3D& lane,
                                                        const std::vector<double>& eta) {
    const std::vector<double> w = laneWeights3D(m, lane);
    double pos=0,neg=0,maxv=-1e300,minv=1e300;
    for (std::size_t i=0;i<eta.size();++i) {
        const double wi = (i < w.size()) ? w[i] : 0.0;
        if (eta[i]>0) pos += eta[i]*wi; else neg += eta[i]*wi;
        maxv=std::max(maxv,eta[i]); minv=std::min(minv,eta[i]);
    }
    const double qMaxArea = lane.q>=0 ? pos : neg;
    const double qMinArea = lane.q>=0 ? neg : pos;
    const double pMaxOrd = lane.p>=0 ? maxv : minv;
    const double pMinOrd = lane.p>=0 ? minv : maxv;
    // 横向分布系数 = 车道数 × 多车道横向折减 × 偏载增大。单梁模型算不出
    // 严格横向影响面（横向只有一个节点），规范上用横向分布系数把多车道
    // 并排布载 + 偏载放大折进纵向影响线包络。单车道默认值 (1,1.0) 退化为
    // 旧的 impact·reduction 乘子。
    const double factor = lane.impact * lane.reduction
                        * static_cast<double>(lane.laneCount) * lane.eccentricity;
    return {factor*(lane.q*qMaxArea + lane.p*pMaxOrd),
            factor*(lane.q*qMinArea + lane.p*pMinOrd)};
}

struct StageLiveEnvelope3D {
    int stage = 0;
    double maxVal = 0.0, minVal = 0.0;
};

// IMPACT 缺省（=1）时由竖向基频自动填 1+μ。显式写成非 1 的倍数则尊重手动值。
// 冲击系数只用于**内力/应力包络**（承载力组合）；位移包络不计冲击，见
// compute3DStageLiveEnvelope。两处用同一套 resolveLaneImpact3D 求 f1/μ，系数一致。
struct LaneImpact3D {
    double f1 = 0.0;
    double mu = 0.0;
    double factor = 1.0;
    bool automatic = false;
};

inline LaneImpact3D resolveLaneImpact3D(const Model3D& mc, const SparseMat3D* Kt = nullptr) {
    LaneImpact3D autoImp;
    try {
        const ModalResult3D modal = computeNaturalFrequencies3D(mc, 8, Kt);
        const DofMap3D dmMod = buildDofMap3D(mc);
        double f1 = verticalBendingFrequency3D(mc, dmMod, modal);
        if (f1 < 0.0 && !modal.frequencies.empty()) f1 = modal.frequencies[0];
        if (f1 > 0.0) {
            autoImp.f1 = f1;
            autoImp.mu = jtgD60ImpactCoefficient(f1);
            autoImp.factor = 1.0 + autoImp.mu;
            autoImp.automatic = true;
        }
    } catch (const std::exception&) {
        autoImp.automatic = false;
    }
    return autoImp;
}

inline std::vector<StageLiveEnvelope3D> compute3DStageLiveEnvelope(
    const Model3D& m, const Model3D::Influence3D& spec) {
    std::vector<StageLiveEnvelope3D> out;
    std::vector<char> active(m.elements.size(),0);
    for (const auto& st : m.stages) {
        for (int id : st.activate)
            for (std::size_t i=0;i<m.elements.size();++i) if(m.elements[i].id==id) active[i]=1;
        for (int id : st.deactivate)
            for (std::size_t i=0;i<m.elements.size();++i) if(m.elements[i].id==id) active[i]=0;
        Model3D mc;
        mc.nodes=m.nodes; mc.constraints=m.constraints;
        mc.creepPhi=m.creepPhi; mc.creepChi=m.creepChi; mc.shrinkageEcs=m.shrinkageEcs;
        for (std::size_t i=0;i<m.elements.size();++i) if(active[i]) mc.elements.push_back(m.elements[i]);
        constrainDisconnected3D(mc);

        double maxAll=-1e300,minAll=1e300;
        // 位移包络不计冲击：JTG D60-2015 表 4.3.2 的挠度限值针对汽车荷载
        // 本身（不计冲击），冲击系数只放大内力（承载力组合）。故这里把
        // lane.impact 压成 1.0，横向折减/偏载增大仍保留。基频也不用算。
        for(const auto& lane:m.lanes){
            Model3D::Lane3D laneUse = lane;
            laneUse.impact = 1.0;
            auto eta=compute3DInfluence(mc,spec,lane);
            auto env=envelope3DFromInfluence(mc,laneUse,eta);
            maxAll=std::max(maxAll,env.first); minAll=std::min(minAll,env.second);
        }
        out.push_back({st.id,maxAll,minAll});
    }
    return out;
}

// 成桥态活载**内力**包络：位移包络拿不到承载力组合要用的 M、V。
//
// 只算最后一个阶段。前几个施工阶段的模型是残缺的（未激活跨的节点没有
// 单元、被 constrainDisconnected3D 约束掉），在其上布载得到的内力没有
// 物理意义——实测阶段 1 的跨中弯矩恒为 0，读起来像"活载不产生内力"。
// 运营期组合本来也只用成桥态。
//
// 为什么批量扫而不是调 compute3DBeamMomentInfluence：那个函数每次都对整条
// 车道扫一遍，逐个单元调用等于把同样的 180 次求解重复 222 遍。
// 这里一次扫描、所有单元同时填影响线，总成本与"一个单元"相同。
//
// 为什么两端分开算：梁端弯矩通常异号，若按位置挑"绝对值大的那端"，不同位置
// 会来自不同端，拼出一条非物理的影响线，再套 q·同号面积 + P·最大竖标
// 会把包络显著放大。这里对两端各求一条包络，最后取外包。
// 成桥态（各施工阶段累加）的梁端内力：运营期荷载组合的基准。
//
// 为什么不能拿 compute3DBeamEnvelope 的阶段包络来组合：逐跨施工的包络里
// 含"第 2 跨还没浇、第 1 跨悬臂端"这种临时构型的内力，它与运营期的活载
// 效应属于两个不同构型，极值相加没有物理意义。施工阶段验算用包络，
// 运营期组合用成桥态——两个问题不能共用一个数。
//
// **为什么必须逐阶段累加**（这条曾经错过）：每个阶段是独立求解的
// （solveModel3DStages 里 mc.loads = st.loads），所以 stages.back().U 只是
// "最后一个阶段那批荷载"在当前构型上的解，不含前面阶段的内力。只取最后一
// 阶段等于把先浇好的几跨恒载丢掉。判决性反例：两阶段模型把全部荷载加在
// 第 1 阶段、第 2 阶段不加载，只取最后一阶段会输出**全零**成桥态。
// 累加的正确性来自施工阶段分析的基本假定：已硬化构件的内力冻结、不因后续
// 构件接入而重分配，故成桥态 = Σ_k（第 k 阶段荷载在第 k 阶段构型上的增量）。
struct BeamEndForces3D {
    int elem = 0, end = 1;
    double N = 0.0, Vy = 0.0, Vz = 0.0, T = 0.0, My = 0.0, Mz = 0.0;
};

inline std::vector<BeamEndForces3D> compute3DFinalStageBeamForces(
    const Model3D& m, const std::vector<StageResult3D>& stages) {
    std::vector<BeamEndForces3D> out;
    if (stages.empty()) return out;
    std::map<std::pair<int, int>, BeamEndForces3D> acc;   // (单元, 端) → 累加值
    for (const auto& sr : stages) {
        const DofMap3D& dm = sr.dm;                       // 每阶段必用自己的映射
        for (const auto& e : m.elements) {
            if (e.type != "BEAM" && e.type != "BEAM7") continue;
            // 该阶段此单元尚未接入（节点自由度不存在）→ 本阶段无贡献
            if (dm.safeDof(e.n1, 5) < 0 || dm.safeDof(e.n2, 5) < 0) continue;
            const bool is7 = (e.type == "BEAM7");
            const int dpe = is7 ? 14 : 12;
            const Eigen::VectorXd fl = withoutCreepRelax(
                is7 ? beam7LocalEndForces(m, dm, e, sr.U)
                    : beam3DLocalEndForces(m, dm, e, sr.U),
                sr.creepRelax, e.id);
            for (int end = 1; end <= 2; ++end) {
                const BeamEndForces f = beamEndForcesAt(fl, dpe, end);
                auto& a = acc[{e.id, end}];
                a.elem = e.id; a.end = end;
                a.N += f.N; a.Vy += f.Vy; a.Vz += f.Vz;
                a.T += f.T; a.My += f.My; a.Mz += f.Mz;
            }
        }
    }
    out.reserve(acc.size());
    for (const auto& [key, v] : acc) out.push_back(v);    // std::map 已按 (单元, 端) 有序
    return out;
}

// 为什么按端输出而不是先外包成"单元包络"：同一单元两端弯矩通常异号，
// 若先把两端包成一对 (min,max) 再与恒载相加，会出现"i 端的恒载 + j 端的
// 活载"——两个不同截面的极值相加，既不是 i 端也不是 j 端。按端给，
// 组合才能落在同一个截面上。
struct LiveBeamEnvelope3D {
    int elem = 0;
    int end = 1;                 // 1 = i 端，2 = j 端
    double MyMin = 0.0, MyMax = 0.0, VzMin = 0.0, VzMax = 0.0;
};

// 成桥态构型：按各 *STAGE 的 activate / deactivate 取最后仍在的单元，
// 断开的孤立节点补约束。活载影响线与模态分析都必须用**同一个构型**——
// 拿"某一跨还没浇完"的临时刚度去算自振频率，和拿全桥刚度算，结果不是
// 一回事；汽车冲击系数的 f1_vertical 用的就是这个 mc。
//
// 无施工阶段时全模型即成桥态（调用方自己决定要不要早退）。
inline Model3D finalStageModel3D(const Model3D& m) {
    if (m.stages.empty()) return m;
    std::vector<char> active(m.elements.size(), 0);
    for (const auto& st : m.stages) {
        for (int id : st.activate)
            for (std::size_t i = 0; i < m.elements.size(); ++i)
                if (m.elements[i].id == id) active[i] = 1;
        for (int id : st.deactivate)
            for (std::size_t i = 0; i < m.elements.size(); ++i)
                if (m.elements[i].id == id) active[i] = 0;
    }
    Model3D mc;
    mc.nodes = m.nodes; mc.constraints = m.constraints;
    mc.creepPhi = m.creepPhi; mc.creepChi = m.creepChi;
    mc.shrinkageEcs = m.shrinkageEcs;
    for (std::size_t i = 0; i < m.elements.size(); ++i)
        if (active[i]) mc.elements.push_back(m.elements[i]);
    constrainDisconnected3D(mc);
    return mc;
}

// 线性屈曲的**参考荷载模型**：成桥态构型 + 成桥态恒载合力。
//
// finalStageModel3D 只给构型，不复制荷载（它是给模态用的：模态只需要
// K/M，不需要 F）。屈曲不一样 —— 几何刚度 Kg 由参考荷载下的轴力决定，
// 没有荷载就没有 Kg，也就没有屈曲模态。
//
// 施工阶段是**增量法**（每阶段 `mc.loads = st.loads`，成桥态 = Σ 各阶段
// 增量），所以这里按 (node,dof) 把各阶段增量累加；无阶段时才用全局
// *CLOAD。**不要把 m.loads 与 stages.loads 相加**：有阶段时阶段求解根本
// 不读 m.loads，两者同时给会重复计入，λ 会偏小一倍。
inline Model3D finalStageBucklingModel3D(const Model3D& m) {
    Model3D ref = finalStageModel3D(m);
    // finalStageModel3D 两条分支口径不一致：无阶段时 `return m` 是连 loads
    // 一起拷的，有阶段时只新建构型、不带 loads。不清空的话无阶段模型会
    // 把同一份 *CLOAD 算两遍 —— 实测悬臂柱 λ 正好减半（10.17 vs Euler
    // 20.33），P_ref 翻倍而 P_cr 不变，是个极隐蔽的错。
    ref.loads.clear();
    std::map<std::pair<int, int>, double> acc;
    if (m.stages.empty()) {
        for (const auto& [node, dof, val] : m.loads) acc[{node, dof}] += val;
    } else {
        for (const auto& st : m.stages)
            for (const auto& [node, dof, val] : st.loads) acc[{node, dof}] += val;
    }
    for (const auto& [key, val] : acc)
        if (val != 0.0) ref.loads.push_back({key.first, key.second, val});
    return ref;
}

inline std::vector<LiveBeamEnvelope3D> compute3DLiveBeamEnvelope(
    const Model3D& m, LaneImpact3D* impactOut = nullptr) {
    std::vector<LiveBeamEnvelope3D> out;
    if (m.lanes.empty() || m.stages.empty()) return out;

    const Model3D mc = finalStageModel3D(m);

    std::vector<const Element3D*> tgt;
    for (const auto& e : mc.elements)
        if (e.type == "BEAM" || e.type == "BEAM7") tgt.push_back(&e);
    if (tgt.empty()) return out;
    const DofMap3D dm = buildDofMap3D(mc);

    std::map<std::pair<int,int>, LiveBeamEnvelope3D> acc;   // (elem, end)
    for (const Element3D* e : tgt)
        for (int end = 1; end <= 2; ++end) {
            auto& slot = acc[{e->id, end}];
            slot.elem = e->id; slot.end = end;
        }

    bool wantAuto = false;
    for (const auto& lane : m.lanes)
        if (lane.impact == 1.0) wantAuto = true;
    const LaneImpact3D autoImp = wantAuto ? resolveLaneImpact3D(mc) : LaneImpact3D{};
    if (impactOut) *impactOut = autoImp;

    for (const auto& lane : m.lanes) {
        const std::size_t npos = lane.nodes.size();
        if (npos == 0) continue;
        Model3D::Lane3D laneUse = lane;
        if (lane.impact == 1.0 && autoImp.automatic)
            laneUse.impact = autoImp.factor;
        // 每单元 × 2 端 × 2 分量
        std::vector<std::vector<double>> eta(tgt.size() * 4,
                                             std::vector<double>(npos, 0.0));
        scanLane3D(mc, lane,
            [&](std::size_t k, const std::vector<double>& U,
                const std::vector<double>& Ub) {
                if (k >= npos) return;
                for (std::size_t t = 0; t < tgt.size(); ++t) {
                    const Element3D& e = *tgt[t];
                    const bool is7 = (e.type == "BEAM7");
                    const int dpe = is7 ? 14 : 12;
                    const Eigen::VectorXd fl =
                        is7 ? beam7LocalEndForces(mc, dm, e, U)
                            : beam3DLocalEndForces(mc, dm, e, U);
                    const Eigen::VectorXd fr =
                        is7 ? beam7LocalEndForces(mc, dm, e, Ub)
                            : beam3DLocalEndForces(mc, dm, e, Ub);
                    for (int end = 1; end <= 2; ++end) {
                        const BeamEndForces fi = beamEndForcesAt(fl, dpe, end);
                        const BeamEndForces fb = beamEndForcesAt(fr, dpe, end);
                        eta[t * 4 + (end - 1) * 2 + 0][k] = fi.My - fb.My;
                        eta[t * 4 + (end - 1) * 2 + 1][k] = fi.Vz - fb.Vz;
                    }
                }
            });

        for (std::size_t t = 0; t < tgt.size(); ++t) {
            for (int end = 1; end <= 2; ++end) {
                const auto eMy = envelope3DFromInfluence(mc, laneUse, eta[t * 4 + (end - 1) * 2 + 0]);
                const auto eVz = envelope3DFromInfluence(mc, laneUse, eta[t * 4 + (end - 1) * 2 + 1]);
                auto& slot = acc[{tgt[t]->id, end}];
                slot.MyMin = std::min(slot.MyMin, eMy.second);
                slot.MyMax = std::max(slot.MyMax, eMy.first);
                slot.VzMin = std::min(slot.VzMin, eVz.second);
                slot.VzMax = std::max(slot.VzMax, eVz.first);
            }
        }
    }
    out.reserve(acc.size());
    for (auto& [key, v] : acc) out.push_back(v);
    return out;
}

// ==================== 切线刚度影响线扫描（几何非线性活载）====================
//
// 几何非线性结构（斜拉桥/悬索桥）的活载影响线必须在**成桥态一致切线刚度**上
// 线性化，而不是线弹性装配：张紧索的切线刚度远大于其弹性刚度（Ernst 模量），
// 用 EA/L 会把索的约束算软、活载内力算得严重失真。
//
// 与 scanLane3D 唯一的差别是"刚度从哪来"：后者内部调用 assembleModel3D（线性
// 弹性装配 + 初应变 base 响应），本函数接受外部已分解的切线刚度（NLTangent3D.K
// 经 Factorization3D 分解）。基态恒为零——切线刚度定义在成桥态平衡位置，单位
// 荷载是纯增量，无需像线性装配那样减去初应变/收缩引起的 base 响应。
//
// fac 必须已分解；isConstrained 与 NLTangent3D.cdofs 对应（约束自由度跳过，
// 否则对角置 1 的行会把单位荷载回代成 1 的位移，物理上错）。
template <class F>
inline void scanTangentLane3D(const Factorization3D& fac,
                              const DofMap3D& dm,
                              const std::vector<char>& isConstrained,
                              const std::vector<int>& laneNodes,
                              int laneDof, F&& f) {
    const int nd = dm.nDofs;
    Eigen::VectorXd Fk = Eigen::VectorXd::Zero(nd);
    for (std::size_t k = 0; k < laneNodes.size(); ++k) {
        Fk.setZero();
        const int d = dm.safeDof(laneNodes[k], laneDof - 1);
        // 加载点若被约束，单位荷载直接由支座承担，对结构无影响。
        if (d >= 0 && d < nd && !(d < static_cast<int>(isConstrained.size())
                                   && isConstrained[d])) Fk(d) += 1.0;
        f(k, toStdVec(fac.solve(Fk)));
    }
}

// 切线刚度下的梁端内力影响线：extract 从 BeamEndForces 取出关注分量（N/My/Vz…）。
// 梁端内力增量用 beam3DLocalEndForces（线性弹性刚度）提取——梁的几何刚度
// （P-Δ）是二阶小量，恒载已由非线性求解精确算入，活载增量按一阶提取是标准做法。
inline std::vector<double> tangentBeamResponseInfluence(
    const Model3D& m, const Factorization3D& fac, const DofMap3D& dm,
    const std::vector<char>& isConstrained,
    int elemId, int end,
    const std::function<double(const BeamEndForces&, const Element3D&)>& extract,
    const std::vector<int>& laneNodes, int laneDof) {
    const Element3D* target = nullptr;
    for (const auto& e : m.elements) if (e.id == elemId) { target = &e; break; }
    if (!target || (target->type != "BEAM" && target->type != "BEAM7"))
        throw std::runtime_error("三维梁单元不存在: " + std::to_string(elemId));
    const bool is7 = (target->type == "BEAM7");
    const int dpe = is7 ? 14 : 12;
    std::vector<double> eta(laneNodes.size(), 0.0);
    scanTangentLane3D(fac, dm, isConstrained, laneNodes, laneDof,
        [&](std::size_t k, const std::vector<double>& U) {
            const Eigen::VectorXd fl = is7 ? beam7LocalEndForces(m, dm, *target, U)
                                           : beam3DLocalEndForces(m, dm, *target, U);
            eta[k] = extract(beamEndForcesAt(fl, dpe, end), *target);
        });
    return eta;
}

// 自动识别加劲梁/主梁的纵向链，作为缺省车道。
//
// 判据：BEAM/BEAM7 单元且方向以水平纵向为主（|Δx| >= max(|Δy|,|Δz|)），
// 把桥塔（竖向）与横梁（横向）排除在外。返回按 x 排序的车道节点序列；
// 空 = 识别不出（调用方应明说"未定义车道"，不静默当作无活载）。
inline Model3D::Lane3D autoIdentifyGirderLane(const Model3D& m) {
    Model3D::Lane3D lane;
    std::set<int> usedNodes;
    bool any = false;
    for (const auto& e : m.elements) {
        if (e.type != "BEAM" && e.type != "BEAM7") continue;
        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
        const double ax = std::abs(d.x()), ay = std::abs(d.y()), az = std::abs(d.z());
        if (ax < std::max(ay, az)) continue;   // 竖向塔 / 横向梁
        any = true;
        usedNodes.insert(e.n1);
        usedNodes.insert(e.n2);
    }
    if (!any) return lane;
    // 主梁是单梁线：每个 x 位置一个节点，按 x 排序即得连续链。
    std::map<double, int> byX;
    for (int nd : usedNodes) byX.emplace(m.nodes.at(nd).x(), nd);
    for (const auto& [x, nd] : byX) lane.nodes.push_back(nd);
    lane.dof = 3;   // 竖向加载
    return lane;
}

inline std::vector<Camber3D> compute3DCamber(const Model3D& m,
                                             const std::vector<StageResult3D>& stages) {
    if (stages.empty()) return {};
    const DofMap3D& dm = stages.back().dm;   // 用最后一个阶段自己的映射
    const std::vector<double>& Ulast = stages.back().U;
    std::vector<Camber3D> out;
    for (const auto& [id, x] : m.nodes) {
        for (int k = 0; k < dm.count(id); ++k) {
            const int d = dm.safeDof(id, k);
            out.push_back({id, k + 1, (d >= 0 && d < static_cast<int>(Ulast.size()))
                                          ? -Ulast[d] : 0.0});
        }
    }
    return out;
}

// JTG D60 自动车道荷载：对标记 CODE=JTG_D60 的车道按跨径填入 qk、Pk，
// 并按车道数自动横向折减（仅影响面）。
inline void applyJtgD60LaneLoads(Model3D& m) {
    for (auto& lane : m.lanes) {
        if (!lane.autoCode || lane.nodes.size() < 2) continue;
        const Vec3 a = m.nodes.at(lane.nodes.front());
        const Vec3 b = m.nodes.at(lane.nodes.back());
        const double L = (b - a).norm();
        lane.q = -jtgD60LaneLoadQk();   // 负号 = 沿 dof 负方向（向下）
        lane.p = -jtgD60LaneLoadPk(L);
    }
    if (m.laneSurface.autoCode && m.laneSurface.xNodes.size() >= 2) {
        const Vec3 a = m.nodes.at(m.laneSurface.xNodes.front());
        const Vec3 b = m.nodes.at(m.laneSurface.xNodes.back());
        const double L = (b - a).norm();
        m.laneSurface.q = -jtgD60LaneLoadQk();
        m.laneSurface.p = -jtgD60LaneLoadPk(L);
        if (m.laneSurface.reduction == 1.0 && m.laneSurface.laneCount > 1)
            m.laneSurface.reduction = jtgD60LaneReduction(m.laneSurface.laneCount);
    }
}

// ==================== 输入文件读取 ====================

// ---- 压缩输入支持：受限四则表达式求值器 ----
// 导出的 .dat 允许把规则排列的单元写成循环（见 BridgeStudio 的
// write_self_input_file / 手册 §7 的“紧凑输出”约定）：
//
//     *DO I 9001 9048
//     *EN,I,2001+I-9000,1001+I-9000
//     *ENDDO
//
// 以及把等间距节点写成起止两点 + *FILL：
//
//     1001 0 0 0
//     1051 200 0 0
//     *FILL 1001 1051
//
// 之所以只做“整数 + 单变量 + 四则运算”的受限求值，而不是塞一个通用脚本
// 引擎：这些表达式完全由导出器生成（形态只有 k*I+c 的几种变体），需要的
// 是可读的报错，不是图灵完备。解析失败必须报出原始文本，绝不静默跳过——
// 静默少一个单元，比报错难查一个数量级。
inline double crEvalSum(const std::string& s, std::size_t& pos,
                        const std::string& var, double varVal, bool& ok);
inline double crEvalProduct(const std::string& s, std::size_t& pos,
                            const std::string& var, double varVal, bool& ok);

inline double crEvalAtom(const std::string& s, std::size_t& pos,
                         const std::string& var, double varVal, bool& ok) {
    while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t')) ++pos;
    if (pos >= s.size()) { ok = false; return 0.0; }
    if (s[pos] == '(') {
        ++pos;
        const double v = crEvalSum(s, pos, var, varVal, ok);
        if (!ok) return 0.0;
        while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t')) ++pos;
        if (pos >= s.size() || s[pos] != ')') { ok = false; return 0.0; }
        ++pos;
        return v;
    }
    if (std::isalpha(static_cast<unsigned char>(s[pos]))) {
        const std::size_t b = pos;
        while (pos < s.size() &&
               (std::isalnum(static_cast<unsigned char>(s[pos])) || s[pos] == '_'))
            ++pos;
        if (s.compare(b, pos - b, var) != 0) { ok = false; return 0.0; }
        return varVal;
    }
    char* endp = nullptr;
    const double v = std::strtod(s.c_str() + pos, &endp);
    if (endp == s.c_str() + pos) { ok = false; return 0.0; }
    pos = static_cast<std::size_t>(endp - s.c_str());
    return v;
}

inline double crEvalProduct(const std::string& s, std::size_t& pos,
                            const std::string& var, double varVal, bool& ok) {
    double v = crEvalAtom(s, pos, var, varVal, ok);
    while (ok) {
        while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t')) ++pos;
        if (pos >= s.size() || (s[pos] != '*' && s[pos] != '/')) break;
        const char op = s[pos++];
        const double r = crEvalAtom(s, pos, var, varVal, ok);
        if (!ok) return 0.0;
        if (op == '*') v *= r;
        else {
            if (r == 0.0) { ok = false; return 0.0; }
            v /= r;
        }
    }
    return v;
}

inline double crEvalSum(const std::string& s, std::size_t& pos,
                        const std::string& var, double varVal, bool& ok) {
    double v = crEvalProduct(s, pos, var, varVal, ok);
    while (ok) {
        while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t')) ++pos;
        if (pos >= s.size() || (s[pos] != '+' && s[pos] != '-')) break;
        const char op = s[pos++];
        const double r = crEvalProduct(s, pos, var, varVal, ok);
        if (!ok) return 0.0;
        v = (op == '+') ? v + r : v - r;
    }
    return v;
}

// 单个整数记号：字面量或含循环变量的表达式。结果的整性不强制（导出器保证），
// 但要求整个 token 被完整消费，避免 "1000abc" 这类被截断后当 1000 用。
inline bool crEvalToken(const std::string& tok, const std::string& var,
                        double varVal, long long& out) {
    std::size_t pos = 0;
    bool ok = true;
    const double v = crEvalSum(tok, pos, var, varVal, ok);
    if (!ok) return false;
    while (pos < tok.size() && (tok[pos] == ' ' || tok[pos] == '\t')) ++pos;
    if (pos != tok.size()) return false;
    out = static_cast<long long>(std::llround(v));
    return true;
}

// 把 LOSSREF=<编号> 解析成有效/张拉比 loss，并记录 σl1~σl6 明细。
//
// 整份输入读完后才调用，因此 *PRESTRESS_LOSS 写在 *STAGE 前后都行。
// 幂等：已解析过的束（hasLossDetail）直接跳过，重复调用不会二次改写。
// 只给 LOSS= 或不写损失的束完全不进这条路径，数值逐位不变。
inline void resolvePrestressLosses3D(Model3D& m) {
    for (auto& st : m.stages) {
        for (auto& ps : st.prestresses) {
            if (ps.lossRef <= 0 || ps.hasLossDetail) continue;
            const std::string ptag = "*STAGE 阶段 " + std::to_string(st.id)
                                   + " 的 PRESTRESS 束 " + std::to_string(ps.id);
            if (ps.lossGiven)
                throw std::runtime_error(ptag + " 的 LOSS 与 LOSSREF 只能指定一个");
            const auto it = m.prestressLosses.find(ps.lossRef);
            if (it == m.prestressLosses.end())
                throw std::runtime_error(ptag + " 的 LOSSREF="
                    + std::to_string(ps.lossRef) + " 引用了不存在的 *PRESTRESS_LOSS 编号");
            const PrestressLossInput& li = it->second;
            if (li.sigmaCon <= 0.0)
                throw std::runtime_error("*PRESTRESS_LOSS 编号 "
                    + std::to_string(ps.lossRef) + " 的 SIGMA_CON 必须 > 0");
            const PrestressLossResult d = computePrestressLosses(li);
            const double ratio = d.sigmaPe / li.sigmaCon;
            if (!(ratio > 0.0 && ratio <= 1.0))
                throw std::runtime_error(ptag + " 分项算得 σpe/σcon="
                    + std::to_string(ratio) + " 不在 (0,1]；检查 σl1/σl2/σl4 是否过大，"
                    + "或该损失表是否适用于这束筋");
            ps.loss = ratio;
            ps.lossDetail = d;
            ps.hasLossDetail = true;
        }
    }
}

inline Model3D readModel3D(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开三维输入文件: " + path);
    Model3D m;
    std::map<int, std::string> types;
    std::map<int, double> mats;
    std::map<int, double> areas;
    std::map<int, double> realIy, realIz, realJ, realG, realHth;
    std::map<int, double> realIstrn;
    std::map<int, double> mpDens, realAddmass;
    std::map<int,double> realWtop, realWbot, realWy, realWt, realIw;
    // 索参数（CABLE 用）：W= 单位长度重量 N/m（垂度），L0= 无应力索长 m
    std::map<int,double> realCableW, realCableL0;
    int activeType = 0, activeMat = 0, activeReal = 0;
    // *DO 循环状态（压缩单元定义）。loopBody 缓存 *DO 与 *ENDDO 之间的原始行，
    // 到 *ENDDO 才对每个循环值重放一遍——不这么做就只能展开第一轮。
    bool loopOn = false;
    std::string loopVar = "I";
    long long loopBegin = 0, loopEnd = 0, loopVal = 0;
    std::vector<std::string> loopBody;
    std::string section;
    std::string line;

    // *EN 解析（含 W=/L0= 覆盖）。直接写的单元与 *DO 展开的单元共用这一份逻辑，
    // 否则"循环展开"和"手写"两条路径迟早行为不一致。
    auto parseElementLine = [&](const std::string& rawArgs) {
        std::string args = rawArgs;
        std::replace(args.begin(), args.end(), ',', ' ');
        std::istringstream ss(args);
        std::string tId, t1, t2;
        if (!(ss >> tId >> t1 >> t2)) {
            std::cerr << "[警告] *EN 需要 id,n1,n2: " << rawArgs << "\n";
            return;
        }
        long long id = 0, n1 = 0, n2 = 0;
        const double iv = static_cast<double>(loopVal);
        if (!crEvalToken(tId, loopVar, iv, id) ||
            !crEvalToken(t1, loopVar, iv, n1) ||
            !crEvalToken(t2, loopVar, iv, n2)) {
            std::cerr << "[警告] *EN 的编号无法解析: " << rawArgs << "\n";
            return;
        }
        Element3D e;
        e.id = static_cast<int>(id);
        e.type = types[activeType]; e.E = mats[activeMat]; e.A = areas[activeReal];
        e.n1 = static_cast<int>(n1); e.n2 = static_cast<int>(n2);
        e.Iy = realIy[activeReal]; e.Iz = realIz[activeReal];
        e.J = realJ[activeReal];   e.G = realG[activeReal];
        e.istrn = realIstrn[activeReal];
        e.dens = mpDens[activeMat]; e.addmass = realAddmass[activeReal];
        e.Wtop = realWtop[activeReal]; e.Wbot = realWbot[activeReal];
        e.Wy = realWy[activeReal];     e.Wt = realWt[activeReal];
        e.Iw = realIw[activeReal];
        e.w = realCableW[activeReal];  e.L0 = realCableL0[activeReal];
        e.hTheory = realHth[activeReal];
        // 单元级覆盖：*EN,id,n1,n2 W=<N/m> L0=<m>
        // 每根索的 W/L0 本来就各不相同（下料长度是逐根量），只靠 *R 共享
        // 一份会逼得同类索拆成多个实常数。
        std::string kv;
        while (ss >> kv) {
            if (kv.rfind("W=", 0) == 0) e.w = std::stod(kv.substr(2));
            else if (kv.rfind("L0=", 0) == 0) e.L0 = std::stod(kv.substr(3));
            else std::cerr << "[警告] *EN 未知关键字: " << kv << "\n";
        }
        m.elements.push_back(e);
    };
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        if (line[0] == '*') {
            const std::size_t sep = line.find_first_of(", \t");
            const std::string kw = (sep == std::string::npos) ? line : line.substr(0, sep);
            if (kw == "*DO") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                std::string var;
                long long a = 0, b = 0;
                if (ss >> var >> a >> b) {
                    loopOn = true; loopVar = var; loopBegin = a; loopEnd = b;
                    loopBody.clear();
                } else {
                    std::cerr << "[警告] *DO 语法应为 *DO <变量> <起> <止>: " << line << "\n";
                }
                continue;
            }
            if (kw == "*ENDDO") {
                if (!loopOn) {
                    std::cerr << "[警告] 多余的 *ENDDO（无匹配的 *DO）\n";
                    continue;
                }
                // 重放循环体。*DO 只压缩"编号按等差排列"的规律段，体内常规只有
                // *EN；出现别的关键字一律报错，不猜测语义——猜错就是静默少单元。
                for (long long v = loopBegin; v <= loopEnd; ++v) {
                    loopVal = v;
                    for (const std::string& body : loopBody) {
                        if (body.empty() || body[0] == '#' || body[0] == '!') continue;
                        if (body[0] != '*') continue;
                        const std::size_t p2 = body.find_first_of(", \t");
                        const std::string kw2 =
                            (p2 == std::string::npos) ? body : body.substr(0, p2);
                        if (kw2 == "*EN")
                            parseElementLine(p2 == std::string::npos ? "" : body.substr(p2 + 1));
                        else
                            std::cerr << "[警告] *DO 循环体只支持 *EN，忽略: " << body << "\n";
                    }
                }
                loopOn = false;
                loopBody.clear();
                continue;
            }
            if (loopOn) {                      // 循环体先缓存，*ENDDO 时统一重放
                loopBody.push_back(line);
                continue;
            }
            if (kw == "*FILL") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                int a = 0, b = 0;
                if (!(ss >> a >> b)) {
                    std::cerr << "[警告] *FILL 语法应为 *FILL <首节点> <末节点>: " << line << "\n";
                } else if (a == b) {
                    // 单点，无需填充
                } else {
                    const auto ia = m.nodes.find(a), ib = m.nodes.find(b);
                    if (ia == m.nodes.end() || ib == m.nodes.end()) {
                        std::cerr << "[警告] *FILL 引用了未定义的端点节点: " << line << "\n";
                    } else if (b - a < 2) {
                        std::cerr << "[警告] *FILL 端点号不连续（需 b > a+1）: " << line << "\n";
                    } else {
                        // 等间距线性插值。导出器只在"三点以上共线且等距"时才写
                        // *FILL，故按等差数列补中间点即可，不必再做曲率判断。
                        const Vec3 pa = ia->second, pb = ib->second;
                        const int nseg = b - a;
                        for (int k = 1; k < nseg; ++k) {
                            if (m.nodes.find(a + k) == m.nodes.end())
                                m.nodes.emplace(a + k, pa + (pb - pa) * (static_cast<double>(k) / nseg));
                        }
                    }
                }
                continue;
            }
            if (kw == "*EN") {
                parseElementLine(sep == std::string::npos ? "" : line.substr(sep + 1));
                continue;
            }
            if (kw == "*TYPE" || kw == "*MAT" || kw == "*REAL") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                if (kw == "*TYPE") ss >> activeType;
                else if (kw == "*MAT") ss >> activeMat;
                else ss >> activeReal;
                continue;
            }
            if (kw == "*CREEP" || kw == "*SHRINKAGE") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                std::string kv;
                while (ss >> kv) {
                    if (kw == "*CREEP") {
                        if (kv.rfind("PHI0=",0)==0) m.creepPhi0=std::stod(kv.substr(5));
                        else if (kv.rfind("PHI=",0)==0) m.creepPhi=std::stod(kv.substr(4));
                        else if (kv.rfind("CHI=",0)==0) m.creepChi=std::stod(kv.substr(4));
                        else if (kv.rfind("T0=",0)==0) m.creepT0=std::stod(kv.substr(3));
                        else if (kv.rfind("BETAH=",0)==0) m.creepBetaH=std::stod(kv.substr(6));
                        else if (kv.rfind("POWER=",0)==0) m.creepPower=std::stod(kv.substr(6));
                        else if (kv.rfind("RH=",0)==0) m.creepRh=std::stod(kv.substr(3));
                        else if (kv.rfind("H=",0)==0) m.creepH=std::stod(kv.substr(2));
                        else if (kv.rfind("FCM=",0)==0) m.creepFcm=std::stod(kv.substr(4));
                        else if (kv.rfind("TS=",0)==0) m.shrinkageTs=std::stod(kv.substr(3));
                    } else {
                        if (kv.rfind("ECS=",0)==0) m.shrinkageEcs=std::stod(kv.substr(4));
                    }
                }
                continue;
            }
            if (kw == "*PRESTRESS_LOSS") {
                // 段头统一写成规范名，避免行尾带注释时 section 与比较串不等。
                section = "*PRESTRESS_LOSS";
                continue;
            }
            if (kw == "*INFLUENCE_3D" || kw == "*LANE_3D" || kw == "*LANE_SURFACE_3D" || kw == "*ELEMENT_INFLUENCE_3D") {
                section = line;
                continue;
            }
            section = line;
            continue;
        }
        std::istringstream ss(line);
        if (section == "*NODE") {
            int id; double x,y,z;
            if (ss >> id >> x >> y >> z) m.nodes[id] = Vec3(x,y,z);
        } else if (section == "*ET") {
            int id; std::string t;
            if (ss >> id >> t) types[id] = t;
        } else if (section == "*MP") {
            int id; std::string kv;
            ss >> id;
            double E=0, dens=0;
            while (ss >> kv) {
                if (kv.rfind("E=",0)==0) E=std::stod(kv.substr(2));
                else if (kv.rfind("DENS=",0)==0) dens=std::stod(kv.substr(5));
                else if (kv.rfind("DENSITY=",0)==0) dens=std::stod(kv.substr(8));   // DENSITY 为 DENS 别名
                else std::cerr << "[警告] *MP 未知关键字: " << kv << "\n";
            }
            mats[id]=E; mpDens[id]=dens;
        } else if (section == "*R") {
            int id; std::string kv;
            ss >> id;
            double A=0, Iy=0, Iz=0, J=0, G=0, addmass=0;
            while (ss >> kv) {
                if (kv.rfind("A=",0)==0) A=std::stod(kv.substr(2));
                else if (kv.rfind("IY=",0)==0) Iy=std::stod(kv.substr(3));
                else if (kv.rfind("IZ=",0)==0) Iz=std::stod(kv.substr(3));
                else if (kv.rfind("J=",0)==0) J=std::stod(kv.substr(2));
                else if (kv.rfind("G=",0)==0) G=std::stod(kv.substr(2));
                else if (kv.rfind("ADDMASS=",0)==0) addmass=std::stod(kv.substr(8));
                else if (kv.rfind("ISTRN=",0)==0) realIstrn[id]=std::stod(kv.substr(6));
                else if (kv.rfind("WT=",0)==0) realWtop[id]=std::stod(kv.substr(3));
                else if (kv.rfind("WB=",0)==0) realWbot[id]=std::stod(kv.substr(3));
                else if (kv.rfind("WY=",0)==0) realWy[id]=std::stod(kv.substr(3));
                else if (kv.rfind("WTORS=",0)==0) realWt[id]=std::stod(kv.substr(7));
                else if (kv.rfind("IW=",0)==0) realIw[id]=std::stod(kv.substr(3));
                // 索参数：W 单位线重 N/m（垂度/Ernst 用），L0 无应力索长 m。
                // 注意 "W=" 与 "WT=/WB=/WY=/WTORS=" 前缀不同，靠精确前两位区分。
                else if (kv.rfind("W=",0)==0) realCableW[id]=std::stod(kv.substr(2));
                else if (kv.rfind("L0=",0)==0) realCableL0[id]=std::stod(kv.substr(3));
                else if (kv.rfind("HTH=",0)==0) realHth[id]=std::stod(kv.substr(4));
                else std::cerr << "[警告] *R 未知关键字: " << kv << "\n";
            }
            areas[id]=A;
            realIy[id]=Iy; realIz[id]=Iz; realJ[id]=J; realG[id]=G;
            realAddmass[id]=addmass;
        } else if (section == "*BOUNDARY") {
            // 格式：<节点> <dof 或 a-b> [指定位移值]
            // 第三列可选，缺省 0；非零即支座沉降 / 强制位移。
            // 第四列可选：K=<弹簧刚度>，>0 即弹性支承。
            int node; std::string ds; double val = 0.0, k = 0.0;
            if (ss >> node >> ds) {
                ss >> val;
                std::string extra;
                if (ss >> extra) {
                    if (extra.rfind("K=",0)==0) extra = extra.substr(2);
                    k = std::stod(extra);
                }
                if (ds.find('-') != std::string::npos) {
                    int a=std::stoi(ds.substr(0,ds.find('-'))); int b=std::stoi(ds.substr(ds.find('-')+1));
                    for (int d=a; d<=b; ++d) m.constraints.push_back({node,d,val,k});
                } else m.constraints.push_back({node,std::stoi(ds),val,k});
            }
        } else if (section == "*CLOAD") {
            int node,dof; double val; if (ss >> node >> dof >> val) m.loads.push_back({node,dof,val});
        } else if (section == "*PRESTRESS_LOSS") {
            // 键名映射与 2D 共用 parsePrestressLossLine（唯一来源）。
            const auto [id, li] = parsePrestressLossLine(line);
            m.prestressLosses[id] = li;
        } else if (section == "*STAGE") {
            int id; if (ss >> id) {
                Model3D::Stage3D st; st.id=id;
                std::string kv;
                while (ss >> kv) {
                    if (kv.rfind("ACTIVE=",0)==0) {
                        std::string list=kv.substr(7);
                        std::stringstream ls(list); std::string tok;
                        while(std::getline(ls,tok,',')) st.activate.push_back(std::stoi(tok));
                    } else if (kv.rfind("DEACTIVE=",0)==0) {
                        std::string list=kv.substr(9);
                        std::stringstream ls(list); std::string tok;
                        while(std::getline(ls,tok,',')) st.deactivate.push_back(std::stoi(tok));
                    } else if (kv.rfind("DPHI=",0)==0) {
                        st.creepDphi = std::stod(kv.substr(5));
                    } else if (kv.rfind("TIME=",0)==0) {
                        st.time = std::stod(kv.substr(5));
                    } else if (kv.rfind("TENSION=",0)==0 || kv.rfind("L0=",0)==0) {
                        // TENSION=12:3.2e6,13:3.1e6
                        // L0=12:10.234,13:10.301
                        // 逐对"单元号:数值"。写成逗号分隔是为了让一行能表达
                        // "本阶段同时张拉这一对索"——斜拉桥对称悬臂施工的标准
                        // 动作，也是这个键存在的理由。
                        const bool isL0 = (kv.rfind("L0=",0)==0);
                        const std::string body = kv.substr(isL0 ? 3 : 8);
                        std::stringstream ls(body); std::string pair;
                        while (std::getline(ls, pair, ',')) {
                            if (pair.empty()) continue;
                            const std::size_t c = pair.find(':');
                            if (c == std::string::npos)
                                throw std::runtime_error(
                                    "三维 *STAGE 的 " + std::string(isL0 ? "L0=" : "TENSION=")
                                    + " 需要写成 <单元号>:<数值>，收到「" + pair + "」");
                            const int eid = std::stoi(pair.substr(0, c));
                            const double v = std::stod(pair.substr(c + 1));
                            if (isL0) st.unstrained.push_back({eid, v});
                            else      st.tensions.push_back({eid, v});
                        }
                    } else if (kv.rfind("PRESTRESS=",0)==0) {
                        // 两种写法（同一束只取其一，由值段是否含 '=' 自动判定）：
                        //   位置模式（直线筋，向后兼容）：
                        //     PRESTRESS=<elems>:<P>:<e>[:<loss>]
                        //   关键字模式（曲线筋用；S 为矢高，S=0 亦退化为直线筋）：
                        //     PRESTRESS=<elems>:P=<P>:E=<e>[:S=<矢高>][:LOSS=<loss>|:LOSSREF=<编号>]
                        //   LOSS= 与 LOSSREF= 互斥：前者把总损失打包成一个比值，
                        //   后者引用 *PRESTRESS_LOSS 表按 JTG 3362-2018 §6.2 分项算。
                        //   elems：a-b 范围或逗号列表（连续梁单元，首尾即锚固端）
                        //   P：张拉控制力 N（正）；e：端部偏心 m（形心上方为正）
                        //   S：矢高 = 跨中偏心 − 端部偏心（向上为正；S<0 → 上拱）
                        //   loss：有效/张拉比，可选，缺省 1.0
                        // 一束一个 PRESTRESS=，束号按本阶段内出现顺序自动编。
                        const std::string body = kv.substr(10);
                        const std::size_t c1 = body.find(':');
                        if (c1 == std::string::npos)
                            throw std::runtime_error(
                                "三维 *STAGE 的 PRESTRESS= 需要写成 "
                                "<单元列表>:<P>:<e>[:<loss>] 或 "
                                "<单元列表>:P=...:E=...[:S=...][:LOSS=...]，收到「"
                                + body + "」");
                        const std::string elemsStr = body.substr(0, c1);
                        const std::string valStr = body.substr(c1 + 1);
                        Model3D::Prestress3D ps;
                        ps.id = static_cast<int>(st.prestresses.size()) + 1;
                        if (valStr.find('=') == std::string::npos) {
                            // 位置模式：<P>:<e>[:<loss>]
                            std::vector<std::string> parts;
                            std::stringstream ps2(valStr); std::string t;
                            while (std::getline(ps2, t, ':')) parts.push_back(t);
                            if (parts.size() < 2 || parts.size() > 3)
                                throw std::runtime_error(
                                    "三维 *STAGE 的 PRESTRESS= 位置模式需要 2~3 个数值 "
                                    "(<P>:<e>[:<loss>])，收到「" + valStr + "」");
                            ps.P = std::stod(parts[0]);
                            ps.e = std::stod(parts[1]);
                            ps.loss = (parts.size() == 3) ? std::stod(parts[2]) : 1.0;
                            ps.lossGiven = (parts.size() == 3);
                        } else {
                            // 关键字模式：P=...:E=...[:S=...][:LOSS=... | :LOSSREF=<编号>]
                            std::stringstream vs(valStr); std::string tok;
                            bool hasP = false, hasE = false, hasLossRef = false;
                            while (std::getline(vs, tok, ':')) {
                                if (tok.empty()) continue;
                                const std::size_t eq = tok.find('=');
                                if (eq == std::string::npos)
                                    throw std::runtime_error(
                                        "三维 *STAGE 的 PRESTRESS= 关键字段需写成 "
                                        "KEY=VALUE，收到「" + tok + "」");
                                const std::string key = tok.substr(0, eq);
                                const std::string val = tok.substr(eq + 1);
                                if (key == "LOSSREF") {
                                    ps.lossRef = std::stoi(val);
                                    hasLossRef = true;
                                    continue;
                                }
                                const double v = std::stod(val);
                                if (key == "P")      { ps.P = v; hasP = true; }
                                else if (key == "E") { ps.e = v; hasE = true; }
                                else if (key == "S") { ps.s = v; }
                                else if (key == "LOSS") { ps.loss = v; ps.lossGiven = true; }
                                else throw std::runtime_error(
                                    "三维 *STAGE 的 PRESTRESS= 未知关键字「" + key
                                    + "」（支持 P/E/S/LOSS/LOSSREF）");
                            }
                            if (!hasP || !hasE)
                                throw std::runtime_error(
                                    "三维 *STAGE 的 PRESTRESS= 关键字模式必须含 P= 与 E=");
                            if (ps.lossGiven && hasLossRef)
                                throw std::runtime_error(
                                    "三维 *STAGE 的 PRESTRESS= 的 LOSS 与 LOSSREF "
                                    "只能指定一个（前者是总损失比，后者是分项损失表编号）");
                            if (hasLossRef && ps.lossRef <= 0)
                                throw std::runtime_error(
                                    "三维 *STAGE 的 PRESTRESS= 的 LOSSREF 必须为正整数");
                        }
                        std::stringstream es(elemsStr); std::string tok;
                        while (std::getline(es, tok, ',')) {
                            if (tok.empty()) continue;
                            const std::size_t d = tok.find('-');
                            if (d != std::string::npos) {
                                const int a = std::stoi(tok.substr(0, d));
                                const int b = std::stoi(tok.substr(d + 1));
                                for (int v = a; v <= b; ++v) ps.elems.push_back(v);
                            } else {
                                ps.elems.push_back(std::stoi(tok));
                            }
                        }
                        st.prestresses.push_back(ps);
                    }
                }
                m.stages.push_back(st);
            }
        } else if (section == "*STAGE_LOAD") {
            int stageId; std::string kv; ss >> stageId;
            int node=0,dof=0; double val=0;
            while (ss >> kv) {
                if (kv.rfind("NODE=",0)==0) node=std::stoi(kv.substr(5));
                else if (kv.rfind("DOF=",0)==0) dof=std::stoi(kv.substr(4));
                else if (kv.rfind("VALUE=",0)==0) val=std::stod(kv.substr(6));
            }
            if (stageId>=1 && stageId<=static_cast<int>(m.stages.size()))
                m.stages[stageId-1].loads.push_back({node,dof,val});
        } else if (section == "*INFLUENCE_3D") {
            int node,dof; if (ss >> node >> dof) m.influences.push_back({node,dof});
        } else if (section == "*ELEMENT_INFLUENCE_3D") {
            int elem; if (ss >> elem) m.elementInfluences.push_back({elem});
        } else if (section == "*LANE_3D") {
            std::string kv;
            std::string nodes;
            Model3D::Lane3D lane;
            while (ss >> kv) {
                if (kv.rfind("NODES=",0)==0) nodes=kv.substr(6);
                else if (kv.rfind("Q=",0)==0) lane.q=std::stod(kv.substr(2));
                else if (kv.rfind("P=",0)==0) lane.p=std::stod(kv.substr(2));
                else if (kv.rfind("DOF=",0)==0) lane.dof=std::stoi(kv.substr(4));
                else if (kv.rfind("REDUCTION=",0)==0) lane.reduction=std::stod(kv.substr(10));
                else if (kv.rfind("IMPACT=",0)==0) lane.impact=std::stod(kv.substr(7));
                else if (kv.rfind("CODE=",0)==0) lane.autoCode = (kv.substr(5)=="JTG_D60");
                else if (kv.rfind("LANE_COUNT=",0)==0) lane.laneCount=std::stoi(kv.substr(11));
                else if (kv.rfind("ECCENTRICITY=",0)==0) lane.eccentricity=std::stod(kv.substr(13));
            }
            std::stringstream ns(nodes); std::string tok;
            while(std::getline(ns,tok,',')) if(!tok.empty()) lane.nodes.push_back(std::stoi(tok));
            m.lanes.push_back(lane);
        } else if (section == "*LANE_SURFACE_3D") {
            std::string kv, xn, yn;
            while (ss >> kv) {
                if (kv.rfind("X_NODES=",0)==0) xn=kv.substr(8);
                else if (kv.rfind("Y_NODES=",0)==0) yn=kv.substr(8);
                else if (kv.rfind("Q=",0)==0) m.laneSurface.q=std::stod(kv.substr(2));
                else if (kv.rfind("P=",0)==0) m.laneSurface.p=std::stod(kv.substr(2));
                else if (kv.rfind("DOF=",0)==0) m.laneSurface.dof=std::stoi(kv.substr(4));
                else if (kv.rfind("REDUCTION=",0)==0) m.laneSurface.reduction=std::stod(kv.substr(10));
                else if (kv.rfind("LANE_COUNT=",0)==0) m.laneSurface.laneCount=std::stoi(kv.substr(11));
                else if (kv.rfind("IMPACT=",0)==0) m.laneSurface.impact=std::stod(kv.substr(7));
                else if (kv.rfind("CODE=",0)==0) m.laneSurface.autoCode = (kv.substr(5)=="JTG_D60");
            }
            auto fill=[&](const std::string& s, std::vector<int>& out){
                std::stringstream ls(s); std::string tok;
                while(std::getline(ls,tok,',')) if(!tok.empty()) out.push_back(std::stoi(tok));
            };
            fill(xn,m.laneSurface.xNodes); fill(yn,m.laneSurface.yNodes);
        }
    }
    // 整份文件读完后才把 LOSSREF 换成有效比，使 *PRESTRESS_LOSS 写在
    // *STAGE 之前或之后都能引用。只给 LOSS= 的束不进这条路径。
    resolvePrestressLosses3D(m);
    applyJtgD60LaneLoads(m);
    return m;
}

} // namespace fem
