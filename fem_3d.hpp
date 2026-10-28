#pragma once

// 复用 fem.hpp 中的稀疏求解器与 Dirichlet 处理（applyDirichlet / solve）。
#include "fem.hpp"
#include "code_loads.hpp"
#include "section.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <map>
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

    struct Stage3D {
        int id = 0;
        std::vector<int> activate;
        std::vector<int> deactivate;
        std::vector<std::tuple<int, int, double>> loads;
    };
    std::vector<Stage3D> stages;
    double creepPhi = 0.0;
    double creepChi = 0.8;
    double shrinkageEcs = 0.0;

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
        if (e.istrn == 0.0 && m.shrinkageEcs == 0.0) continue;
        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
        const double L = std::max(d.norm(), 1e-12);
        const Vec3 c = d / L;
        const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
        const double P = Eeff * e.A * (e.istrn - m.shrinkageEcs);
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

inline ModalResult3D computeNaturalFrequencies3D(const Model3D& m, int nModes = 6) {
    DofMap3D dm;
    const SparseMat3D K = assembleStiffness3D(m, dm);
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
            Kd(a, b) = K.coeff(free[a], free[b]);
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

inline std::vector<double> computeBucklingFactors3D(const Model3D& m, int nModes = 3) {
    const std::vector<double> U = solveModel3D(m);   // 参考荷载下的线性位移
    DofMap3D dm;
    const SparseMat3D K = assembleStiffness3D(m, dm);

    // 几何刚度：由各梁的轴力（受拉正）组装
    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        if (e.type != "BEAM") continue;
        BeamElement b;
        b.E = e.E; b.A = e.A; b.Iy = e.Iy; b.Iz = e.Iz; b.J = e.J; b.G = e.G;
        b.x0 = m.nodes.at(e.n1); b.x1 = m.nodes.at(e.n2);
        const double L = b.length();
        const DenseMat T = b.transformation();
        Eigen::VectorXd ug(12);
        for (int k = 0; k < 6; ++k) {
            const int d1 = dm.safeDof(e.n1, k), d2 = dm.safeDof(e.n2, k);
            ug(k)     = (d1 >= 0 && d1 < static_cast<int>(U.size())) ? U[d1] : 0.0;
            ug(6 + k) = (d2 >= 0 && d2 < static_cast<int>(U.size())) ? U[d2] : 0.0;
        }
        const Eigen::VectorXd ul = T * ug;
        const double N = e.E * e.A / L * (ul(6) - ul(0));
        const DenseMat Kg = T.transpose() * b.geometricStiffness(N) * T;
        const std::vector<int> g = elementGdofs3D(dm, e);
        for (int r = 0; r < 12; ++r)
            for (int c = 0; c < 12; ++c)
                if (Kg(r, c) != 0.0) trips.emplace_back(g[r], g[c], Kg(r, c));
    }
    SparseMat3D Ksg(dm.nDofs, dm.nDofs);
    Ksg.setFromTriplets(trips.begin(), trips.end());
    Ksg.makeCompressed();

    const std::vector<int> free = freeDofs3D(m, dm, K);
    const int nf = static_cast<int>(free.size());
    std::vector<double> out;
    if (nf == 0) return out;

#ifdef FEM_HAS_SPECTRA
    const std::vector<double> sp = solveSparseBuckling(K, Ksg, free, dm.nDofs, nModes);
    if (!sp.empty()) return sp;
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
    if (llt.info() != Eigen::Success) return {};
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
    for (double mval : mu) out.push_back(-1.0 / mval);
    return out;
}

// ==================== 施工阶段 ====================

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

inline std::vector<StageResult3D> solveModel3DStages(const Model3D& m) {
    std::vector<char> active(m.elements.size(), 0);
    std::vector<StageResult3D> out;
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
        mc.creepPhi = m.creepPhi; mc.creepChi = m.creepChi; mc.shrinkageEcs = m.shrinkageEcs;
        for (std::size_t i = 0; i < m.elements.size(); ++i)
            if (active[i]) mc.elements.push_back(m.elements[i]);
        constrainDisconnected3D(mc);

        StageResult3D sr;
        sr.stage = st.id;
        sr.dm = buildDofMap3D(mc);   // 约束补齐后再取，与 solveModel3D 内部一致
        sr.U = solveModel3D(mc);
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
                const Eigen::VectorXd fl = beam7LocalEndForces(m, dm, e, sr.U);
                accumBeamEnv(r, -fl(0),  fl(1),  fl(2), -fl(3), -fl(4), -fl(5), -fl(6));   // i 端
                accumBeamEnv(r,  fl(7), -fl(8), -fl(9),  fl(10), fl(11), fl(12), fl(13)); // j 端
            } else {
                const Eigen::VectorXd fl = beam3DLocalEndForces(m, dm, e, sr.U);
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

inline std::vector<DisplacementEnvelope3D> compute3DDisplacementEnvelope(
    const Model3D& m, const std::vector<StageResult3D>& stages) {
    // 用"最后一个阶段"的映射决定输出哪些自由度；每个阶段再用它自己的映射取值。
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

inline std::pair<double,double> envelope3DFromInfluence(const Model3D::Lane3D& lane,
                                                        const std::vector<double>& eta) {
    double pos=0,neg=0,maxv=-1e300,minv=1e300;
    for (std::size_t i=0;i<eta.size();++i) {
        if (eta[i]>0) pos += eta[i]; else neg += eta[i];
        maxv=std::max(maxv,eta[i]); minv=std::min(minv,eta[i]);
    }
    const double qMaxArea = lane.q>=0 ? pos : neg;
    const double qMinArea = lane.q>=0 ? neg : pos;
    const double pMaxOrd = lane.p>=0 ? maxv : minv;
    const double pMinOrd = lane.p>=0 ? minv : maxv;
    return {lane.impact*lane.reduction*(lane.q*qMaxArea + lane.p*pMaxOrd),
            lane.impact*lane.reduction*(lane.q*qMinArea + lane.p*pMinOrd)};
}

struct StageLiveEnvelope3D {
    int stage = 0;
    double maxVal = 0.0, minVal = 0.0;
};

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
        for(const auto& lane:m.lanes){
            auto eta=compute3DInfluence(mc,spec,lane);
            auto env=envelope3DFromInfluence(lane,eta);
            maxAll=std::max(maxAll,env.first); minAll=std::min(minAll,env.second);
        }
        out.push_back({st.id,maxAll,minAll});
    }
    return out;
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

inline Model3D readModel3D(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开三维输入文件: " + path);
    Model3D m;
    std::map<int, std::string> types;
    std::map<int, double> mats;
    std::map<int, double> areas;
    std::map<int, double> realIy, realIz, realJ, realG;
    std::map<int, double> realIstrn;
    std::map<int, double> mpDens, realAddmass;
    std::map<int,double> realWtop, realWbot, realWy, realWt, realIw;
    int activeType = 0, activeMat = 0, activeReal = 0;
    std::string section;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        if (line[0] == '*') {
            const std::size_t sep = line.find_first_of(", \t");
            const std::string kw = (sep == std::string::npos) ? line : line.substr(0, sep);
            if (kw == "*TYPE" || kw == "*MAT" || kw == "*REAL" || kw == "*EN") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                if (kw == "*TYPE") ss >> activeType;
                else if (kw == "*MAT") ss >> activeMat;
                else if (kw == "*REAL") ss >> activeReal;
                else {
                    int id,n1,n2; if (ss >> id >> n1 >> n2) {
                        Element3D e; e.id=id; e.type=types[activeType]; e.E=mats[activeMat]; e.A=areas[activeReal]; e.n1=n1; e.n2=n2;
                        e.Iy=realIy[activeReal]; e.Iz=realIz[activeReal]; e.J=realJ[activeReal]; e.G=realG[activeReal];
                        e.istrn=realIstrn[activeReal];
                        e.dens=mpDens[activeMat]; e.addmass=realAddmass[activeReal];
                        e.Wtop=realWtop[activeReal]; e.Wbot=realWbot[activeReal]; e.Wy=realWy[activeReal]; e.Wt=realWt[activeReal];
                        e.Iw=realIw[activeReal];
                        m.elements.push_back(e);
                    }
                }
                continue;
            }
            if (kw == "*CREEP" || kw == "*SHRINKAGE") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                std::string kv;
                while (ss >> kv) {
                    if (kw == "*CREEP") {
                        if (kv.rfind("PHI=",0)==0) m.creepPhi=std::stod(kv.substr(4));
                        else if (kv.rfind("CHI=",0)==0) m.creepChi=std::stod(kv.substr(4));
                    } else {
                        if (kv.rfind("ECS=",0)==0) m.shrinkageEcs=std::stod(kv.substr(4));
                    }
                }
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
    applyJtgD60LaneLoads(m);
    return m;
}

} // namespace fem
