#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

// ---- 可选的高性能后端 AMGCL（代数多重网格）----
#if defined(__has_include)
#if __has_include(<amgcl/amg.hpp>) && __has_include(<amgcl/make_solver.hpp>)
#define FEM_HAS_AMGCL 1
#endif
#endif

#ifdef FEM_HAS_AMGCL
#ifndef AMGCL_NO_BOOST
#define AMGCL_NO_BOOST
#endif
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#endif

// ---- 可选的高性能后端 CHOLMOD（SuiteSparse 直接法）----
#if defined(__has_include) && __has_include(<cholmod.h>)
#define FEM_HAS_CHOLMOD 1
#endif

#ifdef FEM_HAS_CHOLMOD
#include <cholmod.h>
#endif

// 轻量有限元框架：支持一维杆、二维/三维桁架、二维平面应力三角形、
// 三维四面体，以及每节点多自由度组装与 Dirichlet 边界条件。
namespace fem {

using SparseMat = Eigen::SparseMatrix<double>;
using DenseMat = Eigen::MatrixXd;
using DenseVec = Eigen::VectorXd;

// 让 Windows 控制台按 UTF-8 解释程序输出，避免中文乱码。
inline void setupConsoleUTF8() {
#if defined(_WIN32)
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
}

// ---- 自由度与组装辅助函数 ----

// 节点 node 的第 localDof 个自由度对应的全局自由度编号。
inline int dofIndex(int node, int localDof, int dofsPerNode) {
    return node * dofsPerNode + localDof;
}

// 一组节点各自的所有自由度编号（按节点顺序展开）。
inline std::vector<int> nodeDofs(const std::vector<int>& nodes, int dofsPerNode) {
    std::vector<int> g;
    g.reserve(nodes.size() * dofsPerNode);
    for (int n : nodes)
        for (int d = 0; d < dofsPerNode; ++d)
            g.push_back(n * dofsPerNode + d);
    return g;
}

// 把单元刚度矩阵按自由度映射散入 Triplet 列表。
inline void scatter(const DenseMat& Ke, const std::vector<int>& gdofs,
                    std::vector<Eigen::Triplet<double>>& triplets) {
    assert(static_cast<int>(gdofs.size()) == Ke.rows());
    for (int i = 0; i < Ke.rows(); ++i)
        for (int j = 0; j < Ke.cols(); ++j)
            if (Ke(i, j) != 0.0)
                triplets.emplace_back(gdofs[i], gdofs[j], Ke(i, j));
}

// 由若干单元刚度矩阵组装全局稀疏刚度矩阵。
inline SparseMat assemble(const std::vector<DenseMat>& kes,
                          const std::vector<std::vector<int>>& gdofsList,
                          int nDofs) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (std::size_t e = 0; e < kes.size(); ++e)
        scatter(kes[e], gdofsList[e], triplets);
    SparseMat K(nDofs, nDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();
    return K;
}

// Dirichlet 边界条件：对指定自由度施加指定位移（行/列消去法）。
inline void applyDirichlet(SparseMat& K, DenseVec& F,
                           const std::vector<int>& constrainedDofs,
                           const std::vector<double>& prescribedValues) {
    assert(constrainedDofs.size() == prescribedValues.size());
    const int n = static_cast<int>(K.rows());
    std::vector<char> isConstrained(n, 0);
    for (int d : constrainedDofs) isConstrained[d] = 1;

    // 1) 把被约束自由度所在列的已知贡献移到右端。
    for (std::size_t idx = 0; idx < constrainedDofs.size(); ++idx) {
        const int d = constrainedDofs[idx];
        const double g = prescribedValues[idx];
        for (int col = 0; col < K.outerSize(); ++col) {
            for (SparseMat::InnerIterator it(K, col); it; ++it) {
                if (it.col() == d && it.row() != d)
                    F(it.row()) -= it.value() * g;
            }
        }
    }

    // 2) 清零被约束自由度所在的行与列（保持对称性）。
    for (int col = 0; col < K.outerSize(); ++col) {
        for (SparseMat::InnerIterator it(K, col); it; ++it) {
            if (isConstrained[it.row()] || isConstrained[it.col()])
                it.valueRef() = 0.0;
        }
    }

    // 3) 对角线置 1，右端置为指定位移。
    for (std::size_t idx = 0; idx < constrainedDofs.size(); ++idx) {
        K.coeffRef(constrainedDofs[idx], constrainedDofs[idx]) = 1.0;
        F(constrainedDofs[idx]) = prescribedValues[idx];
    }
}

// 可切换的求解后端。
enum class SolverType {
    SimplicialLDLT,     // 默认：对称正定直接法（LDL^T）
    SimplicialLLT,      // 对称正定直接法（Cholesky LL^T）
    SparseLU,           // 一般矩阵直接法（SuperLU）
    SparseQR,           // QR 分解（最小二乘 / 欠定 / 满秩均可）
    ConjugateGradient,  // 对称正定迭代法（Jacobi 预条件）
    BiCGSTAB,           // 一般矩阵迭代法（Jacobi 预条件）
    AMGCL,              // 代数多重网格 + BiCGStab（高性能，需 AMGCL 头文件）
    CHOLMOD             // SuiteSparse CHOLMOD 直接法（Cholesky）
};

namespace detail {

template <typename Solver>
inline DenseVec runSolver(Solver& solver, const SparseMat& K, const DenseVec& F) {
    solver.compute(K);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("求解器：分解/预处理失败");
    DenseVec u = solver.solve(F);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("求解器：求解失败（迭代法可能未收敛）");
    return u;
}

template <typename IterativeSolver>
inline void setIterativeDefaults(IterativeSolver& solver, const SparseMat& K) {
    solver.setTolerance(1e-12);
    int maxIter = 100 * static_cast<int>(K.rows());
    if (maxIter < 1000) maxIter = 1000;
    solver.setMaxIterations(maxIter);
}

} // namespace detail

#ifdef FEM_HAS_AMGCL
// AMGCL：代数多重网格预条件 + BiCGStab。
inline DenseVec solveAMGCL(const SparseMat& K, const DenseVec& F) {
    // 转成 AMGCL 需要的 CRS，并去掉显式零。
    SparseMat Kp = K;
    Kp.prune(0.0);
    Kp.makeCompressed();
    Eigen::SparseMatrix<double, Eigen::RowMajor> Kcsr(Kp);

    std::ptrdiff_t n = Kcsr.rows();
    std::vector<std::ptrdiff_t> ptr(Kcsr.outerIndexPtr(),
                                    Kcsr.outerIndexPtr() + n + 1);
    std::vector<std::ptrdiff_t> col(Kcsr.innerIndexPtr(),
                                    Kcsr.innerIndexPtr() + Kcsr.nonZeros());
    std::vector<double> val(Kcsr.valuePtr(),
                            Kcsr.valuePtr() + Kcsr.nonZeros());
    std::vector<double> rhs(F.data(), F.data() + n);
    std::vector<double> x(n, 0.0);

    auto A = std::tie(n, ptr, col, val);
    using Backend = amgcl::backend::builtin<double>;
    using Solver = amgcl::make_solver<
        amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation,
                   amgcl::relaxation::spai0>,
        amgcl::solver::bicgstab<Backend>>;
    Solver::params P;
    P.solver.maxiter = 1000;
    P.solver.tol = 1e-10;
    Solver solver(A, P);

    int iters = 0;
    double error = 0.0;
    std::tie(iters, error) = solver(A, rhs, x);
    std::cerr << "[AMGCL] 迭代次数 = " << iters
              << "，相对残差 = " << error << "\n";

    return Eigen::Map<DenseVec>(x.data(), n);
}
#endif

#ifdef FEM_HAS_CHOLMOD
// SuiteSparse CHOLMOD 封装：支持符号分解(analyze)复用、数值分解(factorize)复用，
// 以及单/多右端项求解。
//   CholmodSolver solver;
//   solver.analyzePattern(K0);   // 符号分解：同一结构只做一次
//   solver.factorize(K1);        // 数值分解：结构不变时复用符号分解
//   DenseVec u = solver.solve(f);
//   DenseMat U = solver.solve(B); // B 每列一个右端项
class CholmodSolver {
public:
    CholmodSolver() {
        cholmod_start(&common_);
        common_.itype = CHOLMOD_INT;
    }

    ~CholmodSolver() {
        clear();
        cholmod_finish(&common_);
    }

    CholmodSolver(const CholmodSolver&) = delete;
    CholmodSolver& operator=(const CholmodSolver&) = delete;

    // 符号分解：仅依赖矩阵非零结构，同一结构只调用一次。
    void analyzePattern(const SparseMat& K) {
        clear();
        n_ = static_cast<int>(K.rows());
        A_ = buildMatrix(K);
        L_ = cholmod_analyze(A_, &common_);
        if (!L_) {
            freeMatrix();
            throw std::runtime_error("CHOLMOD: 符号分解失败");
        }
    }

    // 数值分解：复用已完成的符号分解（结构不变时无需重新 analyze）。
    void factorize(const SparseMat& K) {
        if (!L_)
            throw std::runtime_error("CHOLMOD: 请先调用 analyzePattern()");
        freeMatrix();
        n_ = static_cast<int>(K.rows());
        A_ = buildMatrix(K);
        if (!cholmod_factorize(A_, L_, &common_)) {
            freeMatrix();
            freeFactor();
            throw std::runtime_error("CHOLMOD: 数值分解失败（矩阵可能非正定）");
        }
    }

    // 单右端项。
    DenseVec solve(const DenseVec& F) {
        requireReady();
        cholmod_dense* b = cholmod_allocate_dense(
            static_cast<size_t>(n_), 1, static_cast<size_t>(n_), CHOLMOD_REAL,
            &common_);
        std::copy(F.data(), F.data() + n_, static_cast<double*>(b->x));
        cholmod_dense* x = cholmod_solve(CHOLMOD_A, L_, b, &common_);
        cholmod_free_dense(&b, &common_);
        if (!x) throw std::runtime_error("CHOLMOD: 求解失败");
        DenseVec U(n_);
        std::copy(static_cast<double*>(x->x),
                  static_cast<double*>(x->x) + n_, U.data());
        cholmod_free_dense(&x, &common_);
        return U;
    }

    // 多右端项：B 为 n x k，返回 n x k 的解（每列对应一个右端项）。
    DenseMat solve(const DenseMat& B) {
        requireReady();
        const int k = static_cast<int>(B.cols());
        if (B.rows() != n_)
            throw std::runtime_error("CHOLMOD: 右端项行数与矩阵不匹配");
        cholmod_dense* b = cholmod_allocate_dense(
            static_cast<size_t>(n_), static_cast<size_t>(k),
            static_cast<size_t>(n_), CHOLMOD_REAL, &common_);
        std::copy(B.data(), B.data() + static_cast<long long>(n_) * k,
                  static_cast<double*>(b->x));
        cholmod_dense* x = cholmod_solve(CHOLMOD_A, L_, b, &common_);
        cholmod_free_dense(&b, &common_);
        if (!x) throw std::runtime_error("CHOLMOD: 求解失败");
        DenseMat X(n_, k);
        std::copy(static_cast<double*>(x->x),
                  static_cast<double*>(x->x) + static_cast<long long>(n_) * k,
                  X.data());
        cholmod_free_dense(&x, &common_);
        return X;
    }

    // 释放矩阵与因子（保留 cholmod_common，可再次 analyzePattern 重用）。
    void clear() {
        freeFactor();
        freeMatrix();
        n_ = 0;
    }

private:
    void requireReady() const {
        if (!A_ || !L_)
            throw std::runtime_error("CHOLMOD: 尚未完成分解");
    }

    cholmod_sparse* buildMatrix(const SparseMat& K) {
        SparseMat Kp = K;
        Kp.prune(0.0);
        Kp.makeCompressed();
        const int n = static_cast<int>(Kp.rows());

        std::vector<Eigen::Triplet<double>> lower;
        lower.reserve(Kp.nonZeros());
        for (int k = 0; k < Kp.outerSize(); ++k)
            for (SparseMat::InnerIterator it(Kp, k); it; ++it)
                if (it.row() >= it.col())
                    lower.emplace_back(static_cast<int>(it.row()),
                                       static_cast<int>(it.col()), it.value());
        SparseMat Kl(n, n);
        Kl.setFromTriplets(lower.begin(), lower.end());
        Kl.makeCompressed();

        cholmod_sparse* A = cholmod_allocate_sparse(
            static_cast<size_t>(Kl.rows()), static_cast<size_t>(Kl.cols()),
            static_cast<size_t>(Kl.nonZeros()), true, true, -1, CHOLMOD_REAL,
            &common_);
        if (!A) throw std::runtime_error("CHOLMOD: 分配矩阵失败");
        std::copy(Kl.outerIndexPtr(), Kl.outerIndexPtr() + n + 1,
                  static_cast<int*>(A->p));
        std::copy(Kl.innerIndexPtr(), Kl.innerIndexPtr() + Kl.nonZeros(),
                  static_cast<int*>(A->i));
        std::copy(Kl.valuePtr(), Kl.valuePtr() + Kl.nonZeros(),
                  static_cast<double*>(A->x));
        return A;
    }

    void freeMatrix() {
        if (A_) {
            cholmod_free_sparse(&A_, &common_);
            A_ = nullptr;
        }
    }

    void freeFactor() {
        if (L_) {
            cholmod_free_factor(&L_, &common_);
            L_ = nullptr;
        }
    }

    cholmod_common common_;
    cholmod_sparse* A_ = nullptr;
    cholmod_factor* L_ = nullptr;
    int n_ = 0;
};

// 便捷的一次性求解：内部即 analyzePattern -> factorize -> solve。
inline DenseVec solveCHOLMOD(const SparseMat& K, const DenseVec& F) {
    CholmodSolver solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    return solver.solve(F);
}
#endif

// 通用求解接口：type 默认使用 SimplicialLDLT。
inline DenseVec solve(const SparseMat& K, const DenseVec& F,
                      SolverType type = SolverType::SimplicialLDLT) {
    switch (type) {
        case SolverType::SimplicialLDLT: {
            Eigen::SimplicialLDLT<SparseMat> s;
            return detail::runSolver(s, K, F);
        }
        case SolverType::SimplicialLLT: {
            Eigen::SimplicialLLT<SparseMat> s;
            return detail::runSolver(s, K, F);
        }
        case SolverType::SparseLU: {
            Eigen::SparseLU<SparseMat> s;
            return detail::runSolver(s, K, F);
        }
        case SolverType::SparseQR: {
            Eigen::SparseQR<SparseMat, Eigen::COLAMDOrdering<int>> s;
            return detail::runSolver(s, K, F);
        }
        case SolverType::ConjugateGradient: {
            Eigen::ConjugateGradient<SparseMat, Eigen::Lower,
                                     Eigen::DiagonalPreconditioner<double>> s;
            detail::setIterativeDefaults(s, K);
            return detail::runSolver(s, K, F);
        }
        case SolverType::BiCGSTAB: {
            Eigen::BiCGSTAB<SparseMat, Eigen::DiagonalPreconditioner<double>> s;
            detail::setIterativeDefaults(s, K);
            return detail::runSolver(s, K, F);
        }
        case SolverType::AMGCL:
#ifdef FEM_HAS_AMGCL
            return solveAMGCL(K, F);
#else
            throw std::runtime_error("AMGCL 后端未启用（缺少 AMGCL 头文件）");
#endif
        case SolverType::CHOLMOD:
#ifdef FEM_HAS_CHOLMOD
            return solveCHOLMOD(K, F);
#else
            throw std::runtime_error("CHOLMOD 后端未启用（缺少 SuiteSparse 头文件）");
#endif
    }
    throw std::runtime_error("未知的求解器类型");
}

// 向后兼容：默认仍使用 SimplicialLDLT。
inline DenseVec solveSPD(const SparseMat& K, const DenseVec& F) {
    return solve(K, F, SolverType::SimplicialLDLT);
}

// ==================== 可复用的求解器句柄 ====================
//
// 直接法的成本分成两半：符号分解（analyzePattern：AMD/COLAMD 排序 + 消去树 +
// 符号模式，只依赖非零结构）和数值分解（factorize：真正算因子，依赖数值）。
// 影响线扫描、Newton 迭代、施工阶段推进这三处，矩阵**结构完全不变、只有数值变**
// （或干脆数值也不变），却每次都走一遍 `compute()` = analyze + factorize。
//
// 这个类把两者拆开：结构指纹不变时只做 factorize，并支持多右端项一次回代。
// CHOLMOD 封装本来就有 analyzePattern/factorize 分离，但从来没有接到求解路径上
// （只在 example_cholmod.cpp 演示过一次），这里让它真正可用。
class ReusableSolver {
public:
    explicit ReusableSolver(SolverType type = SolverType::SimplicialLDLT)
        : type_(type),
          reusable_(type == SolverType::SimplicialLDLT ||
                    type == SolverType::SimplicialLLT ||
                    type == SolverType::SparseLU ||
                    type == SolverType::CHOLMOD)
#ifdef FEM_HAS_CHOLMOD
          , cholmod_(nullptr)
#endif
    {}

    SolverType type() const { return type_; }
    // 该后端是否支持 analyze/factorize 分离（迭代法与 AMGCL 不支持）。
    bool reusable() const { return reusable_; }
    int analyzeCount() const { return analyzeCount_; }
    int factorizeCount() const { return factorizeCount_; }
    // 结构未变却只做了 1 次符号分解 → 说明复用生效了（用于自检/统计）。
    bool patternReused() const { return factorizeCount_ > analyzeCount_; }

    // 结构指纹：行数 + 非零元数 + 列偏移 + 行索引。四项都相同才认为结构不变。
    bool samePattern(const SparseMat& K) const {
        if (!patternValid_ || !reusable_) return false;
        if (K.rows() != n_ || K.cols() != n_) return false;
        if (!K.isCompressed()) return false;
        if (K.nonZeros() != nnz_) return false;
        const int* outer = K.outerIndexPtr();
        const int* inner = K.innerIndexPtr();
        for (int i = 0; i <= n_; ++i)
            if (outer[i] != outer_[static_cast<std::size_t>(i)]) return false;
        for (int i = 0; i < nnz_; ++i)
            if (inner[i] != inner_[static_cast<std::size_t>(i)]) return false;
        return true;
    }

    // 符号分解：只在结构变化时需要。
    void analyze(const SparseMat& K) {
        if (!reusable_)
            throw std::runtime_error("该求解后端不支持符号分解复用");
        switch (type_) {
            case SolverType::SimplicialLDLT: ldlt_.analyzePattern(K); break;
            case SolverType::SimplicialLLT:  llt_.analyzePattern(K);  break;
            case SolverType::SparseLU:       lu_.analyzePattern(K);   break;
            case SolverType::CHOLMOD:
#ifdef FEM_HAS_CHOLMOD
                cholmod().analyzePattern(K);
#else
                throw std::runtime_error("CHOLMOD 后端未启用（缺少 SuiteSparse 头文件）");
#endif
                break;
            default:
                throw std::runtime_error("该求解后端不支持符号分解复用");
        }
        capturePattern(K);
        ++analyzeCount_;
    }

    // 数值分解：可反复调用，复用上一次的符号分解结果。
    void factorize(const SparseMat& K) {
        if (!reusable_) {
            Klast_ = K;                       // 迭代法/AMGCL 只能整体重来
            ++factorizeCount_;
            return;
        }
        if (factorizeCount_ == 0 && !patternValid_)
            analyze(K);
        switch (type_) {
            case SolverType::SimplicialLDLT:
                ldlt_.factorize(K);
                if (ldlt_.info() != Eigen::Success)
                    throw std::runtime_error("求解器：LDLT 数值分解失败");
                break;
            case SolverType::SimplicialLLT:
                llt_.factorize(K);
                if (llt_.info() != Eigen::Success)
                    throw std::runtime_error("求解器：LLT 数值分解失败");
                break;
            case SolverType::SparseLU:
                lu_.factorize(K);
                if (lu_.info() != Eigen::Success)
                    throw std::runtime_error("求解器：LU 数值分解失败");
                break;
            case SolverType::CHOLMOD:
#ifdef FEM_HAS_CHOLMOD
                cholmod().factorize(K);
#else
                throw std::runtime_error("CHOLMOD 后端未启用（缺少 SuiteSparse 头文件）");
#endif
                break;
            default:
                throw std::runtime_error("该求解后端不支持分解复用");
        }
        ++factorizeCount_;
    }

    // 便捷入口：结构不变则只 factorize，否则 analyze + factorize。
    void factorizeReuse(const SparseMat& K) {
        if (!samePattern(K)) analyze(K);
        factorize(K);
    }

    DenseVec solve(const DenseVec& F) const {
        if (!reusable_) {
            if (Klast_.rows() == 0)
                throw std::runtime_error("ReusableSolver：尚未 factorize");
            return fem::solve(Klast_, F, type_);
        }
        if (factorizeCount_ == 0)
            throw std::runtime_error("ReusableSolver：尚未 factorize");
        DenseVec u;
        switch (type_) {
            case SolverType::SimplicialLDLT: u = ldlt_.solve(F); break;
            case SolverType::SimplicialLLT:  u = llt_.solve(F);  break;
            case SolverType::SparseLU:       u = lu_.solve(F);   break;
            case SolverType::CHOLMOD:
#ifdef FEM_HAS_CHOLMOD
                if (!cholmod_) throw std::runtime_error("ReusableSolver：CHOLMOD 未初始化");
                u = cholmod_->solve(F);
#else
                throw std::runtime_error("CHOLMOD 后端未启用（缺少 SuiteSparse 头文件）");
#endif
                break;
            default:
                throw std::runtime_error("该求解后端不支持分解复用");
        }
        if (!u.allFinite())
            throw std::runtime_error("求解器：结果含 NaN/Inf");
        return u;
    }

    // 多右端项：B 的每一列一个右端项。CHOLMOD 原生支持；Eigen 逐列回代。
    DenseMat solveMulti(const DenseMat& B) const {
        const int k = static_cast<int>(B.cols());
        if (B.rows() != n_)
            throw std::runtime_error("ReusableSolver：右端项行数与矩阵不匹配");
#ifdef FEM_HAS_CHOLMOD
        if (reusable_ && type_ == SolverType::CHOLMOD && cholmod_)
            return cholmod_->solve(B);
#endif
        DenseMat X(n_, k);
        for (int c = 0; c < k; ++c)
            X.col(c) = solve(DenseVec(B.col(c)));
        return X;
    }

private:
    void capturePattern(const SparseMat& K) {
        n_ = static_cast<int>(K.rows());
        nnz_ = static_cast<int>(K.nonZeros());
        outer_.assign(K.outerIndexPtr(), K.outerIndexPtr() + n_ + 1);
        inner_.assign(K.innerIndexPtr(), K.innerIndexPtr() + nnz_);
        patternValid_ = K.isCompressed();
    }

#ifdef FEM_HAS_CHOLMOD
    CholmodSolver& cholmod() {
        if (!cholmod_) cholmod_ = std::make_unique<CholmodSolver>();
        return *cholmod_;
    }
#endif

    SolverType type_;
    bool reusable_;
    Eigen::SimplicialLDLT<SparseMat> ldlt_;
    Eigen::SimplicialLLT<SparseMat> llt_;
    Eigen::SparseLU<SparseMat> lu_;
#ifdef FEM_HAS_CHOLMOD
    std::unique_ptr<CholmodSolver> cholmod_;
#endif
    SparseMat Klast_;               // 仅非可复用后端需要
    std::vector<int> outer_, inner_;
    int n_ = 0, nnz_ = 0;
    bool patternValid_ = false;
    int analyzeCount_ = 0;
    int factorizeCount_ = 0;
};

// ==================== 单元 ====================

// 两节点杆/桁架单元，适用于 1D / 2D / 3D。
struct TrussElement {
    double E, A;
    int dim; // 1、2 或 3
    Eigen::Vector3d x0, x1;

    int dofsPerNode() const { return dim; }
    int numNodes() const { return 2; }

    DenseMat stiffness() const {
        Eigen::Vector3d d = x1 - x0;
        const double L = d.norm();
        if (L <= 0.0) throw std::runtime_error("桁架单元长度为零");
        const Eigen::Vector3d c = d / L;
        const Eigen::VectorXd dir = c.head(dim);
        const DenseMat T = dir * dir.transpose(); // 方向余弦外积
        DenseMat Ke = DenseMat::Zero(2 * dim, 2 * dim);
        Ke.topLeftCorner(dim, dim) = T;
        Ke.topRightCorner(dim, dim) = -T;
        Ke.bottomLeftCorner(dim, dim) = -T;
        Ke.bottomRightCorner(dim, dim) = T;
        return (E * A / L) * Ke;
    }
};

// 三节点常应变三角形（CST），二维平面应力。
struct TriangleCST {
    double E, nu, t; // t 为厚度
    Eigen::Vector2d x0, x1, x2;

    int dofsPerNode() const { return 2; }
    int numNodes() const { return 3; }

    DenseMat D() const { // 平面应力本构矩阵
        const double c = E / (1.0 - nu * nu);
        DenseMat Dm(3, 3);
        Dm << 1.0, nu, 0.0,
              nu, 1.0, 0.0,
              0.0, 0.0, (1.0 - nu) / 2.0;
        return c * Dm;
    }

    DenseMat B() const { // 应变-位移矩阵 (3 x 6)
        Eigen::Matrix2d J;
        J.col(0) = x0 - x2;
        J.col(1) = x1 - x2;
        if (J.determinant() == 0.0) throw std::runtime_error("退化三角形单元");
        const Eigen::Matrix2d JinvT = J.inverse().transpose();

        // 参考三角形：节点0=(1,0)，节点1=(0,1)，节点2=(0,0)
        const Eigen::Vector2d g[3] = {
            Eigen::Vector2d(1.0, 0.0),
            Eigen::Vector2d(0.0, 1.0),
            Eigen::Vector2d(-1.0, -1.0)};

        DenseMat B = DenseMat::Zero(3, 6);
        for (int a = 0; a < 3; ++a) {
            const Eigen::Vector2d p = JinvT * g[a];
            const double bx = p(0), by = p(1);
            B(0, 2 * a) = bx;
            B(1, 2 * a + 1) = by;
            B(2, 2 * a) = by;
            B(2, 2 * a + 1) = bx;
        }
        return B;
    }

    double area() const {
        Eigen::Matrix2d J;
        J.col(0) = x0 - x2;
        J.col(1) = x1 - x2;
        return 0.5 * std::abs(J.determinant());
    }

    DenseMat stiffness() const {
        const DenseMat Bm = B();
        return t * area() * Bm.transpose() * D() * Bm;
    }
};

// 四节点常应变四面体，三维各向同性实体。
struct Tetra4 {
    double E, nu;
    Eigen::Vector3d x0, x1, x2, x3;

    int dofsPerNode() const { return 3; }
    int numNodes() const { return 4; }

    DenseMat D() const { // 各向同性三维本构矩阵
        const double a = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double b = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double g = E / (2.0 * (1.0 + nu));
        DenseMat Dm(6, 6);
        Dm << a, b, b, 0.0, 0.0, 0.0,
              b, a, b, 0.0, 0.0, 0.0,
              b, b, a, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, g, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, g, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, g;
        return Dm;
    }

    DenseMat B() const { // 应变-位移矩阵 (6 x 12)
        Eigen::Matrix3d J;
        J.col(0) = x0 - x3;
        J.col(1) = x1 - x3;
        J.col(2) = x2 - x3;
        if (J.determinant() == 0.0) throw std::runtime_error("退化四面体单元");
        const Eigen::Matrix3d JinvT = J.inverse().transpose();

        // 参考四面体：节点0=(1,0,0)，节点1=(0,1,0)，节点2=(0,0,1)，节点3=(0,0,0)
        const Eigen::Vector3d g[4] = {
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(0.0, 1.0, 0.0),
            Eigen::Vector3d(0.0, 0.0, 1.0),
            Eigen::Vector3d(-1.0, -1.0, -1.0)};

        DenseMat B = DenseMat::Zero(6, 12);
        for (int a = 0; a < 4; ++a) {
            const Eigen::Vector3d p = JinvT * g[a];
            const double bx = p(0), by = p(1), bz = p(2);
            const int c = 3 * a;
            B(0, c + 0) = bx;
            B(1, c + 1) = by;
            B(2, c + 2) = bz;
            B(3, c + 0) = by; B(3, c + 1) = bx; // gamma_xy
            B(4, c + 1) = bz; B(4, c + 2) = by; // gamma_yz
            B(5, c + 0) = bz; B(5, c + 2) = bx; // gamma_xz
        }
        return B;
    }

    double volume() const {
        Eigen::Matrix3d J;
        J.col(0) = x0 - x3;
        J.col(1) = x1 - x3;
        J.col(2) = x2 - x3;
        return std::abs(J.determinant()) / 6.0;
    }

    DenseMat stiffness() const {
        const DenseMat Bm = B();
        return volume() * Bm.transpose() * D() * Bm;
    }
};

// 两节点二维欧拉-伯努利梁（框架）单元。
// 局部自由度顺序：每节点 [u(轴向), v(横向), θ(转角)]，共 6 个。
struct EulerBeam2D {
    double E, A, I; // 弹性模量, 截面积, 截面惯性矩
    Eigen::Vector2d x0, x1;

    int dofsPerNode() const { return 3; }
    int numNodes() const { return 2; }

    // 局部刚度矩阵（轴向 + 弯曲），局部自由度 [u1,v1,θ1,u2,v2,θ2]。
    DenseMat localStiffness() const {
        const Eigen::Vector2d d = x1 - x0;
        const double L = d.norm();
        if (L <= 0.0) throw std::runtime_error("梁单元长度为零");
        const double k = E * I / (L * L * L); // 弯曲主系数 EI/L^3

        DenseMat Kl = DenseMat::Zero(6, 6);
        Kl << E * A / L, 0.0, 0.0, -E * A / L, 0.0, 0.0,
              0.0, 12.0 * k, 6.0 * k * L, 0.0, -12.0 * k, 6.0 * k * L,
              0.0, 6.0 * k * L, 4.0 * k * L * L, 0.0, -6.0 * k * L, 2.0 * k * L * L,
              -E * A / L, 0.0, 0.0, E * A / L, 0.0, 0.0,
              0.0, -12.0 * k, -6.0 * k * L, 0.0, 12.0 * k, -6.0 * k * L,
              0.0, 6.0 * k * L, 2.0 * k * L * L, 0.0, -6.0 * k * L, 4.0 * k * L * L;
        return Kl;
    }

    // 局部坐标 -> 全局坐标的变换矩阵：T（d_local = T * d_global）。
    DenseMat transformation() const {
        const Eigen::Vector2d d = x1 - x0;
        const double L = d.norm();
        if (L <= 0.0) throw std::runtime_error("梁单元长度为零");
        const double c = d(0) / L;
        const double s = d(1) / L;

        DenseMat T = DenseMat::Zero(6, 6);
        T << c, s, 0.0, 0.0, 0.0, 0.0,
             -s, c, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, c, s, 0.0,
             0.0, 0.0, 0.0, -s, c, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
        return T;
    }

    DenseMat stiffness() const {
        const DenseMat T = transformation();
        const DenseMat Kl = localStiffness();
        return T.transpose() * Kl * T;
    }
};

} // namespace fem
