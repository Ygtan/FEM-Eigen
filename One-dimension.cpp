#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

// 定义单元数据结构
struct Element {
    int node1, node2;
    double E, A, L; // 弹性模量, 面积, 长度
};

int main() {
#if defined(_WIN32)
    SetConsoleOutputCP(CP_UTF8);
#endif
    // ----- 1. 定义问题参数和网格 -----
    const double E = 200e9; // 弹性模量 (Pa)
    const double A = 0.01;  // 截面积 (m^2)
    const double L = 1.0;   // 单元长度 (m)
    const int num_elements = 3; // 单元数量
    const int num_nodes = num_elements + 1; // 节点数量
    const int dof_per_node = 1; // 每个节点1个自由度
    const int total_dof = num_nodes * dof_per_node;

    // 创建单元列表
    std::vector<Element> elements;
    for (int i = 0; i < num_elements; ++i) {
        elements.push_back({i, i+1, E, A, L});
    }

    // ----- 2. 组装全局稀疏矩阵 K 和载荷向量 F -----
    // 使用 Triplet 列表高效组装
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(total_dof);

    // 遍历每个单元
    for (const auto& elem : elements) {
        // 计算单元刚度矩阵 Ke (对于1D杆, 是一个2x2矩阵)
        double k = elem.E * elem.A / elem.L;
        Eigen::Matrix2d Ke;
        Ke <<  k, -k,
              -k,  k;

        // 自由度映射
        int dof1 = elem.node1 * dof_per_node;
        int dof2 = elem.node2 * dof_per_node;
        std::vector<int> elemDofs = {dof1, dof2};

        // 将Ke的每个元素作为Triplet加入列表
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                tripletList.push_back(Eigen::Triplet<double>(elemDofs[i], elemDofs[j], Ke(i, j)));
            }
        }
    }

    // 在节点2施加一个向右的集中力 (索引为2的节点，如果从0开始计数)
    F(2) = 1000.0; // 1000 N

    // 根据Triplet列表构建稀疏矩阵K
    Eigen::SparseMatrix<double> K(total_dof, total_dof);
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    K.makeCompressed(); // 压缩存储格式，优化后续求解

    // ----- 3. 施加边界条件 (置1消行/消列法) -----
    // 固定节点0 (索引0) 的位移 u0 = 0
    // 注意：必须同时清零该自由度所在的行和列，否则会破坏矩阵的对称性，
    // 导致对称求解器 (SimplicialLDLT) 得到错误结果。
    int fixed_dof = 0;
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            if (it.row() == fixed_dof || it.col() == fixed_dof) {
                it.valueRef() = 0.0;
            }
        }
    }
    K.coeffRef(fixed_dof, fixed_dof) = 1.0;
    F(fixed_dof) = 0.0;

    // ----- 4. 求解线性方程组 K * U = F -----
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    if(solver.info()!=Eigen::Success) {
        std::cerr << "分解失败!" << std::endl;
        return 1;
    }
    Eigen::VectorXd U = solver.solve(F);
    if (solver.info() != Eigen::Success) {
        std::cerr << "求解失败!" << std::endl;
        return 1;
    }

    // ----- 5. 输出结果 -----
    std::cout << "节点位移 (m):" << std::endl;
    for (int i = 0; i < total_dof; ++i) {
        std::cout << "节点 " << i << ": " << U(i) << std::endl;
    }

    return 0;
}
