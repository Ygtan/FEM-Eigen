#pragma once
// modal_output_3d.hpp
// ==================
// 把内核**已经算出来**的模态结果写成文本段，供 BridgeStudio 解析。
//
// 这里不做任何特征值计算：K φ = λ M φ 由 fem_3d.hpp 的
// computeNaturalFrequencies3D / cable_solver.cpp 的 modalFromTangent 求解，
// 本文件只负责把同一份 ModalResult3D 格式化输出。这样报告里的各阶频率、
// 振型与既有 f1_vertical（汽车冲击系数依据）必然同源，不会出现"两套频率"。
//
// 段格式（每行定长字段，便于正则解析；数值 %.10g）：
//   <header>                      "modal"（summary.txt）或 "[MODAL]"（cable.txt）
//   NMODES=<n> NDOF=<nd> BASIS=xyz
//   MASS X=<kg> Y=<kg> Z=<kg>       各方向平动总质量 r_d^T M r_d
//   F <i> <f_Hz> <T_s> <mx> <my> <mz> <mnorm>
//                                  i 从 1 起；mx/my/mz 为该阶在 X/Y/Z 方向的
//                                  有效参与质量 Γ_d²/(φ^T M φ)（kg）；
//                                  mnorm = φ^T M φ，M-正交归一时 ≈ 1。
//   S <i> <node> <ux> <uy> <uz>      该阶在该节点的三个平动分量（全局坐标系）
//
// 参与质量的定义：Γ_d = φ^T M r_d，r_d 为方向 d 的刚体平动向量（该方向所有
// 平动自由度置 1，其余为 0）。有效参与质量 = Γ_d² / (φ^T M φ)。Σ 各阶
// 有效参与质量趋近于 r_d^T M r_d，是模态截断是否足够的直接判据。
//
// 自由度顺序与 fem.hpp 的 DofMap 一致：[u, v, w, θx, θy, θz]，故
// local 2 = w = 竖向（Z）。纯轴向节点只有 3 个平动自由度，safeDof 自动处理。

#include "fem_3d.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

inline std::string formatModalSection3D(const DofMap3D& dm,
                                        const SparseMat3D& M,
                                        const ModalResult3D& modal,
                                        const std::string& header) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << header << "\n";

    const int nd = dm.nDofs;
    const int nm = static_cast<int>(modal.frequencies.size());
    ss << "NMODES=" << nm << " NDOF=" << nd << " BASIS=xyz\n";

    if (nm == 0) {
        ss << "MASS X=0 Y=0 Z=0\n";
        return ss.str();
    }

    Eigen::VectorXd r[3], Mr[3];
    double mtot[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < 3; ++d) {
        r[d] = Eigen::VectorXd::Zero(nd);
        for (int nid : dm.nodeIds) {
            const int gd = dm.safeDof(nid, d);
            if (gd >= 0 && gd < nd) r[d][gd] = 1.0;
        }
        Mr[d] = M * r[d];
        mtot[d] = r[d].dot(Mr[d]);
    }
    ss << "MASS X=" << mtot[0] << " Y=" << mtot[1] << " Z=" << mtot[2] << "\n";

    for (int i = 0; i < nm; ++i) {
        const double f = modal.frequencies[i];
        const double period = (f > 0.0) ? 1.0 / f : 0.0;
        const std::vector<double>& mode = modal.modes[i];
        Eigen::VectorXd phi = Eigen::VectorXd::Zero(nd);
        const int nv = std::min(nd, static_cast<int>(mode.size()));
        for (int k = 0; k < nv; ++k) phi[k] = mode[k];
        const Eigen::VectorXd Mphi = M * phi;
        const double mnorm = phi.dot(Mphi);
        double eff[3] = {0.0, 0.0, 0.0};
        if (mnorm > 0.0) {
            for (int d = 0; d < 3; ++d) {
                const double gamma = phi.dot(Mr[d]);
                eff[d] = gamma * gamma / mnorm;
            }
        }
        ss << "F " << (i + 1) << " " << f << " " << period << " "
           << eff[0] << " " << eff[1] << " " << eff[2] << " " << mnorm << "\n";
        for (int nid : dm.nodeIds) {
            double u[3] = {0.0, 0.0, 0.0};
            for (int d = 0; d < 3; ++d) {
                const int gd = dm.safeDof(nid, d);
                if (gd >= 0 && gd < static_cast<int>(mode.size()))
                    u[d] = mode[gd];
            }
            ss << "S " << (i + 1) << " " << nid << " "
               << u[0] << " " << u[1] << " " << u[2] << "\n";
        }
    }
    return ss.str();
}

// 计算不出模态时也要写一段，不能让下游"没解析到就当没这回事"——不静默。
inline std::string modalErrorSection3D(const std::string& header,
                                       const std::string& code) {
    std::ostringstream ss;
    ss << header << "\n"
       << "NMODES=0 NDOF=0 BASIS=xyz\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

} // namespace fem
