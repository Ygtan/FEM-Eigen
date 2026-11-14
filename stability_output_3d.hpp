#pragma once
// stability_output_3d.hpp
// ========================
// 把 stability_path_3d.hpp 已经算出来的球面弧长法平衡路径写成文本段，
// 供 BridgeStudio 解析。这里**不求解、不重算、不补零**：λ、控制点位移、
// 残差、切线奇异指标、收敛标志全部来自内核的逐平衡点记录。
//
// B 类后续项接入规格 §2.3 第十七次：STABILITY=1 是 cable_solver.exe 的
// 独立 opt-in 模式（默认全关），段置于结果文件 [END] 前。既有成桥态闭环
// 的 [BUCKLING] ERROR=nonlinear_state_not_supported 原样保留——线性屈曲对
// 共转成桥态不适用是独立的物理结论，不能与全过程极值点混为一谈。
//
// 段格式（每行定长字段，数值 %.14g，与 CB3D17 snap_probe 同源）：
//   [STABILITY]
//   KIND=spherical_arc_length STATUS=<ok|not_converged>
//       BASIS=reference_configuration PATHMETHOD=SPHERE
//   NSTEP=<n> NPOINT=<n> CONVERGED=<0|1> CTRLNODE=<id> CTRLDOF=<1..6>
//       STRAIN=<0|1> PSI=<w>
//   LAMBDA_FINAL=<lambda> ERROR=<none|code>
//   I STEP=<i> LAMBDA=<λ> CTRL=<u_ctrl> RES=<|R|> ITER=<n>
//       DS=<ds> PIVOT=<min|diag(U)|/max|diag(U)|> CONVERGED=<0|1>
//
// 口径（必须写进报告，不能只给一个 λ 数字）：
//   · λ 是参考荷载 Fref 的整体放大倍数，不是规范稳定安全系数；
//     JTG 3362 无 λ 合格限值，段里不给「通过 / 不通过」判定。
//   · PIVOT 是**全矩阵** LU 的 min|diag(U)| / max|diag(U)|：约束自由度
//     被置 1、自由自由度对角量级约 EA/L，故绝对值本来就很小
//     （两杆桁架基准题约 5e-9 且沿路径近似恒定）。它是量级比、
//     **不是奇异判据**；极限点由 λ 转折判定，不看 PIVOT 绝对值。
//   · STRAIN=0 工程应变 N=EA(L-L0)/L0（fem_3d_nl.hpp 口径），
//     STRAIN=1 对数应变 N=EA*ln(L/L0)（ANSYS LINK180 KEYOPT(2)=1 口径）。
//   · 未收敛 / 缺参数 / 无参考荷载一律走 stabilityErrorSection3D 的
//     ERROR=<code>，下游按中文说明写「做不了」，不摆空表、不补 0。

#include "stability_path_3d.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

inline std::string formatStabilitySection3D(const ArcLengthResult3D& r,
                                            const ArcLengthOptions3D& opt,
                                            const std::string& header) {
    std::ostringstream ss;
    ss << std::setprecision(14);
    ss << header << "\n";
    ss << "KIND=spherical_arc_length STATUS="
       << (r.converged ? "ok" : "not_converged")
       << " BASIS=reference_configuration PATHMETHOD=SPHERE\n";
    ss << "NSTEP=" << r.steps << " NPOINT=" << r.points.size()
       << " CONVERGED=" << (r.converged ? 1 : 0)
       << " CTRLNODE=" << opt.ctrlNode << " CTRLDOF=" << opt.ctrlDof
       << " STRAIN=" << opt.strainMeasure << " PSI=" << opt.psi << "\n";
    ss << "LAMBDA_FINAL=" << r.lambda
       << " ERROR=" << (r.error.empty() ? "none" : r.error) << "\n";
    for (const auto& p : r.points) {
        ss << "I STEP=" << p.step << " LAMBDA=" << p.lambda
           << " CTRL=" << p.ctrlDisp << " RES=" << p.residual
           << " ITER=" << p.iterations << " DS=" << p.ds
           << " PIVOT=" << p.pivotRatio
           << " CONVERGED=" << (p.converged ? 1 : 0) << "\n";
    }
    return ss.str();
}

inline std::string stabilityErrorSection3D(const std::string& header,
                                           const std::string& code) {
    std::ostringstream ss;
    ss << header << "\n"
       << "KIND=none STATUS=not_applicable BASIS=reference_configuration\n"
       << "NSTEP=0 NPOINT=0 CONVERGED=0\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

}  // namespace fem
