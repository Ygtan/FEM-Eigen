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
// 第十九次 ①：新增 BASIS=completion_state（默认仍是
// reference_configuration）。两种 BASIS 的 λ **不是同一个量**：
//   · reference_configuration：λ 放大总参考荷载，u(0)=0；
//   · completion_state：λ 只放大 *STABILITY_LOAD 增量荷载，
//     u(0)=u0 为成桥态收敛位移、F0 为不缩放的恒定荷载。
// completion_state 段额外写一行 STATE CTRL0= / F0_NORM= / FINC_NORM= /
// RES0=，供下游与 [U] 段自洽对拍；缺初态不补 0、不摆空表。
//
// 段格式（每行定长字段，数值 %.14g，与 CB3D17 snap_probe 同源）：
//   [STABILITY]
//   KIND=spherical_arc_length STATUS=<ok|not_converged>
//       BASIS=<reference_configuration|completion_state> PATHMETHOD=SPHERE
//   NSTEP=<n> NPOINT=<n> CONVERGED=<0|1> CTRLNODE=<id> CTRLDOF=<1..6>
//       STRAIN=<0|1> PSI=<w>
//   STATE CTRL0=<u0_ctrl> F0_NORM=<|F0|> FINC_NORM=<|Finc|> RES0=<|R(0)|>
//                                       ^ 仅 completion_state 基准有此行
//   LAMBDA_FINAL=<lambda> ERROR=<none|code>
//   I STEP=<i> LAMBDA=<λ> CTRL=<u_ctrl-u0_ctrl> RES=<|R|> ITER=<n>
//       DS=<ds> PIVOT=<min|diag(U)|/max|diag(U)|> CONVERGED=<0|1>
//   （CTRL 是**相对起点的增量**；参考构形起点为 0，与既有输出逐位一致）
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

// 第十九次 ①：completion_state 基准的段头补充信息。默认 hasState=false，
// 既有调用（3 参数 formatStabilitySection3D）逐位不变。
struct StabilityStateHeader3D {
    bool hasState = false;
    double ctrl0 = 0.0;     // u0(ctrl)：成桥态起点的绝对控制位移
    double f0Norm = 0.0;    // |F0|
    double fincNorm = 0.0;  // |Finc|
    // 第二十次 ③b：施工阶段中间态起点。stageId < 0 时**不输出任何行**，
    // 既有 completion_state / reference_configuration 调用逐位不变。
    int stageId = -1;
    int stageElems = 0;     // 该阶段激活单元数（不是全桥单元数）
    int stageDofs = 0;      // 该阶段自由度数
};

inline std::string formatStabilitySection3D(
    const ArcLengthResult3D& r,
    const ArcLengthOptions3D& opt,
    const std::string& header,
    const std::string& basis = "reference_configuration",
    const StabilityStateHeader3D& sh = {}) {
    std::ostringstream ss;
    ss << std::setprecision(14);
    ss << header << "\n";
    // Wempner–Riks 正交轨迹只在 STAGE 起点 opt-in；SPHERE 字符串逐位不变。
    const bool isRiks = (opt.pathMethod == 1);
    ss << "KIND="
       << (isRiks ? "wempner_riks_orthogonal" : "spherical_arc_length")
       << " STATUS=" << (r.converged ? "ok" : "not_converged")
       << " BASIS=" << basis
       << " PATHMETHOD=" << (isRiks ? "RIKS" : "SPHERE") << "\n";
    ss << "NSTEP=" << r.steps << " NPOINT=" << r.points.size()
       << " CONVERGED=" << (r.converged ? 1 : 0)
       << " CTRLNODE=" << opt.ctrlNode << " CTRLDOF=" << opt.ctrlDof
       << " STRAIN=" << opt.strainMeasure << " PSI=" << opt.psi << "\n";
    if (sh.hasState)
        ss << "STATE CTRL0=" << sh.ctrl0
           << " F0_NORM=" << sh.f0Norm
           << " FINC_NORM=" << sh.fincNorm
           << " RES0=" << r.res0 << "\n";
    // 第二十次 ③b：只有 STABASIS=STAGE 才写这行，证明起点是某个施工
    // 阶段中间态（已激活单元数 / 自由度数），不是成桥态、也不是参考构形。
    if (sh.stageId >= 0)
        ss << "STAGE=" << sh.stageId
           << " NELEM=" << sh.stageElems
           << " NDOF=" << sh.stageDofs << "\n";
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

inline std::string stabilityErrorSection3D(
    const std::string& header,
    const std::string& code,
    const std::string& basis = "reference_configuration") {
    std::ostringstream ss;
    ss << header << "\n"
       << "KIND=none STATUS=not_applicable BASIS=" << basis << "\n"
       << "NSTEP=0 NPOINT=0 CONVERGED=0\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

}  // namespace fem
