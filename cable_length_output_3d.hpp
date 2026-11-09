#pragma once
// cable_length_output_3d.hpp
// ==========================
// 把内核**已经算出来**的索长迭代 / 制造线形几何闭合的逐轮过程写成文本段，
// 供 BridgeStudio 解析。这里不做任何力学计算：影响矩阵 Newton、Woodbury
// 低秩 Newton、回溯线搜索全在 fem_3d_nl.hpp 里，本文件只格式化同一份数据。
// 与 modal_output_3d.hpp / buckling_output_3d.hpp 是同一种"暴露内核已有能力"
// 的做法；IterationTrace3D 只记录、不参与判据或步长，数值口径不变。
//
// 段格式（每行定长字段，数值 %.10g）：
//   <header>                       "[LENGTH]"（cable.txt）
//   KIND=<cable_length|erection_geometry|none> STATUS=<ok|not_converged|no_iteration|not_applicable>
//        ERRKIND=<tension_rel|closure_m|->
//        cable_length      = 索长迭代（MODE=TENSION）：影响矩阵 J = M·D 的 Newton，
//                           未知量是各索无应力长度 L0，控制误差是 max 相对索力偏差。
//        erection_geometry = 制造线形几何闭合（MODE=COMPLETION）：Woodbury 低秩
//                           Newton + 回溯线搜索 + 不动点兜底，未知量是节点坐标
//                           偏移 ΔX，控制误差是 max|ΔX + u|（m，几何闭合）。
//   NCABLES=<n> NITER=<轮数> MAXITER=<上限> TOL=<容差> CONVERGED=<0|1>
//   SENS JAVAIL=<0|1> JREBUILD=<重建次数> JDIAGMIN=<..> JDIAGMAX=<..>
//        JAVAIL=1 才有完整影响矩阵 J；几何闭合只有敏感度算子的对角范围。
//   I <it> <err> <geomErr> <stepMax> <stepRelMax> <alpha> <lineSearch> <flags>
//        it 从 1 起；err 的控制量由 ERRKIND 决定（相对偏差 / m，两者不可混用）；
//        geomErr = max|成桥弦长 − 设计弦长| (m)，只是报告量、不作收敛判据；
//        stepMax = 本轮最大修正量 |ΔL0| 或 |ΔX| (m)；stepRelMax 仅索长迭代有；
//        alpha = 线搜索实际步长因子（索长迭代恒 1，靠 stepCap 限幅）；
//        lineSearch = 本轮线搜索评估次数；
//        flags = REBUILD / CAPPED / FALLBACK / CONVERGED 的逗号串，无则 "-"。
//   E <elem> <Ttarget_N> <Tfinal_N> <L0_m> <chord_m> <Ltar_m>
//        每索终值；Ttarget 为 0 表示该模式下不设目标（几何闭合仍给出值）。
//   NJ=<n>                         仅 JAVAIL=1：影响矩阵阶数
//   J <i> <j> <value>             完整 J（i/j 从 0 起），n×n 行
//   ERROR=<code>                   仅在无迭代时（STATUS=not_applicable）出现
//
// 口径（必须写进报告，不能只给"收敛了"三个字）：
//   · 索长迭代收敛判据只有**索力相对偏差** tolTension；几何偏差 geomErr
//     是成桥变形量级，不参与判据。两者单位不同，图上不能画同一根轴。
//   · 几何闭合的收敛判据是 max|ΔX + u| ≤ tolClosure（m），它决定制造线形
//     预抬高；索力偏差在 COMPLETION 里是**结果**而不是判据。
//   · 影响矩阵 J = ∂T_i/∂L0_j 由 M·D 闭式给出（一次分解 + nc 次回代），
//     冻结复用直到残差停滞才重建；JREBUILD 记录重建次数。
//   · 迭代过程量只说明"数值怎么收敛"，不构成规范意义上的合格判定。
//     未收敛必须如实标 STATUS=not_converged，不摆空表、不冒充通过。

#include "fem_3d_nl.hpp"

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace fem {

// STATUS 的判定只有这一处，Python 侧只读不算，避免两边各判一遍漂移。
inline std::string cableLengthStatus3D(const IterationTrace3D& tr) {
    if (tr.steps.empty()) return "no_iteration";
    if (!tr.converged) return "not_converged";
    return "ok";
}

inline std::string formatCableLengthSection3D(
        const IterationTrace3D& tr,
        const std::vector<double>& Tfinal,
        const std::vector<double>& L0,
        const std::vector<double>& chord,
        const std::vector<double>& Ltarget,
        const std::string& header) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << header << "\n";
    const std::string kind = tr.kind.empty() ? "none" : tr.kind;
    const std::string errKind = tr.errKind.empty() ? "-" : tr.errKind;
    ss << "KIND=" << kind << " STATUS=" << cableLengthStatus3D(tr)
       << " ERRKIND=" << errKind << "\n";
    ss << "NCABLES=" << tr.elems.size()
       << " NITER=" << tr.steps.size()
       << " MAXITER=" << tr.maxIter
       << " TOL=" << tr.tol
       << " CONVERGED=" << (tr.converged ? 1 : 0) << "\n";
    ss << "SENS JAVAIL=" << (tr.jAvailable ? 1 : 0)
       << " JREBUILD=" << tr.jRebuilds
       << " JDIAGMIN=" << tr.sensDiagMin
       << " JDIAGMAX=" << tr.sensDiagMax << "\n";

    for (const auto& s : tr.steps) {
        std::string flags;
        auto add = [&](bool b, const char* tag) {
            if (!b) return;
            if (!flags.empty()) flags += ",";
            flags += tag;
        };
        add(s.rebuilt, "REBUILD");
        add(s.capped, "CAPPED");
        add(s.fallback, "FALLBACK");
        add(s.converged, "CONVERGED");
        if (flags.empty()) flags = "-";
        ss << "I " << s.iter << " " << s.err << " " << s.geomErr << " "
           << s.stepMax << " " << s.stepRelMax << " " << s.alpha << " "
           << s.lineSearch << " " << flags << "\n";
    }

    const std::size_t nc = tr.elems.size();
    for (std::size_t k = 0; k < nc; ++k)
        ss << "E " << tr.elems[k] << " "
           << (k < tr.targetTension.size() ? tr.targetTension[k] : 0.0) << " "
           << (k < Tfinal.size() ? Tfinal[k] : 0.0) << " "
           << (k < L0.size() ? L0[k] : 0.0) << " "
           << (k < chord.size() ? chord[k] : 0.0) << " "
           << (k < Ltarget.size() ? Ltarget[k] : 0.0) << "\n";

    if (tr.jAvailable && tr.J.rows() > 0) {
        ss << "NJ=" << tr.J.rows() << "\n";
        for (int i = 0; i < tr.J.rows(); ++i)
            for (int j = 0; j < tr.J.cols(); ++j)
                ss << "J " << i << " " << j << " " << tr.J(i, j) << "\n";
    }
    return ss.str();
}

// 没有迭代（STATIC / 悬索桥找形是闭式反算 L0，一次解完）也要写一段，
// 不能让下游"没解析到就当没这回事"——不静默。
inline std::string cableLengthErrorSection3D(const std::string& header,
                                             const std::string& code) {
    std::ostringstream ss;
    ss << header << "\n"
       << "KIND=none STATUS=not_applicable ERRKIND=-\n"
       << "NCABLES=0 NITER=0 MAXITER=0 TOL=0 CONVERGED=0\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

} // namespace fem
