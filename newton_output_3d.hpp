#pragma once
// newton_output_3d.hpp
// ====================
// 把内核**已经算出来**的内层 Newton 过程写成文本段，供 BridgeStudio 解析。
// 这里不做任何力学计算：共转列式一致切线、荷载步、残差收敛判据全在
// fem_3d_nl.hpp 的 solveModel3DNonlinear 里，本文件只格式化同一份数据。
//
// 与 cable_length_output_3d.hpp 的分工（两层迭代不能混为一谈）：
//   · [LENGTH] = **外层**迭代：改模型参数（索无应力长度 L0，或制造线形 ΔX），
//     每轮都要重跑一次完整非线性平衡；控制误差是索力相对偏差 / 几何闭合 m。
//   · [NEWTON] = **内层**迭代：给定模型上的共转列式 Newton + 荷载步平衡，
//     未知量是节点位移 u；控制误差是残差 |R| 相对 |Ftarget| 与位移增量 |du|。
// 一个 [LENGTH] 轮次对应一次 [NEWTON] 全过程；成桥态只暴露**收敛那一轮**的
// 内层过程，施工阶段则每阶段各暴露一段，不把几十次中间求解全堆进结果文件。
//
// 段格式（每行定长字段，数值 %.10g）：
//   <header>                       "[NEWTON]"（cable.txt / .nlstage.txt 通用）
//   KIND=newton_cr|none STATUS=<ok|not_converged|no_iteration|not_applicable>
//   STAGE=<id>                     仅非线性施工阶段算例；成桥态没有这一行
//   NSTEP=<荷载步数> NITER=<实际 Newton 迭代总数> MAXITER=<每步上限>
//        TOLF=<tolForce> TOLD=<tolDisp> CONVERGED=<0|1>
//        REFACT=<符号分解重建次数> RESFIRST=<首个 |R|> RESLAST=<末个 |R|>
//   S <loadStep> <iter> <lambda> <|R|> <|du|> <|Ftarget|> <refactor> <converged>
//        loadStep 从 1 起；iter 是该荷载步内的迭代序号（从 1 起）；
//        lambda = loadStep/NSTEP；|R| / |du| / |Ftarget| 只统计自由自由度；
//        refactor / converged 是 0|1；converged=1 表示该步之后 break。
//   ERROR=<code>                   仅在无内层过程时（STATUS=not_applicable）出现
//
// 口径（必须写进报告，不能只给"收敛了"三个字）：
//   · 收敛判据是 |R| ≤ tolForce·max(|Ftarget|, 1e-12) **或**
//     |du| ≤ tolDisp·max(1, |u|)；两者是"或"，不是"与"。
//   · REFACT 只统计符号分解重建（首次 / 稀疏模式变化），数值 refactorize
//     不算；它是性能诊断量，不是力学量。
//   · 这些过程量说明"数值怎么收敛"，不构成规范意义上的合格判定；
//     规范也没有 Newton 迭代轮数的限值。未收敛必须如实标 not_converged。

#include "fem_3d_nl.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

// STATUS 的判定只有这一处，Python 侧只读不算，避免两边各判一遍漂移。
inline std::string newtonStatus3D(const NewtonTrace3D& tr) {
    if (tr.steps.empty()) return "no_iteration";
    if (!tr.converged) return "not_converged";
    return "ok";
}

inline std::string formatNewtonSection3D(const NewtonTrace3D& tr,
                                         const std::string& header,
                                         int stageId = -1) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << header << "\n";
    ss << "KIND=" << (tr.steps.empty() ? "none" : "newton_cr")
       << " STATUS=" << newtonStatus3D(tr) << "\n";
    if (stageId >= 0) ss << "STAGE=" << stageId << "\n";
    ss << "NSTEP=" << tr.loadSteps
       << " NITER=" << static_cast<int>(tr.steps.size())
       << " MAXITER=" << tr.maxIter
       << " TOLF=" << tr.tolForce
       << " TOLD=" << tr.tolDisp
       << " CONVERGED=" << (tr.converged ? 1 : 0)
       << " REFACT=" << tr.refactorizations
       << " RESFIRST=" << tr.residualFirst
       << " RESLAST=" << tr.residualLast << "\n";
    for (const auto& s : tr.steps)
        ss << "S " << s.loadStep << " " << s.iter << " " << s.lambda << " "
           << s.residual << " " << s.du << " " << s.ref << " "
           << (s.refactor ? 1 : 0) << " " << (s.converged ? 1 : 0) << "\n";
    return ss.str();
}

// 无内层过程也要写一段，不能让下游"没解析到就当没这回事"——不静默。
inline std::string newtonErrorSection3D(const std::string& header,
                                        const std::string& code) {
    std::ostringstream ss;
    ss << header << "\n"
       << "KIND=none STATUS=not_applicable\n"
       << "NSTEP=0 NITER=0 MAXITER=0 TOLF=0 TOLD=0 CONVERGED=0 "
          "REFACT=0 RESFIRST=0 RESLAST=0\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

} // namespace fem
