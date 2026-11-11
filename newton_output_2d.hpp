#pragma once

#include "fem_2d.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

// 2D opt-in 内层 Newton 过程量结果段 [NEWTON2D]。
//
// 与 3D 的 [NEWTON] 段同源：段头 / KIND+STATUS / 头部键名 / 逐迭代 S 行 /
// setprecision(10) / 单位（|R|、|Fext| 为 N，|du| 为 m，lambda 无量纲），
// 唯一数据源是 fem_2d.cpp::solveModel2DNonlinear 循环里 **已经算出** 的量，
// 本函数只格式化，不重算力学、不求解、不改判据与数值路径。
//
// 字段按 2D 事实取舍，不把 3D 的列照抄成 0：
//   * 2D 用稠密/直接求解，没有「符号分解重建」概念，故不写 REFACT；
//   * 2D 只有一个容差 tol（同时管残差判据 tol*max(1,|Fext|) 与位移增量
//     判据 tol*max(1,|U|)），故写单个 TOL，不拆成 TOLF / TOLD；
//   * 2D 主算例没有阶段维度，故不设阶段占位、不凭空造阶段号；
//   * S 行的 |du|=0 只出现在「残差已达标、本次未求解增量」的迭代，是事实值，
//     不是用 0 冒充「算了但没有位移」。
//
// 过程量只说明数值如何收敛，不是规范意义上的设计合格判定；规范没有
// Newton 迭代轮数的合格限值。段只在 example_mixed.exe --newton-detail /
// --newton-detail-only 时写出，不带开关时 .result 与控制台逐位不变。
inline std::string formatNewtonSection2D(const NewtonTrace2D& tr) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << "[NEWTON2D]\n";
    ss << "BASIS=2d_global\n";
    if (tr.steps.empty()) {
        ss << "KIND=none STATUS=not_applicable\n";
        ss << "NSTEP=0 NITER=0 MAXITER=0 TOL=0 CONVERGED=0 "
              "RESFIRST=0 RESLAST=0 LINECUTS=0 STEPREFINE=0 LOADRATIO=0\n";
        ss << "ERROR=no_newton_solution\n";
        ss << "NOTE=本次模型不是 NONLINEAR 分析或没有产生 Newton 迭代；内核没有"
              "算过逐迭代过程量，不输出零行冒充收敛，也不摆空表。\n";
        return ss.str();
    }
    ss << "KIND=newton_cr STATUS=" << (tr.converged ? "ok" : "not_converged")
       << "\n";
    ss << "NSTEP=" << tr.loadSteps
       << " NITER=" << tr.steps.size()
       << " MAXITER=" << tr.maxIter
       << " TOL=" << tr.tol
       << " CONVERGED=" << (tr.converged ? 1 : 0)
       << " RESFIRST=" << tr.initialResidual
       << " RESLAST=" << tr.lastResidual
       << " LINECUTS=" << tr.lineSearchCuts
       << " STEPREFINE=" << tr.stepRefinements
       << " LOADRATIO=" << tr.loadRatio << "\n";
    // S <loadStep> <iter> <lambda> <|R|> <|du|> <|Fext|> <alpha> <lineCuts> <converged>
    for (const auto& s : tr.steps) {
        ss << "S " << s.loadStep << " " << s.iter << " " << s.lambda << " "
           << s.residual << " " << s.du << " " << s.ref << " "
           << s.alpha << " " << s.lineSearchCuts << " "
           << (s.converged ? 1 : 0) << "\n";
    }
    ss << "NOTE=2D 方案C 几何非线性 Newton 逐迭代留痕；字段名与单位同 3D "
          "[NEWTON] 段，10 位有效数字。2D 没有符号分解重建，故不写 REFACT；"
          "只有一个 tol 同时管残差与位移增量判据，故写 TOL 不拆 TOLF/TOLD；"
          "没有阶段维度，不设阶段占位。|du|=0 表示该迭代残差已达标、未求解"
          "增量，不是补零。过程量只说明数值如何收敛，不是规范合格判定。\n";
    if (!tr.failureReason.empty())
        ss << "FAILREASON=" << tr.failureReason << "\n";
    return ss.str();
}

} // namespace fem
