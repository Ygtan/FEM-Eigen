#pragma once

#include "fem_2d.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

// 2D opt-in 平面索长闭合结果段 [LENGTH2D]。
//
// 与 3D 的 [LENGTH] 段只是「段骨架同构」，力学不是同一套东西：
//   * 3D 的 L0 命中 T* 用垂度缩短 s(T) + Ernst 切线 + 闭式影响矩阵 J=M·D；
//     2D 是无垂度线性杆，闭式初值 Lu0=L*/(1+T*/EA)，用中心有限差分
//     Jacobian + LM 阻尼最小二乘，唯一判据是索力相对偏差 tension_rel。
//   * 3D 另一种 KIND=erection_geometry 控制误差是几何闭合量 closure_m；
//     2D 没有三维制造线形 ΔX+u，不写 geomErr、不写 SENS/NJ/J 行。
//   * 2D 主算例没有阶段维度，不设 STAGE 占位。
// 段名带 2D 后缀是为了不撞 3D 的 [LENGTH] 解析链（fem_eigen_cable.py
// 的段分派会把同名段的行送进 3D 解析器，造成 S/I/E 行串段）。
//
// 本函数只格式化 iterateCableLength2D() 已经算出的过程量：不重算力学、
// 不求解、不补零、不写不适用列冒充 3D 字段。过程量只说明数值如何收敛，
// 不是规范意义上的设计合格判定。段只在 example_mixed.exe
// --length-detail / --length-detail-only 时写出，不带开关时 .result 与
// 控制台逐位不变、黄金基线不刷新。
inline std::string formatCableLengthSection2D(const LengthIterationResult2D& r) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << "[LENGTH2D]\n";
    ss << "BASIS=2d_global\n";

    if (!r.enabled) {
        // 没有 *CABLE_LENGTH2D 契约：内核根本没有算过索长闭合。
        // 不输出零行冒充收敛、不摆空表、不写 3D 的字段占位。
        ss << "KIND=none STATUS=not_applicable\n";
        ss << "NCABLES=0 NITER=0 MAXITER=" << r.maxIter
           << " TOLT=" << r.tolT << " CONVERGED=0\n";
        ss << "ERROR=" << (r.error.empty() ? "no_length_request" : r.error) << "\n";
        ss << "NOTE=本次模型没有 *CABLE_LENGTH2D 目标索力契约；2D 内核没有"
              "算过无应力下料长度闭合，不输出零行冒充收敛，也不摆空表。"
              "要出真实过程量必须显式给出每根 CABLE 的目标索力 T*。\n";
        return ss.str();
    }

    ss << "KIND=plane_cable_length STATUS=" << r.status << "\n";
    ss << "ERRKIND=tension_rel ANALYSIS=" << r.analysis
       << " JMETHOD=FINITE_DIFFERENCE\n";
    ss << "NCABLES=" << r.ncables
       << " NITER=" << r.iterations
       << " MAXITER=" << r.maxIter
       << " TOLT=" << r.tolT
       << " CONVERGED=" << (r.converged ? 1 : 0) << "\n";
    ss << "ERRFIRST=" << r.errFirst << " ERRLAST=" << r.errLast << "\n";

    if (r.steps.empty()) {
        // 第一次完整求解都没跑完：既没有索力也没有过程量可写。
        ss << "ERROR=" << (r.error.empty() ? "solve_failed" : r.error) << "\n";
        if (!r.failureReason.empty())
            ss << "FAILREASON=" << r.failureReason << "\n";
        ss << "NOTE=内核在第一次完整求解中失败，没有逐轮残差与成桥索力；"
              "不写 I/E 零行冒充过程量。\n";
        return ss.str();
    }

    // I <it> <errMax> <errNorm> <stepMax> <alpha>
    for (const auto& s : r.steps) {
        ss << "I " << s.iter << " " << s.errMax << " " << s.errNorm
           << " " << s.stepMax << " " << s.alpha << "\n";
    }
    // E <elem> <Ttarget_N> <Tfinal_N> <Lu0_m> <Lu_m> <Lgeom_m> <EA_N>
    for (const auto& c : r.cables) {
        ss << "E " << c.elem << " " << c.target << " " << c.finalT
           << " " << c.Lu0 << " " << c.Lu << " " << c.Lgeom
           << " " << c.EA << "\n";
    }

    ss << "NOTE=2D 平面索长闭合（首次真正新造的 2D 力学）：无垂度线性杆"
          "本构 N=EA·((L-Lu)/Lu+istrn)，闭式初值 Lu0=L*/(1+T*/EA)（L* 为"
          "几何初始弦长）；中心有限差分 Jacobian + LM 阻尼最小二乘，不是"
          "3D 的闭式影响矩阵 J=M·D。唯一判据是索力相对偏差"
          "max|N/T*-1|<=TOLT（ERRKIND=tension_rel，无量纲），10 位有效数字。"
          "2D 没有垂度 s(T) / Ernst 等效弹模 / 共转成桥态 / 塔顶 ΣT·cos"
          "平衡 / DECKLOAD / SYSGUARD / 三维制造线形，也不等于 3D 的"
          "[OPT]→[LENGTH]→[NEWTON] 三层迭代；没有阶段维度，不设 STAGE 占位。"
          "I 行 stepMax=0 是收敛轮未解增量的事实值，不是补零；E 行 Lu 在"
          "未收敛时是最后一轮的未采纳值，不是成桥下料长度。过程量只说明"
          "数值如何收敛，不是规范合格判定。\n";
    if (!r.converged) {
        ss << "ERROR=" << (r.error.empty() ? "not_converged" : r.error) << "\n";
        if (!r.failureReason.empty())
            ss << "FAILREASON=" << r.failureReason << "\n";
    }
    return ss.str();
}

} // namespace fem
