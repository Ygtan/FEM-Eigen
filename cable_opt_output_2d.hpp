#pragma once

#include "fem_2d.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

// 2D opt-in 平面索力优化结果段 [OPT2D]（第十四次新造力学）。
//
// 与 3D 的 [OPT] 段只是「段骨架同构」，力学不是同一套东西：
//   * 3D 有 J=M·D 闭式影响矩阵、DECKLOAD、塔顶 ΣT·cos 平衡、SYSGUARD
//     三判据与 MINENERGY；2D 只有中心有限差分 Jacobian + IRLS + 投影梯度，
//     目标只做 MINMAXM，约束只做显式索力箱，不写任何 3D 字段占位。
//   * 段名带 2D 后缀，避免与 3D 的 [OPT] 解析链串段。
//
// 本函数只格式化 iterateCableOpt2D() 已经算出的量：不重算力学、不求解、
// 不补零、不写不适用列冒充 3D 字段。段只在 example_mixed.exe
// --opt-detail / --opt-detail-only 时写出，不带开关时 .result 与控制台
// 逐位不变、黄金基线不刷新。
inline std::string formatCableOptSection2D(const CableOptResult2D& r) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << "[OPT2D]\n";
    ss << "BASIS=2d_global\n";

    if (!r.enabled) {
        // 没有 *CABLE_OPT2D 契约：内核根本没有算过索力优化。不输出零行
        // 冒充收敛、不摆空表、不写 3D 的字段占位。
        ss << "KIND=none STATUS=not_applicable\n";
        ss << "NCABLES=0 NMOMENTS=0 NITER=0\n";
        ss << "ERROR=" << (r.error.empty() ? "no_opt_request" : r.error) << "\n";
        ss << "NOTE=本次模型没有 *CABLE_OPT2D 契约；2D 内核没有算过平面索力"
              "优化，不输出零行冒充收敛，也不摆空表。要出真实 [OPT2D] 段必须"
              "显式给出 TMIN/TMAX、逐索 T0 与 M 登记。不是 3D 的 J=M·D 闭式"
              "影响矩阵；无 DECKLOAD / 塔顶 ΣT·cos / SYSGUARD，不造占位；"
              "不等于 3D 的 [OPT]→[LENGTH]→[NEWTON] 三层迭代；命中 T* 属 "
              "[LENGTH2D] 职责；不写回 .dat / .result。\n";
        return ss.str();
    }

    // 契约非法 / 求解失败：写结构化 ERROR=，不摆 I/T/M 空表、不补零、
    // 不静默兜底，也不把半成品过程量拼凑成成功结果。
    if (!r.error.empty()) {
        ss << "KIND=plane_cable_force_opt STATUS=" << r.status << "\n";
        ss << "NCABLES=" << r.ncables << " NMOMENTS=" << r.nmoments
           << " NITER=" << r.niter << "\n";
        ss << "ERROR=" << r.error << "\n";
        ss << "NOTE=契约非法或求解失败，内核没有算出可采纳的 [OPT2D] 过程量，"
              "写 ERROR= 不摆空表、不补零、不静默兜底。不是 3D 的 J=M·D 闭式"
              "影响矩阵；无 DECKLOAD / 塔顶 ΣT·cos / SYSGUARD，不造占位；"
              "不等于 3D 的 [OPT]→[LENGTH]→[NEWTON] 三层迭代；命中 T* 属 "
              "[LENGTH2D] 职责；不写回 .dat / .result。\n";
        return ss.str();
    }

    // 合法段：只格式化 iterateCableOpt2D() 已算出的量，10 位有效数字。
    ss << "KIND=plane_cable_force_opt STATUS=" << r.status << "\n";
    ss << "OBJECTIVE=MINMAXM OBJECTIVE_KIND=NORMALIZED JMETHOD=FINITE_DIFFERENCE\n";
    ss << "NCABLES=" << r.ncables
       << " NMOMENTS=" << r.nmoments
       << " NITER=" << r.niter
       << " IRLS=" << r.irls
       << " MAXITER=" << r.maxIter
       << " IRLSMAX=" << r.irlsMax
       << " TOL=" << r.tol
       << " ROUNDTOL=" << r.roundTol << "\n";
    ss << "TSCALE=" << r.Tscale
       << " MSCALE=" << r.Mscale
       << " CONVERGED=" << (r.converged ? 1 : 0)
       << " CONVBY=" << r.convBy
       << " APPLIED=" << (r.applied ? 1 : 0)
       << " REJECT=" << r.reject << "\n";
    ss << "TMIN=" << r.tmin << " TMAX=" << r.tmax << " BND_SRC=USER\n";
    ss << "MPEAK_KIND=REGISTERED_BEAM_END\n";
    ss << "MPEAK0=" << r.mpeak0 << " MPEAK=" << r.mpeak << "\n";

    // I <irls> <iter> <objective> <mpeak_linearized> <step>
    for (const auto& s : r.steps) {
        ss << "I " << s.irls << " " << s.iter << " " << s.objective
           << " " << s.mpeakLinearized << " " << s.step << "\n";
    }
    // T <k> <elem> <T0_N> <Tstar_N> <Nfinal_N> <Lu0_m> <Lu_m> <Lgeom_m> <EA_N>
    for (std::size_t k = 0; k < r.cables.size(); ++k) {
        const auto& c = r.cables[k];
        ss << "T " << (k + 1) << " " << c.elem << " " << c.t0 << " " << c.tstar
           << " " << c.nfinal << " " << c.Lu0 << " " << c.Lu
           << " " << c.Lgeom << " " << c.EA << "\n";
    }
    // M <i> <elem> <end> <M0_Nm> <Mfinal_Nm>
    for (std::size_t i = 0; i < r.moments.size(); ++i) {
        const auto& mm = r.moments[i];
        ss << "M " << (i + 1) << " " << mm.elem << " " << mm.end
           << " " << mm.m0 << " " << mm.mfinal << "\n";
    }

    // NOTE 是固定中文串：含 5 条边界措辞与「T 行只是过程量」的口径说明。
    ss << "NOTE=2D 平面索力优化（第十四次新造力学）：变量是各索无量纲名义索力 "
          "x=T/Tscale，经下料长度 Lu(x) 落到第十三次本构 "
          "N=EA·((L-Lu)/Lu+istrn)；中心有限差分 Jacobian + "
          "IRLS(p=2/4/8/16) + 投影梯度最小化用户显式登记的 BEAM 端弯矩 "
          "max|M|，不是 3D 的 J=M·D 闭式影响矩阵。2D 无 DECKLOAD / "
          "塔顶 ΣT·cos / SYSGUARD，不造占位；只做 MINMAXM 与显式索力箱 "
          "TMIN/TMAX（BND_SRC=USER），不做 MINENERGY、不识别主梁、不补零、"
          "不静默兜底。本段不等于 3D 的 [OPT]→[LENGTH]→[NEWTON] 三层迭代；"
          "命中 T* 属 [LENGTH2D] 职责，[OPT2D] 只负责定目标。APPLIED=1 仅表示 "
          "T* 被采纳为建议值，T 行只是过程量，不能当设计索力；不写回 .dat / "
          ".result，也不写进 .result。MPEAK0 是用户显式给的 T0 起点实测 "
          "max|M|（守卫基准），优化可行起点 x0=clamp(T0/Tscale,TMIN/Tscale,"
          "TMAX/Tscale) 保证落入显式索力箱；MPEAK 是最终实测 max|Mfinal|"
          "（非线性化预测），10 位有效数字。\n";
    return ss.str();
}

} // namespace fem
