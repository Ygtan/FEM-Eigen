#pragma once

#include "fem_2d.hpp"

#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

// 2D opt-in 分项预应力损失结果段 [PRESTRESSLOSS2D]。
//
// 与 3D 的 stage prestress loss 段同源：SIGMA_CON / SL1..SL6 / SPE / RATIO
// 字段名、顺序、含义与单位（Pa，RATIO 无量纲）完全一致，同为 12 位有效数字，
// 公式唯一来源仍是 prestress_loss.hpp::computePrestressLosses。
//
// 两处差别都是 2D 的事实边界，不是能力补齐：
//   * 2D 没有阶段维度：非阶段 *PRESTRESS 束不属于任何施工阶段，行首用
//     P<编号>（输入行首的显式束号），不写 STAGE<n> 占位，不凭空造阶段号；
//   * 只收录 LOSSREF 路径的束。只给 LOSS= 的束内核没有算过分项，
//     不输出一堆 0 冒充「各项损失为 0」。
//
// 本函数只读 m.prestressLossDetails，不重算力学量、不求解、不改数值路径。
inline std::string formatPrestressLossSection2D(const Model2D& m) {
    std::ostringstream ss;
    ss << "[PRESTRESSLOSS2D]\n";
    if (m.prestressLossDetails.empty()) {
        ss << "ERROR=no_lossref_prestress\n";
        ss << "NOTE=本次没有引用 *PRESTRESS_LOSS 表的 2D 预应力束；内核没有算过"
              "分项损失，不输出零值冒充各项损失为 0。\n";
        return ss.str();
    }
    ss << "NPRESTRESS=" << m.prestressLossDetails.size()
       << " BASIS=2d_global\n";
    ss << "NOTE=2D 非阶段 *PRESTRESS LOSSREF 束的分项损失明细；字段名与单位"
          "（Pa，RATIO 无量纲）同 3D stage prestress loss 段，12 位有效数字。"
          "2D 没有阶段维度，不设 STAGE 占位、不造阶段号；P<id> 是输入行首的显式束号。"
          "σpe 后张=σcon-σl1-σl2-σl4、先张=σcon-σl2；"
          "σl3/σl5/σl6 是后期损失，不进 σpe 也不进等效荷载 Pe=P*(σpe/σcon)。"
          "σl2 只按直线筋 ΣΔl/l*Ep（曲线筋反摩擦未实现）；"
          "σl6 的 εcs/φ 需按附录 C 查表输入，程序不按龄期自动发展。\n";
    ss << std::setprecision(12);
    for (const auto& [id, d, sc] : m.prestressLossDetails) {
        ss << "P" << id
           << " SIGMA_CON=" << sc
           << " SL1=" << d.sigmaL1
           << " SL2=" << d.sigmaL2
           << " SL3=" << d.sigmaL3
           << " SL4=" << d.sigmaL4
           << " SL5=" << d.sigmaL5
           << " SL6=" << d.sigmaL6
           << " SPE=" << d.sigmaPe
           << " RATIO=" << (sc > 0.0 ? d.sigmaPe / sc : 0.0)
           << "\n";
    }
    return ss.str();
}

} // namespace fem
