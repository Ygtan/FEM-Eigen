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
          "σl2 默认按直线筋 ΣΔl/l*Ep；opt-in L2ANTI=1 按附录 G 反摩擦、"
          "代表值 x=0（张拉端保守），BRANCH=1 为 lf≤l、2 为 lf>l、0 为未启用或无摩擦退化；"
          "σl6 默认由 ECS/PHI 用户直给；opt-in C_SKIP=1 显式声明不参与（C_SOURCE=2），"
          "C2D=1 用附录 C 显式表值龄期（C_SOURCE=1），表值由用户核对录入、"
          "程序不调用规范拟合函数兜底，缺任一表值报错。\n";
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
           << " RATIO=" << (sc > 0.0 ? d.sigmaPe / sc : 0.0);
        // 附录 G 反摩擦 opt-in：只在启用且 lf>0 的两支追加中间量，
        // 未启用 / 无摩擦退化的既有行字段与顺序逐位不变。
        if (d.l2Branch >= 1) {
            ss << " LF=" << d.lf
               << " DSL2D=" << d.deltaSigmaD
               << " DSL2PEAK=" << d.deltaSigmaPeak
               << " BRANCH=" << d.l2Branch;
        }
        // 附录 C 收缩徐变 opt-in：cSource>=1 时追加来源与中间量，
        // 默认路径（cSource=0）既有行字段与顺序逐位不变。
        if (d.cSource >= 1) {
            ss << " C_SOURCE=" << d.cSource;
            if (d.cSource == 1) {
                ss << " ECS_DEV=" << d.epsCsDev
                   << " PHI_DEV=" << d.phiDev
                   << " BETAS=" << d.betaS
                   << " BETAS0=" << d.betaS0
                   << " BETAC=" << d.betaC
                   << " BETAH=" << d.betaH
                   << " RHREP=" << d.cHRepresentative;
            }
        }
        ss << "\n";
    }
    return ss.str();
}

// 第二十五次 opt-in：2D **阶段**分项预应力损失结果段。
//
// 与上面的全局段同用段名 [PRESTRESSLOSS2D]，但头行 BASIS=2d_stage、
// 行首 STAGE<n> P<id>；由 example_stage.exe --stage-loss-detail 写到独立
// <stem>.stage.loss.txt，与全局 <stem>.loss.txt 不同文件、互不混行。
// 字段顺序、单位（Pa / 无量纲）、12 位有效数字与 3D
// example_3d_grillage.cpp 的 stage prestress loss 段完全一致；公式唯一来源
// 仍是 prestress_loss.hpp::computePrestressLosses，Python 侧只解析不重算。
inline std::string formatStagePrestressLossSection2D(
    const std::vector<StageResult>& stages) {
    std::ostringstream ss;
    ss << "[PRESTRESSLOSS2D]\n";
    std::size_t nRows = 0;
    for (const auto& sr : stages)
        nRows += sr.prestressLossDetails.size();
    if (nRows == 0) {
        ss << "ERROR=no_lossref_stage_prestress\n";
        ss << "NOTE=本次没有引用 *PRESTRESS_LOSS 表的 2D 阶段预应力束；"
              "内核没有算过阶段分项损失，不输出零值冒充各项损失为 0。\n";
        return ss.str();
    }
    ss << "NPRESTRESS=" << nRows << " BASIS=2d_stage\n";
    ss << "NOTE=2D 阶段 *STAGE_PRESTRESS LOSSREF 束的分项损失明细；字段名与"
          "单位（Pa，RATIO 无量纲）同 3D stage prestress loss 段，12 位有效数字。"
          "STAGE<n> 是显式施工阶段号，P<id> 是输入行首的显式束号；"
          "组织键 (阶段号, 束号) 与 3D 一致。σpe 后张=σcon-σl1-σl2-σl4、"
          "先张=σcon-σl2；σl3/σl5/σl6 是后期损失，不进 σpe 也不进等效荷载 "
          "Pe=P*(σpe/σcon)。σl2 默认按直线筋 ΣΔl/l*Ep；opt-in L2ANTI=1 按"
          "附录 G 反摩擦、代表值 x=0（张拉端保守），BRANCH=1 为 lf≤l、2 为 lf>l、"
          "0 为未启用或无摩擦退化；σl6 默认由 ECS/PHI 用户直给；opt-in "
          "C_SKIP=1 显式声明不参与（C_SOURCE=2），C2D=1 用附录 C 显式表值龄期"
          "（C_SOURCE=1），表值由用户核对录入、程序不调用规范拟合函数兜底，"
          "缺任一表值报错。\n";
    ss << std::setprecision(12);
    for (const auto& sr : stages) {
        for (const auto& [id, d, sc] : sr.prestressLossDetails) {
            ss << "STAGE" << sr.stage << " P" << id
               << " SIGMA_CON=" << sc
               << " SL1=" << d.sigmaL1
               << " SL2=" << d.sigmaL2
               << " SL3=" << d.sigmaL3
               << " SL4=" << d.sigmaL4
               << " SL5=" << d.sigmaL5
               << " SL6=" << d.sigmaL6
               << " SPE=" << d.sigmaPe
               << " RATIO=" << (sc > 0.0 ? d.sigmaPe / sc : 0.0);
            // 与全局段同源的可选后缀：未启用时字段与顺序逐位不变。
            if (d.l2Branch >= 1) {
                ss << " LF=" << d.lf
                   << " DSL2D=" << d.deltaSigmaD
                   << " DSL2PEAK=" << d.deltaSigmaPeak
                   << " BRANCH=" << d.l2Branch;
            }
            if (d.cSource >= 1) {
                ss << " C_SOURCE=" << d.cSource;
                if (d.cSource == 1) {
                    ss << " ECS_DEV=" << d.epsCsDev
                       << " PHI_DEV=" << d.phiDev
                       << " BETAS=" << d.betaS
                       << " BETAS0=" << d.betaS0
                       << " BETAC=" << d.betaC
                       << " BETAH=" << d.betaH
                       << " RHREP=" << d.cHRepresentative;
                }
            }
            ss << "\n";
        }
    }
    return ss.str();
}

// 第二十七次 opt-in：2D **阶段有效预加力**结果段 [PRESTRESSFORCE2D]。
//
// 补的是 3D 早已有的 `stage prestress force` 段（`STAGE<n> P<id> Pe=`）在 2D
// 侧缺失的输出能力——力学上 2D 一直在算 Pe=P·loss（processPrestress），
// 只是从未写进 StageResult，前端拿不到「分阶段看哪批束何时张拉、有效预应力
// 相对张拉控制应力还剩多少」。
//
// 口径：数值唯一来源是 processPrestress 里那一次 `Pe = ps.P * loss`，
// 本函数只读 StageResult::prestressForces，不重算、不反推、不补零；
// loss 已含 LOSSREF 的 σpe/σcon 解析结果，与 3D 的 ps.P*ps.loss 同一公式。
// 字段顺序、单位（N / Pa / 无量纲）、12 位有效数字固定；
// 与 3D 段对应关系：2D 的 PE= 即 3D 的 Pe=，同参同值。
// 只给 LOSS= 的束 SIGMA_CON=0（未登记控制应力），不拿 P/Ap 反造一个。
inline std::string formatStagePrestressForceSection2D(
    const std::vector<StageResult>& stages) {
    std::ostringstream ss;
    ss << "[PRESTRESSFORCE2D]\n";
    std::size_t nRows = 0;
    for (const auto& sr : stages)
        nRows += sr.prestressForces.size();
    if (nRows == 0) {
        ss << "ERROR=no_stage_prestress\n";
        ss << "NOTE=本次施工阶段没有 *STAGE_PRESTRESS 束；没有张拉的束就没有"
              "有效预加力，不输出零值冒充已张拉。\n";
        return ss.str();
    }
    ss << "NPRESTRESS=" << nRows << " BASIS=2d_stage\n";
    ss << "NOTE=2D 施工阶段每束有效预加力明细；组织键 (阶段号, 束号) 与 3D "
          "stage prestress force 段一致，2D 的 PE= 即 3D 的 Pe=，同参同值。"
          "P 是张拉控制力（N）、SIGMA_CON 是张拉控制应力（Pa，只给 LOSS= 时为 0）、"
          "LOSS 是有效预应力系数（LOSS= 原值或 LOSSREF 解析出的 σpe/σcon，无量纲）、"
          "PE=P*LOSS（N）。不输出 PE/SIGMA_CON 之类量纲不成立的比值（那是 m²），"
          "有效比就是 LOSS 本身。本段只报过程量，不参与任何求解；默认不写文件。\n";
    ss << std::setprecision(12);
    for (const auto& sr : stages) {
        for (const auto& pf : sr.prestressForces) {
            ss << "STAGE" << sr.stage << " P" << pf.id
               << " SIGMA_CON=" << pf.sigmaCon
               << " P=" << pf.P
               << " LOSS=" << pf.loss
               << " PE=" << pf.Pe
               << "\n";
        }
    }
    return ss.str();
}

} // namespace fem
