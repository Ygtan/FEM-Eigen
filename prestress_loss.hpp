#pragma once

#include <cmath>
#include <stdexcept>

namespace fem {

// 预应力损失分项（JTG 3362-2018 第 6.2 条）。单位统一 SI：Pa / m / rad。
// 损失编号与规范一致：
//   σl1 摩擦（6.2.2）、σl2 锚具变形/回缩/接缝压缩（6.2.3）、σl3 温差（6.2.4）、
//   σl4 弹性压缩（6.2.5）、σl5 松弛（6.2.6）、σl6 收缩徐变（6.2.7）。
// 注意：与旧版 JTG D62-2004 相比，σl1/σl2 的物理含义互换（摩擦为 σl1）。
struct PrestressLossInput {
    // 材料 / 几何
    double sigmaCon = 0.0;    // 张拉控制应力 (Pa)
    double Ep = 0.0;          // 预应力筋弹性模量 (Pa)
    double Ec = 0.0;          // 混凝土弹性模量 (Pa)
    double fpk = 0.0;         // 预应力筋抗拉强度标准值 (Pa)

    // σl1 摩擦（6.2.2）：σl1 = σcon·[1 - e^-(μθ + kx)]
    double mu = 0.0;          // 摩擦系数 μ（表 6.2.2）
    double kappa = 0.0;       // 每米局部偏差影响系数 k（表 6.2.2，1/m）
    double theta = 0.0;       // 曲线管道切线夹角之和 θ（rad）
    double x = 0.0;           // 张拉端至计算截面管道长度 x（m）

    // σl2 锚具变形、钢筋回缩、接缝压缩（6.2.3）：σl2 = ΣΔl/l·Ep
    double sumDeltaL = 0.0;   // ΣΔl（表 6.2.3，m）
    double tendonLength = 0.0;// l 张拉端至锚固端距离（m）

    // σl3 温差（6.2.4）：σl3 = 2(t2 - t1)（输入 °C，公式内换算为 Pa）
    //   混凝土加热养护（t2 > t1）时为正损失。
    double t1 = 0.0;          // 张拉钢筋时制造场地温度（°C）
    double t2 = 0.0;          // 混凝土加热养护时受拉钢筋最高温度（°C）

    // σl4 弹性压缩（6.2.5）：σl4 = αEP·σpc4
    //   后张分批张拉：σpc4 = ΣΔσpc；先张法：σpc4 = σpc
    double sigmaPc4 = 0.0;

    // σl5 松弛（6.2.6）：钢丝、钢绞线用 6.2.6-1；预应力螺纹钢筋用 6.2.6-2/3。
    //   钢丝、钢绞线：σl5 = Ψ·ζ·(0.52·σpe/fpk - 0.26)·σpe
    //   预应力螺纹钢筋：一次张拉 σl5 = 0.05·σcon；超张拉 σl5 = 0.035·σcon
    bool threadedBar = false; // true=预应力螺纹钢筋（6.2.6-2/3），false=钢丝/钢绞线（6.2.6-1）
    double psi = 1.0;         // 张拉系数 Ψ：一次张拉 1.0，超张拉 0.9（螺纹钢筋据此取 0.05/0.035）
    double zeta = 0.3;        // 松弛系数 ζ（仅钢丝/钢绞线）：Ⅰ级普通松弛 1.0，Ⅱ级低松弛 0.3

    // 张拉方式：true=后张法（默认），false=先张法。
    //   后张法 σpe = σcon - σl1 - σl2 - σl4；
    //   先张法 σpe = σcon - σl2（无管道摩擦，扣锚具变形/回缩）。
    bool postTensioned = true;

    // σl6 收缩徐变（6.2.7-1）：
    //   σl6 = 0.9[Ep·εcs(t,t0) + αEP·σpc6·φ(t,t0)] / (1 + 15ρ·ρps)
    double sigmaPc6 = 0.0;    // 传力锚固时全部预加力在钢筋重心处产生的混凝土法向压应力 σpc
    double epsCs = 0.0;       // εcs(t,t0) 收缩应变
    double phi = 0.0;         // φ(t,t0) 徐变系数
    double rho = 0.0;         // ρ = (Ap+As)/A 配筋率
    double rhoPs = 1.0;       // ρps = 1 + eps²/i²
};

struct PrestressLossResult {
    double sigmaL1 = 0.0;   // 摩擦
    double sigmaL2 = 0.0;   // 锚具变形/回缩/接缝压缩
    double sigmaL3 = 0.0;   // 温差
    double sigmaL4 = 0.0;   // 弹性压缩
    double sigmaL5 = 0.0;   // 松弛
    double sigmaL6 = 0.0;   // 收缩徐变
    double sigmaPe = 0.0;   // 传力锚固应力（后张法：σcon - σl1 - σl2 - σl4）
};

inline PrestressLossResult computePrestressLosses(const PrestressLossInput& in) {
    if (in.tendonLength <= 0.0)
        throw std::runtime_error("预应力筋长度 l 必须 > 0");
    if (!in.threadedBar && in.fpk <= 0.0)
        throw std::runtime_error("σl5 需要预应力筋抗拉强度标准值 fpk > 0");
    const double alphaEp = (in.Ec > 0.0) ? in.Ep / in.Ec : 0.0;

    PrestressLossResult r;
    r.sigmaL1 = in.sigmaCon * (1.0 - std::exp(-(in.mu * in.theta + in.kappa * in.x)));
    r.sigmaL2 = in.sumDeltaL / in.tendonLength * in.Ep;
    r.sigmaL3 = 2.0 * (in.t2 - in.t1) * 1.0e6;   // °C -> Pa
    r.sigmaL4 = alphaEp * in.sigmaPc4;

    // 传力锚固应力：后张法 σpe = σcon - σl1 - σl2 - σl4；先张法 σpe = σcon - σl2。
    r.sigmaPe = in.postTensioned
        ? in.sigmaCon - r.sigmaL1 - r.sigmaL2 - r.sigmaL4
        : in.sigmaCon - r.sigmaL2;
    if (in.threadedBar) {
        // 6.2.6-2 / 6.2.6-3：预应力螺纹钢筋
        r.sigmaL5 = (in.psi < 1.0 ? 0.035 : 0.05) * in.sigmaCon;
    } else {
        // 6.2.6-1：预应力钢丝、钢绞线
        r.sigmaL5 = in.psi * in.zeta * (0.52 * r.sigmaPe / in.fpk - 0.26) * r.sigmaPe;
    }
    r.sigmaL6 = 0.9 * (in.Ep * in.epsCs + alphaEp * in.sigmaPc6 * in.phi)
              / (1.0 + 15.0 * in.rho * in.rhoPs);
    return r;
}

// 有效预应力比 σpe/σcon（供 *PRESTRESS 的 LOSSREF 用，等价于 1 - Σσl/σcon）。
inline double effectivePrestressRatio(const PrestressLossInput& in) {
    if (in.sigmaCon <= 0.0)
        throw std::runtime_error("有效预应力比需要张拉控制应力 sigmaCon > 0");
    return computePrestressLosses(in).sigmaPe / in.sigmaCon;
}

} // namespace fem
