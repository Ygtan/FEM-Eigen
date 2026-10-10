#pragma once

#include <cmath>
#include <stdexcept>

namespace fem {

// 预应力损失分项（JTG D62 常用简化式，单位统一为 SI：Pa / m）。
// σl1 锚具变形回缩、σl2 摩擦、σl4 混凝土弹性压缩、σl5 钢筋松弛、
// σl6 混凝土收缩徐变；σpe = σcon - Σσl。
struct PrestressLossInput {
    double sigmaCon = 0.0;       // 张拉控制应力 (Pa)
    double Ep = 0.0;             // 预应力筋弹性模量 (Pa)
    double Ec = 0.0;             // 混凝土弹性模量 (Pa)
    double anchorageSlip = 0.0;  // 锚具变形 + 回缩 (m)
    double tendonLength = 0.0;   // 预应力筋长度 (m)
    double mu = 0.0;             // 孔道摩擦系数 μ
    double kappa = 0.0;          // 孔道偏差系数 κ (1/m)
    double theta = 0.0;          // 累计转角 (rad)
    double x = 0.0;              // 计算点距张拉端距离 (m)
    double sigmaPc = 0.0;        // 预应力筋重心处混凝土压应力 (Pa)
    double creepPhi = 0.0;       // 徐变系数 φ
    double epsCs = 0.0;          // 混凝土收缩应变 εcs
    double relaxation = 0.0;     // 松弛系数（低松弛钢可取 0.025）
};

struct PrestressLossResult {
    double sigmaL1 = 0.0;   // 锚具回缩损失
    double sigmaL2 = 0.0;   // 摩擦损失
    double sigmaL4 = 0.0;   // 混凝土弹性压缩损失
    double sigmaL5 = 0.0;   // 钢筋松弛损失
    double sigmaL6 = 0.0;   // 收缩徐变损失
    double sigmaPe = 0.0;   // 有效预应力
};

inline PrestressLossResult computePrestressLosses(const PrestressLossInput& in) {
    if (in.tendonLength <= 0.0)
        throw std::runtime_error("预应力筋长度必须 > 0");
    const double alphaEp = (in.Ec > 0.0) ? in.Ep / in.Ec : 0.0;

    PrestressLossResult r;
    r.sigmaL1 = in.anchorageSlip / in.tendonLength * in.Ep;
    r.sigmaL2 = in.sigmaCon * (1.0 - std::exp(-(in.mu * in.theta + in.kappa * in.x)));
    r.sigmaL4 = alphaEp * in.sigmaPc;
    r.sigmaL5 = in.relaxation * in.sigmaCon;
    r.sigmaL6 = alphaEp * in.sigmaPc * in.creepPhi + in.Ep * in.epsCs;
    r.sigmaPe = in.sigmaCon - r.sigmaL1 - r.sigmaL2 - r.sigmaL4 - r.sigmaL5 - r.sigmaL6;
    return r;
}

} // namespace fem
