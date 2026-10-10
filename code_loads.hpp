#pragma once

#include <algorithm>

namespace fem {

// ==================== JTG D60-2015 规范荷载库 ====================
// 返回 SI 单位（N / N/m / N/m²），与求解器荷载单位一致。

// 车道荷载均布力 qk = 10.5 kN/m（返回 N/m）。
inline double jtgD60LaneLoadQk() { return 10500.0; }

// 车道荷载集中力 Pk（返回 N）：
// L≤5m：270 kN；L≥50m：360 kN；中间线性内插 Pk = 2(L+130) kN。
inline double jtgD60LaneLoadPk(double spanMeters) {
    const double L = std::max(0.0, spanMeters);
    if (L <= 5.0) return 270000.0;
    if (L >= 50.0) return 360000.0;
    return 270000.0 + (360000.0 - 270000.0) / (50.0 - 5.0) * (L - 5.0);
}

// 人群荷载 3.0 kN/m²（返回 N/m²）。
inline double jtgD60PedestrianLoad() { return 3000.0; }

// 多车道横向折减系数。
inline double jtgD60LaneReduction(int nLanes) {
    switch (nLanes) {
        case 1:
        case 2:  return 1.0;
        case 3:  return 0.78;
        case 4:  return 0.67;
        case 5:  return 0.60;
        case 6:  return 0.55;
        case 7:  return 0.52;
        default: return 0.50;   // ≥8 车道
    }
}

} // namespace fem
