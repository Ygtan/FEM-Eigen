#pragma once

#include <algorithm>
#include <cmath>

namespace fem {

// ==================== JTG D60-2015 规范荷载库 ====================
// 返回 SI 单位（N / N/m / N/m²），与求解器荷载单位一致。

// 汽车冲击系数 μ（JTG D60）。f 为结构竖向基频，单位 Hz。
// f<1.5 Hz 取 0.05；1.5–14 Hz 为 0.1767·ln(f)−0.0157；f>14 Hz 取 0.45。
inline double jtgD60ImpactCoefficient(double fHz) {
    if (!(fHz > 0.0)) return 0.0;
    if (fHz < 1.5) return 0.05;
    if (fHz <= 14.0) return 0.1767 * std::log(fHz) - 0.0157;
    return 0.45;
}

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

// JTG 3362-2018 附录 C.2.1 的名义徐变系数和 βH。
// RH 为相对湿度百分数，h 为理论厚度（mm），fcm 为 28 天圆柱体平均强度（MPa），
// t0 为加载龄期（天）。φ(t,t0)=φ0·[(t−t0)/(βH+t−t0)]^0.3，发展项仍由调用方计算。
// 表 C.2.2 对 C50 及以上的 32.4/fck 折减是查表用的，这里不重复乘：β(fcm) 已随强度下降。
// 公式适用 40%≤RH≤99%。写出这个范围以外的湿度时夹到边界。
struct Jtg3362Creep {
    double phi0 = 0.0;
    double betaH = 700.0;
};

inline Jtg3362Creep jtg3362NominalCreep(double rhPercent, double hMm,
                                        double fcmMPa, double t0Days) {
    const double rh = std::min(99.0, std::max(40.0, rhPercent));
    const double h = std::max(hMm, 1.0);
    const double fcm = std::max(fcmMPa, 1.0);
    const double t0 = std::max(t0Days, 0.1);
    const double phiRH = 1.0 + (1.0 - rh / 100.0) / (0.46 * std::pow(h / 100.0, 1.0 / 3.0));
    const double betaF = 5.3 / std::sqrt(fcm / 10.0);
    const double betaT0 = 1.0 / (0.1 + std::pow(t0, 0.2));
    double betaH = 150.0 * (1.0 + std::pow(1.2 * rh / 100.0, 18.0)) * (h / 100.0) + 250.0;
    if (betaH > 1500.0) betaH = 1500.0;
    return {phiRH * betaF * betaT0, betaH};
}

// JTG 3362-2018 附录 C.1 的名义收缩系数 εcso（正值表示缩短）。
// 表 C.1.2 的 C50 及以上折减是查表用的，这里不重复乘。βsc 取一般硅酸盐水泥 5.0。
inline double jtg3362ShrinkageNominal(double rhPercent, double fcmMPa) {
    const double rh = std::min(99.0, std::max(40.0, rhPercent));
    const double fcm = std::max(fcmMPa, 1.0);
    const double epsS = (160.0 + 50.0 * (9.0 - fcm / 10.0)) * 1e-6;
    const double betaRH = 1.55 * (1.0 - std::pow(rh / 100.0, 3.0));
    return epsS * betaRH;
}

// 从干燥开始龄期 ts 到 t 的收缩发展系数 βs。h 为理论厚度 mm。
inline double jtg3362ShrinkageBeta(double tDays, double tsDays, double hMm) {
    const double dt = std::max(0.0, tDays - tsDays);
    const double h = std::max(hMm, 1.0);
    const double den = 350.0 * std::pow(h / 100.0, 2.0) + dt;
    if (den <= 0.0) return 0.0;
    return std::sqrt(dt / den);
}

// 从 tRef 到 t 新发展的收缩应变，正值缩短。ts 为干燥开始龄期。
// 构件在 tRef 才参与受力时，更早的收缩不放进平衡。
inline double jtg3362ShrinkageStrain(double tDays, double tRefDays, double tsDays,
                                     double rhPercent, double hMm, double fcmMPa) {
    const double e0 = jtg3362ShrinkageNominal(rhPercent, fcmMPa);
    return e0 * (jtg3362ShrinkageBeta(tDays, tsDays, hMm)
                 - jtg3362ShrinkageBeta(tRefDays, tsDays, hMm));
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
