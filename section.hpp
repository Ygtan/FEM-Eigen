#pragma once
// 纤维截面模型：把截面离散为若干纤维（y,z,面积,材料），
// 截面属性与温度/徐变/非线性响应都退化为"对纤维循环积分"。
// 这是变截面、组合梁、温度按纤维积分、材料非线性的共同地基。

#include <algorithm>
#include <cmath>
#include <vector>

namespace fem {

struct Fiber {
    double y = 0.0;    // 竖向坐标（相对参考点，向上为正）
    double z = 0.0;    // 横向坐标（相对参考点）
    double area = 0.0;
    double E = 0.0;    // 弹性模量（组合截面时各纤维可不同）
    double hy = 0.0;   // 纤维半高（用于精确的截面上下边缘）
    double hz = 0.0;   // 纤维半宽（用于精确的截面左右边缘）
};

struct SectionProperties {
    double A = 0.0;              // 面积
    double cy = 0.0, cz = 0.0;   // 形心（相对参考点）
    double Iy = 0.0, Iz = 0.0;   // 惯性矩（相对形心；Iy=∫z²dA, Iz=∫y²dA）
    double Iyz = 0.0;            // 惯性积（相对形心）
    double yTop = 0.0, yBot = 0.0;   // 上下缘坐标（相对形心）
    double zLeft = 0.0, zRight = 0.0; // 左右缘坐标（相对形心）

    // 截面模量
    double Wtop() const { return std::abs(yTop) > 0 ? Iz / std::abs(yTop) : 0.0; }
    double Wbot() const { return std::abs(yBot) > 0 ? Iz / std::abs(yBot) : 0.0; }
    double WyLeft() const { return std::abs(zLeft) > 0 ? Iy / std::abs(zLeft) : 0.0; }
    double WyRight() const { return std::abs(zRight) > 0 ? Iy / std::abs(zRight) : 0.0; }
};

// 由纤维积分截面属性。
inline SectionProperties computeSectionProperties(const std::vector<Fiber>& fibers) {
    SectionProperties p;
    for (const auto& f : fibers) {
        p.A += f.area;
        p.cy += f.y * f.area;
        p.cz += f.z * f.area;
    }
    if (p.A <= 0.0) return p;
    p.cy /= p.A;
    p.cz /= p.A;

    p.yTop = -1e300; p.yBot = 1e300;
    p.zLeft = 1e300; p.zRight = -1e300;
    for (const auto& f : fibers) {
        const double dy = f.y - p.cy;
        const double dz = f.z - p.cz;
        p.Iy += dz * dz * f.area;
        p.Iz += dy * dy * f.area;
        p.Iyz += dy * dz * f.area;
        p.yTop = std::max(p.yTop, dy + f.hy);
        p.yBot = std::min(p.yBot, dy - f.hy);
        p.zLeft = std::min(p.zLeft, dz - f.hz);
        p.zRight = std::max(p.zRight, dz + f.hz);
    }
    return p;
}

// 生成矩形截面纤维网格：宽 b（z 向）、高 h（y 向），形心在原点。
inline std::vector<Fiber> makeRectangularFibers(double b, double h,
                                                int nz, int ny, double E) {
    std::vector<Fiber> out;
    out.reserve(static_cast<std::size_t>(nz) * ny);
    const double dz = b / nz, dy = h / ny;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nz; ++i) {
            Fiber f;
            f.z = -b / 2.0 + (i + 0.5) * dz;
            f.y = -h / 2.0 + (j + 0.5) * dy;
            f.area = dz * dy;
            f.E = E;
            f.hy = dy / 2.0;
            f.hz = dz / 2.0;
            out.push_back(f);
        }
    return out;
}

// 生成 I 型截面纤维（翼缘 bf×tf，腹板 tw×hw，形心在原点）。
inline std::vector<Fiber> makeIFibers(double bf, double tf, double tw, double hw,
                                      int n, double E) {
    std::vector<Fiber> out;
    auto addRect = [&](double y0, double z0, double w, double hgt, int nz, int ny) {
        const double dz = w / nz, dy = hgt / ny;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nz; ++i) {
                Fiber f;
                f.z = z0 + (i + 0.5) * dz;
                f.y = y0 + (j + 0.5) * dy;
                f.area = dz * dy;
                f.E = E;
                f.hy = dy / 2.0;
                f.hz = dz / 2.0;
                out.push_back(f);
            }
    };
    addRect(-hw / 2.0 - tf, -bf / 2.0, bf, tf, n, n);             // 下翼缘
    addRect(hw / 2.0, -bf / 2.0, bf, tf, n, n);                   // 上翼缘
    addRect(-hw / 2.0, -tw / 2.0, tw, hw, std::max(1, n / 4), n); // 腹板
    return out;
}

struct SectionResultants {
    double N = 0.0;   // 轴力
    double My = 0.0;  // 绕 y 轴弯矩
    double Mz = 0.0;  // 绕 z 轴弯矩
};

// 线性温度梯度（上下缘温差 dT，跨截面高 h）产生的截面合力：
//   N = ∫E·α·T dA，Mz = ∫E·α·T·y dA，My = ∫E·α·T·z dA（相对形心）。
inline SectionResultants sectionTemperatureGradient(const std::vector<Fiber>& fibers,
                                                    double alpha, double dT, double h) {
    const SectionProperties p = computeSectionProperties(fibers);
    SectionResultants r;
    for (const auto& f : fibers) {
        const double T = dT * f.y / h;   // 线性梯度，T(y)=dT·y/h
        r.N += f.E * alpha * T * f.area;
        r.Mz += f.E * alpha * T * (f.y - p.cy) * f.area;
        r.My += f.E * alpha * T * (f.z - p.cz) * f.area;
    }
    return r;
}

} // namespace fem
