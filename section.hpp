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
    int mat = -1;      // 材料索引（-1 = 用纤维自身 E 的弹性）
};

struct SectionProperties {
    double A = 0.0;              // 面积
    double cy = 0.0, cz = 0.0;   // 形心（相对参考点）
    double Iy = 0.0, Iz = 0.0;   // 惯性矩（相对形心；Iy=∫z²dA, Iz=∫y²dA）
    double Iyz = 0.0;            // 惯性积（相对形心）
    double J = 0.0;              // 扭转常数（自由扭转；纤维积分不自动求，需显式给出）
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

// ==================== 组合截面（等效模量法） ====================

struct Section {
    std::vector<Fiber> fibers;
    double Eref = 0.0;   // 基准弹性模量（换算截面用，通常取混凝土 E）
};

// 添加一个矩形子区域（y0,z0 为左下角，高 h、宽 b）到纤维集合。
inline void addRectFibers(std::vector<Fiber>& out,
                          double y0, double z0, double h, double b,
                          int ny, int nz, double E) {
    const double dz = b / nz, dy = h / ny;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nz; ++i) {
            Fiber f;
            f.y = y0 + (j + 0.5) * dy;
            f.z = z0 + (i + 0.5) * dz;
            f.area = dz * dy;
            f.E = E;
            f.hy = dy / 2.0;
            f.hz = dz / 2.0;
            out.push_back(f);
        }
}

// 换算截面属性：纤维面积按 E/Eref 放大后积分（等效模量法）。
inline SectionProperties computeTransformedProperties(const Section& sec) {
    std::vector<Fiber> scaled = sec.fibers;
    for (auto& f : scaled)
        if (sec.Eref > 0.0) f.area *= f.E / sec.Eref;
    return computeSectionProperties(scaled);
}

// 线性截面应力（平截面假设）：σ(y,z) = N/A - Mz·y/Iz + My·z/Iy（换算截面），
// 再按各纤维 E/Eref 恢复真实应力。返回每个纤维的真实应力。
inline std::vector<double> sectionFiberStresses(const Section& sec,
                                                double N, double My, double Mz) {
    const SectionProperties p = computeTransformedProperties(sec);
    std::vector<double> sig;
    sig.reserve(sec.fibers.size());
    const double c = (sec.Eref > 0.0) ? 1.0 : 0.0;
    for (const auto& f : sec.fibers) {
        const double dy = f.y - p.cy;
        const double dz = f.z - p.cz;
        const double s = N / p.A - Mz * dy / p.Iz + My * dz / p.Iy;
        sig.push_back(c ? s * (f.E / sec.Eref) : s);
    }
    return sig;
}

// ==================== 变截面（沿梁长插值截面属性） ====================

struct VariableSection {
    std::vector<std::pair<double, SectionProperties>> control;  // (x, 属性)，按 x 升序
};

// 线性插值截面属性（在控制点间）。
inline SectionProperties sectionAt(const VariableSection& vs, double x) {
    if (vs.control.empty()) return SectionProperties{};
    if (x <= vs.control.front().first) return vs.control.front().second;
    if (x >= vs.control.back().first) return vs.control.back().second;
    for (std::size_t i = 1; i < vs.control.size(); ++i) {
        if (x <= vs.control[i].first) {
            const double x0 = vs.control[i - 1].first;
            const double x1 = vs.control[i].first;
            const double t = (x - x0) / (x1 - x0);
            const SectionProperties& a = vs.control[i - 1].second;
            const SectionProperties& b = vs.control[i].second;
            SectionProperties p;
            p.A = a.A + (b.A - a.A) * t;
            p.cy = a.cy + (b.cy - a.cy) * t;
            p.cz = a.cz + (b.cz - a.cz) * t;
            p.Iy = a.Iy + (b.Iy - a.Iy) * t;
            p.Iz = a.Iz + (b.Iz - a.Iz) * t;
            p.Iyz = a.Iyz + (b.Iyz - a.Iyz) * t;
            p.J = a.J + (b.J - a.J) * t;
            p.yTop = a.yTop + (b.yTop - a.yTop) * t;
            p.yBot = a.yBot + (b.yBot - a.yBot) * t;
            p.zLeft = a.zLeft + (b.zLeft - a.zLeft) * t;
            p.zRight = a.zRight + (b.zRight - a.zRight) * t;
            return p;
        }
    }
    return vs.control.back().second;
}

// 单元级平均截面（变截面近似）：取单元中点处的截面。
inline SectionProperties elementMidSection(const VariableSection& vs,
                                           double x1, double x2) {
    return sectionAt(vs, 0.5 * (x1 + x2));
}

// ==================== 材料非线性（纤维本构 + 截面切线） ====================

struct MaterialModel {
    double E = 0.0;          // 弹性模量
    double fy = 1e300;       // 屈服强度（钢，拉压对称）
    double ft = 0.0;         // 抗拉强度（混凝土开裂阈值）
    std::string kind = "ELASTIC";  // ELASTIC / STEEL / CONCRETE
};

// 材料应力-应变：返回 (应力, 切线模量)。
inline std::pair<double, double> materialStress(const MaterialModel& m, double eps) {
    if (m.kind == "STEEL") {
        const double s = m.E * eps;
        if (std::abs(s) <= m.fy) return {s, m.E};
        return {std::copysign(m.fy, eps), 0.0};   // 理想弹塑性
    }
    if (m.kind == "CONCRETE") {
        if (eps >= 0.0 && m.ft > 0.0) {           // 受拉开裂退出
            const double s = m.E * eps;
            if (s <= m.ft) return {s, m.E};
            return {0.0, 0.0};
        }
        return {m.E * eps, m.E};                  // 受压线弹性（简化）
    }
    return {m.E * eps, m.E};
}

struct SectionResponse {
    double N = 0.0, My = 0.0, Mz = 0.0;
    // 切线矩阵（3×3 对称）：[ε0, κz, κy] → [N, Mz, My]
    double D11 = 0.0, D12 = 0.0, D13 = 0.0;
    double D22 = 0.0, D23 = 0.0, D33 = 0.0;
};

// 截面响应：给定截面变形 (ε0, κz, κy)，平截面假设 ε=ε0-κz·y+κy·z，
// 逐纤维按本构取 σ 与切线，积分合力与切线矩阵。
inline SectionResponse sectionResponse(const Section& sec,
                                       const std::vector<MaterialModel>& mats,
                                       double eps0, double kz, double ky) {
    const SectionProperties p = computeSectionProperties(sec.fibers);
    SectionResponse r;
    for (const auto& f : sec.fibers) {
        const MaterialModel* mm = (f.mat >= 0 && f.mat < static_cast<int>(mats.size()))
                                      ? &mats[f.mat] : nullptr;
        const double E = mm ? mm->E : f.E;
        const MaterialModel elastic{E, 1e300, 0.0, "ELASTIC"};
        const MaterialModel& m = mm ? *mm : elastic;

        const double dy = f.y - p.cy;
        const double dz = f.z - p.cz;
        const double eps = eps0 - kz * dy + ky * dz;
        const auto [sig, Et] = materialStress(m, eps);
        r.N += sig * f.area;
        r.Mz += -sig * dy * f.area;
        r.My += sig * dz * f.area;
        r.D11 += Et * f.area;
        r.D12 += -Et * dy * f.area;
        r.D13 += Et * dz * f.area;
        r.D22 += Et * dy * dy * f.area;
        r.D23 += -Et * dy * dz * f.area;
        r.D33 += Et * dz * dz * f.area;
    }
    return r;
}

} // namespace fem
