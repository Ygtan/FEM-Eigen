#pragma once

#include <cctype>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

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

    // σl2 曲线筋反摩擦（JTG 3362-2018 附录 G，opt-in，默认全关）。
    //   l2Anti=0：沿用原式 σl2 = ΣΔl/l·Ep（22 键既有行为逐位不变）。
    //   l2Anti=1：按附录 G 考虑反向摩擦，代表截面由 l2ReprX 指定。
    // 代表值口径（用户 2026-11-16 拍板）：l2ReprX=0 取张拉端 x=0 峰值损失，
    //   与 G.0.2-4 `Δσ=2Δσd·lf` 闭式一致，且 σpe 用于索力 / 承载力时偏保守。
    //   规范对代表截面无统一规定；非 0 值一律显式拒绝，不静默套用。
    int l2Anti = 0;
    int l2ReprX = 0;

    // σl6 收缩徐变龄期口径（JTG 3362-2018 附录 C，opt-in，默认全关）。
    // 用户 2026-11-16 拍板：表值由用户在模型里显式录入经核对的规范表值，
    //   程序只按解析发展公式组装 εcs(t,t0)/φ(t,t0)，绝不调用
    //   code_loads.hpp::jtg3362NominalCreep()/jtg3362ShrinkageNominal() 拟合值兜底。
    // c2d=0（默认）：沿用原路径 σl6 = 0.9(Ep·εcs + αEP·σpc6·φ)/(1+15ρρps)，
    //   εcs/φ 即上方 ECS/PHI 用户直接给定的已发展值，22 键行为逐位不变。
    // c2d=1：按附录 C 显式表值龄期口径重算 εcs/φ（见下），缺任一表值 → throw。
    // cSkip=1：显式声明 σl6 不参与（σl6=0，SOURCE=SKIP），与「缺值报错」严格区分；
    //   仅在 c2d=0 时有效，c2d=1 且 cSkip=1 视为语义冲突显式拒绝。
    int c2d = 0;
    int cSkip = 0;
    // —— 附录 C 表值（c2d=1 时全部必填，缺一 throw）——
    // eps0：名义收缩系数 εcs0（表 C.1.2，40~70% 取 0.529e-3、70~99% 取 0.310e-3；
    //   C50 及以上应乘以 32.4/fck 折减后的值），用户填最终核对值。
    // phi0：名义徐变系数 φ0（表 C.2.2 查表值，C50 及以上已乘 32.4/fck），用户填最终核对值。
    double eps0 = 0.0;
    double phi0 = 0.0;
    double cRh = 0.0;     // 环境年平均相对湿度 RH（%），用于 βRH 与 βH 代表值判区间
    double cH = 0.0;      // 构件理论厚度 h（mm），h=2A/u
    double cT0 = 0.0;     // 加载 / 传力锚固龄期 t0（d）
    double cTs = 0.0;     // 收缩开始龄期 ts（d），设计取 3~7d
    double cT = 0.0;      // 计算考虑龄期 t（d）
};

struct PrestressLossResult {
    double sigmaL1 = 0.0;   // 摩擦
    double sigmaL2 = 0.0;   // 锚具变形/回缩/接缝压缩
    double sigmaL3 = 0.0;   // 温差
    double sigmaL4 = 0.0;   // 弹性压缩
    double sigmaL5 = 0.0;   // 松弛
    double sigmaL6 = 0.0;   // 收缩徐变
    double sigmaPe = 0.0;   // 传力锚固应力（后张法：σcon - σl1 - σl2 - σl4）
    // 附录 G 反摩擦中间量（仅 l2Anti=1 填写，默认 0，旧路径不填）：
    double lf = 0.0;          // 反向摩擦影响长度 lf (m)
    double deltaSigmaD = 0.0; // Δσd = σl1(l)/l (Pa/m)，锚固端摩擦梯度
    double deltaSigmaPeak = 0.0; // 代表截面 x=0 的 σl2 峰值 Δσ (Pa)
    int l2Branch = 0;         // 0=未启用/无摩擦退化, 1=lf≤l, 2=lf>l
    // 附录 C 收缩徐变中间量（仅 c2d=1 填写，默认 0，旧路径不填）：
    double epsCsDev = 0.0;    // εcs(t,t0) = εcs0[βs(t,ts)-βs(t0,ts)]，已发展值
    double phiDev = 0.0;      // φ(t,t0) = φ0·βc(t,t0)，已发展值
    double betaS = 0.0;       // βs(t,ts)
    double betaS0 = 0.0;      // βs(t0,ts)
    double betaC = 0.0;       // βc(t,t0)
    double betaH = 0.0;       // βH（βc 用）
    double cHRepresentative = 0.0; // βH 计算用的 RH 代表值（55% 或 80%）
    int cSource = 0;          // 0=用户直给 ECS/PHI（旧路径）, 1=附录 C 表值龄期, 2=显式 SKIP
};

inline PrestressLossResult computePrestressLosses(const PrestressLossInput& in) {
    if (in.tendonLength <= 0.0)
        throw std::runtime_error("预应力筋长度 l 必须 > 0");
    if (!in.threadedBar && in.fpk <= 0.0)
        throw std::runtime_error("σl5 需要预应力筋抗拉强度标准值 fpk > 0");
    const double alphaEp = (in.Ec > 0.0) ? in.Ep / in.Ec : 0.0;

    PrestressLossResult r;
    r.sigmaL1 = in.sigmaCon * (1.0 - std::exp(-(in.mu * in.theta + in.kappa * in.x)));
    // σl2：默认直线式；l2Anti=1 走附录 G 反摩擦（opt-in）。
    if (in.l2Anti == 0) {
        r.sigmaL2 = in.sumDeltaL / in.tendonLength * in.Ep;
    } else {
        if (in.l2ReprX != 0)
            throw std::runtime_error("附录 G 反摩擦当前仅支持代表截面 x=0（张拉端），l2ReprX 必须为 0");
        if (in.sumDeltaL <= 0.0)
            throw std::runtime_error("附录 G 反摩擦需要 ΣΔl（DELTA_L）> 0");
        if (in.Ep <= 0.0)
            throw std::runtime_error("附录 G 反摩擦需要 Ep > 0");
        // 锚固端摩擦损失 σl1(l)：用全长 l，不是计算截面 x。
        const double sigmaL1Anchor = in.sigmaCon *
            (1.0 - std::exp(-(in.mu * in.theta + in.kappa * in.tendonLength)));
        r.deltaSigmaD = sigmaL1Anchor / in.tendonLength;   // Pa/m
        if (!(r.deltaSigmaD > 0.0)) {
            // 无摩擦时不存在反向摩擦效应，回缩损失沿程均匀 —— 退化为原式，
            // 这是物理边界而非「缺参兜底」（附录 G 的 Δσd 分母为 0 无定义）。
            r.sigmaL2 = in.sumDeltaL / in.tendonLength * in.Ep;
            r.l2Branch = 0;
        } else {
            // lf = sqrt(ΣΔl·Ep / Δσd)，SI 单位 (m)；规范原文用 mm/MPa，
            // 本内核 m/Pa 同量纲，无需额外换算（sqrt(m·Pa/(Pa/m))=m）。
            r.lf = std::sqrt(in.sumDeltaL * in.Ep / r.deltaSigmaD);
            if (r.lf <= in.tendonLength) {
                // G.0.2-4 / G.0.2-3：Δσ=2Δσd·lf，代表截面 x=0 即峰值。
                r.deltaSigmaPeak = 2.0 * r.deltaSigmaD * r.lf;
                r.sigmaL2 = r.deltaSigmaPeak;
                r.l2Branch = 1;
            } else {
                // G.0.2-5 等腰梯形：σl2(x)=Δσ′−2x·Δσd，面积协调
                //   ΣΔl·Ep = ∫₀ˡσl2 dx = Δσd·lf²  ⇒  Δσ′=Δσd·(lf²/l + l)
                // 该式在 lf=l 时连续于上一支（Δσ′=2Δσd·l）。
                r.deltaSigmaPeak = r.deltaSigmaD *
                    (r.lf * r.lf / in.tendonLength + in.tendonLength);
                r.sigmaL2 = r.deltaSigmaPeak;
                r.l2Branch = 2;
            }
        }
    }
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
    // σl6 收缩徐变（6.2.7-1）：默认用户直给；opt-in 显式 SKIP 或附录 C 表值龄期。
    // cSkip 与 c2d 语义互斥，避免「显式声明不参与」与「启用表值计算」混淆。
    if (in.cSkip != 0 && in.c2d != 0)
        throw std::runtime_error("σl6：C_SKIP=1（显式不参与）与 C2D=1（表值计算）不能同时启用");
    if (in.cSkip != 0) {
        r.sigmaL6 = 0.0;
        r.cSource = 2;
        return r;                       // 显式声明不参与，不进 σpe，也不报错
    }
    if (in.c2d == 0) {
        // 旧路径：εcs/φ 即 ECS/PHI 用户直接给定的已发展值，逐位不变。
        r.sigmaL6 = 0.9 * (in.Ep * in.epsCs + alphaEp * in.sigmaPc6 * in.phi)
                  / (1.0 + 15.0 * in.rho * in.rhoPs);
        r.cSource = 0;
        return r;
    }
    // c2d=1：JTG 3362-2018 附录 C 显式表值龄期口径。
    // 全部表值必须显式给出，缺一 throw —— 绝不调用 code_loads 拟合值兜底。
    if (in.eps0 <= 0.0)
        throw std::runtime_error("σl6 附录 C：缺 EPS0（表 C.1.2 名义收缩系数 εcs0，C50 以上需已乘 32.4/fck）");
    if (in.phi0 <= 0.0)
        throw std::runtime_error("σl6 附录 C：缺 PHI0（表 C.2.2 名义徐变系数 φ0，C50 以上需已乘 32.4/fck）");
    if (!(in.cRh >= 40.0 && in.cRh < 99.0))
        throw std::runtime_error("σl6 附录 C：RH 必须在 40%≤RH<99%（C.1.1-4 适用范围）");
    if (!(in.cH > 0.0))
        throw std::runtime_error("σl6 附录 C：缺 H（构件理论厚度 h=2A/u，mm）");
    if (!(in.cT0 > 0.0) || !(in.cTs > 0.0) || !(in.cT > 0.0))
        throw std::runtime_error("σl6 附录 C：缺 T0 / TS / T（龄期，d）");
    if (!(in.cT > in.cT0) || !(in.cT0 >= in.cTs))
        throw std::runtime_error("σl6 附录 C：龄期须满足 t > t0 ≥ ts");
    const double hRatio = in.cH / 100.0;              // h/h0，h0=100mm
    // C.1.1-5：βs(t,ts) = [(t-ts)/t1 / (350(h/h0)² + (t-ts)/t1)]^0.5，t1=1d。
    const double denS = 350.0 * hRatio * hRatio;
    const double dtS = in.cT - in.cTs;
    const double dtS0 = in.cT0 - in.cTs;
    r.betaS = std::sqrt(dtS / (denS + dtS));
    r.betaS0 = std::sqrt(dtS0 / (denS + dtS0));
    // C.1.3：εcs(t,t0) = εcs0[βs(t,ts)-βs(t0,ts)]。
    r.epsCsDev = in.eps0 * (r.betaS - r.betaS0);
    // C.2.3-1：βH 用 RH 区间代表值（40~70% 取 55%，70~99% 取 80%）。
    r.cHRepresentative = (in.cRh < 70.0) ? 55.0 : 80.0;
    // C.2.1-7：βH = 150[1+(1.2RH/100)^18](h/h0)+250 ≤ 1500。
    double betaH = 150.0 * (1.0 + std::pow(1.2 * r.cHRepresentative / 100.0, 18.0))
                   * hRatio + 250.0;
    if (betaH > 1500.0) betaH = 1500.0;
    r.betaH = betaH;
    // C.2.1-6：βc(t,t0) = [(t-t0)/(βH + (t-t0))]^0.3，t1=1d 已约去。
    const double dtC = in.cT - in.cT0;
    r.betaC = std::pow(dtC / (betaH + dtC), 0.3);
    // C.2.1-1：φ(t,t0) = φ0·βc(t,t0)。
    r.phiDev = in.phi0 * r.betaC;
    r.cSource = 1;
    r.sigmaL6 = 0.9 * (in.Ep * r.epsCsDev + alphaEp * in.sigmaPc6 * r.phiDev)
              / (1.0 + 15.0 * in.rho * in.rhoPs);
    return r;
}

// 有效预应力比 σpe/σcon（供 *PRESTRESS 的 LOSSREF 用，等价于 1 - Σσl/σcon）。
inline double effectivePrestressRatio(const PrestressLossInput& in) {
    if (in.sigmaCon <= 0.0)
        throw std::runtime_error("有效预应力比需要张拉控制应力 sigmaCon > 0");
    return computePrestressLosses(in).sigmaPe / in.sigmaCon;
}

// 解析一行 `<编号> KEY=VALUE …`。
//
// 2D 的 `*PRESTRESS_LOSS` 与 3D 的 `*PRESTRESS_LOSS` 共用这一份键名映射，
// 避免两边各写一套后出现「同一组参数在 2D/3D 算出不同有效比」。键名大小写
// 不敏感；POST/REBAR 用 0/1 表示布尔。行内不认识的键忽略（与 2D 原行为一致）。
inline std::pair<int, PrestressLossInput> parsePrestressLossLine(const std::string& line) {
    std::istringstream ss(line);
    int id = 0;
    if (!(ss >> id))
        throw std::runtime_error("PRESTRESS_LOSS 行格式错误（缺编号）: " + line);
    std::map<std::string, std::string> p;
    std::string tok;
    while (ss >> tok) {
        const std::size_t eq = tok.find('=');
        if (eq == std::string::npos) continue;
        std::string key = tok.substr(0, eq);
        for (char& c : key)
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        p[key] = tok.substr(eq + 1);
    }
    PrestressLossInput li;
    auto g = [&](const char* key, double& dst) {
        auto it = p.find(key);
        if (it != p.end()) dst = std::stod(it->second);
    };
    g("SIGMA_CON", li.sigmaCon);
    g("EP", li.Ep);
    g("EC", li.Ec);
    g("FPK", li.fpk);
    g("MU", li.mu);
    g("K", li.kappa);
    g("THETA", li.theta);
    g("X", li.x);
    g("DELTA_L", li.sumDeltaL);
    g("L", li.tendonLength);
    g("T1", li.t1);
    g("T2", li.t2);
    g("SIGMA_PC4", li.sigmaPc4);
    g("SIGMA_PC6", li.sigmaPc6);
    g("PSI", li.psi);
    g("ZETA", li.zeta);
    g("ECS", li.epsCs);
    g("PHI", li.phi);
    g("RHO", li.rho);
    g("RHOPS", li.rhoPs);
    auto itPost = p.find("POST");
    if (itPost != p.end()) li.postTensioned = (std::stod(itPost->second) != 0.0);
    auto itRebar = p.find("REBAR");
    if (itRebar != p.end()) li.threadedBar = (std::stod(itRebar->second) != 0.0);
    // 附录 G 反摩擦 opt-in 键（默认全关；不认识的行仍忽略，与既有行为一致）。
    auto itAnti = p.find("L2ANTI");
    if (itAnti != p.end()) li.l2Anti = (std::stoi(itAnti->second) != 0);
    auto itRepr = p.find("L2REPX");
    if (itRepr != p.end()) li.l2ReprX = std::stoi(itRepr->second);
    // σl6 附录 C 显式表值龄期口径 opt-in 键（默认全关，不认识的键仍忽略）。
    auto itC2D = p.find("C2D");
    if (itC2D != p.end()) li.c2d = (std::stoi(itC2D->second) != 0);
    auto itCSkip = p.find("C_SKIP");
    if (itCSkip != p.end()) li.cSkip = (std::stoi(itCSkip->second) != 0);
    g("EPS0", li.eps0);
    g("PHI0", li.phi0);
    g("C_RH", li.cRh);
    g("C_H", li.cH);
    g("C_T0", li.cT0);
    g("C_TS", li.cTs);
    g("C_T", li.cT);
    return {id, li};
}

} // namespace fem
