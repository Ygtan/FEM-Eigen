// cable_solver.cpp —— 斜拉桥成桥态求解器的命令行接缝
// =====================================================
// 为什么单独做一个 exe，而不是塞进 BridgeStudio 的 Python：
//   成桥态的三种工况（静力 / 定索力找形 / 成桥态闭环）都要用到 P3 内核的
//   共转列式非线性、索垂度、影响矩阵、最小弯曲能量法。Python 侧只该做两件
//   事——**生成 .dat** 和**解析结果**，任何力学都在这里。
//
// 输入：FEM-Eigen 3D new format .dat（readModel3D 可读）＋ 两个扩展段：
//
//     *CABLE_BRIDGE
//     MODE=COMPLETION          ; COMPLETION | TENSION | STATIC
//     NLSTEPS=3                ; 非线性荷载步
//     NEWTON=30                ; 每个荷载步的 Newton 上限
//     MAXITER=20               ; 外层迭代（索长/几何闭合）上限
//     TOL=1e-8                 ; 几何闭合容差 m
//     TOLT=1e-7                ; 索力相对容差
//     METHOD=MINENERGY|RIGID   ; 自动定索力方法（仅 COMPLETION 且未全给索力时用）
//     MINT=0.02                ; 保力下限 = 该比例 × 最大索力
//     SENSITIVITY=LOWRANK      ; 仅 COMPLETION：LOWRANK(默认) | FULL | FIXEDPOINT
//                              ; FULL = P3-3b 完整形状敏感度 K_shape（诊断/对拍用，
//                              ; 默认路径不变）；FIXEDPOINT = 纯不动点对照
//     OPTIMIZE=0               ; 1 = 先做带约束成桥索力优化，把 T* 写回 *CABLE
//     OPTOBJ=MINMAXM           ; MINMAXM(min max|M主梁|, IRLS) | MINENERGY
//     OPTBALANCE=1             ; 逐塔 ΣT·cos = 0（罚函数并入二次型）
//     TMAX=0                   ; 索力上限 N；0 = 程序默认取 1.5×RIGID 最大索力
//                              ; （不是规范限值，段里写 TMAXSRC=DEFAULT_1.5_RIGID）
//     VERBOSE=0
//
//     *CABLE
//     <elemId> <梁侧锚点nodeId> [目标索力N]     ; 每行一根；索力缺省/为 0 = 自动
//
//   若省略 *CABLE 段，则自动取模型里全部 CABLE 单元，锚点取 z 较低的那端。
//
// ---------------------------------------------------------------------------
// 悬索桥（自锚式 / 地锚式）：另加三段定位信息
// ---------------------------------------------------------------------------
//     *SUSP_BRIDGE
//     SYSTEM=EARTH|SELF    ; EARTH=地锚式，SELF=自锚式
//     TOWER1=<nodeId>      ; 主缆链上的塔顶节点
//     TOWER2=<nodeId>
//     SAG=<m>              ; 主跨矢高（也可给 SAGRATIO=矢跨比）
//     DECKE0=<elemId>      ; 加劲梁单元范围（用于输出梁内力）
//     DECKE1=<elemId>
//
//     *MAINCABLE           ; 主缆单元，严格按"锚→塔1→塔2→锚"的链序排列
//     <elemId> ...
//
//     *HANGER              ; 吊索单元，每根必须一端接主缆、一端接加劲梁
//     <elemId> ...
//
//   与斜拉桥的实质差别只有一条：**主缆线形由恒载决定**，不能人工给定。
//   斜拉桥是"梁面高程由人定、索力去适应"；悬索桥是"索力由恒载定、线形随之"。
//   所以这里先做索多边形找形（水平力 H 对矢高二分，主跨与边跨共用同一个 H，
//   否则塔顶受不平衡水平力），把设计线形改成找形结果，再闭式反算无应力索长
//   L0 解一次。详见 formFindSuspension() 与 main() 里 susp 分支的注释。
//
//   自锚式与地锚式的差别全部在边界条件里，内核不做特殊处理：
//     地锚式——主缆端锚在锚碇（全约束），H 传给地基；梁端纵向固定。
//     自锚式——主缆端锚在加劲梁端节点上，H 压进加劲梁（梁内轴力 ≈ −H）；
//              梁端纵向必须自由，纵向约束放跨中（对称荷载下反力恒为零）。
//
//   两条硬校验（都是"不报错的错"里最贵的两种）：
//     1. TOWER1/TOWER2 上必须有梁单元。主缆若在塔顶旁边另建一个坐标相同的
//        节点，就成了从锚碇到锚碇的自由索，塔不参与受力——求解照样收敛、
//        位移倾角看着都正常，只有索力掉到百分之一。
//     2. 主缆端节点必须"要么全约束（锚碇）、要么接加劲梁（自锚式）"。
//
// 输出：与 .dat 同目录 <stem>_cable.txt，分段纯文本（Python 侧顺序解析）：
//     mode=... nodes=... elements=... cables=... dofs=...
//     converged=... iterations=... elapsed_ms=...
//     closure_err=... tension_err=...
//         （悬索桥：closure_err = 设计构形与成桥构形的最大偏差 max|u|，
//           即需要预变形的量；tension_err = 达到索力相对目标索力的最大偏差）
//     tension_min/max/sum_kN=...
//     vertical_sum_kN=... total_load_kN=... balance_pct=...
//         （斜拉桥：索力竖向分量 vs 恒载；悬索桥：支座反力 vs 竖向荷载）
//     support_v_sum_kN / support_v_max_kN / applied_v_kN
//     susp_system=... H_kN=... sag_m=... sag_target_m=... span_m=... sag_ratio=...
//     hanger_count / hanger_T_min|max / hanger_V_sum / hanger_T_sum
//     girder_N_max_kN / girder_M_max_kNm
//     [CABLE]    elem anchor T L0 chord Ltarget dL dz
//     [MAINCABLE] elem T L0 z chordEnd
//     [HANGER]   elem deckNode T L0 chord
//     [GIRDER]   elem N My Vz
//         （仅悬索桥；每根加劲梁一行，N/My/Vz **各自**取两端中绝对值较大者。
//           它是给曲线与摘要看的，跨端混搭，**不能**拿来做应力验算——
//           max|N| 在 i 端、max|My| 在 j 端时，这一行描述的是不存在的截面。）
//     [GIRDEREND] elem end N My Vz
//         （两种桥型都有；**每一端**一行，end=1 为 i 端、2 为 j 端。
//           应力验算（σ = N/A ± My/W）必须用这一段：同一端才是同一截面。
//           不在 DECKE 范围内的梁也会出现——主梁与塔柱的区分由前端按
//           (材料, 实常数) 分组决定，内核不替使用者判断"哪根是主梁"。）
//     [REACTION] node dof R             支座反力（向上为正；锚碇为负=抗拔）
//     [CAMBER] node dx dy dz            制造线形相对设计线形的偏移（m）
//     [ERECTED] node x y z              制造线形绝对坐标
//     [U]      node dof value           成桥位移（相对制造线形）
//     [MODAL]  ...                      成桥态切线刚度上的模态（modal_output_3d.hpp）
//     [BUCKLING] ...                   本内核只写 ERROR=nonlinear_state_not_supported
//     [LENGTH] KIND/STATUS/…/I/E/J     索长迭代 / 几何闭合逐轮过程
//               （cable_length_output_3d.hpp；无外层迭代时写 ERROR 段）
//     [END]
//
// 编译：随 CMake 一起（目标 cable_solver）。

#include "fem_3d_nl.hpp"
#include "fem_2d.hpp"   // impactCoefficientFromFrequency（JTG D60 冲击系数）
#include "modal_output_3d.hpp"
#include "buckling_output_3d.hpp"
#include "cable_length_output_3d.hpp"
#include "newton_output_3d.hpp"
#include "cable_opt_3d.hpp"      // [OPT]：带约束成桥索力优化（默认关闭）

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace fem;

namespace {

double nowMs() {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// 在成桥态一致切线刚度 K（已施加 Dirichlet，约束自由度由 isC 标记）上求前
// nev 阶自振频率与模态。K、M 同规模（nDofs × nDofs）。索承体系的基频必须
// 用成桥态切线刚度——张紧索的切线刚度远大于弹性刚度，线弹性装配刚度会把
// 索算软、基频算低。返回频率（Hz，升序）与对应全局位移模态。
static ModalResult3D modalFromTangent(const SparseMat3D& K, const SparseMat3D& M,
                                      const std::vector<char>& isC,
                                      int nDofs, int nev) {
    std::vector<int> free;
    for (int d = 0; d < nDofs; ++d) {
        if (isC[d]) continue;
        if (K.coeff(d, d) == 0.0) continue;
        if (M.coeff(d, d) == 0.0) continue;
        free.push_back(d);
    }
    ModalResult3D out;
    const int nf = static_cast<int>(free.size());
    if (nf == 0) return out;
    const int n = std::min(nev, nf);
    constexpr double pi = 3.14159265358979323846;
#ifdef FEM_HAS_SPECTRA
    const SparseModalResult sm = solveSparseModal(K, M, free, nDofs, n);
    if (!sm.eigenvalues.empty()) {
        for (std::size_t i = 0; i < sm.eigenvalues.size(); ++i) {
            const double lam = sm.eigenvalues[i];
            const double f = (lam > 0.0) ? std::sqrt(lam) / (2.0 * pi) : 0.0;
            out.frequencies.push_back(f);
            out.modes.push_back(sm.modes[i]);
        }
        return out;
    }
#endif
    // 稠密回退（Spectra 未收敛或自由度太少）
    Eigen::MatrixXd Kd(nf, nf), Md(nf, nf);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            Kd(a, b) = K.coeff(free[a], free[b]);
            Md(a, b) = M.coeff(free[a], free[b]);
        }
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Kd, Md);
    for (int i = 0; i < n; ++i) {
        const double lam = ges.eigenvalues()(i);
        const double f = (lam > 0.0) ? std::sqrt(lam) / (2.0 * pi) : 0.0;
        out.frequencies.push_back(f);
        std::vector<double> mode(nDofs, 0.0);
        for (int a = 0; a < nf; ++a) mode[free[a]] = ges.eigenvectors()(a, i);
        out.modes.push_back(std::move(mode));
    }
    return out;
}

// 识别竖向一阶弯曲基频：在（频率升序的）模态里，取主梁节点竖向位移 w 主导
//（w² 占主梁平动 u²+v²+w² 比例 > 0.5）的最低阶。JTG D60 汽车冲击系数针对
// 竖向汽车荷载，须用竖向弯曲基频而非任意低阶（横向/纵向）频率。识别不出
// 返回 -1，调用方回退第一阶自振频率。
static double verticalFlexuralFrequency(const Model3D& m, const DofMap3D& dm,
                                        const ModalResult3D& modal) {
    std::set<int> deckNodes;
    for (const auto& e : m.elements) {
        if (e.type != "BEAM" && e.type != "BEAM7") continue;
        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
        if (std::abs(d.x()) < std::max(std::abs(d.y()), std::abs(d.z()))) continue;
        deckNodes.insert(e.n1);
        deckNodes.insert(e.n2);
    }
    if (deckNodes.empty()) return -1.0;

    for (std::size_t i = 0; i < modal.frequencies.size(); ++i) {
        const double f = modal.frequencies[i];
        if (f <= 0.0) continue;
        const auto& mode = modal.modes[i];
        double w2 = 0.0, uv2 = 0.0;
        for (int nid : deckNodes) {
            const int dw = dm.safeDof(nid, 2);   // w（竖向平动，0-based）
            const int du = dm.safeDof(nid, 0);   // u
            const int dv = dm.safeDof(nid, 1);   // v
            const double w = (dw >= 0 && dw < static_cast<int>(mode.size())) ? mode[dw] : 0.0;
            const double u = (du >= 0 && du < static_cast<int>(mode.size())) ? mode[du] : 0.0;
            const double v = (dv >= 0 && dv < static_cast<int>(mode.size())) ? mode[dv] : 0.0;
            w2 += w * w;
            uv2 += u * u + v * v;
        }
        const double tot = w2 + uv2;
        if (tot <= 0.0) continue;
        if (w2 / tot > 0.5) return f;   // 竖向主导的最低阶
    }
    return -1.0;
}

struct CableReq {
    int elem = 0;
    int anchor = 0;        // 梁侧锚点；0 = 由几何自动判定
    double tension = 0.0;  // <=0 = 自动定索力
};

struct Cfg {
    std::string mode = "COMPLETION";
    int nlSteps = 3;
    int newton = 30;
    int maxIter = 20;
    int maxLenIter = 15;
    double tol = 1e-8;
    double tolTension = 1e-7;
    bool minEnergy = true;          // false = 刚性支承连续梁法
    // 仅 MODE=COMPLETION：几何闭合敏感度口径（默认 LOWRANK，既有算例不变）。
    std::string sensitivity = "LOWRANK";   // LOWRANK | FULL | FIXEDPOINT
    double minTension = 0.0;
    double minTensionFraction = 0.02;
    double deckLoad = 0.0;          // 桥面恒载合计 N（仅用于校核报告，0=未知）
    // 带约束成桥索力优化。默认全关，既有算例的索力口径逐位不变。
    bool optimize = false;
    std::string optObjective = "MINMAXM";   // MINMAXM | MINENERGY
    bool optBalance = true;
    double optTMax = 0.0;            // 0 = 程序默认上界 1.5×RIGID
    // 成桥体系守卫：优化目标只有弯矩、约束只有索力箱，不拦的话会把边跨索
    // 压到下限，索承担恒载跌破 90%、相邻索力比暴涨。写回 T* 前用与计算书
    // 工程校核表同一口径预估三条判据，不合格就 APPLIED=0、不写回。
    bool sysGuard = true;            // 默认开启；SYSGUARD=0 仅用于诊断
    double sysDeckMin = 90.0;        // 索承担桥面恒载下限 %（严格 >）
    double sysDeckMax = 100.0;       // 上限 %（超过即支座上拔）
    double sysAdjMax = 1.5;          // 同跨相邻索力比上限
    bool verbose = false;
    std::vector<CableReq> cables;
    // ---- 悬索桥（*SUSP_BRIDGE / *MAINCABLE / *HANGER）----
    std::string system;             // "" = 非悬索桥；"EARTH" 地锚式；"SELF" 自锚式
    double sag = 0.0;               // 主跨矢高 m（>0 时优先于矢跨比）
    double sagRatio = 0.0;          // 矢跨比（0.1 = 1/10）
    int tower1 = 0, tower2 = 0;     // 两个塔顶节点
    std::vector<int> mainCable;     // 主缆单元（按链序：锚→塔1→塔2→锚）
    std::vector<int> hangers;       // 吊索单元
    double deckE0 = 0.0, deckE1 = 1e300;  // 加劲梁单元号区间（输出轴力用）
    // ---- 悬索桥施工阶段（*ERECTION）----
    bool erection = false;          // 是否做施工阶段模拟
    double erectionQ1 = 0.0;        // 一期恒载 N/m（钢箱梁自重，吊装阶段施加）
    double erectionQ2 = 0.0;        // 二期恒载 N/m（铺装+护栏，合拢后施加）
    int erectionPanels = 1;         // 每次对称吊装的吊点对数（每次对称吊 2 段）
};

// 只扫 *CABLE_BRIDGE / *CABLE 两个扩展段；其余交给 readModel3D。
// 单独扫一遍比往 readModel3D 里塞桥型专用逻辑干净——内核的 .dat 格式
// 应该是"几何 + 载荷"，而"要哪根索拉多少力"是求解工况，不是模型。
Cfg readCfg(const std::string& path) {
    Cfg c;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开输入文件: " + path);
    std::string line, sec;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        if (line[0] == '*') {
            std::string kw = line.substr(1);
            const std::size_t sp = kw.find_first_of(" \t,");
            if (sp != std::string::npos) kw = kw.substr(0, sp);
            sec = upper(kw);
            continue;
        }
        if (sec == "CABLE_BRIDGE") {
            std::istringstream ss(line);
            std::string kv;
            while (ss >> kv) {
                const std::size_t p = kv.find('=');
                if (p == std::string::npos) continue;
                const std::string key = upper(kv.substr(0, p));
                const std::string val = kv.substr(p + 1);
                if (key == "MODE") c.mode = upper(val);
                else if (key == "NLSTEPS") c.nlSteps = std::max(1, std::stoi(val));
                else if (key == "NEWTON") c.newton = std::max(1, std::stoi(val));
                else if (key == "MAXITER") c.maxIter = std::max(1, std::stoi(val));
                else if (key == "MAXLENITER") c.maxLenIter = std::max(1, std::stoi(val));
                else if (key == "TOL") c.tol = std::stod(val);
                else if (key == "TOLT") c.tolTension = std::stod(val);
                else if (key == "METHOD") {
                    c.minEnergy = (upper(val) != "RIGID");
                } else if (key == "MINT") c.minTensionFraction = std::stod(val);
                else if (key == "MINABS") c.minTension = std::stod(val);
                else if (key == "SENSITIVITY") c.sensitivity = upper(val);
                else if (key == "DECKLOAD") c.deckLoad = std::stod(val);
                else if (key == "OPTIMIZE") c.optimize = (std::stod(val) != 0.0);
                else if (key == "OPTOBJ") c.optObjective = upper(val);
                else if (key == "OPTBALANCE") c.optBalance = (std::stod(val) != 0.0);
                else if (key == "TMAX") c.optTMax = std::stod(val);
                else if (key == "SYSGUARD") c.sysGuard = (std::stod(val) != 0.0);
                else if (key == "SYSDECKMIN") c.sysDeckMin = std::stod(val);
                else if (key == "SYSDECKMAX") c.sysDeckMax = std::stod(val);
                else if (key == "SYSADJMAX") c.sysAdjMax = std::stod(val);
                else if (key == "VERBOSE") c.verbose = (std::stod(val) != 0.0);
                else std::fprintf(stderr, "[警告] *CABLE_BRIDGE 未知关键字: %s\n", kv.c_str());
            }
        } else if (sec == "CABLE") {
            std::istringstream ss(line);
            CableReq r;
            if (ss >> r.elem >> r.anchor) {
                ss >> r.tension;          // 可缺省
                c.cables.push_back(r);
            }
        } else if (sec == "SUSP_BRIDGE") {
            std::istringstream ss(line);
            std::string kv;
            while (ss >> kv) {
                const std::size_t p = kv.find('=');
                if (p == std::string::npos) continue;
                const std::string key = upper(kv.substr(0, p));
                const std::string val = kv.substr(p + 1);
                if (key == "SYSTEM") c.system = upper(val);
                else if (key == "SAG") c.sag = std::stod(val);
                else if (key == "SAGRATIO") c.sagRatio = std::stod(val);
                else if (key == "TOWER1") c.tower1 = std::stoi(val);
                else if (key == "TOWER2") c.tower2 = std::stoi(val);
                else if (key == "DECKE0") c.deckE0 = std::stod(val);
                else if (key == "DECKE1") c.deckE1 = std::stod(val);
                else std::fprintf(stderr, "[警告] *SUSP_BRIDGE 未知关键字: %s\n",
                                  kv.c_str());
            }
        } else if (sec == "MAINCABLE" || sec == "HANGER") {
            auto& dst = (sec == "MAINCABLE") ? c.mainCable : c.hangers;
            std::istringstream ss(line);
            int id = 0;
            while (ss >> id) dst.push_back(id);   // 一行多个也可
        } else if (sec == "ERECTION") {
            std::istringstream ss(line);
            std::string kv;
            while (ss >> kv) {
                const std::size_t p = kv.find('=');
                if (p == std::string::npos) continue;
                const std::string key = upper(kv.substr(0, p));
                const std::string val = kv.substr(p + 1);
                if (key == "Q1") c.erectionQ1 = std::stod(val);
                else if (key == "Q2") c.erectionQ2 = std::stod(val);
                else if (key == "PANELS") c.erectionPanels = std::max(1, std::stoi(val));
                else std::fprintf(stderr, "[警告] *ERECTION 未知关键字: %s\n", kv.c_str());
            }
            c.erection = true;
        }
    }
    return c;
}

const Element3D* findElem(const Model3D& m, int id) {
    for (const auto& e : m.elements) if (e.id == id) return &e;
    return nullptr;
}

// 梁侧锚点 = 索较低的那端（斜拉索从塔顶下到梁面）。塔顶始终高于梁面，
// 这条规则对竖索/斜索都成立；真要反过来用 *CABLE 段显式给锚点即可。
int autoAnchor(const Model3D& m, const Element3D& e) {
    const double z1 = m.nodes.at(e.n1).z(), z2 = m.nodes.at(e.n2).z();
    return (z1 <= z2) ? e.n1 : e.n2;
}

const Element3D* findElem(const Model3D& m, int id);

Element3D* findElemMut(Model3D& m, int id) {
    for (auto& e : m.elements) if (e.id == id) return &e;
    return nullptr;
}

// ============================================================
// 悬索桥主缆找形
// ============================================================
//
// 为什么必须有这一步，而不是"给个抛物线初形让成桥态去调"
// ------------------------------------------------------
// 悬索桥与斜拉桥的差别不在索，在**设计线形本身由恒载决定**。斜拉桥的设计
// 线形是建筑师给的（梁面高程），索力去适应它；悬索桥的主缆只能垂成一条
// 在恒载下自平衡的索形——给错了，成桥态闭环就无解：
// iterateErectionGeometry 求的是"制造线形 ΔX 使 ΔX + u = 0"，
// 若设计线形不是平衡构形，则不存在任何 u 能同时满足平衡与回归设计线形，
// 迭代的闭合残差会卡在某个非零值上（实测卡在 30 mm 量级，看起来像"没收敛"，
// 其实是物理上无解）。
//
// 索形为什么按**点荷载**精确积分，而不是用抛物线 h = qL²/(8f)
// ----------------------------------------------------------
// 抛物线公式假定荷载沿水平均匀。实际荷载是每根吊索一个集中力，间距 8~15 m
// 时用均布近似会让索形在吊索之间偏出几毫米到几厘米；更麻烦的是它破坏
// 平衡——索形不是真正的索多边形，弯矩会以"梁不得不自己扛"的形式漏出来。
// 所以这里直接按索多边形积分：
//
//     dz/dx = (V0 − S(x)) / H          S(x) = 该截面以左的全部竖向荷载
//
// 其中 H 为水平分力（全桥恒定），V0 为起始端竖向分力。两个未知量由两个
// 条件定：主跨两端等高 ⇒ V0 由竖向力平衡定出；跨中矢高 = f ⇒ 对 H 二分。
// 这是悬索桥手算定线形的标准做法（"索多边形法"），只是把作图换成了积分。
//
// 为什么 H 在主跨与边跨取同一个值
// ------------------------------
// 主缆连续跨过塔顶鞍座。若两侧 H 不等，塔顶会承受不平衡水平力、塔身偏位，
// 设计上正是靠"边跨索形自动调整到 H 相等"来消除这个力（这也是为什么边跨
// 矢高不能随便给——它由主跨 H、边跨荷载和锚碇高程共同决定）。所以：先由
// 主跨定 H，再用同一个 H 积分边跨、反求边跨端部竖向力。
//
// 自锚式与地锚式的差别**只体现在边界条件上**，算法完全相同：
//   · 地锚式：主缆端部锚在锚碇节点，该节点全约束 → H 由锚碇承担；
//   · 自锚式：主缆端部就是加劲梁端节点 → H 压进加劲梁（这就是自锚式的
//     全部特点，也是它必须做"梁轴压"校核的原因）。

struct SuspSpan {
    std::vector<int> nodeIds;      // 该跨主缆节点（含两端）
    std::vector<double> x, q;      // 水平坐标、节点竖向荷载（向下为正）
    std::vector<double> w;         // 各段单位长重量 N/m（段数 = 节点数-1）
    double z0 = 0.0, zEnd = 0.0;
    bool symmetric = false;        // 两端等高（主跨）
};

struct SuspSpanFit {
    double H = 0.0, V0 = 0.0;
    double sag = 0.0;              // 跨中矢高（两端平均高程 − 跨中高程）
    std::vector<double> z;         // 各节点高程
    std::vector<double> segT;      // 各段轴力
};

// 给定 H，积分出索形。V0 由端部高程条件定（见文件头推导）。
SuspSpanFit integrateSpan(const SuspSpan& s, double H) {
    SuspSpanFit f;
    const std::size_t n = s.nodeIds.size();
    f.z.assign(n, 0.0);
    f.segT.assign(n > 0 ? n - 1 : 0, 0.0);
    f.H = H;
    if (n < 2) return f;
    // 段中点的"左側累计荷载"：节点 k 的荷载作用在段 k 的左端，算在段 k 内；
    // 自重按每段中点的半段 + 之前各段全重累计，漏掉累计项会让主缆自重
    // 对 H 完全不起作用（矢高被吊索力一项占满，自重再重也不加 H）。
    std::vector<double> S(n - 1, 0.0), dx(n - 1, 0.0);
    double cum = s.q.empty() ? 0.0 : s.q[0];
    double cumW = 0.0;
    for (std::size_t k = 0; k + 1 < n; ++k) {
        dx[k] = s.x[k + 1] - s.x[k];
        if (k > 0) cum += s.q[k];
        const double wk = (k < s.w.size()) ? s.w[k] : 0.0;
        S[k] = cum + cumW + 0.5 * wk * dx[k];
        cumW += wk * dx[k];
    }
    double L = 0.0, mom = 0.0;
    for (std::size_t k = 0; k + 1 < n; ++k) { L += dx[k]; mom += S[k] * dx[k]; }
    if (L <= 1e-12) return f;
    // 端部条件：z(L) − z(0) = (S·L − V0·L)/H  ⇒  V0 = (mom − H·Δz)/L
    // 主跨两端等高时 Δz=0，退化为 V0 = mom/L（竖向力自平衡）。
    // 符号写成 (H·Δz + mom)/L 会让边跨索"反向垂"——从锚碇往下掉十几米，
    // 初值离平衡态太远，Newton 随后以 |du|~10 m 振荡发散。
    f.V0 = (mom - H * (s.zEnd - s.z0)) / L;
    f.z[0] = s.z0;
    for (std::size_t k = 0; k + 1 < n; ++k) {
        // 斜率取 dz/dx = (S − V0)/H：索受竖向荷载后**向下垂**，左端的竖向上
        // 拔力 V0 使斜率为负。写成 (V0 − S)/H 会得到一条向上拱的"反索"，
        // 矢高变负、H 的二分也随之失去单调性（实测收敛到 H=0.1 kN）。
        f.z[k + 1] = f.z[k] + (S[k] - f.V0) * dx[k] / H;
        const double V = f.V0 - S[k];
        f.segT[k] = std::sqrt(H * H + V * V);
    }
    // 矢高按跨中量取（两端平均高程 − 跨中高程），与矢跨比的定义一致
    double zm = f.z[0];
    if (n == 2) {
        zm = 0.5 * (f.z[0] + f.z[1]);
    } else {
        std::size_t k = 0;
        while (k + 2 < n && s.x[k + 1] <= 0.5 * (s.x.front() + s.x.back())) ++k;
        const double t = (0.5 * (s.x.front() + s.x.back()) - s.x[k]) / dx[k];
        zm = f.z[k] + t * (f.z[k + 1] - f.z[k]);
    }
    f.sag = 0.5 * (f.z.front() + f.z.back()) - zm;
    return f;
}

struct SuspFit {
    bool ok = false;
    double H = 0.0;                // 主缆水平分力 N
    double sag = 0.0, sagTarget = 0.0, span = 0.0;
    std::vector<int> chainNodes;   // 主缆节点链（锚 → 塔1 → 塔2 → 锚）
    std::map<int, double> zNew;    // 找形后的主缆节点高程
    std::map<int, double> T;       // 单元 → 目标索力（主缆分段 + 吊索）
    std::map<int, int> hangerDeck; // 吊索单元 → 梁侧节点
    double V0Main = 0.0, V0Side1 = 0.0, V0Side2 = 0.0;
};

// 某节点的竖向荷载合计（向下为正）
double nodalDown(const Model3D& m, int node) {
    double v = 0.0;
    for (const auto& [n, d, val] : m.loads)
        if (n == node && d == 3) v += val;
    return -v;
}

// 该节点上是否接着梁单元（用来识别"主缆没挂上塔"这类隐蔽错误）
bool hasBeamAt(const Model3D& m, int node) {
    for (const auto& e : m.elements) {
        if (e.type != "BEAM" && e.type != "BEAM7") continue;
        if (e.n1 == node || e.n2 == node) return true;
    }
    return false;
}

bool isFixedNode(const Model3D& m, int node) {
    int cnt = 0;
    for (const auto& c : m.constraints)
        if (c.node == node && c.k <= 0.0 && c.dof >= 1 && c.dof <= 3) ++cnt;
    return cnt >= 3;
}

// 支反力：在受约束自由度上取 R = F_int − F_ext（向上为正）。
//
// 悬索桥的全局平衡**不能**用"索力竖向分量之和"来校核：主缆各段一段向上、
// 一段向下，求和结果取决于每段"锚点"被 autoAnchor 取到哪一端，得到的数纯
// 粹是噪声（实测能到恒载的 6 倍）。真正的平衡式是
//     Σ 支座竖向反力 = Σ 竖向荷载（节点荷载 + 缆索自重）
// 所以这里直接把支反力算出来。它顺带也是计算书必须给的"支座反力"一栏。
struct Reaction3D {
    std::vector<double> R;      // 每个自由度的支反力（N / N·m）
    double sumV = 0.0;          // 竖向支反力之和（向上为正）
    double appliedV = 0.0;      // 竖向荷载合计（向下为正）
    double maxV = 0.0;          // 最大竖向支反力绝对值
    int nd = 0;
    bool ok = false;
};

Reaction3D computeReactions3D(const Model3D& m, const DofMap3D& dm,
                              const std::vector<double>& u) {
    Reaction3D out;
    const int nd = dm.nDofs;
    out.nd = nd;
    if (static_cast<int>(u.size()) < nd) return out;
    out.R.assign(nd, 0.0);

    Eigen::VectorXd Fint = Eigen::VectorXd::Zero(nd);
    Eigen::VectorXd Fext = Eigen::VectorXd::Zero(nd);

    Eigen::Matrix<double, 14, 1> ue, fe;
    for (const auto& e : m.elements) {
        const CRElemData ed = crMakeElemData(m, e, dm);
        ue.setZero();
        for (int k = 0; k < ed.n; ++k) ue[k] = u[ed.gdofs[k]];
        fe.setZero();
        crElementInternal(ed, ue, fe, false);
        for (int k = 0; k < ed.n; ++k) Fint(ed.gdofs[k]) += fe[k];
    }
    for (const auto& [node, dof, val] : m.loads) {
        const int d = dm.safeDof(node, dof - 1);
        if (d >= 0) Fext(d) += val;
    }
    // 缆索自重的等效节点力，与内核保持一致（两端各半）
    for (const auto& e : m.elements) {
        if (e.type != "CABLE" || e.w <= 0.0) continue;
        const double Lc = (e.L0 > 0.0) ? e.L0
                                       : (m.nodes.at(e.n2) - m.nodes.at(e.n1)).norm();
        const double half = -0.5 * e.w * Lc;
        const int d1 = dm.safeDof(e.n1, 2);
        const int d2 = dm.safeDof(e.n2, 2);
        if (d1 >= 0) Fext(d1) += half;
        if (d2 >= 0) Fext(d2) += half;
        out.appliedV += 2.0 * (-half);
    }
    for (const auto& [node, dof, val] : m.loads)
        if (dof == 3) out.appliedV += -val;

    std::vector<char> isC(nd, 0);
    for (const auto& c : m.constraints) {
        if (c.k > 0.0) continue;
        const int d = dm.safeDof(c.node, c.dof - 1);
        if (d >= 0) isC[d] = 1;
    }
    for (int d = 0; d < nd; ++d) {
        if (!isC[d]) continue;
        out.R[d] = Fint(d) - Fext(d);
    }
    for (const auto& [nid, pos] : dm.index) {
        const int s = dm.dofStart[pos];
        if (dm.dofCount[pos] < 3) continue;
        if (!isC[s + 2]) continue;
        out.sumV += out.R[s + 2];
        out.maxV = std::max(out.maxV, std::abs(out.R[s + 2]));
    }
    out.ok = true;
    return out;
}

// 主缆单元列表 → 节点链。方向由前两个单元的公共端点判断，因此
// *MAINCABLE 段既可写成"锚→塔→锚"也可反着写，这点容错是必要的：
// 手写 .dat 排查时很容易把顺序写反，报"不连续"比猜方向有用。
std::vector<int> buildCableChain(const Model3D& m, const std::vector<int>& elems) {
    std::vector<int> chain;
    if (elems.empty()) return chain;
    const Element3D* e0 = findElem(m, elems[0]);
    if (!e0) throw std::runtime_error("主缆单元不存在: " + std::to_string(elems[0]));
    int start = e0->n1;
    if (elems.size() >= 2) {
        const Element3D* e1 = findElem(m, elems[1]);
        if (e1 && (e0->n2 == e1->n1 || e0->n2 == e1->n2)) start = e0->n1;
        else start = e0->n2;
    }
    chain.push_back(start);
    int cur = start;
    for (int id : elems) {
        const Element3D* e = findElem(m, id);
        if (!e) throw std::runtime_error("主缆单元不存在: " + std::to_string(id));
        int nxt = (e->n1 == cur) ? e->n2 : ((e->n2 == cur) ? e->n1 : -1);
        if (nxt < 0)
            throw std::runtime_error(
                "主缆单元 " + std::to_string(id) + " 与链端节点 " + std::to_string(cur)
                + " 不连续。*MAINCABLE 段必须按主缆的实际连接顺序排列。");
        chain.push_back(nxt);
        cur = nxt;
    }
    return chain;
}

SuspFit formFindSuspension(Model3D& m, const Cfg& cfg) {
    SuspFit out;
    if (cfg.mainCable.empty() || cfg.hangers.empty())
        throw std::runtime_error("悬索桥模式要求 *MAINCABLE 与 *HANGER 段都不为空");
    out.chainNodes = buildCableChain(m, cfg.mainCable);

    // ---- 塔顶在主缆链里的位置 ----
    auto posOf = [&](int nid) {
        for (std::size_t i = 0; i < out.chainNodes.size(); ++i)
            if (out.chainNodes[i] == nid) return static_cast<int>(i);
        return -1;
    };
    const int t1 = posOf(cfg.tower1), t2 = posOf(cfg.tower2);
    if (t1 < 0 || t2 < 0)
        throw std::runtime_error("TOWER1/TOWER2 指定的节点不在主缆链上");
    if (t1 >= t2)
        throw std::runtime_error("TOWER1 必须位于 TOWER2 左侧（链序上更靠前）");

    // ---- 校验：主缆必须**共用**塔顶与端部节点，不能只是坐标重合 ----
    //
    // 这是悬索桥建模最隐蔽的一处错：如果在塔顶位置另建一个坐标相同的主缆
    // 节点，主缆就成了一根从锚碇到锚碇自由悬挂的索——塔不参与受力、加劲梁
    // 也拿不到该拿的力。求解照样收敛、位移与倾角都"看着像座桥"，只有索力
    // 掉到设计值的百分之一、吊索全部松弛、桥面被迫自己受弯。不报错的错最贵，
    // 所以这里直接拦。
    for (int tn : {cfg.tower1, cfg.tower2})
        if (!hasBeamAt(m, tn))
            throw std::runtime_error(
                "TOWER1/TOWER2 = " + std::to_string(tn) +
                " 上没有梁单元：主缆必须与桥塔**共用同一个节点**，不能在塔顶旁边"
                "另建一个坐标相同的节点。否则主缆是一根从锚碇到锚碇的自由索，"
                "塔不参与受力。");
    for (int endNode : {out.chainNodes.front(), out.chainNodes.back()})
        if (!isFixedNode(m, endNode) && !hasBeamAt(m, endNode))
            throw std::runtime_error(
                "主缆端节点 " + std::to_string(endNode) +
                " 既不约束、也不接梁单元：地锚式应把端节点做成全约束的锚碇，"
                "自锚式应直接锚在加劲梁端节点上。");

    // ---- 吊索：把梁端竖向荷载搬到主缆节点上 ----
    std::map<int, double> qAtCable;       // 主缆节点 → 该点吊索力（向下为正）
    std::map<int, double> hangerT;        // 吊索单元 → 目标力
    std::map<int, int> hangerDeck;
    for (int hid : cfg.hangers) {
        const Element3D* e = findElem(m, hid);
        if (!e) throw std::runtime_error("吊索单元不存在: " + std::to_string(hid));
        // 哪端在主缆链上，哪端就是梁端
        bool n1Cable = posOf(e->n1) >= 0, n2Cable = posOf(e->n2) >= 0;
        if (n1Cable == n2Cable)
            throw std::runtime_error("吊索 " + std::to_string(hid)
                                     + " 必须一端接主缆、一端接加劲梁");
        const int cableNode = n1Cable ? e->n1 : e->n2;
        const int deckNode = n1Cable ? e->n2 : e->n1;
        const double Q = nodalDown(m, deckNode);
        if (Q <= 0.0)
            throw std::runtime_error(
                "吊索 " + std::to_string(hid) + " 的梁端节点 " + std::to_string(deckNode)
                + " 上没有向下的荷载：悬索桥的吊索力由该节点承担的恒载确定，"
                "取不到荷载就无法定索力。请检查该节点是否漏了桥面荷载。");
        qAtCable[cableNode] += Q;
        hangerT[hid] = Q;
        hangerDeck[hid] = deckNode;
    }

    // ---- 三跨：锚→塔1、塔1→塔2、塔2→锚 ----
    auto makeSpan = [&](int i0, int i1, bool symmetric) {
        SuspSpan s;
        s.symmetric = symmetric;
        for (int i = i0; i <= i1; ++i) {
            const int nid = out.chainNodes[i];
            s.nodeIds.push_back(nid);
            s.x.push_back(m.nodes.at(nid).x());
            auto it = qAtCable.find(nid);
            s.q.push_back(it == qAtCable.end() ? 0.0 : it->second);
        }
        // 各段单位长重量取该段单元自己的 w（主缆截面沿桥可分段不同）
        for (int i = i0; i < i1; ++i) {
            const Element3D* e = findElem(m, cfg.mainCable[i]);
            s.w.push_back(e ? e->w : 0.0);
        }
        s.z0 = m.nodes.at(s.nodeIds.front()).z();
        s.zEnd = m.nodes.at(s.nodeIds.back()).z();
        return s;
    };

    const SuspSpan spMain = makeSpan(t1, t2, std::abs(m.nodes.at(cfg.tower1).z()
                                                      - m.nodes.at(cfg.tower2).z()) < 1e-9);
    const SuspSpan spS1 = makeSpan(0, t1, false);
    const SuspSpan spS2 = makeSpan(t2, static_cast<int>(out.chainNodes.size()) - 1, false);

    // ---- 主跨：对 H 二分，命中目标矢高 ----
    const double Lm = spMain.x.back() - spMain.x.front();
    if (Lm <= 0.0) throw std::runtime_error("主跨水平长度非正，检查 TOWER1/TOWER2");
    out.span = Lm;
    double f = cfg.sag;
    if (f <= 0.0) f = (cfg.sagRatio > 0.0 ? cfg.sagRatio : 0.1) * Lm;
    out.sagTarget = f;
    // 矢高随 H 单调减小，二分必收敛（索形积分对 H 连续）
    double lo = 1e2, hi = 1e12;
    for (int it = 0; it < 200; ++it) {
        const double mid = std::sqrt(lo * hi);
        const SuspSpanFit t = integrateSpan(spMain, mid);
        if (t.sag > f) lo = mid; else hi = mid;
    }
    out.H = std::sqrt(lo * hi);
    const SuspSpanFit fMain = integrateSpan(spMain, out.H);
    out.sag = fMain.sag;
    out.V0Main = fMain.V0;

    // ---- 边跨：同一个 H，端部高程定竖向力 ----
    const SuspSpanFit fS1 = integrateSpan(spS1, out.H);
    const SuspSpanFit fS2 = integrateSpan(spS2, out.H);
    out.V0Side1 = fS1.V0;
    out.V0Side2 = fS2.V0;

    // ---- 落盘：新高程 + 主缆各段目标力 ----
    auto applySpan = [&](const SuspSpan& s, const SuspSpanFit& t, int i0) {
        for (std::size_t k = 0; k < s.nodeIds.size(); ++k) {
            Vec3 p = m.nodes.at(s.nodeIds[k]);
            p.z() = t.z[k];
            m.nodes[s.nodeIds[k]] = p;
            out.zNew[s.nodeIds[k]] = t.z[k];
        }
        for (std::size_t k = 0; k + 1 < s.nodeIds.size(); ++k)
            out.T[cfg.mainCable[i0 + static_cast<int>(k)]] = t.segT[k];
    };
    applySpan(spS1, fS1, 0);
    applySpan(spMain, fMain, t1);
    applySpan(spS2, fS2, t2);

    for (const auto& [hid, T] : hangerT) out.T[hid] = T;
    out.hangerDeck = hangerDeck;
    out.ok = true;
    return out;
}

// ============================================================
// 悬索桥施工阶段：给定主缆下料长度 L0 的逐阶段找形
// ============================================================
//
// 为什么不能复用 formFindSuspension
// ---------------------------------
// formFindSuspension 是"给定设计矢高、反算 H"——那是**设计态**。施工阶段主缆
// 矢高是**结果不是输入**：主缆下料长度 L0 一旦固定，从空缆（只受自重）到成桥
// （全恒载），主缆逐级下沉，矢高从 f0 增大到成桥矢高 f。所以这里二分的不再是
// "矢高命中"，而是"整根主缆的无应力长度 L0 命中"——给定 L0，反算当前荷载下的
// H，矢高与线形随之出来。这正是"从中间吊装主梁"的施工控制要的量。

struct ErectionState {
    double H = 0.0;                 // 主缆水平分力 N
    double sag = 0.0;               // 主跨矢高 m（结果，不是输入）
    double midZ = 0.0;              // 主跨跨中主缆高程 m
    std::map<int, double> zNew;     // 主缆节点高程（找形后，链序）
    std::map<int, double> segT;     // 主缆段轴力（单元 id → T）
    std::map<int, double> hangerT;  // 吊索力（单元 id → T）
    std::vector<double> mainX, mainZ;   // 主跨主缆节点（塔1→塔2）的 x / 高程
    // 逐阶段对标用（ANSYS 侧无需重演 FEM 内部吊装规则）：
    std::vector<int> mainIds;       // 主跨主缆节点 id，顺序同 mainX/mainZ
    std::vector<double> mainQ;      // 吊索搬运到各主缆节点的竖向荷载 N（向下为正）
    std::vector<double> mainT;     // 主跨各段轴力 N（顺序 mainX[k]→mainX[k+1]）
};

ErectionState findErectionState(const Model3D& m, const Cfg& cfg, double L0target) {
    ErectionState st;
    const std::vector<int> chain = buildCableChain(m, cfg.mainCable);
    auto posOf = [&](int nid) {
        for (std::size_t i = 0; i < chain.size(); ++i)
            if (chain[i] == nid) return static_cast<int>(i);
        return -1;
    };
    const int t1 = posOf(cfg.tower1), t2 = posOf(cfg.tower2);
    if (t1 < 0 || t2 < 0) throw std::runtime_error("TOWER1/TOWER2 不在主缆链上");
    if (t1 >= t2) throw std::runtime_error("TOWER1 必须位于 TOWER2 左侧");

    // 吊索荷载搬运（与 formFindSuspension 同一口径）
    std::map<int, double> qAtCable;
    for (int hid : cfg.hangers) {
        const Element3D* e = findElem(m, hid);
        if (!e) throw std::runtime_error("吊索单元不存在: " + std::to_string(hid));
        const bool n1Cable = posOf(e->n1) >= 0, n2Cable = posOf(e->n2) >= 0;
        if (n1Cable == n2Cable)
            throw std::runtime_error("吊索 " + std::to_string(hid)
                                     + " 必须一端接主缆、一端接加劲梁");
        const int cableNode = n1Cable ? e->n1 : e->n2;
        const int deckNode = n1Cable ? e->n2 : e->n1;
        const double Q = std::max(nodalDown(m, deckNode), 0.0);
        qAtCable[cableNode] += Q;
        st.hangerT[hid] = Q;
    }

    // 三跨 span（同 formFindSuspension 的 makeSpan）
    auto makeSpan = [&](int i0, int i1, bool symmetric) {
        SuspSpan s;
        s.symmetric = symmetric;
        for (int i = i0; i <= i1; ++i) {
            const int nid = chain[i];
            s.nodeIds.push_back(nid);
            s.x.push_back(m.nodes.at(nid).x());
            auto it = qAtCable.find(nid);
            s.q.push_back(it == qAtCable.end() ? 0.0 : it->second);
        }
        for (int i = i0; i < i1; ++i) {
            const Element3D* e = findElem(m, cfg.mainCable[i]);
            s.w.push_back(e ? e->w : 0.0);
        }
        s.z0 = m.nodes.at(s.nodeIds.front()).z();
        s.zEnd = m.nodes.at(s.nodeIds.back()).z();
        return s;
    };
    const SuspSpan spMain = makeSpan(t1, t2,
        std::abs(m.nodes.at(cfg.tower1).z() - m.nodes.at(cfg.tower2).z()) < 1e-9);
    const SuspSpan spS1 = makeSpan(0, t1, false);
    const SuspSpan spS2 = makeSpan(t2, static_cast<int>(chain.size()) - 1, false);

    // 主缆轴向刚度（全桥主缆同材质，取第一段）
    double EA = 0.0;
    if (!cfg.mainCable.empty()) {
        const Element3D* e = findElem(m, cfg.mainCable[0]);
        if (e) EA = e->E * e->A;
    }
    if (EA <= 0.0) throw std::runtime_error("主缆单元缺少 E/A，无法反算下料长度");

    // 总无应力长度：三跨索多边形弧长 − 弹性伸长。各段用与 COMPLETION 相同的
    // 闭式口径 cableUnstrainedLengthForTension（弦长 + 垂度，再除 1+T/EA）。
    auto totalL0 = [&](double H) {
        double l0 = 0.0;
        auto acc = [&](const SuspSpan& s) {
            const SuspSpanFit f = integrateSpan(s, H);
            for (std::size_t k = 0; k + 1 < s.nodeIds.size(); ++k) {
                const double dx = s.x[k + 1] - s.x[k];
                const double dz = f.z[k + 1] - f.z[k];
                const double chord = std::sqrt(dx * dx + dz * dz);
                const double w = (k < s.w.size()) ? s.w[k] : 0.0;
                l0 += cableUnstrainedLengthForTension(chord, f.segT[k], EA, w, dx);
            }
        };
        acc(spS1);
        acc(spMain);
        acc(spS2);
        return l0;
    };

    // 二分 H：L0(H) 随 H 单调递减（H 越大主缆越平、弧长越短）
    double lo = 1e2, hi = 1e12;
    for (int it = 0; it < 200; ++it) {
        const double mid = std::sqrt(lo * hi);
        if (totalL0(mid) > L0target) lo = mid; else hi = mid;
    }
    st.H = std::sqrt(lo * hi);

    // 落盘最终状态
    const SuspSpanFit fMain = integrateSpan(spMain, st.H);
    const SuspSpanFit fS1 = integrateSpan(spS1, st.H);
    const SuspSpanFit fS2 = integrateSpan(spS2, st.H);
    st.sag = fMain.sag;
    auto apply = [&](const SuspSpan& s, const SuspSpanFit& f, int i0) {
        for (std::size_t k = 0; k < s.nodeIds.size(); ++k)
            st.zNew[s.nodeIds[k]] = f.z[k];
        for (std::size_t k = 0; k + 1 < s.nodeIds.size(); ++k)
            st.segT[cfg.mainCable[i0 + static_cast<int>(k)]] = f.segT[k];
    };
    apply(spS1, fS1, 0);
    apply(spMain, fMain, t1);
    apply(spS2, fS2, t2);

    // 主跨跨中主缆高程（离跨中最近的节点）
    const double xmid = 0.5 * (spMain.x.front() + spMain.x.back());
    double bd = 1e300;
    for (std::size_t k = 0; k < spMain.nodeIds.size(); ++k) {
        const double d = std::abs(spMain.x[k] - xmid);
        if (d < bd) { bd = d; st.midZ = fMain.z[k]; }
    }
    st.mainX = spMain.x;
    st.mainZ = fMain.z;
    st.mainT = fMain.segT;
    st.mainIds = spMain.nodeIds;
    st.mainQ.assign(spMain.nodeIds.size(), 0.0);
    for (std::size_t k = 0; k < spMain.nodeIds.size(); ++k) {
        auto qit = qAtCable.find(spMain.nodeIds[k]);
        if (qit != qAtCable.end()) st.mainQ[k] = qit->second;
    }
    return st;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("用法: %s <model.dat>\n"
                    "  读 *CABLE_BRIDGE / *CABLE，输出 <stem>_cable.txt\n", argv[0]);
        return 1;
    }
    const std::string dat = argv[1];
    try {
        Cfg cfg = readCfg(dat);
        Model3D m = readModel3D(dat);          // 参考构形（只读副本，后续传值）

        // ---- 悬索桥：先找形，再走成桥态闭环 ----
        // 找形会改写 m.nodes 的主缆高程（设计线形 = 恒载索形），这是悬索桥
        // 与斜拉桥在求解流程上唯一的实质差别。
        const bool susp = !cfg.system.empty();
        SuspFit fit;
        if (susp) {
            if (cfg.mode == "STATIC")
                throw std::runtime_error(
                    "悬索桥不支持 MODE=STATIC：主缆线形必须由找形确定，"
                    "不能用 .dat 里的几何直接做静力分析。请用 COMPLETION。");
            fit = formFindSuspension(m, cfg);
            // 目标索力全部由找形与荷载推出；*CABLE 段里写的力不再采用
            cfg.cables.clear();
            for (int id : cfg.mainCable) cfg.cables.push_back({id, 0, fit.T[id]});
            for (int id : cfg.hangers)   cfg.cables.push_back({id, 0, fit.T[id]});
            std::printf("悬索桥找形：H=%.1f kN  矢高 %.3f m（目标 %.3f m，矢跨比 1/%.1f）\n",
                        fit.H / 1e3, fit.sag, fit.sagTarget,
                        fit.sag > 0 ? fit.span / fit.sag : 0.0);
        }

        if (cfg.cables.empty()) {
            for (const auto& e : m.elements) {
                if (e.type != "CABLE") continue;
                cfg.cables.push_back({e.id, autoAnchor(m, e), 0.0});
            }
        }
        if (cfg.cables.empty())
            throw std::runtime_error("模型中没有 CABLE 单元，也没有 *CABLE 段");

        for (auto& r : cfg.cables) {
            const Element3D* e = findElem(m, r.elem);
            if (!e) throw std::runtime_error("索单元不存在: " + std::to_string(r.elem));
            if (e->type != "CABLE")
                throw std::runtime_error("单元 " + std::to_string(r.elem) + " 不是 CABLE");
            if (r.anchor <= 0) r.anchor = autoAnchor(m, *e);
            if (r.anchor != e->n1 && r.anchor != e->n2)
                throw std::runtime_error("锚点 " + std::to_string(r.anchor)
                                         + " 不是索单元 " + std::to_string(e->id) + " 的端点");
        }

        // ---- 带约束成桥索力优化（OPTIMIZE=1，默认关闭） ------------------
        // 在三种 MODE 之前先求 T* 并写回 cfg.cables，后面 TENSION / COMPLETION
        // 仍走各自既有的「命中给定索力」路径：优化只负责"定目标"，不另造闭环。
        // 失败不静默：写 [OPT] ERROR 段、保留原索力、stderr 给警告，绝不
        // 拿一组半吊子的 T* 冒充优化结果。
        CableOptResult3D optResult;
        bool optRequested = cfg.optimize && !susp;
        // anchors 提到块外：成桥态求解后回填 SYSCHECK_FINAL 还要用同一组
        // 索几何与顺序，不能另造一份。
        std::vector<CableAnchor3D> anchors;
        if (cfg.optimize && susp)
            std::fprintf(stderr,
                         "[警告] 悬索桥找形已定主缆/吊索索力，OPTIMIZE=1 不适用，已忽略。\n");
        if (optRequested) {
            anchors.reserve(cfg.cables.size());
            std::vector<double> givenT;
            givenT.reserve(cfg.cables.size());
            for (const auto& c : cfg.cables) {
                anchors.push_back({c.elem, c.anchor});
                givenT.push_back(c.tension);   // <=0 = 该索未给定
            }
            CableOptOptions3D oo;
            oo.objective = (cfg.optObjective == "MINENERGY")
                               ? CableOptOptions3D::Objective::MinBendingEnergy
                               : CableOptOptions3D::Objective::MinMaxMoment;
            oo.towerBalance = cfg.optBalance;
            oo.tMin = cfg.minTension;
            oo.tMax = cfg.optTMax;
            oo.minTensionFraction = cfg.minTensionFraction;
            oo.verbose = cfg.verbose;
            // 用户把 *CABLE 段目标索力全给齐时，优化从这组索力出发
            // （START=GIVEN）；否则退回刚性支承连续梁法解析起点（START=RIGID）。
            optResult = optimizeCableTensions3D(m, anchors, oo, givenT);
            optResult.sysCheckEnabled = true;
            // 体系守卫：写回 T* **之前**，用与计算书工程校核表同一口径预估
            // 三条体系判据（索承担恒载 / 支座不上拔 / 同跨相邻索力比）。
            // 不合格就不写回，STATUS=system_check_rejected、APPLIED=0，
            // 下游不得拿 T 行当设计索力。DECKLOAD 缺失或支座断面不足时
            // sysPred.available=false，段里写 NOT_AVAILABLE，不假判合格。
            optResult.sysPred = evaluateCableOptSystem3D(
                m, anchors, optResult.T, cfg.deckLoad,
                cfg.sysDeckMin, cfg.sysDeckMax, cfg.sysAdjMax, cfg.sysGuard);
            // REJECT 行始终如实列出体系判据的原因码，与 STATUS 独立：数学
            // 守卫 objective_not_improved 优先时 STATUS 不改成体系拒绝，但
            // SYSCHECK_PASS=0 却没有原因码会让下游自检判定自相矛盾。
            if (optResult.sysPred.available && !optResult.sysPred.pass)
                optResult.rejectCodes = optResult.sysPred.reject;
            if (optResult.error.empty()
                && optResult.status != "objective_not_improved"
                && optResult.sysPred.mode == "GUARDED"
                && !optResult.sysPred.pass) {
                optResult.status = "system_check_rejected";
            }
            const bool ok = optResult.error.empty()
                && optResult.status != "objective_not_improved"
                && optResult.status != "system_check_rejected"
                && optResult.T.size() == cfg.cables.size();
            if (ok) {
                for (std::size_t k = 0; k < cfg.cables.size(); ++k)
                    cfg.cables[k].tension = optResult.T[k];
            } else {
                // 不静默：优化失败、线性化目标变差、体系判据不合格都保留
                // 原索力，[OPT] 段写 APPLIED=0 并列明拒绝原因码。
                std::string why;
                if (!optResult.error.empty()) why = optResult.error;
                else if (optResult.status == "objective_not_improved")
                    why = "线性化目标未改善";
                else if (optResult.status == "system_check_rejected") {
                    why = "体系判据不合格：";
                    for (std::size_t i = 0; i < optResult.rejectCodes.size(); ++i) {
                        if (i) why += ",";
                        why += optResult.rejectCodes[i];
                    }
                } else why = "结果长度不符";
                std::fprintf(stderr,
                             "[警告] 成桥索力优化未采纳（%s），保留原索力口径。\n",
                             why.c_str());
            }
        }
        const bool optApplied = optRequested && optResult.error.empty()
            && optResult.status != "objective_not_improved"
            && optResult.status != "system_check_rejected"
            && optResult.T.size() == cfg.cables.size();

        const DofMap3D dm = buildDofMap3D(m);
        const std::size_t nc = cfg.cables.size();

        // 总恒载（竖向向下为正），用于全局平衡校核
        double Wtot = 0.0;
        for (const auto& [node, dof, val] : m.loads) if (dof == 3) Wtot -= val;

        NonlinearOptions3D nl;
        nl.loadSteps = cfg.nlSteps;
        nl.maxIter = cfg.newton;
        nl.tolForce = 1e-8;
        nl.verbose = cfg.verbose;

        bool converged = false;
        int iterations = 0;
        double closureErr = 0.0, tensionErr = 0.0, elapsed = 0.0;
        std::vector<double> T(nc, 0.0), L0(nc, 0.0), chord(nc, 0.0), Ltar(nc, 0.0);
        std::vector<double> u;
        std::map<int, Vec3> erected;
        std::vector<double> camber;
        std::string method = "given";
        // 逐轮过程留痕：MODE=TENSION / COMPLETION 才填，STATIC 与悬索桥找形
        // 是闭式反算 L0 解一次，本来就没有外层迭代（输出 ERROR 段说明）。
        IterationTrace3D lengthTrace;
        // 内层 Newton 过程：只保留**最终成桥态那一次**非线性平衡（外层每轮都
        // 会重跑，把几十次中间求解全写进结果文件没有意义也不好读）。
        NewtonTrace3D finalNewton;

        const double t0 = nowMs();
        if (susp) {
            // 悬索桥成桥态：无应力索长反算 + 一次几何非线性求解。
            //
            // 为什么不用 solveCompletionState3D（斜拉桥那套"制造线形闭环"）
            // -----------------------------------------------------------------
            // 那套方法的未知量是**节点坐标偏移** ΔX（求 ΔX 使 ΔX+u=0）。
            // 对斜拉桥成立：梁面高程由人给定、可预抬。对悬索桥不成立——主缆
            // 的节点位置不是可制造量，可制造的是**无应力索长**（下料长度）；
            // 主缆形状是恒载索形的结果。把主缆节点坐标也当未知量，方程比
            // 未知量少（21 根索要管 100 多个平动自由度），迭代会漂到索力只剩
            // 几十 kN 的"松弛构型"上，还报 converged=0（实测 ΣV 只有恒载的 3%）。
            //
            // 现在这条路是：找形给出设计几何与各索目标力 → 由"设计弦长 + 目标
            // 索长力"闭式反算 L0 → 解一次。因为设计几何**就是**恒载下的索形，
            // 平衡方程在设计几何上近似恒等成立，一次求解的位移就是弹性压缩量
            // 那一级（毫米级）。这个位移取负就是预拱度建议值。
            // 索长迭代：L0 反算用**变形后弦长**，不是设计弦长。
            //
            // 设计弦长与平衡弦长差一个弹性缩短量——自锚式可达 100 mm 以上
            // （梁被 H 压缩，主缆跨度随之缩短）。用它反算 L0，达到的索力会
            // 系统性偏低 1~6%，这点差额最后全部变成加劲梁的恒载弯矩。
            // 每轮只改 L0、不动几何，2~3 轮即收敛；因为弦长对 L0 的敏感度
            // 是 EA/L 量级（10⁸ N/m），本身就是个压缩映射。
            std::vector<double> chordNow(nc, 0.0);
            for (std::size_t k = 0; k < nc; ++k) {
                const Element3D* e = findElem(m, cfg.cables[k].elem);
                chordNow[k] = (m.nodes.at(e->n2) - m.nodes.at(e->n1)).norm();
            }
            NonlinearResult3D r;
            for (int cyc = 0; cyc < 8; ++cyc) {
                for (std::size_t k = 0; k < nc; ++k) {
                    Element3D* e = findElemMut(m, cfg.cables[k].elem);
                    if (!e) throw std::runtime_error("找不到索单元 "
                                                     + std::to_string(cfg.cables[k].elem));
                    const Vec3 d3 = m.nodes.at(e->n2) - m.nodes.at(e->n1);
                    const double Lh = std::sqrt(d3.x() * d3.x() + d3.y() * d3.y());
                    e->L0 = cableUnstrainedLengthForTension(chordNow[k],
                                                            cfg.cables[k].tension,
                                                            e->E * e->A, e->w, Lh);
                }
                r = solveModel3DNonlinear(m, nl);
                finalNewton = r.trace;
                if (!r.converged) break;
                double dev = 0.0;
                for (std::size_t k = 0; k < nc; ++k) {
                    std::size_t idx = 0;
                    while (idx < r.elemIds.size() && r.elemIds[idx] != cfg.cables[k].elem)
                        ++idx;
                    if (idx >= r.elemIds.size()) break;
                    const double tt = cfg.cables[k].tension;
                    if (tt > 0.0) dev = std::max(dev, std::abs(r.axial[idx] - tt) / tt);
                    chordNow[k] = r.Ln[idx];
                }
                if (dev < 1e-3) break;         // 索力已达标，不必再轮
            }
            converged = r.converged;
            iterations = r.totalIterations;
            u = r.u;
            T.assign(nc, 0.0); L0.assign(nc, 0.0);
            chord.assign(nc, 0.0); Ltar.assign(nc, 0.0);
            for (std::size_t k = 0; k < nc; ++k) {
                std::size_t idx = 0;
                while (idx < r.elemIds.size() && r.elemIds[idx] != cfg.cables[k].elem) ++idx;
                if (idx >= r.elemIds.size())
                    throw std::runtime_error("结果缺少索单元 "
                                             + std::to_string(cfg.cables[k].elem));
                T[k] = r.axial[idx];
                chord[k] = r.Ln[idx];
                const Element3D* e = findElem(m, cfg.cables[k].elem);
                L0[k] = e->L0;
                Ltar[k] = (m.nodes.at(e->n2) - m.nodes.at(e->n1)).norm();
            }
            // 预拱度 = −u：制造线形 = 设计线形 − 成桥位移
            camber.assign(static_cast<std::size_t>(dm.nDofs), 0.0);
            for (std::size_t i = 0; i < u.size() && i < camber.size(); ++i) camber[i] = -u[i];
            for (const auto& [id, x] : m.nodes) {
                Vec3 p = x;
                const int s = dm.dofStart[dm.index.at(id)];
                const int cnt = std::min(dm.dofCount[dm.index.at(id)], 3);
                for (int a = 0; a < cnt; ++a) p[a] += camber[s + a];
                erected[id] = p;
            }
            method = "FormFinding";
            // 悬索桥的两个诚实指标：
            //   closure_err  = 设计几何与平衡构形的最大偏差（max|u|）。
            //                  不是"弦长与下料长之差"——那对索就是弹性伸长
            //                  （实测 2 cm，看着像没收敛，其实是正常伸长量）。
            //   tension_err  = 达到索力相对目标索力的最大偏差。悬索桥的索力是
            //                  找形给的，能否达到是这条路线的核心校核。
            closureErr = 0.0;
            for (const auto& [id, p] : m.nodes) {
                (void)p;
                auto it = dm.index.find(id);
                if (it == dm.index.end()) continue;
                const int s = dm.dofStart[it->second];
                const int cnt = std::min(dm.dofCount[it->second], 3);
                double du2 = 0.0;
                for (int a = 0; a < cnt; ++a) du2 += u[s + a] * u[s + a];
                closureErr = std::max(closureErr, std::sqrt(du2));
            }
            tensionErr = 0.0;
            for (std::size_t k = 0; k < nc; ++k) {
                const double tt = cfg.cables[k].tension;
                if (tt > 0.0) tensionErr = std::max(tensionErr, std::abs(T[k] - tt) / tt);
            }
        } else if (cfg.mode == "STATIC") {
            // 只算：索的 L0 来自 .dat（每根索的 *R W=/L0= 或 *EN 覆盖）。
            const NonlinearResult3D r = solveModel3DNonlinear(m, nl);
            finalNewton = r.trace;
            converged = r.converged;
            iterations = r.totalIterations;
            u = r.u;
            for (std::size_t k = 0; k < nc; ++k) {
                std::size_t idx = 0;
                while (idx < r.elemIds.size() && r.elemIds[idx] != cfg.cables[k].elem) ++idx;
                if (idx >= r.elemIds.size())
                    throw std::runtime_error("结果缺少索单元 " +
                                             std::to_string(cfg.cables[k].elem));
                T[k] = r.axial[idx];
                chord[k] = r.Ln[idx];
                const Element3D* e = findElem(m, cfg.cables[k].elem);
                L0[k] = (e->L0 > 0.0) ? e->L0
                                      : (m.nodes.at(e->n2) - m.nodes.at(e->n1)).norm();
                Ltar[k] = (m.nodes.at(e->n2) - m.nodes.at(e->n1)).norm();
            }
        } else if (cfg.mode == "TENSION") {
            std::vector<CableTarget3D> tg;
            for (const auto& c : cfg.cables) {
                if (c.tension <= 0.0)
                    throw std::runtime_error(
                        "MODE=TENSION 要求 *CABLE 段给全目标索力（单元 "
                        + std::to_string(c.elem) + " 未给）");
                tg.push_back({c.elem, c.tension});
            }
            CableLengthOptions3D co;
            co.maxIter = cfg.maxLenIter;
            co.tolTension = cfg.tolTension;
            co.nl = nl;
            const CableLengthResult3D r = iterateCableUnstrainedLengths(m, tg, co);
            lengthTrace = r.trace;
            finalNewton = r.lastNewton;
            converged = r.converged;
            iterations = r.iterations;
            T = r.tension;
            L0 = r.L0;
            chord = r.chord;
            Ltar = r.Ltarget;
            u = r.u;
            tensionErr = r.maxTensionErr;
            closureErr = r.maxGeomErr;
            method = optApplied ? "ConstrainedOpt" : "target";
        } else if (cfg.mode == "COMPLETION") {
            std::vector<CompletionCable3D> cc;
            for (const auto& c : cfg.cables) cc.push_back({c.elem, c.anchor, c.tension});
            CompletionStateOptions3D co;
            co.maxIter = cfg.maxIter;
            co.tolClosure = cfg.tol;
            co.verbose = cfg.verbose;
            co.cableForce = cfg.minEnergy
                                ? CompletionStateOptions3D::CableForceMethod::MinBendingEnergy
                                : CompletionStateOptions3D::CableForceMethod::RigidSupport;
            co.minEnergy.minTension = cfg.minTension;
            co.minEnergy.minTensionFraction = cfg.minTensionFraction;
            co.minEnergy.verbose = cfg.verbose;
            if (cfg.sensitivity == "FULL")
                co.sensitivity = SensitivityMode3D::Full;
            else if (cfg.sensitivity == "FIXEDPOINT")
                co.sensitivity = SensitivityMode3D::FixedPoint;
            else if (cfg.sensitivity == "LOWRANK")
                co.sensitivity = SensitivityMode3D::LowRank;
            else
                throw std::runtime_error("未知 SENSITIVITY: " + cfg.sensitivity
                    + "（应为 LOWRANK / FULL / FIXEDPOINT）");
            co.nl = nl;
            bool anyAuto = false, allAuto = true;
            for (const auto& c : cfg.cables) {
                if (c.tension > 0.0) allAuto = false; else anyAuto = true;
            }
            const CompletionStateResult3D r = solveCompletionState3D(m, cc, co);
            lengthTrace = r.trace;
            finalNewton = r.lastNewton;
            converged = r.converged;
            iterations = r.iterations;
            T = r.tension;
            L0 = r.L0;
            chord = r.chord;
            Ltar = r.Ltarget;
            u = r.u;
            erected = r.erected;
            camber = r.camber;
            closureErr = r.maxClosureErr;
            tensionErr = r.maxTensionErr;
            method = optApplied ? "ConstrainedOpt"
                                : (allAuto ? (cfg.minEnergy ? "MinBendingEnergy"
                                                            : "RigidSupport")
                                           : (anyAuto ? "Mixed" : "given"));
        } else {
            throw std::runtime_error("未知 MODE: " + cfg.mode
                                     + "（应为 COMPLETION / TENSION / STATIC）");
        }
        elapsed = nowMs() - t0;

        std::set<int> hangerSet(cfg.hangers.begin(), cfg.hangers.end());
        double sumT = 0.0, sumV = 0.0, tmin = 1e300, tmax = -1e300;
        double hangerTmin = 1e300, hangerTmax = -1e300, hangerTsum = 0.0, hangerVsum = 0.0;
        // 汽车冲击系数（JTG D60）：基于成桥态切线刚度的竖向一阶弯曲基频。
        // 仅在活载包络段内填充；未收敛/无车道时保持 0/0/1.0。
        double f1Vertical = 0.0;    // 竖向一阶弯曲基频 Hz（0 = 未计算）
        double impactMu = 0.0;      // 冲击系数 μ
        double impactFactor = 1.0;  // 活载冲击放大 (1+μ)
        for (std::size_t k = 0; k < nc; ++k) {
            const Element3D* e = findElem(m, cfg.cables[k].elem);
            const int A = cfg.cables[k].anchor;
            const int B = (e->n1 == A) ? e->n2 : e->n1;
            const Vec3 ev = (m.nodes.at(B) - m.nodes.at(A)).normalized();
            sumT += T[k];
            sumV += T[k] * ev.z();
            tmin = std::min(tmin, T[k]);
            tmax = std::max(tmax, T[k]);
            if (susp && hangerSet.count(cfg.cables[k].elem)) {
                hangerTsum += T[k];
                hangerVsum += T[k] * ev.z();
                hangerTmin = std::min(hangerTmin, T[k]);
                hangerTmax = std::max(hangerTmax, T[k]);
            }
        }

        // 体系守卫实测回填：用**成桥态索力 T**（不是 T* 过程量）按同一口径
        // 再算一次，写进 [OPT] 让计算书把「写回前的守卫预估」与「本次实际
        // 成桥态」并排看到。APPLIED=0 时这里自然是原方案的实测值。
        if (optRequested && !anchors.empty()) {
            optResult.sysFinal = evaluateCableOptSystem3D(
                m, anchors, T, cfg.deckLoad,
                cfg.sysDeckMin, cfg.sysDeckMax, cfg.sysAdjMax, cfg.sysGuard);
            optResult.sysFinalFilled = optResult.sysFinal.available;
        }

        // 全局竖向平衡
        //
        // 斜拉桥：ΣT·ez 有意义——每根索都锚在梁上，其竖向分量之和就是"索承担
        //         的桥面恒载"，差额由支座承担。
        // 悬索桥：ΣT·ez 是噪声（主缆各段上下交替，求和取决于 autoAnchor 取到
        //         哪一端，实测能到恒载的 6 倍）。改用真正的平衡式：
        //             Σ 支座竖向反力 = Σ 竖向荷载
        const Reaction3D rea = computeReactions3D(m, dm, u);
        const bool useReact = susp && rea.ok;
        const double balSum = useReact ? rea.sumV : sumV;
        const double balTot = useReact ? rea.appliedV : Wtot;

        // ---- 活载梁内力包络（成桥态切线刚度影响线）----
        //
        // 斜拉桥/悬索桥的加劲梁是压弯构件，活载包络必须基于成桥态一致切线
        // 刚度（张紧索的切线刚度远大于弹性刚度），线弹性影响线会严重失真。
        // 恒载已由非线性求解精确算入；活载增量按一阶在切线刚度上提取——
        // 梁的几何刚度（P-Δ）是二阶小量，属标准做法，文档如实标注该近似。
        struct GirderLiveEnd {
            int elem = 0, end = 0;
            double Nmin = 0, Nmax = 0, MyMin = 0, MyMax = 0, VzMin = 0, VzMax = 0;
        };
        std::vector<GirderLiveEnd> girderLive;
        // 成桥态切线刚度上的完整模态（频率 / 周期 / 参与质量 / 节点振型）。
        // 在下面算 f1Vertical 的同一次 modalFromTangent 里就地格式化成文本，
        // 避免把 K/M 的生命周期拖到输出段；没有就是没有，输出 ERROR 段。
        std::string modalSection;
        if (converged) {
            // 车道：显式 *LANE_3D 优先；没写则自动识别主梁纵向链。
            // 识别不出时 girderLive 留空，前端按"未定义车道"提示，不静默。
            if (m.lanes.empty()) {
                Model3D::Lane3D lane = autoIdentifyGirderLane(m);
                if (!lane.nodes.empty()) { lane.autoCode = true; m.lanes.push_back(lane); }
            }
            applyJtgD60LaneLoads(m);   // 自动识别的新车道补填 qk/Pk

            if (!m.lanes.empty()) {
                // 成桥态后 m 的 L0 已是最终下料长度，再解一次收敛到同一平衡态
                // 并返回一致切线刚度（用于活载线性化）。
                NLTangent3D tangent;
                const NonlinearResult3D rLive = solveModel3DNonlinear(m, nl, &tangent);
                if (rLive.converged) {
                    const DofMap3D dmLive = buildDofMap3D(m);
                    std::vector<char> isC(tangent.K.rows(), 0);
                    for (int d : tangent.cdofs)
                        if (d >= 0 && d < static_cast<int>(isC.size())) isC[d] = 1;
                    const Factorization3D fac(tangent.K);

                    // ---- 汽车冲击系数 μ（JTG D60，竖向一阶弯曲基频）----
                    // 索承体系基频必须用成桥态一致切线刚度（张紧索切线刚度
                    // >> 弹性刚度），线弹性装配刚度会把索算软、基频算低、μ
                    // 失真。质量矩阵用材料 DENS + ADDMASS 的一致质量；竖向
                    // 一阶弯曲频率 = 主梁节点竖向位移 w 主导的最低阶。
                    {
                        const SparseMat3D Mlive = assembleMassMatrix3D(m, dmLive);
                        if (Mlive.nonZeros() > 0) {
                            const ModalResult3D modal = modalFromTangent(
                                tangent.K, Mlive, isC, dmLive.nDofs, 12);
                            if (!modal.frequencies.empty()) {
                                f1Vertical = verticalFlexuralFrequency(m, dmLive, modal);
                                if (f1Vertical < 0.0) f1Vertical = modal.frequencies[0];
                                impactMu = impactCoefficientFromFrequency(f1Vertical);
                                impactFactor = 1.0 + impactMu;
                                // 同一份 ModalResult3D：报告逐阶数值与
                                // f1_vertical / μ 必然出自同一次特征求解。
                                modalSection = formatModalSection3D(
                                    dmLive, Mlive, modal, "[MODAL]");
                            }
                        }
                    }

                    // 主梁（水平纵向梁）：活载包络只对桥面纵梁有意义。
                    std::vector<const Element3D*> tgt;
                    for (const auto& e : m.elements) {
                        if (e.type != "BEAM" && e.type != "BEAM7") continue;
                        const Vec3 d = m.nodes.at(e.n2) - m.nodes.at(e.n1);
                        const double ax = std::abs(d.x()), ay = std::abs(d.y()),
                                     az = std::abs(d.z());
                        if (ax < std::max(ay, az)) continue;
                        tgt.push_back(&e);
                    }

                    for (const auto& lane : m.lanes) {
                        const std::size_t npos = lane.nodes.size();
                        if (npos == 0) continue;
                        // 汽车冲击：默认 IMPACT=1.0 时用自动 (1+μ)；用户显式
                        // 写了非 1.0 的 IMPACT 则尊重手动值（视为额外放大）。
                        Model3D::Lane3D laneEff = lane;
                        laneEff.impact = (lane.impact != 1.0) ? lane.impact : impactFactor;
                        // 每单元 × 2 端 × 3 分量（N/My/Vz）
                        std::vector<std::vector<double>> eta(tgt.size() * 6,
                                                             std::vector<double>(npos, 0.0));
                        scanTangentLane3D(fac, dmLive, isC, lane.nodes, lane.dof,
                            [&](std::size_t k, const std::vector<double>& U) {
                                if (k >= npos) return;
                                for (std::size_t t = 0; t < tgt.size(); ++t) {
                                    const Element3D& e = *tgt[t];
                                    const bool is7 = (e.type == "BEAM7");
                                    const int dpe = is7 ? 14 : 12;
                                    const Eigen::VectorXd fl =
                                        is7 ? beam7LocalEndForces(m, dmLive, e, U)
                                            : beam3DLocalEndForces(m, dmLive, e, U);
                                    for (int end = 1; end <= 2; ++end) {
                                        const BeamEndForces bf = beamEndForcesAt(fl, dpe, end);
                                        eta[t * 6 + (end - 1) * 3 + 0][k] = bf.N;
                                        eta[t * 6 + (end - 1) * 3 + 1][k] = bf.My;
                                        eta[t * 6 + (end - 1) * 3 + 2][k] = bf.Vz;
                                    }
                                }
                            });
                        for (std::size_t t = 0; t < tgt.size(); ++t)
                            for (int end = 1; end <= 2; ++end) {
                                const auto eN  = envelope3DFromInfluence(
                                    m, laneEff, eta[t * 6 + (end - 1) * 3 + 0]);
                                const auto eMy = envelope3DFromInfluence(
                                    m, laneEff, eta[t * 6 + (end - 1) * 3 + 1]);
                                const auto eVz = envelope3DFromInfluence(
                                    m, laneEff, eta[t * 6 + (end - 1) * 3 + 2]);
                                girderLive.push_back({tgt[t]->id, end,
                                                      eN.second, eN.first,
                                                      eMy.second, eMy.first,
                                                      eVz.second, eVz.first});
                            }
                    }
                }
            }
        }

        // ---- 悬索桥施工阶段（分级找形：空缆 → 对称吊装 → 合拢 → 铺装）----
        //
        // 主缆下料长度 L0 固定（成桥态反算得到），从空缆态起逐级加载：主缆
        // 随吊装逐级下沉，矢高从 f0 增大到成桥矢高 f。加劲梁自重（一期 Q1）在
        // 吊装阶段铰接上梁、不产生加劲梁弯矩；二期铺装 Q2 在合拢后施加。
        struct ErectionRow {
            int stage = 0;
            std::string name;
            double H = 0.0, sag = 0.0, midZ = 0.0;
            double hangerTmin = 0.0, hangerTmax = 0.0;
            std::vector<double> mainX, mainZ;
            std::vector<int> chainIds;        // 主跨主缆节点 id
            std::vector<double> chainX, qLoad; // 节点 x / 该阶段竖向荷载 N
            std::vector<double> segT;         // 主跨各段轴力 N
        };
        std::vector<ErectionRow> erectionRows;
        if (susp && cfg.erection && converged) {
            // 反算成桥主缆下料长度 L0（**纯找形口径**，不用 COMPLETION 里那套
            // solveModel3DNonlinear 迭代出的 L0）。原因：施工阶段模拟用的是
            // 同一套"索多边形找形 + 闭式弹性伸长"口径，若拿非线性解出的 L0 来
            // 二分，铺装态（= 全恒载）不能精确重现成桥矢高——会系统性偏大
            // 2~3%（弹性伸长在 1000 m 主缆上就是数米的量级）。
            double L0total = 0.0;
            for (std::size_t k = 0; k + 1 < fit.chainNodes.size(); ++k) {
                const int n1 = fit.chainNodes[k], n2 = fit.chainNodes[k + 1];
                const double x1 = m.nodes.at(n1).x(), x2 = m.nodes.at(n2).x();
                auto z1 = fit.zNew.find(n1), z2 = fit.zNew.find(n2);
                if (z1 == fit.zNew.end() || z2 == fit.zNew.end()) continue;
                const double dx = x2 - x1, dz = z2->second - z1->second;
                const double chord = std::sqrt(dx * dx + dz * dz);
                if (k >= cfg.mainCable.size()) continue;
                const int id = cfg.mainCable[k];
                auto tit = fit.T.find(id);
                const Element3D* e = findElem(m, id);
                if (tit == fit.T.end() || !e) continue;
                const double T = tit->second, EAe = e->E * e->A;
                L0total += cableUnstrainedLengthForTension(chord, T, EAe, e->w, dx);
            }
            if (L0total <= 0.0) {
                std::fprintf(stderr, "[警告] 施工阶段模拟跳过：成桥态未给出主缆下料长度\n");
            } else {
                struct DeckPt { int deck; int hanger; double x; };
                std::vector<DeckPt> pts;
                const double xt1 = m.nodes.at(cfg.tower1).x();
                const double xt2 = m.nodes.at(cfg.tower2).x();
                for (int hid : cfg.hangers) {
                    auto it = fit.hangerDeck.find(hid);
                    if (it == fit.hangerDeck.end()) continue;
                    const double x = m.nodes.at(it->second).x();
                    if (x > xt1 + 1e-6 && x < xt2 - 1e-6)
                        pts.push_back({it->second, hid, x});
                }
                const double xmid = 0.5 * (xt1 + xt2);
                std::sort(pts.begin(), pts.end(), [&](const DeckPt& a, const DeckPt& b) {
                    return std::abs(a.x - xmid) < std::abs(b.x - xmid);
                });

                const double qsum = cfg.erectionQ1 + cfg.erectionQ2;
                const double q1frac = qsum > 0.0 ? cfg.erectionQ1 / qsum : 0.0;

                std::set<int> deckSet;
                for (int hid : cfg.hangers)
                    if (fit.hangerDeck.count(hid)) deckSet.insert(fit.hangerDeck.at(hid));

                auto solveStage = [&](const std::set<int>& active, double frac,
                                      const std::string& name) {
                    Model3D ms = m;
                    ms.loads.clear();
                    for (const auto& [node, dof, val] : m.loads) {
                        if (dof == 3 && deckSet.count(node))
                            ms.loads.push_back({node, dof,
                                                active.count(node) ? val * frac : 0.0});
                        else
                            ms.loads.push_back({node, dof, val});
                    }
                    const ErectionState es = findErectionState(ms, cfg, L0total);
                    ErectionRow row;
                    row.stage = static_cast<int>(erectionRows.size());
                    row.name = name;
                    row.H = es.H;
                    row.sag = es.sag;
                    row.midZ = es.midZ;
                    row.hangerTmin = 1e300;
                    row.hangerTmax = -1e300;
                    for (const auto& [hid, T] : es.hangerT) {
                        row.hangerTmin = std::min(row.hangerTmin, T);
                        row.hangerTmax = std::max(row.hangerTmax, T);
                    }
                    if (row.hangerTmin > row.hangerTmax) {
                        row.hangerTmin = row.hangerTmax = 0.0;
                    }
                    row.mainX = es.mainX;
                    row.mainZ = es.mainZ;
                    row.chainIds = es.mainIds;
                    row.chainX = es.mainX;
                    row.qLoad = es.mainQ;
                    row.segT = es.mainT;
                    erectionRows.push_back(row);
                };

                // 空缆态
                solveStage({}, 0.0, "空缆");
                // 对称吊装：从跨中向两侧，每步 cfg.erectionPanels 对（2 段）
                const int step = 2 * cfg.erectionPanels;
                for (std::size_t done = 1; done <= pts.size(); done += step) {
                    std::set<int> act;
                    for (std::size_t i = 0; i < done && i < pts.size(); ++i)
                        act.insert(pts[i].deck);
                    solveStage(act, q1frac, "吊装");
                }
                // 合拢：一期恒载全上（主跨加劲梁合拢；边跨加劲梁同步就位，全桥 Q1）
                // 注意：必须用 deckSet（全部吊索梁端节点，含边跨），不能用 pts（仅主跨）。
                // 否则边跨加劲梁的恒载被清零，铺装态少算边跨恒载、H 偏低、主跨矢高偏大，
                // 铺装态无法精确重现成桥态。
                solveStage(deckSet, q1frac, "合拢");
                // 铺装：二期恒载全上（成桥态，全桥 Q1+Q2）
                solveStage(deckSet, 1.0, "铺装");
            }
        }

        // ---- 输出 ----

        // ---- 梁端内力（成桥态）----
        // 自锚式悬索桥的全部特点都在这里：主缆水平力 H 直接压进加劲梁，
        // 梁内轴力 ≈ −H。不做这一项，自锚式就退化成"少了个锚碇的悬索桥"，
        // 而它真正要控制的（梁的轴压稳定、截面压应力）就被漏掉了。
        // 位移 u 是相对**制造线形**的，所以这里要用 erected 重建参考构形，
        // 否则算出来的梁内力会把预拱度当成变形。
        //
        // 斜拉桥同样需要：主梁在恒载下是压弯构件，索力把梁压住、桥面荷载把它
        // 弯下去，两者在截面上互相抵消。没有梁内力，索力再准也答不出"梁扛不
        // 扛得住"。早期这一整段被 `if (susp)` 挡住，斜拉桥的应力验算因此一直
        // 只报"做不了"——现在两种桥型都算。
        struct GirderEnd { int elem; int end; double N, My, Vz; };
        std::vector<GirderEnd> girderEnds;
        if (!erected.empty() && !u.empty()) {
            Model3D m2 = m;
            for (const auto& [id, p] : erected)
                if (m2.nodes.count(id)) m2.nodes[id] = p;
            const DofMap3D dm2 = buildDofMap3D(m2);
            for (const auto& e : m2.elements) {
                if (e.type != "BEAM" && e.type != "BEAM7") continue;
                const bool is7 = (e.type == "BEAM7");
                const int dpe = is7 ? 14 : 12;
                const Eigen::VectorXd fl = is7
                    ? beam7LocalEndForces(m2, dm2, e, u)
                    : beam3DLocalEndForces(m2, dm2, e, u);
                for (int end = 1; end <= 2; ++end) {
                    const BeamEndForces f = beamEndForcesAt(fl, dpe, end);
                    girderEnds.push_back({e.id, end, f.N, f.My, f.Vz});
                }
            }
        }

        // [GIRDER]：悬索桥加劲梁的一行式摘要（跨端取各分量绝对值最大者）。
        // 只用于曲线与 girder_N_max / girder_M_max 这两个摘要量。
        struct GirderForce { int elem; double N, My, Vz; };
        std::vector<GirderForce> girder;
        if (susp) {
            std::map<int, GirderForce> byElem;
            std::vector<int> order;
            for (const auto& g : girderEnds) {
                if (g.elem < cfg.deckE0 || g.elem > cfg.deckE1) continue;
                auto it = byElem.find(g.elem);
                if (it == byElem.end()) {
                    byElem[g.elem] = {g.elem, g.N, g.My, g.Vz};
                    order.push_back(g.elem);
                    continue;
                }
                if (std::abs(g.N) > std::abs(it->second.N)) it->second.N = g.N;
                if (std::abs(g.My) > std::abs(it->second.My)) it->second.My = g.My;
                if (std::abs(g.Vz) > std::abs(it->second.Vz)) it->second.Vz = g.Vz;
            }
            for (int id : order) girder.push_back(byElem[id]);
        }
        double girderNmax = 0.0, girderMmax = 0.0;
        for (const auto& g : girder) {
            girderNmax = std::max(girderNmax, std::abs(g.N));
            girderMmax = std::max(girderMmax, std::abs(g.My));
        }

        // ---- 输出 ----
        std::string outPath = dat;
        const std::size_t dot = outPath.find_last_of('.');
        const std::size_t sep = outPath.find_last_of("/\\");
        if (dot != std::string::npos && (sep == std::string::npos || dot > sep))
            outPath = outPath.substr(0, dot);
        outPath += "_cable.txt";

        std::FILE* f = std::fopen(outPath.c_str(), "w");
        if (!f) throw std::runtime_error("无法写入 " + outPath);
        std::fprintf(f, "cable_solver v1\n");
        std::fprintf(f, "mode=%s\n", cfg.mode.c_str());
        std::fprintf(f, "tension_method=%s\n", method.c_str());
        std::fprintf(f, "model nodes=%d elements=%d cables=%d dofs=%d\n",
                     static_cast<int>(m.nodes.size()),
                     static_cast<int>(m.elements.size()),
                     static_cast<int>(nc), dm.nDofs);
        std::fprintf(f, "converged=%d iterations=%d elapsed_ms=%.1f\n",
                     converged ? 1 : 0, iterations, elapsed);
        std::fprintf(f, "closure_err=%.6e tension_err=%.6e\n", closureErr, tensionErr);
        std::fprintf(f, "tension_min_kN=%.6g tension_max_kN=%.6g tension_sum_kN=%.6g\n",
                     tmin / 1e3, tmax / 1e3, sumT / 1e3);
        std::fprintf(f, "vertical_sum_kN=%.6g total_load_kN=%.6g balance_pct=%.4f\n",
                     balSum / 1e3, balTot / 1e3,
                     balTot != 0.0 ? 100.0 * (balSum - balTot) / std::abs(balTot) : 0.0);
        std::fprintf(f, "f1_vertical=%.6g impact_mu=%.6g impact_factor=%.6g\n",
                     f1Vertical, impactMu, impactFactor);
        if (rea.ok)
            std::fprintf(f, "support_v_sum_kN=%.6g support_v_max_kN=%.6g applied_v_kN=%.6g\n",
                         rea.sumV / 1e3, rea.maxV / 1e3, rea.appliedV / 1e3);
        // 桥面恒载校核：索的竖向分量应当 ≤ 桥面恒载，差额由支座承担。
        // 用 total_load（含塔自重）算比例会得到"索只承担一半"的假象——塔自重
        // 由塔底直接传走，永远不经过索。
        // 悬索桥不用这一项：它该看的是"吊索竖向力之和 = 桥面恒载"，见下面 *SUSP。
        if (cfg.deckLoad > 0.0 && !susp)
            std::fprintf(f, "deck_load_kN=%.6g deck_share_pct=%.2f support_share_kN=%.6g\n",
                         cfg.deckLoad / 1e3, 100.0 * sumV / cfg.deckLoad,
                         (cfg.deckLoad - sumV) / 1e3);

        if (susp) {
            std::fprintf(f, "susp_system=%s H_kN=%.6g sag_m=%.6g sag_target_m=%.6g\n",
                         cfg.system.c_str(), fit.H / 1e3, fit.sag, fit.sagTarget);
            std::fprintf(f, "span_m=%.6g sag_ratio=%.6f main_cable_segments=%d\n",
                         fit.span, fit.sag > 0 ? fit.span / fit.sag : 0.0,
                         static_cast<int>(cfg.mainCable.size()));
            std::fprintf(f, "hanger_count=%d hanger_T_min_kN=%.6g hanger_T_max_kN=%.6g\n",
                         static_cast<int>(cfg.hangers.size()),
                         hangerTmin / 1e3, hangerTmax / 1e3);
            std::fprintf(f, "hanger_V_sum_kN=%.6g hanger_T_sum_kN=%.6g\n",
                         hangerVsum / 1e3, hangerTsum / 1e3);
            std::fprintf(f, "girder_N_max_kN=%.6g girder_M_max_kNm=%.6g\n",
                         girderNmax / 1e3, girderMmax / 1e3);
        }

        std::fprintf(f, "[CABLE]\n");
        for (std::size_t k = 0; k < nc; ++k)
            std::fprintf(f, "%d %d %.10g %.10g %.10g %.10g %.10g\n",
                         cfg.cables[k].elem, cfg.cables[k].anchor,
                         T[k], L0[k], chord[k], Ltar[k], chord[k] - Ltar[k]);

        if (susp) {
            // 主缆分段：给定 H 与 V0 的索多边形解，T = √(H²+V²)
            std::fprintf(f, "[MAINCABLE]\n");
            for (std::size_t k = 0; k < cfg.mainCable.size(); ++k) {
                const int id = cfg.mainCable[k];
                std::size_t idx = 0;
                while (idx < nc && cfg.cables[idx].elem != id) ++idx;
                if (idx >= nc) continue;
                // 第 4 列是找形给出的该段远端主缆高程（zNew），第 5 列是
                // **平衡态弦长**。弦长必须印出来：L0 是对着平衡弦长反算的，
                // 拿设计弦长去校核 L0 会得到 1~6% 的假偏差——而这点偏差正好
                // 是把索力校核判成"不合格"的量级。
                std::fprintf(f, "%d %.10g %.10g %.10g %.10g\n", id, T[idx], L0[idx],
                             fit.zNew.count(cfg.cables[idx].anchor)
                                 ? fit.zNew[cfg.cables[idx].anchor] : 0.0,
                             chord[idx]);
            }
            // 找形后的主缆节点坐标。绘图与"设计线形"陈述只能用这一份，
            // **不能**拿 .dat 里的初始抛物线初值当设计线形：边跨主缆在锚点
            // 处是水平的（V=0），与"锚点到塔顶的直线"差得能有几米。
            std::fprintf(f, "[CABLESHAPE]\n");
            for (const auto& kv : fit.zNew) {
                const auto it = m.nodes.find(kv.first);
                if (it == m.nodes.end()) continue;
                std::fprintf(f, "%d %.10g %.10g\n", kv.first,
                             it->second.x(), kv.second);
            }
            std::fprintf(f, "[HANGER]\n");
            for (std::size_t k = 0; k < cfg.hangers.size(); ++k) {
                const int id = cfg.hangers[k];
                std::size_t idx = 0;
                while (idx < nc && cfg.cables[idx].elem != id) ++idx;
                if (idx >= nc) continue;
                std::fprintf(f, "%d %d %.10g %.10g %.10g\n", id,
                             fit.hangerDeck.count(id) ? fit.hangerDeck.at(id) : 0,
                             T[idx], L0[idx], chord[idx]);
            }
            std::fprintf(f, "[GIRDER]\n");
            for (const auto& g : girder)
                std::fprintf(f, "%d %.10g %.10g %.10g\n", g.elem, g.N, g.My, g.Vz);

            // 施工阶段（分级找形）：每阶段一行 H/矢高/跨中高程/吊索力。
            if (!erectionRows.empty()) {
                std::fprintf(f, "[ERECTION]\n");
                for (const auto& r : erectionRows)
                    std::fprintf(f, "%d %s H=%.10g sag=%.10g midz=%.10g "
                                 "hmin=%.10g hmax=%.10g\n",
                                 r.stage, r.name.c_str(), r.H, r.sag, r.midZ,
                                 r.hangerTmin, r.hangerTmax);
                // 每阶段主跨主缆线形（塔1→塔2 的高程序列），画"主缆随吊装下沉"用。
                std::fprintf(f, "[ERECTION_SHAPE]\n");
                for (const auto& r : erectionRows) {
                    std::fprintf(f, "%d", r.stage);
                    for (std::size_t k = 0; k < r.mainZ.size(); ++k)
                        std::fprintf(f, " %.10g", r.mainZ[k]);
                    std::fprintf(f, "\n");
                }
                // 逐阶段 ANSYS 独立复算所需的同源输入。这两段只在 *ERECTION
                // 启用时输出，既有算例的 _cable.txt 逐字节不变。
                //   [ERECTION_CHAIN]  <node id> <x>   主跨主缆节点链（各阶段相同）
                //   [ERECTION_LOAD]   <stage> <q...>  吊索搬运到各主缆节点的
                //                                    竖向荷载 N（向下为正）
                //   [ERECTION_FORCE]  <stage> <T...>  主跨各段轴力 N
                const auto& ids0 = erectionRows.front().chainIds;
                const auto& xs0 = erectionRows.front().chainX;
                std::fprintf(f, "[ERECTION_CHAIN]\n");
                for (std::size_t k = 0; k < ids0.size(); ++k)
                    std::fprintf(f, "%d %.10g\n", ids0[k], xs0[k]);
                std::fprintf(f, "[ERECTION_LOAD]\n");
                for (const auto& r : erectionRows) {
                    std::fprintf(f, "%d", r.stage);
                    for (double q : r.qLoad) std::fprintf(f, " %.10g", q);
                    std::fprintf(f, "\n");
                }
                std::fprintf(f, "[ERECTION_FORCE]\n");
                for (const auto& r : erectionRows) {
                    std::fprintf(f, "%d", r.stage);
                    for (double T : r.segT) std::fprintf(f, " %.10g", T);
                    std::fprintf(f, "\n");
                }
            }
        }

        // 梁端内力：每端一行，两种桥型都输出。应力验算只能用这一段——同一端
        // 才是同一截面，[GIRDER] 那种"跨端各取最大"混搭出来的 (N, My) 组合
        // 在梁上并不存在。
        if (!girderEnds.empty()) {
            std::fprintf(f, "[GIRDEREND]\n");
            for (const auto& g : girderEnds)
                std::fprintf(f, "%d %d %.10g %.10g %.10g\n",
                             g.elem, g.end, g.N, g.My, g.Vz);
        }

        // 活载梁内力包络（成桥态切线刚度影响线）。每端一行，min/max 成对。
        // 与 [GIRDEREND] 恒载标准值配对使用：荷载组合里活载项取自这里。
        if (!girderLive.empty()) {
            std::fprintf(f, "[GIRDER_LIVE]\n");
            for (const auto& g : girderLive)
                std::fprintf(f, "%d %d NMin=%.10g NMax=%.10g MyMin=%.10g "
                             "MyMax=%.10g VzMin=%.10g VzMax=%.10g\n",
                             g.elem, g.end, g.Nmin, g.Nmax,
                             g.MyMin, g.MyMax, g.VzMin, g.VzMax);
        }

        // 支座反力（平动）。计算书的"支座反力"一览直接取这里；
        // 也是自锚式"主缆水平力 H 压进梁端"的量化证据：梁端节点的 R_x 应 ≈ −H。
        if (rea.ok) {
            std::fprintf(f, "[REACTION]\n");
            for (const auto& kv : m.nodes) {
                const int nid = kv.first;
                auto it = dm.index.find(nid);
                if (it == dm.index.end()) continue;
                const int cnt = std::min(dm.dofCount[it->second], 3);
                for (int a = 0; a < cnt; ++a) {
                    const int d = dm.dofStart[it->second] + a;
                    if (std::abs(rea.R[d]) < 1.0) continue;
                    std::fprintf(f, "%d %d %.10g\n", nid, a + 1, rea.R[d]);
                }
            }
        }

        // 预抬高只对节点平动有意义：camber/u 是按自由度排的向量，靠 DofMap 取 dz
        auto dofOf = [&](int node, int a) {
            auto it = dm.index.find(node);
            if (it == dm.index.end()) return -1;
            const int cnt = dm.dofCount[it->second];
            if (a >= cnt) return -1;
            return dm.dofStart[it->second] + a;
        };
        if (!camber.empty()) {
            std::fprintf(f, "[CAMBER]\n");
            for (const auto& [id, x] : m.nodes) {
                (void)x;
                const int dx = dofOf(id, 0), dy = dofOf(id, 1), dz = dofOf(id, 2);
                std::fprintf(f, "%d %.10g %.10g %.10g\n", id,
                             dx >= 0 ? camber[dx] : 0.0,
                             dy >= 0 ? camber[dy] : 0.0,
                             dz >= 0 ? camber[dz] : 0.0);
            }
        }
        if (!erected.empty()) {
            std::fprintf(f, "[ERECTED]\n");
            for (const auto& [id, p] : erected)
                std::fprintf(f, "%d %.10g %.10g %.10g\n", id, p.x(), p.y(), p.z());
        }
        if (!u.empty()) {
            std::fprintf(f, "[U]\n");
            for (const auto& [id, x] : m.nodes) {
                (void)x;
                for (int a = 0; a < 6; ++a) {
                    const int d = dofOf(id, a);
                    if (d >= 0) std::fprintf(f, "%d %d %.10g\n", id, a + 1, u[d]);
                }
            }
        }
        // 模态段放在 [END] 前。没算出模态（未收敛 / 无车道 / 无质量）时
        // 也要写一段 ERROR=…，下游据此明确标注「本次无模态结果」。
        if (!modalSection.empty())
            std::fputs(modalSection.c_str(), f);
        else
            std::fputs(modalErrorSection3D(
                "[MODAL]", converged ? "no_lane_or_no_mass" : "not_converged").c_str(), f);
        // ---- 线性屈曲：本内核不适用，明确写出原因，不静默 -----------------
        // computeBucklingFactors3D 是「未变形构型 + 线弹性 K + 参考荷载 Kg」
        // 的一次特征值问题。这里算的是共转列式几何非线性成桥态，一致切线
        // 刚度 Kt 已经包含成桥态几何刚度；再把线性 λ 套上去会得到一个既
        // 不是"荷载放大倍数"也不是"稳定安全系数"的数。真实需求是极值点/
        // 荷载-位移全过程的非线性稳定分析，本版本未覆盖，写 ERROR 段让
        // 报告层如实标注「本次未产出」。
        std::fputs(bucklingErrorSection3D(
            "[BUCKLING]", "nonlinear_state_not_supported").c_str(), f);
        // ---- 带约束成桥索力优化（[OPT]） --------------------------------
        // 与 [LENGTH] 不是同一层：这里是"定目标索力 T*"的外层优化过程，
        // [LENGTH] 是"命中 T*"的索长/几何闭合迭代。没请求优化时明确写
        // ERROR 段，下游按「不适用」解释，不摆空表也不冒充已优化。
        if (optApplied) {
            // 成桥态实际弯矩峰值：与 [GIRDEREND] 同源，回填进 [OPT] 让下游
            // 并排看到「优化用的线性化目标」与「成桥态实际效应」。
            // 实测二者可能反向变化，报告必须写清不能互相冒充。
            double mp = 0.0;
            for (const auto& ge : girderEnds) mp = std::max(mp, std::abs(ge.My));
            optResult.momentPeakFinal = mp;
        }
        if (optRequested)
            std::fputs(formatCableOptSection3D(optResult, "[OPT]",
                                               optApplied).c_str(), f);
        else
            std::fputs(cableOptErrorSection3D(
                "[OPT]", susp && cfg.optimize
                             ? "not_applicable_suspension_formfinding"
                             : "no_optimization_requested").c_str(), f);
        // ---- 索长迭代 / 制造线形几何闭合的逐轮过程 ----------------------
        // 内核本来就把每轮误差打印到控制台，但结果文件里只有 converged /
        // iterations / closure_err / tension_err 四个汇总标量。这里把同一份
        // IterationTrace3D 写进段，Python 侧不再另算一遍。没有外层迭代的
        // 模式（STATIC、悬索桥找形）明确写 ERROR 段，不摆空表。
        if (!lengthTrace.steps.empty())
            std::fputs(formatCableLengthSection3D(lengthTrace, T, L0, chord, Ltar,
                                                  "[LENGTH]").c_str(), f);
        else
            std::fputs(cableLengthErrorSection3D(
                "[LENGTH]", susp ? "no_iteration_suspension_formfinding"
                                 : "no_iteration_static_mode").c_str(), f);
        // ---- 内层 Newton（共转列式）过程 --------------------------------
        // 与 [LENGTH] 是两层不同的迭代：这里记录给定模型上的荷载步 / 残差
        // 收敛，只暴露最终成桥态那一次平衡（外层每轮都会重跑，把几十次中间
        // 求解全写进结果文件没有意义）。旧内核无段时下游写 no_newton_section。
        if (!finalNewton.steps.empty())
            std::fputs(formatNewtonSection3D(finalNewton, "[NEWTON]").c_str(), f);
        else
            std::fputs(newtonErrorSection3D(
                "[NEWTON]", "no_newton_solution").c_str(), f);
        std::fprintf(f, "[END]\n");
        std::fclose(f);

        std::printf("cable_solver: mode=%s cables=%d converged=%d iterations=%d "
                    "%.1f ms\n", cfg.mode.c_str(), static_cast<int>(nc),
                    converged ? 1 : 0, iterations, elapsed);
        std::printf("  索力 %.4g ~ %.4g kN（Σ=%.4g kN）；竖向平衡 %.4g vs %.4g kN，偏差 %.2f%%\n",
                    tmin / 1e3, tmax / 1e3, sumT / 1e3, balSum / 1e3, balTot / 1e3,
                    balTot != 0.0 ? 100.0 * (balSum - balTot) / std::abs(balTot) : 0.0);
        std::printf("  结果: %s\n", outPath.c_str());
        return converged ? 0 : 3;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "cable_solver 失败: %s\n", e.what());
        return 2;
    }
}
