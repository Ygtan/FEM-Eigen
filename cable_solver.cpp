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
//     [REACTION] node dof R             支座反力（向上为正；锚碇为负=抗拔）
//     [CAMBER] node dx dy dz            制造线形相对设计线形的偏移（m）
//     [ERECTED] node x y z              制造线形绝对坐标
//     [U]      node dof value           成桥位移（相对制造线形）
//     [END]
//
// 编译：随 CMake 一起（目标 cable_solver）。

#include "fem_3d_nl.hpp"

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
    double minTension = 0.0;
    double minTensionFraction = 0.02;
    double deckLoad = 0.0;          // 桥面恒载合计 N（仅用于校核报告，0=未知）
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
                else if (key == "DECKLOAD") c.deckLoad = std::stod(val);
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
            converged = r.converged;
            iterations = r.iterations;
            T = r.tension;
            L0 = r.L0;
            chord = r.chord;
            Ltar = r.Ltarget;
            u = r.u;
            tensionErr = r.maxTensionErr;
            closureErr = r.maxGeomErr;
            method = "target";
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
            co.nl = nl;
            bool anyAuto = false, allAuto = true;
            for (const auto& c : cfg.cables) {
                if (c.tension > 0.0) allAuto = false; else anyAuto = true;
            }
            const CompletionStateResult3D r = solveCompletionState3D(m, cc, co);
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
            method = allAuto ? (cfg.minEnergy ? "MinBendingEnergy" : "RigidSupport")
                             : (anyAuto ? "Mixed" : "given");
        } else {
            throw std::runtime_error("未知 MODE: " + cfg.mode
                                     + "（应为 COMPLETION / TENSION / STATIC）");
        }
        elapsed = nowMs() - t0;

        std::set<int> hangerSet(cfg.hangers.begin(), cfg.hangers.end());
        double sumT = 0.0, sumV = 0.0, tmin = 1e300, tmax = -1e300;
        double hangerTmin = 1e300, hangerTmax = -1e300, hangerTsum = 0.0, hangerVsum = 0.0;
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

        // ---- 输出 ----

        // ---- 加劲梁内力（成桥态）----
        // 自锚式悬索桥的全部特点都在这里：主缆水平力 H 直接压进加劲梁，
        // 梁内轴力 ≈ −H。不做这一项，自锚式就退化成"少了个锚碇的悬索桥"，
        // 而它真正要控制的（梁的轴压稳定、截面压应力）就被漏掉了。
        // 位移 u 是相对**制造线形**的，所以这里要用 erected 重建参考构形，
        // 否则算出来的梁内力会把预拱度当成变形。
        struct GirderForce { int elem; double N, My, Vz; };
        std::vector<GirderForce> girder;
        if (susp && !erected.empty() && !u.empty()) {
            Model3D m2 = m;
            for (const auto& [id, p] : erected)
                if (m2.nodes.count(id)) m2.nodes[id] = p;
            const DofMap3D dm2 = buildDofMap3D(m2);
            for (const auto& e : m2.elements) {
                if (e.type != "BEAM" && e.type != "BEAM7") continue;
                if (e.id < cfg.deckE0 || e.id > cfg.deckE1) continue;
                const Eigen::VectorXd fl = (e.type == "BEAM7")
                    ? beam7LocalEndForces(m2, dm2, e, u)
                    : beam3DLocalEndForces(m2, dm2, e, u);
                const BeamEndForces fi = beamEndForcesAt(fl, 12, 1);
                const BeamEndForces fj = beamEndForcesAt(fl, 12, 2);
                // 取两端中 |N| 较大、|My| 较大的组合，一行报完一根梁
                const double N = std::abs(fi.N) >= std::abs(fj.N) ? fi.N : fj.N;
                const double My = std::abs(fi.My) >= std::abs(fj.My) ? fi.My : fj.My;
                const double Vz = std::abs(fi.Vz) >= std::abs(fj.Vz) ? fi.Vz : fj.Vz;
                girder.push_back({e.id, N, My, Vz});
            }
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
