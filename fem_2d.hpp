#pragma once

#include "fem.hpp"
#include "section.hpp"
#include "code_loads.hpp"
#include "prestress_loss.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// 二维结构层（fem_2d）：输入文件读取与组装。
//
// 输入文件约定（关键字风格，大小写不敏感，`#`/`!` 开头为注释）：
//
//   *TITLE
//   <标题文本>
//
//   *NODE
//   <节点号> <x> <y>
//
//   *ELEMENT, TYPE=TRUSS, E=<值>, A=<值>
//   <单元号> <节点i> <节点j>
//
//   *ELEMENT, TYPE=BEAM, E=<值>, A=<值>, I=<值>
//   <单元号> <节点i> <节点j>
//
//   *BOUNDARY
//   <节点号> <dof> <值>     # dof: 1=u(x), 2=v(y), 3=转角θ
//
//   *CLOAD
//   <节点号> <dof> <值>
//
//   *END
//
// 二维结构统一按“每节点 3 个自由度 [u, v, θ]”编号；桁架只贡献 u/v，
// 梁贡献 u/v/θ，因此可以混合使用。
namespace fem {

// 由纤维截面构造 2D 梁单元（换算截面：E=Eref，A=A*，I=Iz*）。
// 组合截面时 Eref·A* = Σ E_i·A_i、Eref·Iz* = Σ E_i·(y_i-cy)²·A_i，即真实轴向/抗弯刚度。
inline EulerBeam2D beamFromSection2D(const Section& sec,
                                     const Eigen::Vector2d& x0, const Eigen::Vector2d& x1) {
    const SectionProperties p = computeTransformedProperties(sec);
    const double E = (sec.Eref > 0.0) ? sec.Eref : 1.0;
    return EulerBeam2D{E, p.A, p.Iz, x0, x1};
}

struct Model2D {
    std::string title;
    std::map<int, Eigen::Vector2d> nodes; // 节点号 -> 坐标

    struct Element {
        int id = 0;
        std::string type;   // "TRUSS" 或 "BEAM"
        double E = 0.0;
        double A = 0.0;
        double I = 0.0;     // 仅梁使用
        double alpha = 0.0; // 热膨胀系数（温度荷载）
        double h = 0.0;     // 梁截面高度（温度梯度）
        double Wt = 0.0;    // 梁上翼缘抗弯模量（应力）
        double Wb = 0.0;    // 梁下翼缘抗弯模量（应力）
        double dens = 0.0;  // 材料密度（自动自重）
        double addmass = 0.0; // 附加质量 kg/m（自动自重）
        double istrn = 0.0; // 初应变（cable/strut/truss 预应力）
        double nu = 0.0;    // 泊松比（>0 时梁按 Timoshenko 考虑剪切变形）
        double K = 5.0 / 6.0; // 剪切修正系数（矩形默认 5/6）
        std::vector<int> conn; // 两个节点号
        // 端部释放：被释放的局部自由度索引（0=u_i,1=v_i,2=θ_i,3=u_j,4=v_j,5=θ_j）。
        // 仅对 BEAM 有意义；M 释放=铰（θ 释放），V=剪力释放，N=轴力释放。
        std::vector<int> releases;
    };
    std::vector<Element> elements;

    struct Constraint {
        int node;
        int dof;    // 1/2/3
        double value;
        double k = 0.0;   // 弹性支承刚度（>0 为弹簧，0 为刚性约束）
    };
    std::vector<Constraint> constraints;

    struct Load {
        int node;
        int dof;    // 1/2/3
        double value;
    };
    std::vector<Load> loads;

    // ---- 新增荷载类型 ----
    struct ElemLoad {         // 单元均布/梯形荷载
        int elem = 0;
        std::string dir;      // GX/GY/LX/LY 或 G<角度>/L<角度>
        double q1 = 0.0, q2 = 0.0;
    };
    std::vector<ElemLoad> elemLoads;

    double gx = 0.0, gy = -9.81; // 重力加速度（全局，+Y 向上，默认向下）

    struct SelfWeight {
        bool all = false;
        int elem = 0;
        std::string rule; // "RHO"(面积×密度×g) 或 "MASS"(线质量×g)
        double value = 0.0;
    };
    std::vector<SelfWeight> selfWeights;

    double uniformTemp = 0.0;  // 体系温度升降（0 表示未设置）

    struct ElemTemp {
        int elem = 0;
        std::string kind; // "UNIFORM"(整体升温) 或 "GRADIENT"(温度梯度)
        double dT = 0.0;
    };
    std::vector<ElemTemp> elemTemps;

    struct Prestress {
        int id = 0;
        std::vector<int> elems;     // 施加预应力的梁单元号（按 i→j 顺序连接）
        double P = 0.0;             // 张拉控制力（N）
        double F = 0.0;             // 跨中相对两端连线的下垂量（m，向下为正）
        double E1 = 0.0;            // 起点偏心距（m，向下为正）
        double E2 = 0.0;            // 终点偏心距（m，向下为正）
        double loss = 1.0;          // 总损失系数，有效预应力 = P*loss
        int lossRef = 0;            // 引用 *PRESTRESS_LOSS 的 id；>0 时用 σpe/σcon 代替 loss
        double Ep = 195e9;          // 预应力筋弹性模量（Pa）
        double Ap = 0.0;            // 预应力筋面积（m²），TENDON 方式必填
        std::string method = "AUTO"; // AUTO/EQUIV/TENDON
    };
    std::vector<Prestress> prestresses;
    std::map<int, PrestressLossInput> prestressLosses;   // *PRESTRESS_LOSS 表

    struct InfluenceRequest {
        int id = 0;
        std::string kind;      // DISP / REACTION / FORCE / STRESS
        int node = 0;          // DISP/REACTION 用
        int dof = 0;           // DISP/REACTION 用
        int elem = 0;          // FORCE/STRESS 用
        int end = 1;           // 1 或 2
        std::string component; // FORCE: N/V/M；STRESS: AXIAL/TOP/BOT
    };
    std::vector<InfluenceRequest> influences;

    struct LaneLoad {
        int id = 0;
        std::vector<int> nodes; // 车道纵向节点，按顺序排列
        double q = 0.0;         // 均布荷载强度，沿 DIR 方向，可为负
        double P = 0.0;         // 集中荷载，沿 DIR 方向，可为负
        std::string dir = "GY"; // GX 或 GY
        std::string method = "LINEAR"; // LINEAR=影响线叠加；NONLINEAR=直接非线性扫描
        bool autoCode = false;  // CODE=JTG_D60：按跨径自动生成车道荷载
        // 1 表示按竖向基频自动取 1+μ。写成别的数则当作放大倍数本身。
        double impact = 1.0;
    };
    std::vector<LaneLoad> laneLoads;

    struct ConstraintEquationTerm {
        int node = 0;
        int dof = 0;
        double coef = 0.0;
    };

    struct ConstraintEquation {
        int id = 0;
        double constant = 0.0;
        std::vector<ConstraintEquationTerm> terms;
    };
    std::vector<ConstraintEquation> ces;
    double cePenalty = 1e12; // CE 罚系数

    struct CreepSpec {
        bool all = false;
        std::vector<int> elems;
        double phi = 0.0;
        double chi = 0.8;
    };
    std::vector<CreepSpec> creepSpecs;

    // 徐变折减基准：单元 id -> **未折减**的初始弹性模量。
    // 解析完成后由 snapshotBaseModulus() 自动填充。
    // AEMM 的每次折减都从这张表重新计算（而不是就地 e.E *= factor），
    // 因此处理徐变是幂等的、可重入的。
    std::map<int, double> baseModulus;

    struct ShrinkageSpec {
        bool all = false;
        std::vector<int> elems;
        double ecs = 0.0; // 收缩应变，正值为缩短
    };
    std::vector<ShrinkageSpec> shrinkageSpecs;

    struct TempProfileSpec {
        bool all = false;
        std::vector<int> elems;
        double T1 = 0.0; // 梁顶温度
        double T2 = 0.0; // 距梁顶 A 处温度
        double A = 0.10; // T2 控制点距梁顶距离
        double B = 0.0;  // 温度降为 0 的深度，0 表示取梁高
        double H = 0.0;  // 梁高，0 表示取实常数 H
    };
    std::vector<TempProfileSpec> tempProfiles;

    struct Stage {
        int id = 0;
        double time = 0.0; // 施工阶段龄期（天），用于时变曲线
        std::vector<int> activate;
        std::vector<int> deactivate;
        std::vector<Load> loads;
        std::vector<ElemLoad> elemLoads;
        std::vector<Prestress> prestresses;
        std::vector<Constraint> addedConstraints;
        std::vector<std::pair<int, int>> removedConstraints;
        std::vector<CreepSpec> creepSpecs;
        std::vector<ShrinkageSpec> shrinkageSpecs;
        std::vector<TempProfileSpec> tempProfiles;
        std::vector<std::pair<int, int>> unloadNodeDofs;
        std::vector<int> unloadElems;
    };
    std::vector<Stage> stages;

    struct TimeCurve {
        std::vector<std::pair<double, double>> points;
    };
    TimeCurve creepCurve;
    TimeCurve shrinkageCurve;

    std::string creepModel = "CURVE";      // CURVE / CEB90
    double creepPhi0 = 2.0;
    double creepT0 = 7.0;
    double creepBetaH = 700.0;
    double creepPower = 0.3;
    // 三者都大于 0 时按 JTG 3362 附录 C.2.1 取名义徐变系数和 βH。
    double creepRh = 0.0;
    double creepH = 0.0;
    double creepFcm = 0.0;

    std::string shrinkageModel = "CURVE";  // CURVE / CEB90
    double shrinkageEcsInf = 3.0e-4;
    double shrinkageTs = 7.0;
    double shrinkageAlphaDs = 35.0;
    double shrinkagePower = 0.5;

    struct StageTarget {
        int node = 0;
        int dof = 0;
        double value = 0.0;
        double weight = 1.0;
        std::string kind = "DISP"; // DISP / REACTION / FORCE / STRESS
        int elem = 0;
        int end = 1;
        std::string component;
    };
    std::vector<StageTarget> stageTargets;

    struct StageControl {
        int stage = 0;          // 阶段号
        int prestressIndex = 0; // 该阶段 *STAGE_PRESTRESS 的序号，1-based
        double lower = 0.0;     // 缩放系数下限
        double upper = 1e9;     // 缩放系数上限
    };
    std::vector<StageControl> stageControls;

    struct StageConstraint {
        StageTarget target;
        double lower = -1e30;
        double upper = 1e30;
    };
    std::vector<StageConstraint> stageConstraints;

    struct CamberSpec {
        int node = 0;
        int dof = 0;
        double value = 0.0; // 目标最终位移，默认 0
    };
    std::vector<CamberSpec> cambers;

    struct StageFeedback {
        int node = 0;
        int dof = 0;
        double measured = 0.0;
    };
    std::vector<StageFeedback> feedbacks;

    // ---- 新格式：单元类型 / 材料 / 实常数表（声明一次，多处复用）----
    std::map<int, std::string> elementTypes; // 类型号 -> "BEAM"/"TRUSS"

    struct Material {
        double E = 0.0;      // 弹性模量
        double alpha = 0.0;  // 热膨胀系数
        double dens = 0.0;   // 密度（自重备用）
        double nu = 0.0;     // 泊松比（Timoshenko 剪切用）
    };
    std::map<int, Material> materials;

    struct RealConst {
        double A = 0.0;   // 截面积
        double I = 0.0;   // 惯性矩（梁）
        double Wt = 0.0;  // 上翼缘抗弯模量（梁应力）
        double Wb = 0.0;  // 下翼缘抗弯模量（梁应力）
        double h = 0.0;   // 截面高度（梁温度梯度）
        double addmass = 0.0; // 附加质量 kg/m
        double istrn = 0.0;   // 初应变（cable/strut 预应力）
        double K = 5.0 / 6.0; // 剪切修正系数
    };
    std::map<int, RealConst> realConsts;

    std::string analysisType = "LINEAR"; // "LINEAR"(方案A) 或 "NONLINEAR"(方案C)
    double nonlinearTol = 1e-6;          // 方案C Newton 相对收敛容差
    int nonlinearMaxIter = 50;           // 方案C Newton 最大迭代次数
    // 荷载分级步数。1 = 不分級（与旧行为一致）；大位移/索松弛工况建议 4~10。
    // 某一步不收敛时会自动细分步长重来。
    int nonlinearLoadSteps = 1;
    // 是否启用 dU 阻尼线搜索（残差不下降时把步长折半）。
    // 正常情况下 α=1 直接被接受，结果与不启用时完全一致。
    bool nonlinearLineSearch = true;
    // 不收敛时是否抛异常。设为 false 则返回"部分成功"的结果，
    // 由 ModelResult::nonlinearConverged 告知调用方未收敛。
    bool nonlinearThrowOnFailure = true;
};

// 2D 模型的统一自由度映射：所有节点固定 3 个自由度 [u,v,θ]，
// 与 solveModel2D 内部“3*idx + (dof-1)”的编号完全一致。
inline DofMap buildDofMap2D(const Model2D& m) {
    DofMap dm;
    for (const auto& [id, x] : m.nodes) {
        dm.index[id] = static_cast<int>(dm.nodeIds.size());
        dm.nodeIds.push_back(id);
        dm.dofStart.push_back(3 * static_cast<int>(dm.nodeIds.size() - 1));
        dm.dofCount.push_back(3);
    }
    dm.nDofs = 3 * static_cast<int>(m.nodes.size());
    return dm;
}

// JTG D60 自动车道荷载（二维）：CODE=JTG_D60 时按跨径填入 Q、P（负号 = 向下）。
inline void applyJtgD60LaneLoads2D(Model2D& m) {
    for (auto& lane : m.laneLoads) {
        if (!lane.autoCode || lane.nodes.size() < 2) continue;
        const Eigen::Vector2d a = m.nodes.at(lane.nodes.front());
        const Eigen::Vector2d b = m.nodes.at(lane.nodes.back());
        const double L = (b - a).norm();
        lane.q = -jtgD60LaneLoadQk();
        lane.P = -jtgD60LaneLoadPk(L);
    }
}

namespace detail {

inline std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

inline std::string toUpper(std::string s) {
    for (char& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    return s;
}

inline std::map<std::string, std::string> parseParams(const std::string& s) {
    std::map<std::string, std::string> p;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim(tok);
        const std::size_t eq = tok.find('=');
        if (eq == std::string::npos) continue;
        p[toUpper(trim(tok.substr(0, eq)))] = trim(tok.substr(eq + 1));
    }
    return p;
}

inline double toDouble(const std::map<std::string, std::string>& p, const std::string& key) {
    auto it = p.find(toUpper(key));
    if (it == p.end())
        throw std::runtime_error("输入文件缺少参数: " + key);
    return std::stod(it->second);
}

inline std::string toLower(std::string s) {
    for (char& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

// 解析空格分隔的 key=value（如 "E=2e11 ALPX=1.2e-5"）。
inline std::map<std::string, std::string> parseKeyValues(const std::string& s) {
    std::map<std::string, std::string> p;
    std::istringstream ss(s);
    std::string tok;
    while (ss >> tok) {
        const std::size_t eq = tok.find('=');
        if (eq == std::string::npos) continue;
        p[toUpper(trim(tok.substr(0, eq)))] = trim(tok.substr(eq + 1));
    }
    return p;
}

// 解析 "5" 或 "2-5" 为整数列表。
inline std::vector<int> parseRange(const std::string& s) {
    std::vector<int> ids;
    const std::size_t dash = s.find('-');
    if (dash != std::string::npos && dash > 0 && dash + 1 < s.size()) {
        const int a = std::stoi(s.substr(0, dash));
        const int b = std::stoi(s.substr(dash + 1));
        if (b < a) throw std::runtime_error("区间起始大于结束: " + s);
        for (int v = a; v <= b; ++v) ids.push_back(v);
    } else {
        ids.push_back(std::stoi(s));
    }
    return ids;
}

// 解析 "3"、"1-5"、"1-5,8,10-12" 为整数列表。
inline std::vector<int> parseIntList(const std::string& s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const std::string t = trim(tok);
        if (t.empty()) continue;
        const auto part = parseRange(t);
        out.insert(out.end(), part.begin(), part.end());
    }
    return out;
}

inline bool isIdentChar(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

// 把一行里的循环变量替换成数值，支持 VAR、VAR±m、k*VAR、k*VAR±m（大小写不敏感）。
inline std::string substituteLoopVar(const std::string& line, const std::string& var, int value) {
    const std::string low = toLower(line);
    const std::string vlow = toLower(var);
    std::string out;
    for (std::size_t i = 0; i < line.size(); ) {
        bool replaced = false;
        if (i + vlow.size() <= line.size() &&
            low.compare(i, vlow.size(), vlow) == 0) {
            const bool varLeftOk = (i == 0) || !isIdentChar(line[i - 1]);
            const bool varRightOk = (i + vlow.size() == line.size()) || !isIdentChar(line[i + vlow.size()]);
            if (varLeftOk && varRightOk) {
                // 向前找可选系数 k*（正整数）。
                int coeff = 1;
                std::size_t start = i;
                if (i >= 2 && line[i - 1] == '*') {
                    std::size_t d = i - 1; // d 指向 '*'
                    while (d > 0 && std::isdigit(static_cast<unsigned char>(line[d - 1]))) --d;
                    if (d < i - 1) { // 存在数字
                        if (d == 0 || !isIdentChar(line[d - 1])) {
                            int num = 0;
                            for (std::size_t q = d; q < i - 1; ++q)
                                num = num * 10 + (line[q] - '0');
                            coeff = num;
                            start = d;
                        }
                    }
                }
                // 向后找可选偏移 ±m。
                int offset = 0;
                std::size_t j = i + vlow.size();
                if (j < line.size() && (line[j] == '+' || line[j] == '-')) {
                    const char sign = line[j];
                    std::size_t q = j + 1;
                    if (q < line.size() && std::isdigit(static_cast<unsigned char>(line[q]))) {
                        int num = 0;
                        while (q < line.size() && std::isdigit(static_cast<unsigned char>(line[q]))) {
                            num = num * 10 + (line[q] - '0');
                            ++q;
                        }
                        offset = (sign == '+') ? num : -num;
                        j = q;
                    }
                }
                // 系数 k* 位于变量之前，且已被逐字符输出，需先撤销。
                if (start < i) out.erase(out.size() - (i - start));
                out += std::to_string(coeff * value + offset);
                i = j;
                replaced = true;
            }
        }
        if (!replaced) {
            out += line[i];
            ++i;
        }
    }
    return out;
}

// 递归展开 *FILL 与 *DO...*ENDDO。
inline void expandLines(const std::vector<std::string>& lines, std::size_t& idx,
                        std::size_t end, std::vector<std::string>& out,
                        std::map<int, std::pair<double, double>>& nodeCoord,
                        bool& inNode) {
    while (idx < end) {
        const std::string raw = lines[idx];
        const std::string t = trim(raw);
        const std::string up = toUpper(t);

        if (up.rfind("*DO", 0) == 0 &&
            (t.size() == 3 || std::isspace(static_cast<unsigned char>(t[3])))) {
            std::istringstream ss(t);
            std::string kw, var;
            int start = 0, endv = 0, step = 1;
            ss >> kw >> var >> start >> endv;
            if (!(ss >> step)) step = 1;
            if (step == 0) throw std::runtime_error("*DO 步长不能为 0");

            std::size_t bodyStart = idx + 1;
            std::size_t bodyEnd = bodyStart;
            int depth = 1;
            while (bodyEnd < end) {
                const std::string bu = toUpper(trim(lines[bodyEnd]));
                if (bu.rfind("*DO", 0) == 0) ++depth;
                else if (bu == "*ENDDO") --depth;
                if (depth == 0) break;
                ++bodyEnd;
            }
            if (depth != 0) throw std::runtime_error("*DO 缺少匹配的 *ENDDO");

            for (int v = start; (step > 0 ? v <= endv : v >= endv); v += step) {
                for (std::size_t k = bodyStart; k < bodyEnd; ++k)
                    out.push_back(substituteLoopVar(lines[k], var, v));
            }
            idx = bodyEnd + 1; // 跳过 *ENDDO
        } else if (up.rfind("*FILL", 0) == 0 &&
                   (t.size() == 5 || std::isspace(static_cast<unsigned char>(t[5])))) {
            std::istringstream ss(t);
            std::string kw;
            int a = 0, b = 0;
            ss >> kw >> a >> b;
            auto ita = nodeCoord.find(a);
            auto itb = nodeCoord.find(b);
            if (ita == nodeCoord.end() || itb == nodeCoord.end())
                throw std::runtime_error("*FILL 端点节点不存在");
            if (b <= a) throw std::runtime_error("*FILL 要求第二个节点号大于第一个");
            const double x1 = ita->second.first, y1 = ita->second.second;
            const double x2 = itb->second.first, y2 = itb->second.second;
            for (int id = a + 1; id < b; ++id) {
                const double f = static_cast<double>(id - a) / static_cast<double>(b - a);
                std::ostringstream oss;
                oss << std::setprecision(15) << id << ' ' << (x1 + f * (x2 - x1)) << ' ' << (y1 + f * (y2 - y1));
                out.push_back(oss.str());
            }
            ++idx;
        } else {
            if (!t.empty() && t[0] == '*') {
                inNode = (up == "*NODE");
            } else if (inNode && !t.empty() && t[0] != '#' && t[0] != '!') {
                std::istringstream ss(t);
                int id = 0;
                double x = 0.0, y = 0.0;
                if (ss >> id >> x >> y) nodeCoord[id] = {x, y};
            }
            out.push_back(raw);
            ++idx;
        }
    }
}

// 预处理：把 *FILL / *DO...*ENDDO 展开成普通行。
inline std::vector<std::string> preprocess(std::istream& in) {
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);
    std::vector<std::string> out;
    std::map<int, std::pair<double, double>> nodeCoord;
    bool inNode = false;
    std::size_t idx = 0;
    expandLines(lines, idx, lines.size(), out, nodeCoord, inNode);
    return out;
}

} // namespace detail

// ==================== 荷载等效节点荷载 ====================

inline constexpr double kPi = 3.14159265358979323846;

// 判断是否为“只含轴向刚度”的单元：普通桁架、只受拉索、只受压杆。
bool isAxial(const std::string& type) ;

bool isRigid(const std::string& type) ;

// 销接刚臂：只约束梁节点与筋节点之间的平动刚体关系，不向筋节点传递弯矩。
// 用于无粘结预应力筋，避免在转向点产生虚假附加弯矩。
DenseMat rigidPinStiffness(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                  double penalty = 1e11) ;

// 计算轴向单元的机械轴力（受拉为正）。
// 约定：N = E·A/L·ΔL + E·A·ISTRN。
double axialForce(const Model2D::Element& e, const DenseVec& U,
                         int i1, int i2,
                         const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 二维单元全局自由度编号（桁架 4 个，梁 6 个）。
std::vector<int> elementGdofs2D(const std::string& type, int i1, int i2) ;

// 把方向字段解析成全局单位方向向量。
// GX/GY：全局轴；LX/LY：局部轴/横向；G<θ>：与全局X成θ角；L<θ>：与杆轴(i→j)成θ角。
Eigen::Vector2d resolveDir(const std::string& dir,
                                  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 分布荷载等效节点荷载（全局坐标）。q1/q2 为 i/j 端沿方向的强度。
DenseVec distributedLoadEquivalent(const std::string& type,
                                          const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                          const Eigen::Vector2d& gdir,
                                          double q1, double q2) ;

// 梁单元局部横向分布荷载（向下为正），取 i/j 两端平均值；用于弯矩图的抛物线修正。
// 仅含分布荷载与自重（不含温度梯度——温度梯度产生的是常弯矩，不产生抛物线项）。
double localTransverseUDL(const Model2D& m, const Model2D::Element& e,
                          const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 梁剪切参数 Φ = 12EI/(G·As·L²) = 24(1+ν)·I/(K·A·L²)。
// ν<=0 或 K<=0 或 A<=0 时返回 0（退化为欧拉-伯努利梁）。
double beamShearParameter(const Model2D::Element& e, double L) ;

// 梁端部释放：静态凝聚。输入局部刚度 Kl(6x6)、局部等效荷载 Fq(6)、释放的局部自由度。
// 释放的局部自由度内力被强制为 0（M 释放=铰、V=剪力释放、N=轴力释放）。
struct BeamRelease2D {
    DenseMat Kl;                 // 凝聚后局部刚度
    std::vector<int> freeLocal;  // 保留的局部自由度索引
    DenseVec Fq;                 // 凝聚后局部等效荷载
    DenseMat KrrInvKrf;          // K_rr^{-1} K_rf（内力反算用）
    DenseMat KrrInv;             // K_rr^{-1}
    std::vector<int> relLocal;   // 释放的局部自由度索引
};
BeamRelease2D releaseBeam2D(const DenseMat& Kl, const DenseVec& Fq,
                            const std::vector<int>& releasedLocal);

// 均匀温度变化等效节点荷载（轴力）。
DenseVec uniformTempEquivalent(const std::string& type, double E, double A, double alpha,
                                      double dT, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 梁温度梯度等效节点荷载（弯曲）。
DenseVec gradientTempEquivalent(double E, double I, double alpha, double h, double dT,
                                       const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 规范竖向温度曲线 -> 截面温度轴力/弯矩。
std::pair<double, double> tempProfileResultants(double E, double A, double H,
                                                       double alpha,
                                                       double T1, double T2,
                                                       double A0, double B) ;

DenseVec tempProfileEquivalent(const Model2D::Element& e,
                                      const Eigen::Vector2d& p1,
                                      const Eigen::Vector2d& p2,
                                      double N, double M) ;

// 把单元等效荷载散入全局 F。
void addElementLoad(const std::vector<int>& gdofs, const DenseVec& fe, DenseVec& F) ;

// 单个单元上所有分布/自重/温度荷载的等效节点荷载（全局坐标，按单元 gdofs 顺序）。
DenseVec elementEquivalentLoad(const Model2D& m, const Model2D::Element& e,
                                      const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// ==================== 方案 C：几何非线性轴向单元 ====================

// 当前构形下的机械轴力（受拉为正），工程应变 + 初应变。
double axialForceNonlinear(const Model2D::Element& e, const DenseVec& U,
                                  int i1, int i2,
                                  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 轴向单元当前构形下的内部节点力向量：f_int = [-N c, +N c]。
DenseVec axialInternalForceNonlinear(const Model2D::Element& e, const DenseVec& U,
                                             int i1, int i2,
                                             const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 轴向单元一致切线刚度矩阵（对应上述 f_int）。
DenseMat axialTangentNonlinear(const Model2D::Element& e, const DenseVec& U,
                                      int i1, int i2,
                                      const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 方案 C 外部等效节点荷载：与线性路径一致，但不包含初应变 ISTRN。
// ISTRN 已计入轴向单元内部力。
DenseVec elementExternalEquivalentLoad(const Model2D& m, const Model2D::Element& e,
                                              const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) ;

// 处理 *PRESTRESS：等效荷载法或自动生成预应力筋单元。
void processPrestress(Model2D& m) ;

// 记录徐变折减基准（未折减的初始弹性模量）。
// 必须在所有会新增单元的步骤（如 processPrestress 生成预应力筋）之后调用，
// 这样新生成的钢束单元也有正确的基准值。
void snapshotBaseModulus(Model2D& m) ;

// 取单元的徐变折减基准模量；没有记录时退化为单元自身的 E。
double baseModulusOf(const Model2D& m, const Model2D::Element& e) ;

// AEMM（龄期调整有效模量法）徐变折减系数。
//   Eeff = E0 / (1 + chi * Δφ)，  Δφ = φ(t_now) - φ(t_load) >= 0
double creepReductionFactor(double chi, double dphi) ;

double interpolateTimeCurve(const Model2D::TimeCurve& curve, double time) ;

double creepCoefficientAt(const Model2D& m, double time) ;

double shrinkageStrainAt(const Model2D& m, double time) ;

// 处理混凝土徐变的模量折减（AEMM 龄期调整有效模量法）。
//
// 与旧实现的三点关键区别：
//   1) **可重入**：不再 `e.E *= factor` 就地累乘，而是每次都从 baseModulus
//      （未折减基准）重新算一遍。重复调用同一阶段不会把模量折减两次。
//   2) **多 spec 不再连乘**：同一单元命中多条 *CREEP 时，取折减最强的一条
//      （max chi*Δφ），语义明确且与书写顺序无关。
//   3) **逐阶段时效**：给出 now > 0 与 loadAge 时，按
//      Δφ = φ(now) - φ(t_load) 折减；不同单元因加载龄期不同而得到不同的有效
//      模量，这才是 AEMM 的核心特征。不给 now 时退化为旧的"直接给 phi"。
//
// 注意：本函数只改模量。收缩按等效温度应变在单元荷载向量里处理，不在这里。
void processTimeEffects(Model2D& m,
                               const std::map<int, double>* loadAge = nullptr,
                               double now = 0.0) ;

// 读取二维结构输入文件。
// 支持两种风格：
//   1) 新风格（表声明 + *EN 命令 + *DO/*FILL 预处理）；
//   2) 旧风格（*ELEMENT, TYPE=..., E=..., A=... 内联模板），保持兼容。
Model2D readModel2D(const std::string& path) ;

struct ModelResult {
    DenseVec U;                        // 位移（按 [u,v,θ] 每节点排列）
    DenseVec R;                        // 反力 = K0*U - F0
    std::vector<int> nodeIds;          // 0-based 索引 -> 原节点号
    std::vector<Eigen::Vector2d> coords; // 对应坐标
    int nonlinearIterations = 0;       // 方案C Newton 迭代次数（线性分析为 0）
    double nonlinearResidual = 0.0;    // 方案C 收敛时自由残差范数（线性分析为 0）
    bool nonlinearConverged = true;    // 方案C 是否收敛（关闭 throwOnFailure 时可能为 false）
};

// ==================== 方案C Newton 的诊断信息 ====================
//
// 旧的失败信息只有一句"未收敛（迭代 50 次，最后残差=xxx）"，用户既不知道
// 卡在哪个自由度、也不知道残差是发散还是在震荡、更拿不到部分成功的结果。
// 这里把迭代过程的量记下来，失败时拼成可定位的诊断文本。
struct NonlinearDiag {
    int iterations = 0;              // 累计 Newton 迭代次数
    int loadSteps = 0;               // 已完成的荷载分级步数
    int lineSearchCuts = 0;          // 线搜索把步长折半的次数
    int stepRefinements = 0;         // 荷载步因不收敛被细分的次数
    double loadRatio = 0.0;          // 失败时已成功施加的荷载比例
    double initialResidual = 0.0;
    double lastResidual = 0.0;
    double residualRef = 0.0;        // 收敛判据的参考量 tol*max(1,|Fext|)
    std::vector<double> residualHistory;
    int maxResidualNode = -1;        // 最后一步残差最大的节点号
    int maxResidualDof = -1;         // 1..3
    double maxResidualValue = 0.0;
    // 条件数的廉价代理：max|Kii| / min|Kii|（自由自由度）。
    // 真正的条件数估计要做极值特征值求解，代价太高，这里只用它做量级提示。
    double diagCondEst = 0.0;
    bool converged = false;
    std::string failureReason;

    std::string describe() const {
        std::ostringstream os;
        os << "方案C：Newton 迭代未收敛\n";
        if (!failureReason.empty()) os << "  失败原因      : " << failureReason << "\n";
        os << "  已施加荷载比例: " << std::fixed << std::setprecision(1)
           << loadRatio * 100.0 << "%（分 " << loadSteps << " 步完成）\n";
        os << "  累计迭代次数  : " << iterations << "\n";
        os << "  残差          : 最后=" << std::scientific << std::setprecision(3)
           << lastResidual << "，初始=" << initialResidual
           << "，判据阈值=" << residualRef << "\n";
        if (maxResidualNode >= 0)
            os << "  残差最大自由度: 节点 " << maxResidualNode
               << " dof " << maxResidualDof
               << "（|r|=" << std::scientific << std::setprecision(3)
               << maxResidualValue << "）\n";
        if (diagCondEst > 0.0)
            os << "  对角条件数估计: " << std::scientific << std::setprecision(3)
               << diagCondEst << "（max|Kii|/min|Kii|，仅作量级参考）\n";
        os << "  线搜索折减    : " << lineSearchCuts
           << " 次；荷载步细分 " << stepRefinements << " 次\n";
        if (!residualHistory.empty()) {
            os << "  残差历史（最后 " << std::min<std::size_t>(10, residualHistory.size())
               << " 步）: ";
            const std::size_t from = residualHistory.size() > 10
                                         ? residualHistory.size() - 10 : 0;
            for (std::size_t i = from; i < residualHistory.size(); ++i)
                os << std::scientific << std::setprecision(2)
                   << residualHistory[i] << (i + 1 < residualHistory.size() ? " " : "");
            os << "\n";
        }
        return os.str();
    }
};

// ==================== CE 约束方程（罚函数法） ====================

void appendConstraintPenaltyTriplets(const Model2D& m,
                                            const std::map<int, int>& index,
                                            int nDofs,
                                            double penalty,
                                            std::vector<Eigen::Triplet<double>>& triplets) ;

// 若用户未显式指定罚系数（<=0），按单元刚度最大对角元的 1e8 倍自适应取值，
// 避免固定 1e12 在细长杆/小刚度模型中条件数恶化。
double effectiveCePenalty(const std::vector<Eigen::Triplet<double>>& triplets,
                                 double userPenalty) ;

DenseVec constraintPenaltyConstantForce(const Model2D& m,
                                               const std::map<int, int>& index,
                                               int nDofs,
                                               double penalty) ;

DenseVec constraintPenaltyResidualForce(const Model2D& m,
                                               const std::map<int, int>& index,
                                               int nDofs,
                                               double penalty,
                                               const DenseVec& U) ;

// 自动约束"无转动刚度的转角自由度"。
// 二维统一按每节点 3 个自由度 [u, v, θ] 编号，但 TRUSS/CABLE/STRUT 只贡献
// 平动刚度。纯桁架节点的 θ 对角元恒为 0，若不由用户显式约束，总体矩阵必然
// 奇异，且报错信息对使用者很不友好。这里自动补零约束（等价于 ANSYS 等软件
// 对未用自由度的处理），平移自由度仍然保持严格的零刚度检查。
void autoConstrainFreeRotations(const SparseMat& K, int nDofs,
                                       const std::vector<int>& baseCdofs,
                                       const std::vector<double>& baseCvals,
                                       std::vector<int>& cdofs,
                                       std::vector<double>& cvals) ;

// 方案 C：几何非线性轴向单元 + 线性欧拉梁的 Newton-Raphson 求解。
// 当前版本非线性部分覆盖 TRUSS/CABLE/STRUT；BEAM 仍使用线性切线。
ModelResult solveModel2DNonlinear(const Model2D& m,
                                         SolverType type = SolverType::SimplicialLDLT,
                                         const DenseVec* U0 = nullptr,
                                         const std::vector<char>* active0 = nullptr,
                                         NonlinearDiag* diagOut = nullptr,
                                         ReusableSolver* reuse = nullptr) ;

// ==================== 线性系统的装配（供多次求解复用）====================

// 给定激活状态，装配总体刚度矩阵 K 与荷载向量 F（尚未施加约束）。
// solveModel2D 与影响线扫描共用这一段，避免两处逻辑走偏。
void assembleLinear2D(const Model2D& m,
                             const std::vector<char>& active,
                             const std::map<int, int>& index,
                             int nDofs,
                             SparseMat& K, DenseVec& F) ;

// 装配结果：影响线/车道扫描可复用同一份 K 的分解，避免逐位置重新装配+重新分解。
struct Assembled2D {
    SparseMat K0;                 // 未施加约束（反力 = K0*U - F0）
    DenseVec  F0;
    SparseMat K;                  // 已施加 Dirichlet，可直接分解
    DenseVec  F;                  // 已施加 Dirichlet
    std::vector<int> cdofs;
    std::vector<double> cvals;
    std::vector<char> isConstrained;
    std::vector<int> ids;         // 0-based 索引 -> 节点号
    std::map<int, int> index;
    int nDofs = 0;
};

// 一次性装配（全部单元激活）+ 施加 Dirichlet。
// 仅适用于线性模型；含 CABLE/STRUT 时激活状态依赖位移，必须逐次迭代。
Assembled2D assembleModel2DLinear(const Model2D& m) ;

// 一致质量矩阵（全局稀疏，nDofs×nDofs）。线质量密度 = dens·A + addmass。
SparseMat assembleMassMatrix2D(const Model2D& m) ;

// 广义特征值问题 K·φ = λ·M·φ 的自振频率（Hz，升序）与振型。
struct ModalResult2D {
    std::vector<double> frequencies;  // Hz，升序
    std::vector<DenseVec> modes;      // 全局振型（约束自由度处为 0）
};
ModalResult2D computeNaturalFrequencies2D(const Model2D& m, int nModes = 6) ;

// 线性屈曲荷载因子：参考荷载由 m.loads 指定，返回前 n 阶屈曲因子 λ（P_cr = λ·P_ref）。
// 仅考虑 BEAM 的轴力几何刚度（P-Δ）；受压为负、受拉为正。
std::vector<double> computeBucklingFactors2D(const Model2D& m, int nModes = 3) ;

// 冲击系数（按基频 f，Hz）：JTG D60 简支梁/连续梁近似。
double impactCoefficientFromFrequency(double fHz) ;

struct LaneImpact2D {
    double f1 = 0.0;
    double mu = 0.0;
    double factor = 1.0;
    bool automatic = false;
};

// IMPACT=1 时由竖向弯曲基频按 JTG D60 取 1+μ。算不出基频则 automatic 为假。
LaneImpact2D resolveLaneImpact2D(const Model2D& m);
double laneImpactFactor(const Model2D::LaneLoad& lane, const LaneImpact2D& imp);

// 由位移与装配结果构造 ModelResult（含反力），供响应计算复用。
ModelResult makeResult2D(const Model2D& m, const Assembled2D& a,
                                const DenseVec& U, const DenseVec& F0) ;

// 组装并求解二维结构（默认 SimplicialLDLT）。
ModelResult solveModel2D(const Model2D& m,
                                SolverType type = SolverType::SimplicialLDLT,
                                const DenseVec* U0 = nullptr,
                                const std::vector<char>* active0 = nullptr,
                                NonlinearDiag* diagOut = nullptr,
                                ReusableSolver* reuse = nullptr) ;

// 单元内力：轴力 N（受拉为正）；梁单元还给出两端剪力 V 与弯矩 M。
struct ElementForce {
    int id = 0;
    std::string type;
    double N1 = 0.0;  // 节点1端轴力（受拉为正）
    double N2 = 0.0;  // 节点2端轴力（受拉为正）
    double V1 = 0.0, M1 = 0.0; // 梁：节点 1 剪力/弯矩
    double V2 = 0.0, M2 = 0.0; // 梁：节点 2 剪力/弯矩
};

// 由节点位移计算各单元内力。
std::vector<ElementForce> elementForces2D(const Model2D& m, const DenseVec& U) ;

// 单元应力：桁架轴向应力；梁上下翼缘正应力。
struct ElementStress {
    int id = 0;
    std::string type;
    double s1 = 0.0, s2 = 0.0;       // 桁架：两端轴向应力
    double stop1 = 0.0, sbot1 = 0.0; // 梁：节点1端上/下翼缘应力
    double stop2 = 0.0, sbot2 = 0.0; // 梁：节点2端上/下翼缘应力
};

// 由内力计算应力（受拉为正）。
std::vector<ElementStress> elementStress2D(const Model2D& m,
                                                  const std::vector<ElementForce>& forces) ;

// ==================== 施工阶段正装分析 ====================

struct StageResult {
    int stage = 0;
    DenseVec U;
    DenseVec R;
    std::vector<ElementForce> forces;
    std::vector<ElementStress> stresses;
};

std::vector<StageResult> solveModel2DStages(const Model2D& m,
                                                   SolverType type = SolverType::SimplicialLDLT) ;

struct CamberResult {
    int node = 0;
    int dof = 0;
    double camber = 0.0;      // 预拱度 = -最终位移
    double finalDisplacement = 0.0;
    double original = 0.0;
    double installation = 0.0;
};

std::vector<CamberResult> computeCamber(const Model2D& m,
                                               SolverType type = SolverType::SimplicialLDLT) ;

struct CamberIterationResult {
    DenseVec cambers;
    DenseVec finalDisplacements;
    int iterations = 0;
    std::vector<StageResult> stages;
};

CamberIterationResult iterateCamber(const Model2D& m,
                                           SolverType type = SolverType::SimplicialLDLT,
                                           double tol = 1e-8,
                                           int maxIter = 30) ;

struct FeedbackCorrection {
    double current = 0.0;
    double derivative = 0.0;
    double correction = 0.0;
};

FeedbackCorrection computeFeedbackCorrection(const Model2D& m,
                                                     SolverType type = SolverType::SimplicialLDLT) ;

std::vector<StageResult> solveModel2DStagesBackward(const Model2D& m,
                                                           SolverType type = SolverType::SimplicialLDLT) ;

double forwardBackwardStageResidual(const Model2D& m,
                                           SolverType type = SolverType::SimplicialLDLT) ;

struct StageIterationResult {
    DenseVec scales;
    DenseVec finalValues;
    int iterations = 0;
    std::vector<StageResult> stages;
};

StageIterationResult iterateStagePrestress(const Model2D& m,
                                                  SolverType type = SolverType::SimplicialLDLT,
                                                  double tol = 1e-8,
                                                  int maxIter = 30) ;

// ==================== 影响线与车道荷载包络 ====================

struct InfluenceCurve {
    std::string label;
    std::vector<double> s;       // 车道纵向位置
    std::vector<double> eta;     // 单位荷载影响值
    double vmax = 0.0, vmin = 0.0;
    double smax = 0.0, smin = 0.0;
};

struct EnvelopeValue {
    std::string label;
    double maxVal = 0.0, minVal = 0.0;
    double smax = 0.0, smin = 0.0;
};

struct EnvelopeSpec {
    Model2D::InfluenceRequest req;
    std::string label;
};

double evaluateInfluenceResponse(const Model2D& m, const ModelResult& r,
                                        const Model2D::InfluenceRequest& spec) ;

InfluenceCurve computeInfluenceCurve(const Model2D& m,
                                             const Model2D::LaneLoad& lane,
                                             const Model2D::InfluenceRequest& spec,
                                             SolverType type = SolverType::SimplicialLDLT) ;

std::vector<double> laneTributaryWeights(const Model2D& m,
                                                const Model2D::LaneLoad& lane) ;

EnvelopeValue envelopeFromInfluence(const Model2D& m,
                                           const Model2D::LaneLoad& lane,
                                           const InfluenceCurve& curve) ;

EnvelopeValue nonlinearEnvelopeFromLane(const Model2D& m,
                                               const Model2D::LaneLoad& lane,
                                               const Model2D::InfluenceRequest& spec,
                                               SolverType type = SolverType::SimplicialLDLT) ;

std::vector<EnvelopeSpec> defaultEnvelopeSpecs(const Model2D& m) ;

// 把节点位移、单元内力、支座反力写入结果文件。
void writeResult(const std::string& path, const Model2D& model,
                        const ModelResult& result,
                        const std::vector<ElementForce>& forces) ;

} // namespace fem
