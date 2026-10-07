#pragma once

#include "fem.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// 二维结构输入文件读取与组装。
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
        std::vector<int> conn; // 两个节点号
    };
    std::vector<Element> elements;

    struct Constraint {
        int node;
        int dof;    // 1/2/3
        double value;
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
        double Ep = 195e9;          // 预应力筋弹性模量（Pa）
        double Ap = 0.0;            // 预应力筋面积（m²），TENDON 方式必填
        std::string method = "AUTO"; // AUTO/EQUIV/TENDON
    };
    std::vector<Prestress> prestresses;

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
inline bool isAxial(const std::string& type) {
    return type == "TRUSS" || type == "CABLE" || type == "STRUT";
}

inline bool isRigid(const std::string& type) {
    return type == "RIGID" || type == "RIGID_PIN";
}

// 销接刚臂：只约束梁节点与筋节点之间的平动刚体关系，不向筋节点传递弯矩。
// 用于无粘结预应力筋，避免在转向点产生虚假附加弯矩。
inline DenseMat rigidPinStiffness(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                  double penalty = 1e11) {
    const Eigen::Vector2d r = p2 - p1; // 梁节点 -> 筋节点

    // DOF: [u_beam, v_beam, θ_beam, u_tendon, v_tendon, θ_tendon]
    // 前两个约束为：u_t - u_b + r_y·θ_b = 0
    //               v_t - v_b - r_x·θ_b = 0
    Eigen::Matrix<double, 2, 6> C = Eigen::Matrix<double, 2, 6>::Zero();
    C(0, 0) = -1.0; C(0, 2) =  r(1); C(0, 3) = 1.0;
    C(1, 1) = -1.0; C(1, 2) = -r(0); C(1, 4) = 1.0;

    return penalty * (C.transpose() * C);
}

// 计算轴向单元的机械轴力（受拉为正）。
// 约定：N = E·A/L·ΔL + E·A·ISTRN。
inline double axialForce(const Model2D::Element& e, const DenseVec& U,
                         int i1, int i2,
                         const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const Eigen::Vector2d c = d / L;
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const double dL = (u2 - u1).dot(c);
    return e.E * e.A / L * dL + e.E * e.A * e.istrn;
}

// 二维单元全局自由度编号（桁架 4 个，梁 6 个）。
inline std::vector<int> elementGdofs2D(const std::string& type, int i1, int i2) {
    if (isAxial(type))
        return {3 * i1 + 0, 3 * i1 + 1, 3 * i2 + 0, 3 * i2 + 1};
    return {3 * i1 + 0, 3 * i1 + 1, 3 * i1 + 2,
            3 * i2 + 0, 3 * i2 + 1, 3 * i2 + 2};
}

// 把方向字段解析成全局单位方向向量。
// GX/GY：全局轴；LX/LY：局部轴/横向；G<θ>：与全局X成θ角；L<θ>：与杆轴(i→j)成θ角。
inline Eigen::Vector2d resolveDir(const std::string& dir,
                                  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const std::string t = detail::toUpper(detail::trim(dir));
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const double c = d(0) / L, s = d(1) / L;
    if (t == "GX") return {1.0, 0.0};
    if (t == "GY") return {0.0, 1.0};
    if (t == "LX") return {c, s};
    if (t == "LY") return {-s, c};
    if (t.size() >= 2 && (t[0] == 'G' || t[0] == 'L')) {
        const double th = std::stod(t.substr(1)) * kPi / 180.0;
        const double co = std::cos(th), si = std::sin(th);
        if (t[0] == 'G') return {co, si};
        return {c * co - s * si, s * co + c * si}; // 局部 -> 全局
    }
    throw std::runtime_error("未知荷载方向: " + dir);
}

// 分布荷载等效节点荷载（全局坐标）。q1/q2 为 i/j 端沿方向的强度。
inline DenseVec distributedLoadEquivalent(const std::string& type,
                                          const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                          const Eigen::Vector2d& gdir,
                                          double q1, double q2) {
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const double c = d(0) / L, s = d(1) / L;

    if (isAxial(type)) {
        DenseVec fe(4);
        fe << L / 6.0 * (2.0 * q1 + q2) * gdir(0),
              L / 6.0 * (2.0 * q1 + q2) * gdir(1),
              L / 6.0 * (q1 + 2.0 * q2) * gdir(0),
              L / 6.0 * (q1 + 2.0 * q2) * gdir(1);
        return fe;
    }

    // BEAM：全局方向分解到局部轴向/横向。
    const double lx =  c * gdir(0) + s * gdir(1);
    const double ly = -s * gdir(0) + c * gdir(1);
    const double qx1 = q1 * lx, qx2 = q2 * lx;
    const double qy1 = q1 * ly, qy2 = q2 * ly;

    DenseVec fl(6);
    fl << L / 6.0 * (2.0 * qx1 + qx2),
          L / 20.0 * (7.0 * qy1 + 3.0 * qy2),
          L * L / 60.0 * (3.0 * qy1 + 2.0 * qy2),
          L / 6.0 * (qx1 + 2.0 * qx2),
          L / 20.0 * (3.0 * qy1 + 7.0 * qy2),
          -L * L / 60.0 * (2.0 * qy1 + 3.0 * qy2);

    EulerBeam2D be{1.0, 1.0, 1.0, p1, p2}; // 只需其坐标变换
    return be.transformation().transpose() * fl;
}

// 均匀温度变化等效节点荷载（轴力）。
inline DenseVec uniformTempEquivalent(const std::string& type, double E, double A, double alpha,
                                      double dT, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const double NT = E * A * alpha * dT;
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const Eigen::Vector2d c = d / L;
    if (isAxial(type)) {
        DenseVec fe(4);
        fe << -NT * c(0), -NT * c(1), NT * c(0), NT * c(1);
        return fe;
    }
    DenseVec fl = DenseVec::Zero(6);
    fl(0) = -NT;
    fl(3) = NT;
    EulerBeam2D be{1.0, 1.0, 1.0, p1, p2};
    return be.transformation().transpose() * fl;
}

// 梁温度梯度等效节点荷载（弯曲）。
inline DenseVec gradientTempEquivalent(double E, double I, double alpha, double h, double dT,
                                       const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    if (h <= 0.0)
        throw std::runtime_error("梁温度梯度需要 *ELEMENT 里的 H");
    const double MT = E * I * alpha * dT / h;
    DenseVec fl = DenseVec::Zero(6);
    fl(2) = -MT;
    fl(5) = MT;
    EulerBeam2D be{1.0, 1.0, 1.0, p1, p2};
    return be.transformation().transpose() * fl;
}

// 规范竖向温度曲线 -> 截面温度轴力/弯矩。
inline std::pair<double, double> tempProfileResultants(double E, double A, double H,
                                                       double alpha,
                                                       double T1, double T2,
                                                       double A0, double B) {
    if (H <= 0.0) throw std::runtime_error("TEMP_PROFILE 需要梁高 H");
    if (A0 <= 0.0) throw std::runtime_error("TEMP_PROFILE A 必须大于 0");
    const double b = A / H;
    const double B0 = (B <= 0.0 || B > H) ? H : B;
    if (B0 < A0) throw std::runtime_error("TEMP_PROFILE B 不能小于 A");

    auto T = [&](double z) -> double {
        if (z <= 0.0) return T1;
        if (z <= A0)
            return T1 + (T2 - T1) * (z / A0);
        if (z <= B0)
            return T2 + (0.0 - T2) * ((z - A0) / (B0 - A0));
        return 0.0;
    };

    const int n = 400;
    double N = 0.0, M = 0.0;
    for (int i = 0; i < n; ++i) {
        const double z0 = i * H / n;
        const double z1 = (i + 1) * H / n;
        const double zm = 0.5 * (z0 + z1);
        const double y = H * 0.5 - zm; // 向上为正
        const double dz = z1 - z0;
        N += b * T(zm) * dz;
        M += b * T(zm) * y * dz;
    }
    return {E * alpha * N, E * alpha * M};
}

inline DenseVec tempProfileEquivalent(const Model2D::Element& e,
                                      const Eigen::Vector2d& p1,
                                      const Eigen::Vector2d& p2,
                                      double N, double M) {
    DenseVec fl = DenseVec::Zero(6);
    fl(0) = -N;
    fl(3) = N;
    fl(2) = -M;
    fl(5) = M;
    EulerBeam2D be{e.E, e.A, e.I, p1, p2};
    return be.transformation().transpose() * fl;
}

// 把单元等效荷载散入全局 F。
inline void addElementLoad(const std::vector<int>& gdofs, const DenseVec& fe, DenseVec& F) {
    for (std::size_t i = 0; i < gdofs.size(); ++i)
        F(gdofs[i]) += fe(i);
}

// 单个单元上所有分布/自重/温度荷载的等效节点荷载（全局坐标，按单元 gdofs 顺序）。
inline DenseVec elementEquivalentLoad(const Model2D& m, const Model2D::Element& e,
                                      const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    if (isRigid(e.type))
        return DenseVec::Zero(6);
    const int nd = isAxial(e.type) ? 4 : 6;
    DenseVec fe = DenseVec::Zero(nd);

    for (const auto& el : m.elemLoads) {
        if (el.elem != e.id) continue;
        const Eigen::Vector2d dir = resolveDir(el.dir, p1, p2);
        fe += distributedLoadEquivalent(e.type, p1, p2, dir, el.q1, el.q2);
    }

    const double gmag = std::sqrt(m.gx * m.gx + m.gy * m.gy);
    Eigen::Vector2d gdir(0.0, 0.0);
    if (gmag > 0.0) gdir = Eigen::Vector2d(m.gx / gmag, m.gy / gmag);

    // 自动自重：材料密度×截面积 + 实常数附加质量（线密度 kg/m）。
    // e.addmass 既参与静力自重，也保留在 Element 里供后续动力分析使用。
    const double autoLineMass = e.dens * e.A + e.addmass;
    if (autoLineMass > 0.0) {
        const double w = autoLineMass * gmag;
        fe += distributedLoadEquivalent(e.type, p1, p2, gdir, w, w);
    }

    // 兼容旧格式的显式 *SELF_WEIGHT。
    for (const auto& sw : m.selfWeights) {
        if (!sw.all && sw.elem != e.id) continue;
        double w = 0.0;
        if (sw.rule == "RHO") w = sw.value * gmag * e.A;
        else if (sw.rule == "MASS") w = sw.value * gmag;
        else throw std::runtime_error("未知自重规则: " + sw.rule);
        fe += distributedLoadEquivalent(e.type, p1, p2, gdir, w, w);
    }

    if (m.uniformTemp != 0.0)
        fe += uniformTempEquivalent(e.type, e.E, e.A, e.alpha, m.uniformTemp, p1, p2);

    for (const auto& et : m.elemTemps) {
        if (et.elem != e.id) continue;
        if (et.kind == "UNIFORM")
            fe += uniformTempEquivalent(e.type, e.E, e.A, e.alpha, et.dT, p1, p2);
        else if (et.kind == "GRADIENT")
            fe += gradientTempEquivalent(e.E, e.I, e.alpha, e.h, et.dT, p1, p2);
    }

    // 收缩应变：正值为缩短，按负温度应变等效为节点荷载。
    for (const auto& sh : m.shrinkageSpecs) {
        if (!sh.all && std::find(sh.elems.begin(), sh.elems.end(), e.id) == sh.elems.end())
            continue;
        fe += uniformTempEquivalent(e.type, e.E, e.A, 1.0, -sh.ecs, p1, p2);
    }

    // 规范竖向温度梯度：仅用于梁单元。
    for (const auto& tp : m.tempProfiles) {
        if (e.type != "BEAM") continue;
        if (!tp.all && std::find(tp.elems.begin(), tp.elems.end(), e.id) == tp.elems.end())
            continue;
        const double h = (tp.H > 0.0) ? tp.H : e.h;
        const auto [N, M] = tempProfileResultants(e.E, e.A, h, e.alpha,
                                                  tp.T1, tp.T2, tp.A, tp.B);
        fe += tempProfileEquivalent(e, p1, p2, N, M);
    }

    // 初应变 ISTRN（cable/strut/truss 预应力）：等效节点力 ±EA·ISTRN 沿杆轴。
    if (e.istrn != 0.0 && isAxial(e.type)) {
        const Eigen::Vector2d d = p2 - p1;
        const double L = d.norm();
        const Eigen::Vector2d c = d / L;
        const double P = e.E * e.A * e.istrn;
        fe(0) +=  P * c(0); fe(1) +=  P * c(1);   // 节点 i 沿 +c 方向
        fe(2) += -P * c(0); fe(3) += -P * c(1);   // 节点 j 沿 -c 方向
    }
    return fe;
}

// ==================== 方案 C：几何非线性轴向单元 ====================

// 当前构形下的机械轴力（受拉为正），工程应变 + 初应变。
inline double axialForceNonlinear(const Model2D::Element& e, const DenseVec& U,
                                  int i1, int i2,
                                  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d0 = p2 - p1;
    const double L0 = d0.norm();
    if (L0 <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const Eigen::Vector2d d = d0 + u2 - u1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
    return e.E * e.A * ((L - L0) / L0 + e.istrn);
}

// 轴向单元当前构形下的内部节点力向量：f_int = [-N c, +N c]。
inline DenseVec axialInternalForceNonlinear(const Model2D::Element& e, const DenseVec& U,
                                             int i1, int i2,
                                             const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d0 = p2 - p1;
    const double L0 = d0.norm();
    if (L0 <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const Eigen::Vector2d d = d0 + u2 - u1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
    const Eigen::Vector2d c = d / L;
    const double N = e.E * e.A * ((L - L0) / L0 + e.istrn);

    DenseVec fe(4);
    fe << -N * c(0), -N * c(1), N * c(0), N * c(1);
    return fe;
}

// 轴向单元一致切线刚度矩阵（对应上述 f_int）。
inline DenseMat axialTangentNonlinear(const Model2D::Element& e, const DenseVec& U,
                                      int i1, int i2,
                                      const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d0 = p2 - p1;
    const double L0 = d0.norm();
    if (L0 <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const Eigen::Vector2d d = d0 + u2 - u1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
    const Eigen::Vector2d c = d / L;
    const double N = e.E * e.A * ((L - L0) / L0 + e.istrn);

    const double kMat = e.E * e.A / L0;
    const double kGeo = N / L;
    const Eigen::Matrix2d A =
        (kMat - kGeo) * (c * c.transpose()) + kGeo * Eigen::Matrix2d::Identity();

    DenseMat Kt = DenseMat::Zero(4, 4);
    Kt.block<2, 2>(0, 0) = A;
    Kt.block<2, 2>(0, 2) = -A;
    Kt.block<2, 2>(2, 0) = -A;
    Kt.block<2, 2>(2, 2) = A;
    return Kt;
}

// 方案 C 外部等效节点荷载：与线性路径一致，但不包含初应变 ISTRN。
// ISTRN 已计入轴向单元内部力。
inline DenseVec elementExternalEquivalentLoad(const Model2D& m, const Model2D::Element& e,
                                              const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    DenseVec fe = elementEquivalentLoad(m, e, p1, p2);
    if (e.istrn != 0.0 && isAxial(e.type)) {
        const Eigen::Vector2d d = p2 - p1;
        const double L = d.norm();
        if (L <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
        const Eigen::Vector2d c = d / L;
        const double P = e.E * e.A * e.istrn;
        fe(0) -= P * c(0);
        fe(1) -= P * c(1);
        fe(2) += P * c(0);
        fe(3) += P * c(1);
    }
    return fe;
}

// 处理 *PRESTRESS：等效荷载法或自动生成预应力筋单元。
inline void processPrestress(Model2D& m) {
    if (m.prestresses.empty()) return;

    auto getElem = [&](int id) -> Model2D::Element* {
        for (auto& e : m.elements)
            if (e.id == id) return &e;
        return nullptr;
    };

    auto isConstrained = [&](int node, int dof) {
        for (const auto& c : m.constraints)
            if (c.node == node && c.dof == dof) return true;
        return false;
    };

    int nextNode = 0, nextElem = 0, nextType = 0, nextMat = 0, nextReal = 0;
    for (const auto& [id, c] : m.nodes) nextNode = std::max(nextNode, id);
    for (const auto& e : m.elements) nextElem = std::max(nextElem, e.id);
    for (const auto& [id, v] : m.elementTypes) nextType = std::max(nextType, id);
    for (const auto& [id, v] : m.materials) nextMat = std::max(nextMat, id);
    for (const auto& [id, v] : m.realConsts) nextReal = std::max(nextReal, id);

    for (const auto& ps : m.prestresses) {
        if (ps.elems.empty())
            throw std::runtime_error("*PRESTRESS 的 ELEM 列表为空");
        const double Pe = ps.P * ps.loss;

        // 按输入顺序建立梁节点链。
        std::vector<int> beamNodes;
        int prevEnd = 0;
        for (std::size_t i = 0; i < ps.elems.size(); ++i) {
            Model2D::Element* e = getElem(ps.elems[i]);
            if (!e) throw std::runtime_error("*PRESTRESS 引用了不存在的单元: " + std::to_string(ps.elems[i]));
            if (e->type != "BEAM")
                throw std::runtime_error("*PRESTRESS 目前只支持施加到 BEAM 单元");
            if (i == 0) {
                beamNodes.push_back(e->conn[0]);
                beamNodes.push_back(e->conn[1]);
            } else {
                if (e->conn[0] != prevEnd)
                    throw std::runtime_error("*PRESTRESS 的 ELEM 必须按节点顺序连续");
                beamNodes.push_back(e->conn[1]);
            }
            prevEnd = e->conn[1];
        }

        const int startNode = beamNodes.front();
        const int endNode = beamNodes.back();
        const Eigen::Vector2d pstart = m.nodes.at(startNode);
        const Eigen::Vector2d pend = m.nodes.at(endNode);
        const Eigen::Vector2d dspan = pend - pstart;
        const double L = dspan.norm();
        if (L <= 0.0) throw std::runtime_error("*PRESTRESS 梁组长度为零");
        const Eigen::Vector2d t = dspan / L;
        const Eigen::Vector2d n(-t(1), t(0)); // 局部 +y

        std::string method = ps.method;
        if (method == "AUTO")
            method = (m.analysisType == "NONLINEAR") ? "TENDON" : "EQUIV";

        if (method == "EQUIV") {
            if (ps.F != 0.0) {
                const double w = 8.0 * Pe * ps.F / (L * L);
                for (int elem : ps.elems)
                    m.elemLoads.push_back({elem, "LY", w, w});
            }

            if (ps.E1 != 0.0)
                m.loads.push_back({startNode, 3, Pe * ps.E1});
            if (ps.E2 != 0.0)
                m.loads.push_back({endNode, 3, -Pe * ps.E2});

            const bool startAxialFixed = isConstrained(startNode, 1);
            const bool endAxialFixed = isConstrained(endNode, 1);
            if (startAxialFixed && !endAxialFixed)
                m.loads.push_back({endNode, 1, -Pe});
            else if (!startAxialFixed && endAxialFixed)
                m.loads.push_back({startNode, 1, Pe});
            else if (startAxialFixed && endAxialFixed)
                throw std::runtime_error("*PRESTRESS EQUIV：两端 dof1 均固定时无法用外荷载模拟弹性缩短，请使用 METHOD=TENDON");
            else
                throw std::runtime_error("*PRESTRESS EQUIV：梁组至少需要一端约束 dof1");
            continue;
        }

        if (method != "TENDON")
            throw std::runtime_error("*PRESTRESS 未知 METHOD: " + ps.method);
        if (ps.Ap <= 0.0)
            throw std::runtime_error("*PRESTRESS TENDON 需要 APS/AP 预应力筋面积");

        const double istrn = Pe / (ps.Ep * ps.Ap);

        const int tendonType = ++nextType;
        const int linkType = ++nextType;
        const int tendonMat = ++nextMat;
        const int linkMat = ++nextMat;
        const int tendonReal = ++nextReal;
        const int linkReal = ++nextReal;

        m.elementTypes[tendonType] = "CABLE";
        m.elementTypes[linkType] = "RIGID_PIN";

        Model2D::Material tm;
        tm.E = ps.Ep;
        tm.alpha = 0.0;
        tm.dens = 0.0;
        m.materials[tendonMat] = tm;

        Model2D::Material lm;
        lm.E = 1e11;
        lm.alpha = 0.0;
        lm.dens = 0.0;
        m.materials[linkMat] = lm;

        Model2D::RealConst tr;
        tr.A = ps.Ap;
        tr.istrn = istrn;
        m.realConsts[tendonReal] = tr;

        Model2D::RealConst lr;
        lr.A = 1e-2;
        lr.I = 1e-5;
        m.realConsts[linkReal] = lr;

        std::vector<int> tendonNodes;
        tendonNodes.reserve(beamNodes.size());
        for (std::size_t k = 0; k < beamNodes.size(); ++k) {
            const int bn = beamNodes[k];
            const Eigen::Vector2d pb = m.nodes.at(bn);
            const double s = std::max(0.0, std::min(L, (pb - pstart).dot(t)));
            const double xi = L > 0.0 ? s / L : 0.0;
            const double e = (1.0 - xi) * ps.E1 + xi * ps.E2 + 4.0 * ps.F * xi * (1.0 - xi);
            const Eigen::Vector2d pt = pb - e * n;

            if ((pt - pb).norm() < 1e-8) {
                tendonNodes.push_back(bn); // 锚固点或偏心为 0 时直接共用梁节点
                continue;
            }

            const int tn = ++nextNode;
            m.nodes[tn] = pt;
            tendonNodes.push_back(tn);
            // 筋节点只有 u/v 参与受力；θ 为无刚度自由度，固定为 0。
            m.constraints.push_back({tn, 3, 0.0});

            Model2D::Element link;
            link.id = ++nextElem;
            link.type = "RIGID_PIN";
            link.E = lm.E;
            link.A = lr.A;
            link.I = lr.I;
            link.conn = {bn, tn};
            m.elements.push_back(link);
        }

        for (std::size_t k = 0; k + 1 < tendonNodes.size(); ++k) {
            Model2D::Element cable;
            cable.id = ++nextElem;
            cable.type = "CABLE";
            cable.E = ps.Ep;
            cable.A = ps.Ap;
            cable.istrn = istrn;
            cable.conn = {tendonNodes[k], tendonNodes[k + 1]};
            m.elements.push_back(cable);
        }
    }
}

// 记录徐变折减基准（未折减的初始弹性模量）。
// 必须在所有会新增单元的步骤（如 processPrestress 生成预应力筋）之后调用，
// 这样新生成的钢束单元也有正确的基准值。
inline void snapshotBaseModulus(Model2D& m) {
    for (const auto& e : m.elements)
        m.baseModulus[e.id] = e.E;
}

// 取单元的徐变折减基准模量；没有记录时退化为单元自身的 E。
inline double baseModulusOf(const Model2D& m, const Model2D::Element& e) {
    auto it = m.baseModulus.find(e.id);
    return (it != m.baseModulus.end()) ? it->second : e.E;
}

// AEMM（龄期调整有效模量法）徐变折减系数。
//   Eeff = E0 / (1 + chi * Δφ)，  Δφ = φ(t_now) - φ(t_load) >= 0
inline double creepReductionFactor(double chi, double dphi) {
    const double r = chi * std::max(0.0, dphi);
    return 1.0 / (1.0 + r);
}

inline double interpolateTimeCurve(const Model2D::TimeCurve& curve, double time) {
    if (curve.points.empty())
        throw std::runtime_error("时变曲线没有数据点");
    std::vector<std::pair<double, double>> pts = curve.points;
    std::sort(pts.begin(), pts.end());
    if (time <= pts.front().first) return pts.front().second;
    if (time >= pts.back().first) return pts.back().second;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (time <= pts[i].first) {
            const double t0 = pts[i - 1].first, t1 = pts[i].first;
            const double v0 = pts[i - 1].second, v1 = pts[i].second;
            return v0 + (v1 - v0) * (time - t0) / (t1 - t0);
        }
    }
    return pts.back().second;
}

inline double creepCoefficientAt(const Model2D& m, double time) {
    if (m.creepModel == "CURVE" && !m.creepCurve.points.empty())
        return interpolateTimeCurve(m.creepCurve, time);
    const double dt = std::max(0.0, time - m.creepT0);
    return m.creepPhi0 * std::pow(dt / (m.creepBetaH + dt), m.creepPower);
}

inline double shrinkageStrainAt(const Model2D& m, double time) {
    if (m.shrinkageModel == "CURVE" && !m.shrinkageCurve.points.empty())
        return interpolateTimeCurve(m.shrinkageCurve, time);
    const double dt = std::max(0.0, time - m.shrinkageTs);
    return m.shrinkageEcsInf * std::pow(dt / (m.shrinkageAlphaDs + dt), m.shrinkagePower);
}

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
inline void processTimeEffects(Model2D& m,
                               const std::map<int, double>* loadAge = nullptr,
                               double now = 0.0) {
    if (m.creepSpecs.empty()) return;
    const bool useAge = (now > 0.0);

    for (auto& e : m.elements) {
        if (isRigid(e.type)) continue;

        double bestChi = 0.0, bestDphi = 0.0;
        bool hit = false;
        for (const auto& cs : m.creepSpecs) {
            if (!cs.all && std::find(cs.elems.begin(), cs.elems.end(), e.id) == cs.elems.end())
                continue;
            double dphi = cs.phi;
            if (useAge) {
                double t0 = 0.0;
                if (loadAge) {
                    auto it = loadAge->find(e.id);
                    if (it != loadAge->end()) t0 = it->second;
                }
                // AEMM：只计加载龄期之后发展的那部分徐变
                dphi = creepCoefficientAt(m, now) - creepCoefficientAt(m, t0);
            }
            const double red = cs.chi * std::max(0.0, dphi);
            if (!hit || red > bestChi * bestDphi) {
                bestChi = cs.chi;
                bestDphi = std::max(0.0, dphi);
                hit = true;
            }
        }
        if (!hit) continue;

        const double E0 = baseModulusOf(m, e);
        e.E = E0 * creepReductionFactor(bestChi, bestDphi);
    }
}

// 读取二维结构输入文件。
// 支持两种风格：
//   1) 新风格（表声明 + *EN 命令 + *DO/*FILL 预处理）；
//   2) 旧风格（*ELEMENT, TYPE=..., E=..., A=... 内联模板），保持兼容。
inline Model2D readModel2D(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开文件: " + path);

    const std::vector<std::string> lines = detail::preprocess(in);

    Model2D m;
    std::string section;
    Model2D::Element legacyTemplate; // 旧式 *ELEMENT 段的内联模板
    int activeType = 0, activeMat = 0, activeReal = 0;

    for (const std::string& line : lines) {
        const std::string t = detail::trim(line);
        if (t.empty() || t[0] == '#' || t[0] == '!') continue;

        if (t[0] == '*') {
            // 关键字 = 第一个空白或逗号之前的部分；其余为参数。
            const std::size_t sep = t.find_first_of(", \t");
            const std::string kw = detail::toUpper(sep == std::string::npos ? t : t.substr(0, sep));
            const std::string rest = sep == std::string::npos ? std::string() : detail::trim(t.substr(sep + 1));

            // ---- 命令：*TYPE / *MAT / *REAL / *EN ----
            if (kw == "*TYPE" || kw == "*MAT" || kw == "*REAL" || kw == "*EN") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                if (kw == "*TYPE") {
                    ss >> activeType;
                } else if (kw == "*MAT") {
                    ss >> activeMat;
                } else if (kw == "*REAL") {
                    ss >> activeReal;
                } else { // *EN,eid,nodei,nodej
                    int eid = 0, ni = 0, nj = 0;
                    if (!(ss >> eid >> ni >> nj))
                        throw std::runtime_error("*EN 行格式错误: " + t);
                    if (activeType <= 0 || activeMat <= 0 || activeReal <= 0)
                        throw std::runtime_error("*EN 前需先激活 *TYPE/*MAT/*REAL");
                    auto tit = m.elementTypes.find(activeType);
                    auto mit = m.materials.find(activeMat);
                    auto rit = m.realConsts.find(activeReal);
                    if (tit == m.elementTypes.end() || mit == m.materials.end() || rit == m.realConsts.end())
                        throw std::runtime_error("*EN 引用的 TYPE/MAT/REAL 未定义");
                    Model2D::Element e;
                    e.id = eid;
                    e.type = tit->second;
                    e.E = mit->second.E;
                    e.alpha = mit->second.alpha;
                    e.dens = mit->second.dens;
                    e.A = rit->second.A;
                    e.I = rit->second.I;
                    e.Wt = rit->second.Wt;
                    e.Wb = rit->second.Wb;
                    e.h = rit->second.h;
                    e.addmass = rit->second.addmass;
                    e.istrn = rit->second.istrn;
                    e.conn = {ni, nj};
                    if (e.type == "BEAM" && e.I <= 0.0)
                        throw std::runtime_error("梁单元缺少实常数 I");
                    if (e.type == "BEAM" && (e.Wt <= 0.0 || e.Wb <= 0.0))
                        throw std::runtime_error("梁单元缺少 WT/WB（截面抵抗矩，应力计算必需）");
                    if (isAxial(e.type) && e.A <= 0.0)
                        throw std::runtime_error("桁架/索/压杆单元缺少实常数 A");
                    m.elements.push_back(e);
                }
                continue;
            }

            // *TEMPERATURE 既支持命令式（同行带值），也支持段式（下一行给值）。
            if (kw == "*TEMPERATURE" && !rest.empty()) {
                m.uniformTemp = std::stod(rest);
                continue;
            }

            // *ANALYSIS：分析类型。LINEAR=方案A（默认），NONLINEAR=方案C。
            // 可选参数：TOL=<相对容差> MAXITER=<最大迭代次数>。
            if (kw == "*ANALYSIS") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                auto p = detail::parseKeyValues(args);
                std::string v;
                std::istringstream probe(args);
                std::string first;
                probe >> first;
                const std::string firstUp = detail::toUpper(detail::trim(first));
                if (firstUp == "LINEAR" || firstUp == "NONLINEAR") {
                    v = firstUp;
                } else {
                    auto it = p.find("TYPE");
                    if (it != p.end())
                        v = detail::toUpper(detail::trim(it->second));
                }

                if (v == "LINEAR" || v == "NONLINEAR")
                    m.analysisType = v;
                else
                    throw std::runtime_error("*ANALYSIS 未知类型: " + rest);

                auto it = p.find("TOL");
                if (it != p.end()) m.nonlinearTol = std::stod(it->second);
                it = p.find("MAXITER");
                if (it != p.end()) m.nonlinearMaxIter = std::stoi(it->second);
                // STEPS：荷载分级步数（默认 1，即不分级）。
                it = p.find("STEPS");
                if (it != p.end()) m.nonlinearLoadSteps = std::stoi(it->second);
                // LINESEARCH：ON/OFF，ΔU 阻尼回溯（默认 ON）。
                it = p.find("LINESEARCH");
                if (it != p.end()) {
                    const std::string s = detail::toUpper(detail::trim(it->second));
                    if (s == "ON" || s == "1" || s == "TRUE" || s == "YES")
                        m.nonlinearLineSearch = true;
                    else if (s == "OFF" || s == "0" || s == "FALSE" || s == "NO")
                        m.nonlinearLineSearch = false;
                    else
                        throw std::runtime_error("*ANALYSIS LINESEARCH 只能为 ON 或 OFF");
                }
                // THROW：不收敛时是否抛异常（默认 ON）。OFF 则返回部分成功结果。
                it = p.find("THROW");
                if (it != p.end()) {
                    const std::string s = detail::toUpper(detail::trim(it->second));
                    if (s == "ON" || s == "1" || s == "TRUE" || s == "YES")
                        m.nonlinearThrowOnFailure = true;
                    else if (s == "OFF" || s == "0" || s == "FALSE" || s == "NO")
                        m.nonlinearThrowOnFailure = false;
                    else
                        throw std::runtime_error("*ANALYSIS THROW 只能为 ON 或 OFF");
                }
                if (m.nonlinearTol <= 0.0)
                    throw std::runtime_error("*ANALYSIS TOL 必须大于 0");
                if (m.nonlinearMaxIter < 1)
                    throw std::runtime_error("*ANALYSIS MAXITER 必须至少为 1");
                if (m.nonlinearLoadSteps < 1)
                    throw std::runtime_error("*ANALYSIS STEPS 必须至少为 1");
                continue;
            }

            if (kw == "*PRESTRESS") {
                section = kw;
                continue;
            }

            if (kw == "*INFLUENCE" || kw == "*LANE_LOAD") {
                section = kw;
                continue;
            }

            if (kw == "*CE_PARAMETER") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                auto p = detail::parseKeyValues(args);
                auto it = p.find("PENALTY");
                if (it != p.end()) m.cePenalty = std::stod(it->second);
                it = p.find("ALPHA");
                if (it != p.end()) m.cePenalty = std::stod(it->second);
                if (m.cePenalty <= 0.0)
                    throw std::runtime_error("*CE_PARAMETER PENALTY 必须大于 0");
                continue;
            }

            if (kw == "*CE") {
                section = kw;
                continue;
            }

            if (kw == "*CREEP_MODEL" || kw == "*SHRINKAGE_MODEL") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                auto p = detail::parseKeyValues(args);
                if (kw == "*CREEP_MODEL") {
                    auto it = p.find("MODEL"); if (it != p.end()) m.creepModel = detail::toUpper(it->second);
                    it = p.find("PHI0"); if (it != p.end()) m.creepPhi0 = std::stod(it->second);
                    it = p.find("T0"); if (it != p.end()) m.creepT0 = std::stod(it->second);
                    it = p.find("BETA_H"); if (it != p.end()) m.creepBetaH = std::stod(it->second);
                    it = p.find("POWER"); if (it != p.end()) m.creepPower = std::stod(it->second);
                } else {
                    auto it = p.find("MODEL"); if (it != p.end()) m.shrinkageModel = detail::toUpper(it->second);
                    it = p.find("ECS_INF"); if (it != p.end()) m.shrinkageEcsInf = std::stod(it->second);
                    it = p.find("TS"); if (it != p.end()) m.shrinkageTs = std::stod(it->second);
                    it = p.find("ALPHA_DS"); if (it != p.end()) m.shrinkageAlphaDs = std::stod(it->second);
                    it = p.find("POWER"); if (it != p.end()) m.shrinkagePower = std::stod(it->second);
                }
                continue;
            }

            if (kw == "*CREEP" || kw == "*SHRINKAGE") {
                section = kw;
                continue;
            }

            if (kw == "*TEMP_PROFILE") {
                section = kw;
                continue;
            }

            if (kw == "*STAGE" || kw == "*STAGE_LOAD" || kw == "*STAGE_ELEMENT_LOAD" ||
                kw == "*STAGE_PRESTRESS" || kw == "*STAGE_BOUNDARY" || kw == "*STAGE_UNBOUNDARY" ||
                kw == "*STAGE_CREEP" || kw == "*STAGE_SHRINKAGE" || kw == "*STAGE_UNLOAD" ||
                kw == "*STAGE_TARGET" || kw == "*STAGE_CONTROL" || kw == "*STAGE_TEMP_PROFILE" ||
                kw == "*CREEP_CURVE" || kw == "*SHRINKAGE_CURVE" || kw == "*STAGE_CONSTRAINT" ||
                kw == "*CAMBER" || kw == "*FEEDBACK" || kw == "*CREEP_MODEL" || kw == "*SHRINKAGE_MODEL") {
                section = kw;
                continue;
            }

            section = kw;

            if (kw == "*ELEMENT") {
                // 旧式：*ELEMENT, TYPE=..., E=..., A=..., I=..., ALPHA=..., H=..., WT=..., WB=...
                auto p = detail::parseParams(rest);
                legacyTemplate = Model2D::Element();
                legacyTemplate.type = detail::toUpper(p["TYPE"]);
                legacyTemplate.E = detail::toDouble(p, "E");
                legacyTemplate.A = detail::toDouble(p, "A");
                legacyTemplate.I = 0.0;
                legacyTemplate.alpha = 0.0;
                legacyTemplate.h = 0.0;
                legacyTemplate.Wt = 0.0;
                legacyTemplate.Wb = 0.0;
                legacyTemplate.conn.clear();
                auto it = p.find("I");
                if (it != p.end()) legacyTemplate.I = std::stod(it->second);
                it = p.find("ALPHA");
                if (it != p.end()) legacyTemplate.alpha = std::stod(it->second);
                it = p.find("H");
                if (it != p.end()) legacyTemplate.h = std::stod(it->second);
                it = p.find("WT");
                if (it != p.end()) legacyTemplate.Wt = std::stod(it->second);
                it = p.find("WB");
                if (it != p.end()) legacyTemplate.Wb = std::stod(it->second);
            }
            continue;
        }

        // ---- 数据行，按当前 section 解析 ----
        std::istringstream ss(t);
        if (section == "*TITLE") {
            m.title = t;
        } else if (section == "*NODE") {
            int id = 0;
            double x = 0.0, y = 0.0;
            if (!(ss >> id >> x >> y))
                throw std::runtime_error("NODE 行格式错误: " + t);
            m.nodes[id] = {x, y};
        } else if (section == "*ELEMENT") {
            int eid = 0, a = 0, b = 0;
            if (!(ss >> eid >> a >> b))
                throw std::runtime_error("ELEMENT 行格式错误: " + t);
            Model2D::Element e = legacyTemplate;
            e.id = eid;
            e.conn = {a, b};
            if (e.type == "BEAM" && e.I <= 0.0)
                throw std::runtime_error("旧式 *ELEMENT：梁单元缺少 I");
            if (e.type == "BEAM" && (e.Wt <= 0.0 || e.Wb <= 0.0))
                throw std::runtime_error("旧式 *ELEMENT：梁单元缺少 WT/WB（截面抵抗矩，应力计算必需）");
            if (isAxial(e.type) && e.A <= 0.0)
                throw std::runtime_error("旧式 *ELEMENT：桁架/索/压杆单元缺少 A");
            m.elements.push_back(e);
        } else if (section == "*ET") {
            int id = 0;
            std::string kind;
            if (!(ss >> id >> kind))
                throw std::runtime_error("ET 行格式错误: " + t);
            m.elementTypes[id] = detail::toUpper(kind);
        } else if (section == "*MP") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("MP 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Material mat;
            auto it = p.find("E"); if (it != p.end()) mat.E = std::stod(it->second);
            it = p.find("ALPX"); if (it != p.end()) mat.alpha = std::stod(it->second);
            it = p.find("DENS"); if (it != p.end()) mat.dens = std::stod(it->second);
            m.materials[id] = mat;
        } else if (section == "*R") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("R 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::RealConst rc;
            auto it = p.find("A"); if (it != p.end()) rc.A = std::stod(it->second);
            it = p.find("I"); if (it != p.end()) rc.I = std::stod(it->second);
            it = p.find("WT"); if (it != p.end()) rc.Wt = std::stod(it->second);
            it = p.find("WB"); if (it != p.end()) rc.Wb = std::stod(it->second);
            it = p.find("H"); if (it != p.end()) rc.h = std::stod(it->second);
            it = p.find("ADDMASS"); if (it != p.end()) rc.addmass = std::stod(it->second);
            it = p.find("ISTRN"); if (it != p.end()) rc.istrn = std::stod(it->second);
            m.realConsts[id] = rc;
        } else if (section == "*BOUNDARY") {
            std::string ns, ds, vs;
            if (!(ss >> ns >> ds >> vs))
                throw std::runtime_error("BOUNDARY 行格式错误: " + t);
            const double value = std::stod(vs);
            for (int node : detail::parseRange(ns))
                for (int dof : detail::parseRange(ds))
                    m.constraints.push_back({node, dof, value});
        } else if (section == "*CLOAD") {
            std::string ns, ds, vs;
            if (!(ss >> ns >> ds >> vs))
                throw std::runtime_error("CLOAD 行格式错误: " + t);
            const double value = std::stod(vs);
            for (int node : detail::parseRange(ns))
                for (int dof : detail::parseRange(ds))
                    m.loads.push_back({node, dof, value});
        } else if (section == "*ELEMENT_LOAD") {
            std::string es, dir;
            double q1 = 0.0, q2 = 0.0;
            if (!(ss >> es >> dir >> q1 >> q2))
                throw std::runtime_error("ELEMENT_LOAD 行格式错误: " + t);
            for (int elem : detail::parseRange(es))
                m.elemLoads.push_back({elem, dir, q1, q2});
        } else if (section == "*GRAVITY") {
            double gx = 0.0, gy = 0.0;
            if (!(ss >> gx >> gy))
                throw std::runtime_error("GRAVITY 行格式错误: " + t);
            m.gx = gx;
            m.gy = gy;
        } else if (section == "*SELF_WEIGHT") {
            std::string who, rule;
            double value = 0.0;
            if (!(ss >> who >> rule >> value))
                throw std::runtime_error("SELF_WEIGHT 行格式错误: " + t);
            Model2D::SelfWeight sw;
            sw.rule = detail::toUpper(rule);
            sw.value = value;
            if (detail::toUpper(who) == "ALL") sw.all = true;
            else sw.elem = std::stoi(who);
            m.selfWeights.push_back(sw);
        } else if (section == "*TEMPERATURE") {
            double dT = 0.0;
            if (!(ss >> dT))
                throw std::runtime_error("TEMPERATURE 行格式错误: " + t);
            m.uniformTemp = dT;
        } else if (section == "*ELEMENT_TEMPERATURE") {
            std::string es, kind;
            double dT = 0.0;
            if (!(ss >> es >> kind >> dT))
                throw std::runtime_error("ELEMENT_TEMPERATURE 行格式错误: " + t);
            for (int elem : detail::parseRange(es))
                m.elemTemps.push_back({elem, detail::toUpper(kind), dT});
        } else if (section == "*PRESTRESS") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("PRESTRESS 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Prestress ps;
            ps.id = id;
            auto it = p.find("ELEM");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS 缺少 ELEM");
            ps.elems = detail::parseIntList(it->second);
            it = p.find("P");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS 缺少 P");
            ps.P = std::stod(it->second);
            it = p.find("F"); if (it != p.end()) ps.F = std::stod(it->second);
            it = p.find("E1"); if (it != p.end()) ps.E1 = std::stod(it->second);
            it = p.find("E2"); if (it != p.end()) ps.E2 = std::stod(it->second);
            it = p.find("LOSS");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS / *STAGE_PRESTRESS 缺少 LOSS（损失系数，必须显式指定）");
            ps.loss = std::stod(it->second);
            it = p.find("EPS"); if (it != p.end()) ps.Ep = std::stod(it->second);
            it = p.find("AP"); if (it != p.end()) ps.Ap = std::stod(it->second);
            it = p.find("APS"); if (it != p.end()) ps.Ap = std::stod(it->second);
            it = p.find("METHOD"); if (it != p.end()) ps.method = detail::toUpper(it->second);
            if (ps.P <= 0.0 || ps.loss <= 0.0 || ps.Ep <= 0.0)
                throw std::runtime_error("*PRESTRESS 的 P/LOSS/EPS 必须大于 0");
            if (ps.method != "AUTO" && ps.method != "EQUIV" && ps.method != "TENDON")
                throw std::runtime_error("*PRESTRESS 未知 METHOD: " + ps.method);
            m.prestresses.push_back(ps);
        } else if (section == "*INFLUENCE") {
            int id = 0;
            std::string kind;
            if (!(ss >> id >> kind))
                throw std::runtime_error("INFLUENCE 行格式错误: " + t);
            Model2D::InfluenceRequest ir;
            ir.id = id;
            ir.kind = detail::toUpper(kind);
            if (ir.kind == "DISP" || ir.kind == "REACTION") {
                if (!(ss >> ir.node >> ir.dof))
                    throw std::runtime_error("INFLUENCE DISP/REACTION 需要 node dof");
            } else if (ir.kind == "FORCE" || ir.kind == "STRESS") {
                std::string comp;
                if (!(ss >> ir.elem >> comp >> ir.end))
                    throw std::runtime_error("INFLUENCE FORCE/STRESS 需要 elem component end");
                ir.component = detail::toUpper(comp);
            } else {
                throw std::runtime_error("INFLUENCE 未知类型: " + kind);
            }
            if (ir.end != 1 && ir.end != 2)
                throw std::runtime_error("INFLUENCE end 只能为 1 或 2");
            m.influences.push_back(ir);
        } else if (section == "*LANE_LOAD") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("LANE_LOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::LaneLoad ll;
            ll.id = id;
            auto it = p.find("NODES");
            if (it == p.end()) throw std::runtime_error("*LANE_LOAD 缺少 NODES");
            ll.nodes = detail::parseIntList(it->second);
            it = p.find("Q"); if (it != p.end()) ll.q = std::stod(it->second);
            it = p.find("P"); if (it != p.end()) ll.P = std::stod(it->second);
            it = p.find("DIR"); if (it != p.end()) ll.dir = detail::toUpper(it->second);
            it = p.find("METHOD"); if (it != p.end()) ll.method = detail::toUpper(it->second);
            if (ll.dir != "GX" && ll.dir != "GY")
                throw std::runtime_error("*LANE_LOAD DIR 只能为 GX 或 GY");
            if (ll.method != "LINEAR" && ll.method != "NONLINEAR")
                throw std::runtime_error("*LANE_LOAD METHOD 只能为 LINEAR 或 NONLINEAR");
            m.laneLoads.push_back(ll);
        } else if (section == "*CE") {
            Model2D::ConstraintEquation ce;
            if (!(ss >> ce.id >> ce.constant))
                throw std::runtime_error("CE 行格式错误: " + t);
            int node = 0, dof = 0;
            double coef = 0.0;
            while (ss >> node >> dof >> coef) {
                if (dof < 1 || dof > 3)
                    throw std::runtime_error("CE 自由度必须在 1~3 之间");
                ce.terms.push_back({node, dof, coef});
            }
            if (ce.terms.empty())
                throw std::runtime_error("CE 至少需要一个自由度项");
            m.ces.push_back(ce);
        } else if (section == "*CREEP") {
            std::string who;
            if (!(ss >> who))
                throw std::runtime_error("CREEP 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::CreepSpec cs;
            if (detail::toUpper(who) == "ALL") cs.all = true;
            else cs.elems = detail::parseIntList(who);
            auto it = p.find("PHI");
            if (it == p.end()) throw std::runtime_error("*CREEP 缺少 PHI");
            cs.phi = std::stod(it->second);
            it = p.find("CHI"); if (it != p.end()) cs.chi = std::stod(it->second);
            if (cs.phi < 0.0 || cs.chi <= 0.0)
                throw std::runtime_error("*CREEP PHI/CHI 参数错误");
            m.creepSpecs.push_back(cs);
        } else if (section == "*SHRINKAGE") {
            std::string who;
            if (!(ss >> who))
                throw std::runtime_error("SHRINKAGE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::ShrinkageSpec ss_;
            if (detail::toUpper(who) == "ALL") ss_.all = true;
            else ss_.elems = detail::parseIntList(who);
            auto it = p.find("ECS");
            if (it == p.end()) throw std::runtime_error("*SHRINKAGE 缺少 ECS");
            ss_.ecs = std::stod(it->second);
            m.shrinkageSpecs.push_back(ss_);
        } else if (section == "*TEMP_PROFILE") {
            std::string who;
            if (!(ss >> who))
                throw std::runtime_error("TEMP_PROFILE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::TempProfileSpec tp;
            if (detail::toUpper(who) == "ALL") tp.all = true;
            else tp.elems = detail::parseIntList(who);
            auto it = p.find("T1");
            if (it != p.end()) tp.T1 = std::stod(it->second);
            it = p.find("T2");
            if (it != p.end()) tp.T2 = std::stod(it->second);
            it = p.find("A"); if (it != p.end()) tp.A = std::stod(it->second);
            it = p.find("B"); if (it != p.end()) tp.B = std::stod(it->second);
            it = p.find("H"); if (it != p.end()) tp.H = std::stod(it->second);
            if (tp.A <= 0.0)
                throw std::runtime_error("*TEMP_PROFILE A 必须大于 0");
            m.tempProfiles.push_back(tp);
        } else if (section == "*STAGE") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("STAGE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Stage st;
            st.id = id;
            auto it = p.find("ACTIVE");
            if (it != p.end()) st.activate = detail::parseIntList(it->second);
            it = p.find("DEACTIVE");
            if (it != p.end()) st.deactivate = detail::parseIntList(it->second);
            it = p.find("TIME");
            if (it != p.end()) st.time = std::stod(it->second);
            m.stages.push_back(st);
        } else if (section == "*STAGE_LOAD") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_LOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            auto itV = p.find("VALUE");
            if (itN == p.end() || itD == p.end() || itV == p.end())
                throw std::runtime_error("*STAGE_LOAD 需要 NODE/DOF/VALUE");
            for (int node : detail::parseIntList(itN->second))
                for (int dof : detail::parseIntList(itD->second))
                    m.stages[stageId - 1].loads.push_back({node, dof, std::stod(itV->second)});
        } else if (section == "*STAGE_ELEMENT_LOAD") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_ELEMENT_LOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itE = p.find("ELEM");
            auto itDir = p.find("DIR");
            auto itQ1 = p.find("Q1");
            auto itQ2 = p.find("Q2");
            if (itE == p.end() || itDir == p.end() || itQ1 == p.end() || itQ2 == p.end())
                throw std::runtime_error("*STAGE_ELEMENT_LOAD 需要 ELEM/DIR/Q1/Q2");
            for (int elem : detail::parseIntList(itE->second))
                m.stages[stageId - 1].elemLoads.push_back({elem, itDir->second, std::stod(itQ1->second), std::stod(itQ2->second)});
        } else if (section == "*STAGE_PRESTRESS") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_PRESTRESS 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Prestress ps;
            auto it = p.find("ELEM");
            if (it == p.end()) throw std::runtime_error("*STAGE_PRESTRESS 缺少 ELEM");
            ps.elems = detail::parseIntList(it->second);
            it = p.find("P");
            if (it == p.end()) throw std::runtime_error("*STAGE_PRESTRESS 缺少 P");
            ps.P = std::stod(it->second);
            it = p.find("F"); if (it != p.end()) ps.F = std::stod(it->second);
            it = p.find("E1"); if (it != p.end()) ps.E1 = std::stod(it->second);
            it = p.find("E2"); if (it != p.end()) ps.E2 = std::stod(it->second);
            it = p.find("LOSS");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS / *STAGE_PRESTRESS 缺少 LOSS（损失系数，必须显式指定）");
            ps.loss = std::stod(it->second);
            it = p.find("EPS"); if (it != p.end()) ps.Ep = std::stod(it->second);
            it = p.find("APS"); if (it != p.end()) ps.Ap = std::stod(it->second);
            ps.method = "EQUIV"; // 施工阶段第一版只使用等效荷载法
            m.stages[stageId - 1].prestresses.push_back(ps);
        } else if (section == "*STAGE_BOUNDARY" || section == "*STAGE_UNBOUNDARY") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_BOUNDARY/UNBOUNDARY 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            if (itN == p.end() || itD == p.end())
                throw std::runtime_error("*STAGE_BOUNDARY/UNBOUNDARY 需要 NODE/DOF");
            const auto nodes = detail::parseIntList(itN->second);
            const auto dofs = detail::parseIntList(itD->second);
            for (int node : nodes) {
                for (int dof : dofs) {
                    if (section == "*STAGE_BOUNDARY") {
                        double value = 0.0;
                        auto itV = p.find("VALUE");
                        if (itV != p.end()) value = std::stod(itV->second);
                        m.stages[stageId - 1].addedConstraints.push_back({node, dof, value});
                    } else {
                        m.stages[stageId - 1].removedConstraints.push_back({node, dof});
                    }
                }
            }
        } else if (section == "*STAGE_CREEP" || section == "*STAGE_SHRINKAGE") {
            int stageId = 0;
            std::string who;
            if (!(ss >> stageId >> who))
                throw std::runtime_error("STAGE_CREEP/SHRINKAGE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            if (section == "*STAGE_CREEP") {
                Model2D::CreepSpec cs;
                if (detail::toUpper(who) == "ALL") cs.all = true;
                else cs.elems = detail::parseIntList(who);
                auto it = p.find("PHI");
                if (it == p.end()) throw std::runtime_error("*STAGE_CREEP 缺少 PHI");
                cs.phi = std::stod(it->second);
                it = p.find("CHI"); if (it != p.end()) cs.chi = std::stod(it->second);
                m.stages[stageId - 1].creepSpecs.push_back(cs);
            } else {
                Model2D::ShrinkageSpec sh;
                if (detail::toUpper(who) == "ALL") sh.all = true;
                else sh.elems = detail::parseIntList(who);
                auto it = p.find("ECS");
                if (it == p.end()) throw std::runtime_error("*STAGE_SHRINKAGE 缺少 ECS");
                sh.ecs = std::stod(it->second);
                m.stages[stageId - 1].shrinkageSpecs.push_back(sh);
            }
        } else if (section == "*STAGE_TEMP_PROFILE") {
            int stageId = 0;
            std::string who;
            if (!(ss >> stageId >> who))
                throw std::runtime_error("STAGE_TEMP_PROFILE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::TempProfileSpec tp;
            if (detail::toUpper(who) == "ALL") tp.all = true;
            else tp.elems = detail::parseIntList(who);
            auto it = p.find("T1"); if (it != p.end()) tp.T1 = std::stod(it->second);
            it = p.find("T2"); if (it != p.end()) tp.T2 = std::stod(it->second);
            it = p.find("A"); if (it != p.end()) tp.A = std::stod(it->second);
            it = p.find("B"); if (it != p.end()) tp.B = std::stod(it->second);
            it = p.find("H"); if (it != p.end()) tp.H = std::stod(it->second);
            if (tp.A <= 0.0)
                throw std::runtime_error("*STAGE_TEMP_PROFILE A 必须大于 0");
            m.stages[stageId - 1].tempProfiles.push_back(tp);
        } else if (section == "*CREEP_CURVE" || section == "*SHRINKAGE_CURVE") {
            double time = 0.0, value = 0.0;
            if (!(ss >> time >> value))
                throw std::runtime_error("CREEP_CURVE/SHRINKAGE_CURVE 行格式错误: " + t);
            if (section == "*CREEP_CURVE") m.creepCurve.points.push_back({time, value});
            else m.shrinkageCurve.points.push_back({time, value});
        } else if (section == "*CREEP_MODEL" || section == "*SHRINKAGE_MODEL") {
            auto p = detail::parseKeyValues(t);
            if (section == "*CREEP_MODEL") {
                auto it = p.find("MODEL"); if (it != p.end()) m.creepModel = detail::toUpper(it->second);
                it = p.find("PHI0"); if (it != p.end()) m.creepPhi0 = std::stod(it->second);
                it = p.find("T0"); if (it != p.end()) m.creepT0 = std::stod(it->second);
                it = p.find("BETA_H"); if (it != p.end()) m.creepBetaH = std::stod(it->second);
                it = p.find("POWER"); if (it != p.end()) m.creepPower = std::stod(it->second);
            } else {
                auto it = p.find("MODEL"); if (it != p.end()) m.shrinkageModel = detail::toUpper(it->second);
                it = p.find("ECS_INF"); if (it != p.end()) m.shrinkageEcsInf = std::stod(it->second);
                it = p.find("TS"); if (it != p.end()) m.shrinkageTs = std::stod(it->second);
                it = p.find("ALPHA_DS"); if (it != p.end()) m.shrinkageAlphaDs = std::stod(it->second);
                it = p.find("POWER"); if (it != p.end()) m.shrinkagePower = std::stod(it->second);
            }
        } else if (section == "*STAGE_UNLOAD") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_UNLOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            auto itE = p.find("ELEM");
            if (itN != p.end() && itD != p.end()) {
                for (int node : detail::parseIntList(itN->second))
                    for (int dof : detail::parseIntList(itD->second))
                        m.stages[stageId - 1].unloadNodeDofs.push_back({node, dof});
            } else if (itE != p.end()) {
                const auto elems = detail::parseIntList(itE->second);
                m.stages[stageId - 1].unloadElems.insert(m.stages[stageId - 1].unloadElems.end(),
                                                        elems.begin(), elems.end());
            } else {
                throw std::runtime_error("*STAGE_UNLOAD 需要 NODE/DOF 或 ELEM");
            }
        } else if (section == "*STAGE_TARGET") {
            std::string rest2 = t;
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            auto itV = p.find("VALUE");
            Model2D::StageTarget tg;
            if (itV == p.end())
                throw std::runtime_error("*STAGE_TARGET 需要 VALUE");
            tg.value = std::stod(itV->second);
            auto itW = p.find("WEIGHT");
            if (itW != p.end()) tg.weight = std::stod(itW->second);
            auto itK = p.find("TYPE");
            tg.kind = itK != p.end() ? detail::toUpper(itK->second) : "DISP";
            if (tg.kind == "DISP" || tg.kind == "REACTION") {
                if (itN == p.end() || itD == p.end())
                    throw std::runtime_error("*STAGE_TARGET DISP/REACTION 需要 NODE/DOF");
                tg.node = std::stoi(itN->second);
                tg.dof = std::stoi(itD->second);
            } else if (tg.kind == "FORCE" || tg.kind == "STRESS") {
                auto itE = p.find("ELEM");
                auto itC = p.find("COMPONENT");
                auto itEnd = p.find("END");
                if (itE == p.end() || itC == p.end())
                    throw std::runtime_error("*STAGE_TARGET FORCE/STRESS 需要 ELEM/COMPONENT");
                tg.elem = std::stoi(itE->second);
                tg.component = detail::toUpper(itC->second);
                if (itEnd != p.end()) tg.end = std::stoi(itEnd->second);
            } else {
                throw std::runtime_error("*STAGE_TARGET 未知 TYPE: " + tg.kind);
            }
            m.stageTargets.push_back(tg);
        } else if (section == "*STAGE_CONTROL") {
            auto p = detail::parseKeyValues(t);
            auto itS = p.find("STAGE");
            auto itP = p.find("PRESTRESS");
            if (itS == p.end() || itP == p.end())
                throw std::runtime_error("*STAGE_CONTROL 需要 STAGE/PRESTRESS");
            Model2D::StageControl ctrl;
            ctrl.stage = std::stoi(itS->second);
            ctrl.prestressIndex = std::stoi(itP->second);
            auto itL = p.find("LOWER");
            if (itL != p.end()) ctrl.lower = std::stod(itL->second);
            auto itU = p.find("UPPER");
            if (itU != p.end()) ctrl.upper = std::stod(itU->second);
            if (ctrl.lower > ctrl.upper)
                throw std::runtime_error("*STAGE_CONTROL LOWER 不能大于 UPPER");
            m.stageControls.push_back(ctrl);
        } else if (section == "*STAGE_CONSTRAINT") {
            auto p = detail::parseKeyValues(t);
            Model2D::StageConstraint sc;
            auto itV = p.find("VALUE");
            if (itV != p.end()) sc.target.value = std::stod(itV->second);
            auto itW = p.find("WEIGHT");
            if (itW != p.end()) sc.target.weight = std::stod(itW->second);
            auto itK = p.find("TYPE");
            sc.target.kind = itK != p.end() ? detail::toUpper(itK->second) : "DISP";
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            if (sc.target.kind == "DISP" || sc.target.kind == "REACTION") {
                if (itN == p.end() || itD == p.end())
                    throw std::runtime_error("*STAGE_CONSTRAINT DISP/REACTION 需要 NODE/DOF");
                sc.target.node = std::stoi(itN->second);
                sc.target.dof = std::stoi(itD->second);
            } else if (sc.target.kind == "FORCE" || sc.target.kind == "STRESS") {
                auto itE = p.find("ELEM");
                auto itC = p.find("COMPONENT");
                auto itEnd = p.find("END");
                if (itE == p.end() || itC == p.end())
                    throw std::runtime_error("*STAGE_CONSTRAINT FORCE/STRESS 需要 ELEM/COMPONENT");
                sc.target.elem = std::stoi(itE->second);
                sc.target.component = detail::toUpper(itC->second);
                if (itEnd != p.end()) sc.target.end = std::stoi(itEnd->second);
            } else {
                throw std::runtime_error("*STAGE_CONSTRAINT 未知 TYPE");
            }
            auto itL = p.find("LOWER");
            if (itL != p.end()) sc.lower = std::stod(itL->second);
            auto itU = p.find("UPPER");
            if (itU != p.end()) sc.upper = std::stod(itU->second);
            m.stageConstraints.push_back(sc);
        } else if (section == "*CAMBER") {
            int node = 0, dof = 0;
            double value = 0.0;
            if (!(ss >> node >> dof))
                throw std::runtime_error("CAMBER 行格式错误: " + t);
            ss >> value;
            if (dof < 1 || dof > 3)
                throw std::runtime_error("CAMBER dof 必须在 1~3 之间");
            m.cambers.push_back({node, dof, value});
        } else if (section == "*FEEDBACK") {
            int node = 0, dof = 0;
            double measured = 0.0;
            if (!(ss >> node >> dof >> measured))
                throw std::runtime_error("FEEDBACK 行格式错误: " + t);
            m.feedbacks.push_back({node, dof, measured});
        } else if (section == "*END") {
            break;
        }
    }

    if (m.nodes.empty())
        throw std::runtime_error("输入文件没有定义任何节点");
    processPrestress(m);
    // 先记录未折减的基准模量，再做徐变折减。
    // 之后每一次徐变处理都从这个基准重算，保证幂等。
    snapshotBaseModulus(m);
    processTimeEffects(m);

    // 危险默认值检查：若存在重力/自重荷载但单元密度为 0，给出警告。
    const bool hasGravity = (m.gx != 0.0 || m.gy != 0.0 || !m.selfWeights.empty());
    if (hasGravity) {
        for (const auto& e : m.elements) {
            if (e.dens <= 0.0 && e.addmass <= 0.0 && !isRigid(e.type)) {
                std::cerr << "[警告] 单元 " << e.id << " 的密度 dens=0 且 addmass=0，"
                          << "但模型存在 *GRAVITY 或 *SELF_WEIGHT，该单元将不产生自重。"
                          << "若需自重，请在 *MAT 或 *ELEMENT 中指定 DENS。\n";
            }
        }
    }
    return m;
}

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

inline void appendConstraintPenaltyTriplets(const Model2D& m,
                                            const std::map<int, int>& index,
                                            int nDofs,
                                            double penalty,
                                            std::vector<Eigen::Triplet<double>>& triplets) {
    for (const auto& ce : m.ces) {
        std::map<int, double> c;
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            c[d] += term.coef;
        }
        for (const auto& [i, ci] : c) {
            for (const auto& [j, cj] : c) {
                const double v = penalty * ci * cj;
                if (v != 0.0)
                    triplets.emplace_back(i, j, v);
            }
        }
    }
}

// 若用户未显式指定罚系数（<=0），按单元刚度最大对角元的 1e8 倍自适应取值，
// 避免固定 1e12 在细长杆/小刚度模型中条件数恶化。
inline double effectiveCePenalty(const std::vector<Eigen::Triplet<double>>& triplets,
                                 double userPenalty) {
    if (userPenalty > 0.0) return userPenalty;
    double maxDiag = 0.0;
    for (const auto& t : triplets) {
        if (t.row() == t.col())
            maxDiag = std::max(maxDiag, std::abs(t.value()));
    }
    if (maxDiag <= 0.0) maxDiag = 1.0;
    return 1e8 * maxDiag;
}

inline DenseVec constraintPenaltyConstantForce(const Model2D& m,
                                               const std::map<int, int>& index,
                                               int nDofs,
                                               double penalty) {
    DenseVec Fce = DenseVec::Zero(nDofs);
    for (const auto& ce : m.ces) {
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            Fce(d) += penalty * ce.constant * term.coef;
        }
    }
    return Fce;
}

inline DenseVec constraintPenaltyResidualForce(const Model2D& m,
                                               const std::map<int, int>& index,
                                               int nDofs,
                                               double penalty,
                                               const DenseVec& U) {
    DenseVec Fce = DenseVec::Zero(nDofs);
    for (const auto& ce : m.ces) {
        double g = -ce.constant;
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            g += term.coef * U(d);
        }
        const double a = penalty * g;
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            Fce(d) += a * term.coef;
        }
    }
    return Fce;
}

// 自动约束"无转动刚度的转角自由度"。
// 二维统一按每节点 3 个自由度 [u, v, θ] 编号，但 TRUSS/CABLE/STRUT 只贡献
// 平动刚度。纯桁架节点的 θ 对角元恒为 0，若不由用户显式约束，总体矩阵必然
// 奇异，且报错信息对使用者很不友好。这里自动补零约束（等价于 ANSYS 等软件
// 对未用自由度的处理），平移自由度仍然保持严格的零刚度检查。
inline void autoConstrainFreeRotations(const SparseMat& K, int nDofs,
                                       const std::vector<int>& baseCdofs,
                                       const std::vector<double>& baseCvals,
                                       std::vector<int>& cdofs,
                                       std::vector<double>& cvals) {
    cdofs = baseCdofs;
    cvals = baseCvals;
    std::vector<char> isC(nDofs, 0);
    for (int d : cdofs) isC[d] = 1;
    for (int d = 2; d < nDofs; d += 3) {   // 每节点的第 3 个自由度 = 转角 θ
        if (isC[d]) continue;
        if (K.coeff(d, d) == 0.0) {
            cdofs.push_back(d);
            cvals.push_back(0.0);
            isC[d] = 1;
        }
    }
}

// 方案 C：几何非线性轴向单元 + 线性欧拉梁的 Newton-Raphson 求解。
// 当前版本非线性部分覆盖 TRUSS/CABLE/STRUT；BEAM 仍使用线性切线。
inline ModelResult solveModel2DNonlinear(const Model2D& m,
                                         SolverType type = SolverType::SimplicialLDLT,
                                         const DenseVec* U0 = nullptr,
                                         const std::vector<char>* active0 = nullptr,
                                         NonlinearDiag* diagOut = nullptr,
                                         ReusableSolver* reuse = nullptr) {
    std::vector<int> ids;
    std::map<int, int> index;
    for (const auto& [id, c] : m.nodes) {
        index[id] = static_cast<int>(ids.size());
        ids.push_back(id);
    }
    const int nNodes = static_cast<int>(ids.size());
    const int nDofs = 3 * nNodes;

    // 约束自由度与指定位移。
    std::vector<int> cdofs;
    std::vector<double> cvals;
    std::vector<char> isConstrained(nDofs, 0);
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        const int d = 3 * index.at(c.node) + (c.dof - 1);
        cdofs.push_back(d);
        cvals.push_back(c.value);
        isConstrained[d] = 1;
    }

    DenseVec U = DenseVec::Zero(nDofs);
    if (U0) {
        if (U0->size() != nDofs)
            throw std::runtime_error("施工阶段初始位移向量长度不匹配");
        U = *U0;
    }
    for (std::size_t i = 0; i < cdofs.size(); ++i)
        U(cdofs[i]) = cvals[i];

    // cable/strut 初始全部激活，后续根据当前构形切换。
    std::vector<char> active(m.elements.size(), 1);
    if (active0 && active0->size() == active.size())
        active = *active0;

    const double tol = m.nonlinearTol;
    const int maxIter = m.nonlinearMaxIter;

    auto axialActiveTol = [](const Model2D::Element& e) {
        return 1e-12 * std::max(1.0, std::abs(e.E * e.A));
    };

    // 荷载分级比例：assembleState 只把**外荷载**乘以它。
    // 1.0 = 全部荷载（默认，与旧行为一致）。
    double loadScale = 1.0;

    // 把指定位移"钉"回向量里。荷载分级/线搜索期间每一步都要保证约束成立。
    auto applyPrescribed = [&](DenseVec& V) {
        for (std::size_t i = 0; i < cdofs.size(); ++i)
            V(cdofs[i]) = cvals[i];
    };

    auto assembleState = [&](const std::vector<char>& act, const DenseVec& Ucur,
                             DenseVec& Fint, DenseVec& Fext,
                             std::vector<Eigen::Triplet<double>>* triplets) {
        Fint = DenseVec::Zero(nDofs);
        Fext = DenseVec::Zero(nDofs);
        if (triplets) triplets->clear();

        // 节点集中荷载始终施加（按当前荷载分级比例缩放）。
        for (const auto& l : m.loads) {
            if (l.dof < 1 || l.dof > 3)
                throw std::runtime_error("载荷自由度必须在 1~3 之间");
            Fext(3 * index.at(l.node) + (l.dof - 1)) += loadScale * l.value;
        }

        for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
            if (!act[ei]) continue;
            const auto& e = m.elements[ei];
            if (e.conn.size() != 2)
                throw std::runtime_error("二维单元必须有两个节点");
            const int i1 = index.at(e.conn[0]);
            const int i2 = index.at(e.conn[1]);
            const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
            const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

            if (isAxial(e.type)) {
                const std::vector<int> g = elementGdofs2D(e.type, i1, i2);
                DenseVec fext = elementExternalEquivalentLoad(m, e, p1, p2);
                fext *= loadScale;
                const DenseVec fint = axialInternalForceNonlinear(e, Ucur, i1, i2, p1, p2);
                addElementLoad(g, fext, Fext);
                addElementLoad(g, fint, Fint);
                if (triplets) {
                    const DenseMat Ke = axialTangentNonlinear(e, Ucur, i1, i2, p1, p2);
                    scatter(Ke, g, *triplets);
                }
            } else if (e.type == "RIGID_PIN") {
                const std::vector<int> g = elementGdofs2D(e.type, i1, i2);
                Eigen::VectorXd ug(6);
                ug << Ucur(3 * i1 + 0), Ucur(3 * i1 + 1), Ucur(3 * i1 + 2),
                      Ucur(3 * i2 + 0), Ucur(3 * i2 + 1), Ucur(3 * i2 + 2);
                const DenseMat Ke = rigidPinStiffness(p1, p2);
                addElementLoad(g, Ke * ug, Fint);
                addElementLoad(g, DenseVec::Zero(6), Fext);
                if (triplets) scatter(Ke, g, *triplets);
            } else if (e.type == "BEAM" || e.type == "RIGID") {
                EulerBeam2D be{e.E, e.A, e.I, p1, p2};
                const std::vector<int> g = elementGdofs2D(e.type, i1, i2);
                Eigen::VectorXd ug(6);
                ug << Ucur(3 * i1 + 0), Ucur(3 * i1 + 1), Ucur(3 * i1 + 2),
                      Ucur(3 * i2 + 0), Ucur(3 * i2 + 1), Ucur(3 * i2 + 2);
                const DenseMat Ke = be.stiffness();
                const DenseVec fint = Ke * ug;
                DenseVec fext = elementExternalEquivalentLoad(m, e, p1, p2);
                fext *= loadScale;
                addElementLoad(g, fext, Fext);
                addElementLoad(g, fint, Fint);
                if (triplets) scatter(Ke, g, *triplets);
            } else {
                throw std::runtime_error("未知单元类型: " + e.type);
            }
        }
    };

    DenseVec Fint, Fext;
    NonlinearDiag diag;
    std::vector<double> zeroPrescribed(cdofs.size(), 0.0);
    bool haveInitialResidual = false;

    // 给定构形与激活状态，返回自由自由度的残差范数（只组装内力/外力，不算切线）。
    // 线搜索每一步都用它判断"这一步有没有让残差下降"。
    auto residualNorm = [&](const std::vector<char>& act, const DenseVec& Ucur) -> double {
        DenseVec fi, fe;
        assembleState(act, Ucur, fi, fe, nullptr);
        DenseVec rr = fe - fi;
        if (!m.ces.empty())
            rr -= constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, Ucur);
        double s = 0.0;
        for (int d = 0; d < nDofs; ++d)
            if (!isConstrained[d]) s += rr(d) * rr(d);
        return std::sqrt(s);
    };

    // 依据当前构形更新 cable/strut 的单向激活状态。
    auto updateActive = [&](const DenseVec& Ucur) -> std::vector<char> {
        std::vector<char> na(m.elements.size(), 1);
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
            const auto& e = m.elements[ei];
            if (e.type != "CABLE" && e.type != "STRUT") continue;
            const int i1 = index.at(e.conn[0]);
            const int i2 = index.at(e.conn[1]);
            const double N = axialForceNonlinear(e, Ucur, i1, i2,
                                                 m.nodes.at(e.conn[0]),
                                                 m.nodes.at(e.conn[1]));
            const double atol = axialActiveTol(e);
            na[ei] = (e.type == "CABLE") ? (N > atol) : (N < -atol);
        }
        return na;
    };

    // 线性路径已经会为"只连桁架的节点"自动约束多余转角自由度（P1-3 的二维版）。
    // 非线性路径原先没有这一步，纯桁架/索模型必须先做一次装配才知道哪些 θ 没有刚度，
    // 因此这里在初始构形上预装配一次，补上这些约束后再进入迭代。
    {
        std::vector<Eigen::Triplet<double>> trip0;
        DenseVec fi0, fe0;
        assembleState(active, U, fi0, fe0, &trip0);
        if (!m.ces.empty()) {
            const double penalty = effectiveCePenalty(trip0, m.cePenalty);
            appendConstraintPenaltyTriplets(m, index, nDofs, penalty, trip0);
        }
        SparseMat K0(nDofs, nDofs);
        K0.setFromTriplets(trip0.begin(), trip0.end());
        K0.makeCompressed();

        std::vector<int> cd2;
        std::vector<double> cv2;
        autoConstrainFreeRotations(K0, nDofs, cdofs, cvals, cd2, cv2);
        if (cd2.size() != cdofs.size()) {
            cdofs.swap(cd2);
            cvals.swap(cv2);
            isConstrained.assign(nDofs, 0);
            for (int d : cdofs) isConstrained[d] = 1;
            zeroPrescribed.assign(cdofs.size(), 0.0);
            applyPrescribed(U);
        }
    }

    // ---- 外层：荷载分级；内层：Newton-Raphson ----
    const int nSteps = std::max(1, m.nonlinearLoadSteps);
    double lam = 0.0;                           // 已成功施加的荷载比例
    double dlam = 1.0 / static_cast<double>(nSteps);
    const double minLam = 1.0 / 4096.0;         // 再细分下去没有意义
    const int maxLineSearch = 12;               // 步长最多折半 12 次

    int finalIter = 0;
    double finalResidual = 0.0;
    bool converged = false;

    while (lam < 1.0 - 1e-12) {
        const double lamTarget = std::min(1.0, lam + dlam);
        loadScale = lamTarget;

        // 本步起点：不收敛时要回滚到它，再用更小的步长重来。
        const DenseVec Usaved = U;
        const std::vector<char> actSaved = active;

        bool stepOk = false;
        int stepIter = 0;

        for (int iter = 0; iter <= maxIter; ++iter) {
            std::vector<Eigen::Triplet<double>> triplets;
            assembleState(active, U, Fint, Fext, &triplets);
            if (!m.ces.empty()) {
                const double penalty = effectiveCePenalty(triplets, m.cePenalty);
                appendConstraintPenaltyTriplets(m, index, nDofs, penalty, triplets);
            }

            SparseMat Kt(nDofs, nDofs);
            Kt.setFromTriplets(triplets.begin(), triplets.end());
            Kt.makeCompressed();

            // 自由自由度残差范数，用于收敛判断。
            DenseVec r = Fext - Fint;
            if (!m.ces.empty())
                r -= constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, U);
            double rFree2 = 0.0, fFree2 = 0.0, rMaxAbs = 0.0;
            int rMaxDof = -1;
            double kDiagMax = 0.0;
            double kDiagMin = std::numeric_limits<double>::infinity();
            for (int d = 0; d < nDofs; ++d) {
                if (isConstrained[d]) continue;
                rFree2 += r(d) * r(d);
                fFree2 += Fext(d) * Fext(d);
                if (std::abs(r(d)) > rMaxAbs) { rMaxAbs = std::abs(r(d)); rMaxDof = d; }
                const double kd = std::abs(Kt.coeff(d, d));
                if (kd > kDiagMax) kDiagMax = kd;
                if (kd < kDiagMin) kDiagMin = kd;
            }
            const double rFree = std::sqrt(rFree2);
            const double fFree = std::sqrt(fFree2);
            const double resRef = tol * std::max(1.0, fFree);
            stepIter = iter + 1;

            diag.residualHistory.push_back(rFree);
            diag.lastResidual = rFree;
            diag.residualRef = resRef;
            if (!haveInitialResidual) { diag.initialResidual = rFree; haveInitialResidual = true; }
            if (kDiagMin > 0.0) diag.diagCondEst = kDiagMax / kDiagMin;
            if (rMaxDof >= 0) {
                diag.maxResidualNode = ids[rMaxDof / 3];
                diag.maxResidualDof = (rMaxDof % 3) + 1;
                diag.maxResidualValue = rMaxAbs;
            }

            // 自由 DOF 若没有刚度，直接给出可诊断错误。
            for (int d = 0; d < nDofs; ++d) {
                if (!isConstrained[d] && Kt.coeff(d, d) == 0.0) {
                    const int ni = d / 3;
                    const int ld = d % 3;
                    throw std::runtime_error(
                        "方案C：存在零刚度自由度：节点 " + std::to_string(ids[ni]) +
                        " dof " + std::to_string(ld + 1) +
                        "（未约束且未被任何激活单元连接）");
                }
            }

            // 残差判据：当前构形已平衡，无需再解增量。
            if (rFree <= resRef) { stepOk = true; break; }

            // 求解增量：Kt·ΔU = Fext - Fint。
            // Newton 迭代中切线矩阵的**非零结构恒定**，只有数值在变，
            // 因此符号分解只需做一次（reuse 由调用方在施工阶段/扫描间传递）。
            DenseVec rhs = r;
            applyDirichlet(Kt, rhs, cdofs, zeroPrescribed);
            DenseVec dU;
            if (reuse && reuse->reusable()) {
                reuse->factorizeReuse(Kt);
                dU = reuse->solve(rhs);
            } else {
                dU = solve(Kt, rhs, type);
            }
            if (!dU.allFinite()) {
                diag.failureReason = "切线方程解出非有限增量 ΔU（矩阵奇异或严重病态）";
                break;
            }

            // ---- 线搜索（backtracking）：α ∈ {1, 1/2, 1/4, ...} ----
            // 索松弛、跳跃失稳时全步长 ΔU 会让残差不降反升甚至发散。
            // 正常情况下 α=1 立刻被接受，结果与不启用线搜索时一致。
            double alpha = 1.0;
            DenseVec Ubest = U + dU;
            applyPrescribed(Ubest);
            double rBest = residualNorm(active, Ubest);
            if (m.nonlinearLineSearch && !(rBest < rFree)) {
                for (int ls = 0; ls < maxLineSearch; ++ls) {
                    alpha *= 0.5;
                    ++diag.lineSearchCuts;
                    Ubest = U + alpha * dU;
                    applyPrescribed(Ubest);
                    rBest = residualNorm(active, Ubest);
                    if (rBest < rFree) break;
                }
            }

            U = Ubest;

            // 依据新构形更新 cable/strut 激活状态。
            std::vector<char> newActive = updateActive(U);
            const bool statusChanged = (newActive != active);
            active.swap(newActive);

            const double duOk = (alpha * dU.norm()) <= tol * std::max(1.0, U.norm());
            const bool resOk = rBest <= resRef;

            // 发散检测：残差涨到初始值的 1e8 倍以上，再迭代没有意义。
            if (!U.allFinite() ||
                (diag.initialResidual > 0.0 && rBest > 1e8 * diag.initialResidual)) {
                diag.failureReason = "残差发散（超过初始残差的 1e8 倍）";
                break;
            }

            if (!statusChanged && duOk && resOk) { stepOk = true; break; }
        }

        finalIter += stepIter;

        if (stepOk) {
            lam = lamTarget;
            ++diag.loadSteps;
            continue;
        }

        // 本荷载步失败：回滚，步长减半重来（自适应荷载分级）。
        U = Usaved;
        active = actSaved;
        if (diag.failureReason.empty())
            diag.failureReason = "单个荷载步内达到最大迭代次数 " + std::to_string(maxIter);
        dlam *= 0.5;
        ++diag.stepRefinements;
        if (dlam < minLam) break;
    }

    converged = (lam >= 1.0 - 1e-12);
    diag.iterations = finalIter;
    diag.loadRatio = lam;
    diag.converged = converged;

    // 用最终位移和最终激活状态重新组装一次，得到精确内力/反力/残差。
    loadScale = 1.0;
    assembleState(active, U, Fint, Fext, nullptr);
    DenseVec rFinal = Fext - Fint;
    if (!m.ces.empty())
        rFinal -= constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, U);
    double rf2 = 0.0;
    for (int d = 0; d < nDofs; ++d)
        if (!isConstrained[d]) rf2 += rFinal(d) * rFinal(d);
    finalResidual = std::sqrt(rf2);
    diag.lastResidual = finalResidual;

    DenseVec R = Fint - Fext; // 约束反力 = 内部节点力 - 外部节点力
    if (!m.ces.empty())
        R += constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, U);

    if (!converged && m.nonlinearThrowOnFailure)
        throw std::runtime_error(diag.describe());

    ModelResult result;
    result.U = U;
    result.R = R;
    result.nodeIds = ids;
    result.coords.reserve(nNodes);
    for (int id : ids) result.coords.push_back(m.nodes.at(id));
    result.nonlinearIterations = finalIter;
    result.nonlinearResidual = finalResidual;
    result.nonlinearConverged = converged;
    if (diagOut) *diagOut = diag;
    return result;
}

// ==================== 线性系统的装配（供多次求解复用）====================

// 给定激活状态，装配总体刚度矩阵 K 与荷载向量 F（尚未施加约束）。
// solveModel2D 与影响线扫描共用这一段，避免两处逻辑走偏。
inline void assembleLinear2D(const Model2D& m,
                             const std::vector<char>& active,
                             const std::map<int, int>& index,
                             int nDofs,
                             SparseMat& K, DenseVec& F) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
        if (!active[ei]) continue;
        const auto& e = m.elements[ei];
        if (e.conn.size() != 2)
            throw std::runtime_error("二维单元必须有两个节点");
        const int i1 = index.at(e.conn[0]);
        const int i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

        if (isAxial(e.type)) {
            TrussElement te{e.E, e.A, 2,
                            Eigen::Vector3d(p1(0), p1(1), 0.0),
                            Eigen::Vector3d(p2(0), p2(1), 0.0)};
            const DenseMat Ke = te.stiffness(); // 4x4
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i2, 3 * i2 + 1};
            scatter(Ke, g, triplets);
        } else if (e.type == "RIGID_PIN") {
            const DenseMat Ke = rigidPinStiffness(p1, p2);
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i1 + 2,
                                        3 * i2, 3 * i2 + 1, 3 * i2 + 2};
            scatter(Ke, g, triplets);
        } else if (e.type == "BEAM" || e.type == "RIGID") {
            EulerBeam2D be{e.E, e.A, e.I, p1, p2};
            const DenseMat Ke = be.stiffness(); // 6x6
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i1 + 2,
                                        3 * i2, 3 * i2 + 1, 3 * i2 + 2};
            scatter(Ke, g, triplets);
        } else {
            throw std::runtime_error("未知单元类型: " + e.type);
        }
    }
    if (!m.ces.empty()) {
        const double penalty = effectiveCePenalty(triplets, m.cePenalty);
        appendConstraintPenaltyTriplets(m, index, nDofs, penalty, triplets);
    }
    K.resize(nDofs, nDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();

    F = DenseVec::Zero(nDofs);
    for (const auto& l : m.loads) {
        if (l.dof < 1 || l.dof > 3)
            throw std::runtime_error("载荷自由度必须在 1~3 之间");
        F(3 * index.at(l.node) + (l.dof - 1)) += l.value;
    }
    for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
        if (!active[ei]) continue;
        const auto& e = m.elements[ei];
        const int i1 = index.at(e.conn[0]), i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);
        const DenseVec fe = elementEquivalentLoad(m, e, p1, p2);
        addElementLoad(elementGdofs2D(e.type, i1, i2), fe, F);
    }
    if (!m.ces.empty())
        F += constraintPenaltyConstantForce(m, index, nDofs, m.cePenalty);
}

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
inline Assembled2D assembleModel2DLinear(const Model2D& m) {
    Assembled2D a;
    for (const auto& [id, c] : m.nodes) {
        a.index[id] = static_cast<int>(a.ids.size());
        a.ids.push_back(id);
    }
    a.nDofs = 3 * static_cast<int>(a.ids.size());

    const std::vector<char> active(m.elements.size(), 1);
    SparseMat K;
    DenseVec F;
    assembleLinear2D(m, active, a.index, a.nDofs, K, F);
    a.K0 = K;
    a.F0 = F;

    a.isConstrained.assign(a.nDofs, 0);
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        const int d = 3 * a.index.at(c.node) + (c.dof - 1);
        a.cdofs.push_back(d);
        a.cvals.push_back(c.value);
        a.isConstrained[d] = 1;
    }
    std::vector<int> cd;
    std::vector<double> cv;
    autoConstrainFreeRotations(K, a.nDofs, a.cdofs, a.cvals, cd, cv);
    a.cdofs = cd;
    a.cvals = cv;
    a.isConstrained.assign(a.nDofs, 0);
    for (int d : a.cdofs) a.isConstrained[d] = 1;

    a.K = K;
    a.F = F;
    applyDirichlet(a.K, a.F, a.cdofs, a.cvals);
    return a;
}

// 由位移与装配结果构造 ModelResult（含反力），供响应计算复用。
inline ModelResult makeResult2D(const Model2D& m, const Assembled2D& a,
                                const DenseVec& U, const DenseVec& F0) {
    ModelResult r;
    r.U = U;
    r.R = a.K0 * U - F0;
    r.nodeIds = a.ids;
    r.coords.reserve(a.ids.size());
    for (int id : a.ids) r.coords.push_back(m.nodes.at(id));
    return r;
}

// 组装并求解二维结构（默认 SimplicialLDLT）。
inline ModelResult solveModel2D(const Model2D& m,
                                SolverType type = SolverType::SimplicialLDLT,
                                const DenseVec* U0 = nullptr,
                                const std::vector<char>* active0 = nullptr,
                                NonlinearDiag* diagOut = nullptr,
                                ReusableSolver* reuse = nullptr) {
    if (m.analysisType == "NONLINEAR")
        return solveModel2DNonlinear(m, type, U0, active0, diagOut, reuse);

    // 节点号 -> 0-based 索引（std::map 按键有序）。
    std::vector<int> ids;
    std::map<int, int> index;
    for (const auto& [id, c] : m.nodes) {
        index[id] = static_cast<int>(ids.size());
        ids.push_back(id);
    }
    const int nNodes = static_cast<int>(ids.size());
    const int nDofs = 3 * nNodes;

    // cable/strut 激活状态（初值全激活）。
    std::vector<char> active(m.elements.size(), 1);

    // 约束自由度。
    std::vector<int> cdofs;
    std::vector<double> cvals;
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        cdofs.push_back(3 * index.at(c.node) + (c.dof - 1));
        cvals.push_back(c.value);
    }

    SparseMat K0;
    DenseVec F0, U;
    const int maxIter = 100;

    for (int iter = 0; ; ++iter) {
        // 组装 K / F（只含激活单元）。
        SparseMat K;
        DenseVec F;
        assembleLinear2D(m, active, index, nDofs, K, F);

        K0 = K;   // 施加约束前，用于反力计算
        F0 = F;
        std::vector<int> cd;
        std::vector<double> cv;
        autoConstrainFreeRotations(K, nDofs, cdofs, cvals, cd, cv);
        applyDirichlet(K, F, cd, cv);

        // 诊断：零对角线自由度。
        for (int d = 0; d < nDofs; ++d) {
            if (K.coeff(d, d) == 0.0) {
                const int ni = d / 3;
                const int ld = d % 3;
                throw std::runtime_error(
                    "存在零刚度自由度：节点 " + std::to_string(ids[ni]) +
                    " dof " + std::to_string(ld + 1) +
                    "（未约束且未被任何单元连接）");
            }
        }

        // 传入 reuse 时复用符号分解（结构不变则只做数值分解）。
        // 施工阶段推进、影响线扫描这类"结构不变、只有数值变"的场景收益最直接。
        if (reuse && reuse->reusable()) {
            reuse->factorizeReuse(K);
            U = reuse->solve(F);
        } else {
            U = solve(K, F, type);
        }
        validateSolution(K, F, U, "二维求解");

        // 检查 cable/strut 的开关状态。
        bool changed = false;
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
            if (!active[ei]) continue;
            const auto& e = m.elements[ei];
            if (e.type != "CABLE" && e.type != "STRUT") continue;
            const double N = axialForce(e, U, index.at(e.conn[0]), index.at(e.conn[1]),
                                        m.nodes.at(e.conn[0]), m.nodes.at(e.conn[1]));
            if (e.type == "CABLE" && N <= 0.0) { active[ei] = 0; changed = true; }
            else if (e.type == "STRUT" && N >= 0.0) { active[ei] = 0; changed = true; }
        }
        if (!changed) break;
        if (iter + 1 >= maxIter)
            throw std::runtime_error("cable/strut 开关迭代未收敛");
    }

    const DenseVec R = K0 * U - F0;

    ModelResult result;
    result.U = U;
    result.R = R;
    result.nodeIds = ids;
    result.coords.reserve(nNodes);
    for (int id : ids) result.coords.push_back(m.nodes.at(id));
    return result;
}

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
inline std::vector<ElementForce> elementForces2D(const Model2D& m, const DenseVec& U) {
    std::map<int, int> index;
    int n = 0;
    for (const auto& [id, c] : m.nodes) index[id] = n++;

    std::vector<ElementForce> forces;
    forces.reserve(m.elements.size());

    for (const auto& e : m.elements) {
        const int i1 = index.at(e.conn[0]);
        const int i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

        ElementForce f;
        f.id = e.id;
        f.type = e.type;

        if (isRigid(e.type)) {
            forces.push_back(f);
            continue;
        }

        if (isAxial(e.type)) {
            if (m.analysisType == "NONLINEAR") {
                const Eigen::Vector2d d0 = p2 - p1;
                const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
                const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
                const Eigen::Vector2d d = d0 + u2 - u1;
                const double L = d.norm();
                if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
                const Eigen::Vector2d c = d / L;
                const DenseVec fint = axialInternalForceNonlinear(e, U, i1, i2, p1, p2);
                const DenseVec fext = elementExternalEquivalentLoad(m, e, p1, p2);
                const double fi = fint(0) * c(0) + fint(1) * c(1);
                const double fj = fint(2) * c(0) + fint(3) * c(1);
                const double ei = fext(0) * c(0) + fext(1) * c(1);
                const double ej = fext(2) * c(0) + fext(3) * c(1);
                f.N1 = -fi + ei;
                f.N2 = fj - ej;
            } else {
                const Eigen::Vector2d d = p2 - p1;
                const double L = d.norm();
                const Eigen::Vector2d c = d / L;
                const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
                const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
                const double Nrestore = (e.E * e.A / L) * (u2 - u1).dot(c); // 恢复力（受拉）
                const DenseVec fe = elementEquivalentLoad(m, e, p1, p2); // [i,j] 全局
                const double fe_i = fe(0) * c(0) + fe(1) * c(1);
                const double fe_j = fe(2) * c(0) + fe(3) * c(1);
                f.N1 = Nrestore + fe_i;  // 节点1端轴力（受拉为正）
                f.N2 = Nrestore - fe_j;  // 节点2端轴力（受拉为正）
            }
            // 索/压杆松弛时内力按 0 输出。
            if (e.type == "CABLE" && f.N1 <= 0.0) { f.N1 = 0.0; f.N2 = 0.0; }
            else if (e.type == "STRUT" && f.N1 >= 0.0) { f.N1 = 0.0; f.N2 = 0.0; }
        } else if (e.type == "BEAM") {
            EulerBeam2D be{e.E, e.A, e.I, p1, p2};
            const DenseMat T = be.transformation();
            const DenseMat Kl = be.localStiffness();
            Eigen::VectorXd ug(6);
            ug << U(3 * i1 + 0), U(3 * i1 + 1), U(3 * i1 + 2),
                  U(3 * i2 + 0), U(3 * i2 + 1), U(3 * i2 + 2);
            const Eigen::VectorXd ul = T * ug;   // 局部位移
            const Eigen::VectorXd fl = Kl * ul;  // 恢复杆端力（局部）
            const DenseVec fe = elementEquivalentLoad(m, e, p1, p2); // 全局等效荷载
            const Eigen::VectorXd fel = T * fe;   // 全局 -> 局部
            const Eigen::VectorXd fphys = fl - fel; // 物理杆端力
            f.N1 = -fphys(0); // 节点1端轴力（受拉）
            f.N2 = fphys(3);  // 节点2端轴力（受拉）
            f.V1 = fphys(1);   // 截面剪力（左端）
            f.V2 = -fphys(4);  // 截面剪力（右端）
            f.M1 = -fphys(2);  // 截面弯矩（左端，受拉为正/下侧受拉）
            f.M2 = fphys(5);   // 截面弯矩（右端）
        }
        forces.push_back(f);
    }
    return forces;
}

// 单元应力：桁架轴向应力；梁上下翼缘正应力。
struct ElementStress {
    int id = 0;
    std::string type;
    double s1 = 0.0, s2 = 0.0;       // 桁架：两端轴向应力
    double stop1 = 0.0, sbot1 = 0.0; // 梁：节点1端上/下翼缘应力
    double stop2 = 0.0, sbot2 = 0.0; // 梁：节点2端上/下翼缘应力
};

// 由内力计算应力（受拉为正）。
inline std::vector<ElementStress> elementStress2D(const Model2D& m,
                                                  const std::vector<ElementForce>& forces) {
    std::map<int, const Model2D::Element*> byId;
    for (const auto& e : m.elements) byId[e.id] = &e;

    std::vector<ElementStress> out;
    out.reserve(forces.size());
    for (const auto& f : forces) {
        const Model2D::Element* e = byId.at(f.id);
        ElementStress s;
        s.id = f.id;
        s.type = f.type;
        if (isRigid(f.type)) {
            // 刚性连接杆不参与应力输出。
        } else if (isAxial(f.type)) {
            s.s1 = f.N1 / e->A;
            s.s2 = f.N2 / e->A;
        } else { // BEAM
            if (e->Wt <= 0.0 || e->Wb <= 0.0)
                throw std::runtime_error("梁单元应力需要 *ELEMENT 里的 WT 和 WB");
            s.stop1 = f.N1 / e->A - f.M1 / e->Wt;
            s.sbot1 = f.N1 / e->A + f.M1 / e->Wb;
            s.stop2 = f.N2 / e->A - f.M2 / e->Wt;
            s.sbot2 = f.N2 / e->A + f.M2 / e->Wb;
        }
        out.push_back(s);
    }
    return out;
}

// ==================== 施工阶段正装分析 ====================

struct StageResult {
    int stage = 0;
    DenseVec U;
    DenseVec R;
    std::vector<ElementForce> forces;
    std::vector<ElementStress> stresses;
};

inline std::vector<StageResult> solveModel2DStages(const Model2D& m,
                                                   SolverType type = SolverType::SimplicialLDLT) {
    if (m.stages.empty())
        throw std::runtime_error("模型没有定义施工阶段 *STAGE");

    std::vector<Model2D::Stage> stages = m.stages;
    std::sort(stages.begin(), stages.end(),
              [](const Model2D::Stage& a, const Model2D::Stage& b) { return a.id < b.id; });

    std::vector<StageResult> results;
    std::vector<char> active(m.elements.size(), 0);
    std::vector<Model2D::Element> workingElements = m.elements;

    // 徐变折减基准：readModel2D 已填好；但程序化建模（不经文件解析）时是空的，
    // 这里兜底补一次——否则 baseModulusOf 会退化成 e.E，
    // 而 e.E 在上一阶段已被折减过，折减就会被逐级累乘（正是 P1-5 要修的 bug）。
    std::map<int, double> baseE = m.baseModulus;
    if (baseE.empty())
        for (const auto& e : m.elements) baseE[e.id] = e.E;

    // 各阶段的总体矩阵结构只在"单元激活/钝化"时才变，徐变折减与荷载变化都不改结构。
    // 用可复用句柄把符号分解跨阶段保留下来（P1-7）。
    ReusableSolver reuse(type);
    std::vector<Model2D::Load> cumulativeLoads;
    std::vector<Model2D::ElemLoad> cumulativeElemLoads;
    DenseVec prevU;
    bool hasPrevU = false;
    double prevStageTime = 0.0;
    std::map<int, double> prevAxialForce;
    std::map<int, double> prevM1, prevM2;
    std::vector<char> prevActive;
    bool hasPrevActive = false;
    // 单元加载龄期（天）：单元被激活的那个阶段的龄期。
    // 逐阶段 AEMM 用它算 Δφ = φ(t_now) - φ(t_load)。
    // 未记录（例如阶段未给 TIME）的单元按 0 处理，即"从龄期 0 开始持荷"。
    std::map<int, double> loadAge;

    for (const auto& st : stages) {
        for (int id : st.activate) {
            bool found = false;
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) { active[ei] = 1; found = true; break; }
            if (!found) throw std::runtime_error("*STAGE 激活了不存在的单元: " + std::to_string(id));
            // *STAGE 的 TIME 是"该阶段末的龄期"，单元在阶段开始时就安装并受力，
            // 因此加载龄期取**上一阶段末**的龄期（首个阶段为 0）。
            // 若误取 st.time，则"第 1 阶段 TIME=100"会得到 Δφ=φ(100)-φ(100)=0，
            // 徐变被整个抹掉。
            loadAge[id] = prevStageTime;
        }
        for (int id : st.deactivate) {
            bool found = false;
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) { active[ei] = 0; found = true; break; }
            if (!found) throw std::runtime_error("*STAGE 钝化了不存在的单元: " + std::to_string(id));
            // 钝化再激活视为重新加载，龄期重新起算。
            loadAge.erase(id);
        }

        double creepDPhi = 0.0;

        Model2D mc = m;
        mc.baseModulus = baseE;
        mc.elements.clear();
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
            if (active[ei]) mc.elements.push_back(workingElements[ei]);

        // 本阶段先移除指定施工临时荷载。
        cumulativeLoads.erase(
            std::remove_if(cumulativeLoads.begin(), cumulativeLoads.end(),
                           [&](const Model2D::Load& l) {
                               for (const auto& [node, dof] : st.unloadNodeDofs)
                                   if (l.node == node && l.dof == dof) return true;
                               return false;
                           }),
            cumulativeLoads.end());
        cumulativeElemLoads.erase(
            std::remove_if(cumulativeElemLoads.begin(), cumulativeElemLoads.end(),
                           [&](const Model2D::ElemLoad& el) {
                               return std::find(st.unloadElems.begin(), st.unloadElems.end(), el.elem) != st.unloadElems.end();
                           }),
            cumulativeElemLoads.end());

        cumulativeLoads.insert(cumulativeLoads.end(), st.loads.begin(), st.loads.end());
        cumulativeElemLoads.insert(cumulativeElemLoads.end(), st.elemLoads.begin(), st.elemLoads.end());
        mc.loads = cumulativeLoads;
        mc.elemLoads = cumulativeElemLoads;
        mc.prestresses = st.prestresses;
        // 本阶段是否走"按时变曲线逐阶段 AEMM"路径。
        bool creepByAge = false;
        if (st.time > 0.0 && (m.creepModel != "CURVE" || !m.creepCurve.points.empty())) {
            Model2D::CreepSpec cs;
            cs.all = true;
            const double phiNow = creepCoefficientAt(m, st.time);
            const double phiPrev = prevStageTime > 0.0
                                       ? creepCoefficientAt(m, prevStageTime)
                                       : 0.0;
            creepDPhi = (prevStageTime > 0.0) ? (phiNow - phiPrev) : phiNow;
            cs.phi = phiNow;   // 仅作为无龄期信息时的兜底值
            cs.chi = 0.8;
            mc.creepSpecs = {cs};
            creepByAge = true;
        } else {
            mc.creepSpecs = st.creepSpecs;
        }
        if (st.time > 0.0 && (m.shrinkageModel != "CURVE" || !m.shrinkageCurve.points.empty())) {
            Model2D::ShrinkageSpec sh;
            sh.all = true;
            const double ecsNow = shrinkageStrainAt(m, st.time);
            const double ecsPrev = prevStageTime > 0.0
                                       ? shrinkageStrainAt(m, prevStageTime)
                                       : 0.0;
            sh.ecs = (prevStageTime > 0.0) ? (ecsNow - ecsPrev) : ecsNow;
            mc.shrinkageSpecs = {sh};
        } else {
            mc.shrinkageSpecs = st.shrinkageSpecs;
        }
        mc.tempProfiles = st.tempProfiles;
        mc.uniformTemp = 0.0;
        mc.elemTemps.clear();
        mc.selfWeights.clear();

        // 施工阶段支座变化：先移除，再添加。
        auto constraints = mc.constraints;
        constraints.erase(
            std::remove_if(constraints.begin(), constraints.end(),
                           [&](const Model2D::Constraint& c) {
                               for (const auto& [node, dof] : st.removedConstraints)
                                   if (c.node == node && c.dof == dof) return true;
                               return false;
                           }),
            constraints.end());
        constraints.insert(constraints.end(), st.addedConstraints.begin(), st.addedConstraints.end());
        mc.constraints = constraints;

        // 徐变折减：逐阶段 AEMM。
        // 关键：模量始终从 baseModulus 重算，而不是在上一阶段折减后的
        // 基础上再乘一次 —— 旧实现那样做会让折减被逐级累乘
        // （E0/((1+χφ1)(1+χφ2)…)，而不是 AEMM 的 E0/(1+χ(φ(t)-φ(t0)))）。
        if (!mc.creepSpecs.empty())
            processTimeEffects(mc, creepByAge ? &loadAge : nullptr,
                               creepByAge ? st.time : 0.0);

        if (creepDPhi != 0.0 && !prevAxialForce.empty()) {
            for (const auto& e : mc.elements) {
                auto it = prevAxialForce.find(e.id);
                if (it == prevAxialForce.end() || it->second == 0.0) continue;
                const double P = creepDPhi * it->second;
                const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
                const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);
                const Eigen::Vector2d d = p2 - p1;
                const double L = d.norm();
                if (L <= 0.0) continue;
                const Eigen::Vector2d c = d / L;
                if (isAxial(e.type)) {
                    mc.loads.push_back({e.conn[0], 1, -P * c(0)});
                    mc.loads.push_back({e.conn[0], 2, -P * c(1)});
                    mc.loads.push_back({e.conn[1], 1, P * c(0)});
                    mc.loads.push_back({e.conn[1], 2, P * c(1)});
                } else if (e.type == "BEAM") {
                    DenseVec fl = DenseVec::Zero(6);
                    fl(0) = -P;
                    fl(3) = P;
                    if (prevM1.count(e.id)) fl(2) = -creepDPhi * prevM1[e.id];
                    if (prevM2.count(e.id)) fl(5) = creepDPhi * prevM2[e.id];
                    EulerBeam2D be{e.E, e.A, e.I, p1, p2};
                    const DenseVec fe = be.transformation().transpose() * fl;
                    for (int d = 0; d < 3; ++d) {
                        mc.loads.push_back({e.conn[0], d + 1, fe(d)});
                        mc.loads.push_back({e.conn[1], d + 1, fe(d + 3)});
                    }
                }
            }
        }

        if (!mc.prestresses.empty())
            processPrestress(mc);

        if (st.time > 0.0)
            prevStageTime = st.time;

        cumulativeLoads = mc.loads;
        cumulativeElemLoads = mc.elemLoads;
        for (const auto& me : mc.elements) {
            for (auto& we : workingElements)
                if (we.id == me.id) { we.E = me.E; break; }
        }

        // 未激活节点暂按固定处理，避免出现零刚度自由度。
        std::map<int, bool> connected;
        for (const auto& e : mc.elements) {
            connected[e.conn[0]] = true;
            connected[e.conn[1]] = true;
        }
        for (const auto& [id, c] : m.nodes) {
            if (connected[id]) continue;
            for (int dof = 1; dof <= 3; ++dof)
                mc.constraints.push_back({id, dof, 0.0});
        }

        const DenseVec* U0 = (m.analysisType == "NONLINEAR" && hasPrevU) ? &prevU : nullptr;
        const std::vector<char>* A0 = (m.analysisType == "NONLINEAR" && hasPrevActive) ? &prevActive : nullptr;
        const ModelResult r = solveModel2D(mc, type, U0, A0, nullptr, &reuse);
        prevU = r.U;
        hasPrevU = true;
        StageResult sr;
        sr.stage = st.id;
        sr.U = r.U;
        sr.R = r.R;
        sr.forces = elementForces2D(mc, r.U);
        if (m.analysisType == "NONLINEAR") {
            prevActive.assign(mc.elements.size(), 1);
            for (std::size_t ei = 0; ei < mc.elements.size(); ++ei) {
                const auto& e = mc.elements[ei];
                if (e.type == "CABLE") prevActive[ei] = (sr.forces[ei].N1 > 0.0) ? 1 : 0;
                else if (e.type == "STRUT") prevActive[ei] = (sr.forces[ei].N1 < 0.0) ? 1 : 0;
                else prevActive[ei] = 1;
            }
            hasPrevActive = true;
        }
        for (const auto& f : sr.forces) {
            prevAxialForce[f.id] = f.N1;
            if (f.type == "BEAM") {
                prevM1[f.id] = f.M1;
                prevM2[f.id] = f.M2;
            }
        }
        try {
            sr.stresses = elementStress2D(mc, sr.forces);
        } catch (const std::exception&) {
            sr.stresses.clear();
        }
        results.push_back(sr);
    }
    return results;
}

struct CamberResult {
    int node = 0;
    int dof = 0;
    double camber = 0.0;      // 预拱度 = -最终位移
    double finalDisplacement = 0.0;
    double original = 0.0;
    double installation = 0.0;
};

inline std::vector<CamberResult> computeCamber(const Model2D& m,
                                               SolverType type = SolverType::SimplicialLDLT) {
    const auto stages = solveModel2DStages(m, type);
    if (stages.empty()) return {};
    const DenseVec U = stages.back().U;
    std::map<int, int> index;
    int n = 0;
    for (const auto& [id, c] : m.nodes) index[id] = n++;

    std::vector<CamberResult> out;
    out.reserve(m.cambers.size());
    for (const auto& c : m.cambers) {
        const int d = 3 * index.at(c.node) + (c.dof - 1);
        const double disp = U(d);
        const Eigen::Vector2d coord = m.nodes.at(c.node);
        const double orig = (c.dof == 1) ? coord(0) : (c.dof == 2) ? coord(1) : 0.0;
        const double camber = c.value - disp;
        out.push_back({c.node, c.dof, camber, disp, orig, orig + camber});
    }
    return out;
}

struct CamberIterationResult {
    DenseVec cambers;
    DenseVec finalDisplacements;
    int iterations = 0;
    std::vector<StageResult> stages;
};

inline CamberIterationResult iterateCamber(const Model2D& m,
                                           SolverType type = SolverType::SimplicialLDLT,
                                           double tol = 1e-8,
                                           int maxIter = 30) {
    if (m.cambers.empty())
        throw std::runtime_error("非线性预拱度迭代需要 *CAMBER");
    const int n = static_cast<int>(m.cambers.size());
    DenseVec x = DenseVec::Zero(n);

    std::map<int, int> nid;
    int nn = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = nn++;

    auto eval = [&](const DenseVec& s) -> DenseVec {
        Model2D mc = m;
        for (int k = 0; k < n; ++k) {
            const auto& cp = m.cambers[k];
            auto& coord = mc.nodes.at(cp.node);
            if (cp.dof == 1) coord(0) += s(k);
            else if (cp.dof == 2) coord(1) += s(k);
            else throw std::runtime_error("非线性预拱度目前只支持 UX/UY");
        }
        const auto stages = solveModel2DStages(mc, type);
        if (stages.empty()) throw std::runtime_error("施工阶段结果为空");
        const DenseVec U = stages.back().U;
        DenseVec y(n);
        for (int k = 0; k < n; ++k) {
            const int d = 3 * nid.at(m.cambers[k].node) + (m.cambers[k].dof - 1);
            y(k) = s(k) + U(d);
        }
        return y;
    };

    for (int iter = 0; iter < maxIter; ++iter) {
        const DenseVec y = eval(x);
        DenseVec r(n);
        for (int k = 0; k < n; ++k) r(k) = y(k) - m.cambers[k].value;
        if (r.norm() < tol) {
            CamberIterationResult res;
            res.cambers = x;
            res.finalDisplacements = y;
            res.iterations = iter;
            Model2D mc = m;
            for (int k = 0; k < n; ++k) {
                auto& coord = mc.nodes.at(m.cambers[k].node);
                if (m.cambers[k].dof == 1) coord(0) += x(k);
                else coord(1) += x(k);
            }
            res.stages = solveModel2DStages(mc, type);
            return res;
        }

        Eigen::MatrixXd J(n, n);
        for (int j = 0; j < n; ++j) {
            const double h = 1e-6 * std::max(1.0, std::abs(x(j)));
            DenseVec xp = x, xm = x;
            xp(j) += h;
            xm(j) -= h;
            J.col(j) = (eval(xp) - eval(xm)) / (2.0 * h);
        }
        const double lambda = 1e-8;
        Eigen::MatrixXd A = J.transpose() * J;
        for (int j = 0; j < n; ++j) A(j, j) += lambda * std::max(1.0, A(j, j));
        const Eigen::VectorXd dx = A.fullPivLu().solve(-(J.transpose() * r));
        if (!dx.allFinite())
            throw std::runtime_error("非线性预拱度迭代线性求解失败");
        x += dx;
    }
    throw std::runtime_error("非线性预拱度迭代未收敛");
}

struct FeedbackCorrection {
    double current = 0.0;
    double derivative = 0.0;
    double correction = 0.0;
};

inline FeedbackCorrection computeFeedbackCorrection(const Model2D& m,
                                                     SolverType type = SolverType::SimplicialLDLT) {
    if (m.feedbacks.empty())
        throw std::runtime_error("反馈修正需要 *FEEDBACK");
    if (m.stageControls.empty())
        throw std::runtime_error("反馈修正需要 *STAGE_CONTROL");
    const auto& fb = m.feedbacks.front();
    std::map<int, int> nid;
    int n = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = n++;

    auto response = [&](double scale) {
        Model2D mc = m;
        const auto& ctrl = m.stageControls.front();
        for (auto& st : mc.stages) {
            if (st.id == ctrl.stage) {
                st.prestresses[ctrl.prestressIndex - 1].P *= scale;
                break;
            }
        }
        const auto stages = solveModel2DStages(mc, type);
        const int d = 3 * nid.at(fb.node) + (fb.dof - 1);
        return stages.back().U(d);
    };

    const double s0 = 1.0;
    const double s1 = 1.01;
    const double y0 = response(s0);
    const double y1 = response(s1);
    const double deriv = (y1 - y0) / (s1 - s0);
    FeedbackCorrection out;
    out.current = y0;
    out.derivative = deriv;
    out.correction = (std::abs(deriv) > 1e-20) ? (fb.measured - y0) / deriv : 0.0;
    return out;
}

inline std::vector<StageResult> solveModel2DStagesBackward(const Model2D& m,
                                                           SolverType type = SolverType::SimplicialLDLT) {
    if (m.stages.empty())
        throw std::runtime_error("模型没有定义施工阶段 *STAGE");

    std::vector<Model2D::Stage> stages = m.stages;
    std::sort(stages.begin(), stages.end(),
              [](const Model2D::Stage& a, const Model2D::Stage& b) { return a.id < b.id; });

    std::vector<char> active(m.elements.size(), 0);
    std::vector<Model2D::Load> cumLoads;
    std::vector<Model2D::ElemLoad> cumElemLoads;
    for (const auto& st : stages) {
        for (int id : st.activate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 1;
        }
        for (int id : st.deactivate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 0;
        }
        cumLoads.insert(cumLoads.end(), st.loads.begin(), st.loads.end());
        cumElemLoads.insert(cumElemLoads.end(), st.elemLoads.begin(), st.elemLoads.end());
    }

    auto removeLoads = [&](std::vector<Model2D::Load>& vec, const std::vector<Model2D::Load>& rem) {
        vec.erase(std::remove_if(vec.begin(), vec.end(),
                                 [&](const Model2D::Load& l) {
                                     for (const auto& r : rem)
                                         if (l.node == r.node && l.dof == r.dof && l.value == r.value)
                                             return true;
                                     return false;
                                 }),
                  vec.end());
    };
    auto removeElemLoads = [&](std::vector<Model2D::ElemLoad>& vec,
                               const std::vector<Model2D::ElemLoad>& rem) {
        vec.erase(std::remove_if(vec.begin(), vec.end(),
                                 [&](const Model2D::ElemLoad& el) {
                                     for (const auto& r : rem)
                                         if (el.elem == r.elem && el.dir == r.dir &&
                                             el.q1 == r.q1 && el.q2 == r.q2)
                                             return true;
                                     return false;
                                 }),
                  vec.end());
    };

    auto buildActiveModel = [&](Model2D& mc) {
        mc.elements.clear();
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
            if (active[ei]) mc.elements.push_back(m.elements[ei]);
        mc.loads = cumLoads;
        mc.elemLoads = cumElemLoads;
        mc.prestresses.clear();
        mc.creepSpecs.clear();
        mc.shrinkageSpecs.clear();
        mc.tempProfiles.clear();
        mc.uniformTemp = 0.0;
        mc.elemTemps.clear();
        mc.selfWeights.clear();

        std::map<int, bool> connected;
        for (const auto& e : mc.elements) {
            connected[e.conn[0]] = true;
            connected[e.conn[1]] = true;
        }
        for (const auto& [id, c] : m.nodes) {
            if (connected[id]) continue;
            for (int dof = 1; dof <= 3; ++dof)
                mc.constraints.push_back({id, dof, 0.0});
        }
    };

    DenseVec prevU;
    bool hasPrevU = false;
    std::vector<char> prevActive;
    bool hasPrevActive = false;
    std::vector<StageResult> results;
    // 同正装：结构不变时跨阶段复用符号分解（P1-7）。
    ReusableSolver reuse(type);

    for (auto it = stages.rbegin(); it != stages.rend(); ++it) {
        const Model2D::Stage& st = *it;
        Model2D mc = m;
        buildActiveModel(mc);

        const DenseVec* U0 = (m.analysisType == "NONLINEAR" && hasPrevU) ? &prevU : nullptr;
        const std::vector<char>* A0 = (m.analysisType == "NONLINEAR" && hasPrevActive) ? &prevActive : nullptr;
        const ModelResult r = solveModel2D(mc, type, U0, A0, nullptr, &reuse);
        prevU = r.U;
        hasPrevU = true;

        StageResult sr;
        sr.stage = st.id;
        sr.U = r.U;
        sr.R = r.R;
        sr.forces = elementForces2D(mc, r.U);
        if (m.analysisType == "NONLINEAR") {
            prevActive.assign(mc.elements.size(), 1);
            for (std::size_t ei = 0; ei < mc.elements.size(); ++ei) {
                const auto& e = mc.elements[ei];
                if (e.type == "CABLE") prevActive[ei] = (sr.forces[ei].N1 > 0.0) ? 1 : 0;
                else if (e.type == "STRUT") prevActive[ei] = (sr.forces[ei].N1 < 0.0) ? 1 : 0;
                else prevActive[ei] = 1;
            }
            hasPrevActive = true;
        }
        try {
            sr.stresses = elementStress2D(mc, sr.forces);
        } catch (const std::exception&) {
            sr.stresses.clear();
        }
        results.push_back(sr);

        removeLoads(cumLoads, st.loads);
        removeElemLoads(cumElemLoads, st.elemLoads);
        for (int id : st.activate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 0;
        }
        for (int id : st.deactivate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 1;
        }
    }
    return results;
}

inline double forwardBackwardStageResidual(const Model2D& m,
                                           SolverType type = SolverType::SimplicialLDLT) {
    const auto forward = solveModel2DStages(m, type);
    const auto backward = solveModel2DStagesBackward(m, type);
    if (forward.empty() || backward.empty())
        return 0.0;
    return (forward.back().U - backward.front().U).norm();
}

struct StageIterationResult {
    DenseVec scales;
    DenseVec finalValues;
    int iterations = 0;
    std::vector<StageResult> stages;
};

inline StageIterationResult iterateStagePrestress(const Model2D& m,
                                                  SolverType type = SolverType::SimplicialLDLT,
                                                  double tol = 1e-8,
                                                  int maxIter = 30) {
    if (m.stageControls.empty() || m.stageTargets.empty())
        throw std::runtime_error("多变量正装-倒拆迭代需要 *STAGE_CONTROL 和 *STAGE_TARGET");
    const int nv = static_cast<int>(m.stageControls.size());
    const int nt = static_cast<int>(m.stageTargets.size());
    const int nc = static_cast<int>(m.stageConstraints.size());
    DenseVec x = DenseVec::Ones(nv);

    std::map<int, int> nid;
    int nn = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = nn++;

    auto eval = [&](const DenseVec& s) -> DenseVec {
        Model2D mc = m;
        for (int k = 0; k < nv; ++k) {
            const auto& ctrl = m.stageControls[k];
            bool found = false;
            for (auto& st : mc.stages) {
                if (st.id != ctrl.stage) continue;
                if (ctrl.prestressIndex < 1 ||
                    ctrl.prestressIndex > static_cast<int>(st.prestresses.size()))
                    throw std::runtime_error("*STAGE_CONTROL 的 PRESTRESS 序号不存在");
                st.prestresses[ctrl.prestressIndex - 1].P *= s(k);
                found = true;
                break;
            }
            if (!found) throw std::runtime_error("*STAGE_CONTROL 的阶段号不存在");
        }

        const auto stages = solveModel2DStages(mc, type);
        if (stages.empty()) throw std::runtime_error("施工阶段结果为空");
        const DenseVec U = stages.back().U;
        DenseVec y(nt);
        for (int j = 0; j < nt; ++j) {
            const auto& tg = m.stageTargets[j];
            if (tg.kind == "DISP" || tg.kind == "REACTION") {
                const int d = 3 * nid.at(tg.node) + (tg.dof - 1);
                y(j) = (tg.kind == "REACTION") ? stages.back().R(d) : U(d);
            } else if (tg.kind == "FORCE") {
                double val = 0.0;
                bool found = false;
                for (const auto& f : stages.back().forces) {
                    if (f.id != tg.elem) continue;
                    found = true;
                    if (tg.component == "N") val = (tg.end == 1) ? f.N1 : f.N2;
                    else if (tg.component == "V") val = (tg.end == 1) ? f.V1 : f.V2;
                    else if (tg.component == "M") val = (tg.end == 1) ? f.M1 : f.M2;
                    else throw std::runtime_error("*STAGE_TARGET FORCE 未知分量: " + tg.component);
                    break;
                }
                if (!found) throw std::runtime_error("*STAGE_TARGET FORCE 找不到单元");
                y(j) = val;
            } else if (tg.kind == "STRESS") {
                double val = 0.0;
                bool found = false;
                for (const auto& s : stages.back().stresses) {
                    if (s.id != tg.elem) continue;
                    found = true;
                    if (tg.component == "AXIAL") val = (tg.end == 1) ? s.s1 : s.s2;
                    else if (tg.component == "TOP") val = (tg.end == 1) ? s.stop1 : s.stop2;
                    else if (tg.component == "BOT") val = (tg.end == 1) ? s.sbot1 : s.sbot2;
                    else throw std::runtime_error("*STAGE_TARGET STRESS 未知分量: " + tg.component);
                    break;
                }
                if (!found) throw std::runtime_error("*STAGE_TARGET STRESS 找不到单元");
                y(j) = val;
            } else {
                throw std::runtime_error("*STAGE_TARGET 未知类型: " + tg.kind);
            }
        }
        return y;
    };

    auto evalConstraints = [&](const DenseVec& s) -> DenseVec {
        Model2D mc = m;
        for (int k = 0; k < nv; ++k) {
            const auto& ctrl = m.stageControls[k];
            for (auto& st : mc.stages) {
                if (st.id != ctrl.stage) continue;
                st.prestresses[ctrl.prestressIndex - 1].P *= s(k);
                break;
            }
        }
        const auto stages = solveModel2DStages(mc, type);
        DenseVec c(nc);
        for (int k = 0; k < nc; ++k) {
            const auto& sc = m.stageConstraints[k];
            const auto& tg = sc.target;
            double resp = 0.0;
            if (tg.kind == "DISP" || tg.kind == "REACTION") {
                const int d = 3 * nid.at(tg.node) + (tg.dof - 1);
                resp = (tg.kind == "REACTION") ? stages.back().R(d) : stages.back().U(d);
            } else if (tg.kind == "FORCE") {
                for (const auto& f : stages.back().forces) {
                    if (f.id != tg.elem) continue;
                    if (tg.component == "N") resp = (tg.end == 1) ? f.N1 : f.N2;
                    else if (tg.component == "V") resp = (tg.end == 1) ? f.V1 : f.V2;
                    else if (tg.component == "M") resp = (tg.end == 1) ? f.M1 : f.M2;
                    break;
                }
            } else if (tg.kind == "STRESS") {
                for (const auto& sts : stages.back().stresses) {
                    if (sts.id != tg.elem) continue;
                    if (tg.component == "AXIAL") resp = (tg.end == 1) ? sts.s1 : sts.s2;
                    else if (tg.component == "TOP") resp = (tg.end == 1) ? sts.stop1 : sts.stop2;
                    else if (tg.component == "BOT") resp = (tg.end == 1) ? sts.sbot1 : sts.sbot2;
                    break;
                }
            }
            double violation = 0.0;
            if (resp < sc.lower) violation = resp - sc.lower;
            else if (resp > sc.upper) violation = resp - sc.upper;
            c(k) = sc.target.weight * violation;
        }
        return c;
    };

    for (int iter = 0; iter < maxIter; ++iter) {
        const DenseVec y = eval(x);
        DenseVec r(nt + nc);
        for (int j = 0; j < nt; ++j)
            r(j) = m.stageTargets[j].weight * (y(j) - m.stageTargets[j].value);
        if (nc > 0) r.tail(nc) = evalConstraints(x);
        if (r.norm() < tol) {
            StageIterationResult res;
            res.scales = x;
            res.finalValues = y;
            res.iterations = iter;
            Model2D mc = m;
            for (int k = 0; k < nv; ++k) {
                for (auto& st : mc.stages) {
                    if (st.id == m.stageControls[k].stage) {
                        st.prestresses[m.stageControls[k].prestressIndex - 1].P *= x(k);
                        break;
                    }
                }
            }
            res.stages = solveModel2DStages(mc, type);
            return res;
        }

        Eigen::MatrixXd J(nt + nc, nv);
        for (int j = 0; j < nv; ++j) {
            const double h = 1e-6 * std::max(1.0, std::abs(x(j)));
            DenseVec xp = x, xm = x;
            xp(j) += h;
            xm(j) -= h;
            const DenseVec yp = eval(xp);
            const DenseVec ym = eval(xm);
            DenseVec gp = (yp - ym) / (2.0 * h);
            for (int j = 0; j < nt; ++j) gp(j) *= m.stageTargets[j].weight;
            if (nc > 0) {
                const DenseVec cp = evalConstraints(xp);
                const DenseVec cm = evalConstraints(xm);
                const DenseVec gpc = (cp - cm) / (2.0 * h);
                DenseVec gfull(nt + nc);
                gfull.head(nt) = gp;
                gfull.tail(nc) = gpc;
                J.col(j) = gfull;
            } else {
                J.col(j) = gp;
            }
        }
        const double lambda = 1e-8;
        Eigen::MatrixXd A = J.transpose() * J;
        for (int j = 0; j < nv; ++j) A(j, j) += lambda * std::max(1.0, A(j, j));
        const Eigen::VectorXd rhs = -(J.transpose() * r);
        const Eigen::VectorXd dx = A.fullPivLu().solve(rhs);
        if (!dx.allFinite())
            throw std::runtime_error("多变量正装-倒拆迭代线性求解失败");
        x += dx;
        for (int k = 0; k < nv; ++k) {
            x(k) = std::max(m.stageControls[k].lower, std::min(m.stageControls[k].upper, x(k)));
        }
    }
    throw std::runtime_error("多变量正装-倒拆迭代未收敛");
}

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

inline double evaluateInfluenceResponse(const Model2D& m, const ModelResult& r,
                                        const Model2D::InfluenceRequest& spec) {
    std::map<int, int> nid;
    int n = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = n++;

    if (spec.kind == "DISP") {
        const int d = 3 * nid.at(spec.node) + (spec.dof - 1);
        return r.U(d);
    }
    if (spec.kind == "REACTION") {
        const int d = 3 * nid.at(spec.node) + (spec.dof - 1);
        return r.R(d);
    }

    const std::vector<ElementForce> forces = elementForces2D(m, r.U);
    if (spec.kind == "FORCE") {
        for (const auto& f : forces) {
            if (f.id != spec.elem) continue;
            if (spec.component == "N") return spec.end == 1 ? f.N1 : f.N2;
            if (spec.component == "V") return spec.end == 1 ? f.V1 : f.V2;
            if (spec.component == "M") return spec.end == 1 ? f.M1 : f.M2;
            throw std::runtime_error("INFLUENCE FORCE 未知分量: " + spec.component);
        }
        throw std::runtime_error("INFLUENCE FORCE 找不到单元: " + std::to_string(spec.elem));
    }

    if (spec.kind == "STRESS") {
        const std::vector<ElementStress> stresses = elementStress2D(m, forces);
        for (const auto& s : stresses) {
            if (s.id != spec.elem) continue;
            if (spec.component == "AXIAL") return spec.end == 1 ? s.s1 : s.s2;
            if (spec.component == "TOP") return spec.end == 1 ? s.stop1 : s.stop2;
            if (spec.component == "BOT") return spec.end == 1 ? s.sbot1 : s.sbot2;
            throw std::runtime_error("INFLUENCE STRESS 未知分量: " + spec.component);
        }
        throw std::runtime_error("INFLUENCE STRESS 找不到单元: " + std::to_string(spec.elem));
    }

    throw std::runtime_error("INFLUENCE 未知类型: " + spec.kind);
}

inline InfluenceCurve computeInfluenceCurve(const Model2D& m,
                                             const Model2D::LaneLoad& lane,
                                             const Model2D::InfluenceRequest& spec,
                                             SolverType type = SolverType::SimplicialLDLT) {
    if (lane.nodes.empty())
        throw std::runtime_error("LANE_LOAD 节点列表为空");

    InfluenceCurve out;
    out.label = spec.kind + " id=" + std::to_string(spec.id);
    out.s.reserve(lane.nodes.size());
    out.eta.reserve(lane.nodes.size());

    const int loadDof = (lane.dir == "GX") ? 1 : 2;

    // 车道纵向里程（与求解方式无关，先算出来）。
    double s = 0.0;
    for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
        if (k > 0) {
            const Eigen::Vector2d pa = m.nodes.at(lane.nodes[k - 1]);
            const Eigen::Vector2d pb = m.nodes.at(lane.nodes[k]);
            s += (pb - pa).norm();
        }
        out.s.push_back(s);
    }

    // 索/压杆的激活状态依赖位移，且非线性分析每步都要重算切线，
    // 这两种情况只能逐位置重新求解。
    bool needsReSolve = (m.analysisType == "NONLINEAR");
    for (const auto& e : m.elements)
        if (e.type == "CABLE" || e.type == "STRUT") { needsReSolve = true; break; }

    if (needsReSolve) {
        const ModelResult base = solveModel2D(m, type);
        const double baseValue = evaluateInfluenceResponse(m, base, spec);
        for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
            Model2D mc = m;
            mc.loads.push_back({lane.nodes[k], loadDof, 1.0});
            const ModelResult r = solveModel2D(mc, type);
            out.eta.push_back(evaluateInfluenceResponse(m, r, spec) - baseValue);
        }
    } else {
        // 线性模型：K 与加载位置无关，装配一次、分解一次，
        // 之后每个位置只做一次三角回代（旧实现每个位置都重新装配+重新分解）。
        const Assembled2D a = assembleModel2DLinear(m);

        // 旧实现在这里硬编码 SimplicialLDLT，直接无视调用方传入的 type，
        // 于是 CHOLMOD 等后端在影响线扫描里根本用不上。改用可复用句柄：
        // 既尊重 type，也把 analyze/factorize 分离的能力真正用起来。
        ReusableSolver solver(type);
        solver.factorizeReuse(a.K);

        const DenseVec Ubase = solver.solve(a.F);
        const ModelResult base = makeResult2D(m, a, Ubase, a.F0);
        const double baseValue = evaluateInfluenceResponse(m, base, spec);

        DenseVec F0k(a.nDofs), Fk(a.nDofs);
        for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
            F0k = a.F0;
            Fk  = a.F;
            const int d = 3 * a.index.at(lane.nodes[k]) + (loadDof - 1);
            // 加载点若被约束，单位荷载由支座直接承担，对结构无影响。
            if (!a.isConstrained[d]) { F0k(d) += 1.0; Fk(d) += 1.0; }
            const DenseVec U = solver.solve(Fk);
            const ModelResult r = makeResult2D(m, a, U, F0k);
            out.eta.push_back(evaluateInfluenceResponse(m, r, spec) - baseValue);
        }
    }

    out.vmax = out.eta[0];
    out.vmin = out.eta[0];
    out.smax = out.s[0];
    out.smin = out.s[0];
    for (std::size_t k = 0; k < out.eta.size(); ++k) {
        if (out.eta[k] > out.vmax) { out.vmax = out.eta[k]; out.smax = out.s[k]; }
        if (out.eta[k] < out.vmin) { out.vmin = out.eta[k]; out.smin = out.s[k]; }
    }
    return out;
}

inline std::vector<double> laneTributaryWeights(const Model2D& m,
                                                const Model2D::LaneLoad& lane) {
    const std::size_t n = lane.nodes.size();
    std::vector<double> w(n, 0.0);
    if (n == 1) {
        w[0] = 1.0;
        return w;
    }

    std::vector<double> s(n, 0.0);
    for (std::size_t k = 1; k < n; ++k) {
        const Eigen::Vector2d a = m.nodes.at(lane.nodes[k - 1]);
        const Eigen::Vector2d b = m.nodes.at(lane.nodes[k]);
        s[k] = s[k - 1] + (b - a).norm();
    }

    w[0] = 0.5 * (s[1] - s[0]);
    for (std::size_t k = 1; k + 1 < n; ++k)
        w[k] = 0.5 * (s[k + 1] - s[k - 1]);
    w[n - 1] = 0.5 * (s[n - 1] - s[n - 2]);
    return w;
}

inline EnvelopeValue envelopeFromInfluence(const Model2D& m,
                                           const Model2D::LaneLoad& lane,
                                           const InfluenceCurve& curve) {
    const std::vector<double> w = laneTributaryWeights(m, lane);

    double posArea = 0.0, negArea = 0.0;
    for (std::size_t k = 0; k < curve.eta.size(); ++k) {
        if (curve.eta[k] > 0.0) posArea += curve.eta[k] * w[k];
        else negArea += curve.eta[k] * w[k];
    }

    const double qMaxArea = (lane.q >= 0.0) ? posArea : negArea;
    const double qMinArea = (lane.q >= 0.0) ? negArea : posArea;
    const double pMaxOrd = (lane.P >= 0.0) ? curve.vmax : curve.vmin;
    const double pMinOrd = (lane.P >= 0.0) ? curve.vmin : curve.vmax;
    const double pSMax = (lane.P >= 0.0) ? curve.smax : curve.smin;
    const double pSMin = (lane.P >= 0.0) ? curve.smin : curve.smax;

    EnvelopeValue ev;
    ev.label = curve.label;
    ev.maxVal = lane.q * qMaxArea + lane.P * pMaxOrd;
    ev.minVal = lane.q * qMinArea + lane.P * pMinOrd;
    ev.smax = pSMax;
    ev.smin = pSMin;
    return ev;
}

inline EnvelopeValue nonlinearEnvelopeFromLane(const Model2D& m,
                                               const Model2D::LaneLoad& lane,
                                               const Model2D::InfluenceRequest& spec,
                                               SolverType type = SolverType::SimplicialLDLT) {
    // 先用线性影响线确定 q 的正负加载区，然后用真实 q、P 做非线性工况扫描。
    const InfluenceCurve lin = computeInfluenceCurve(m, lane, spec, type);
    const std::vector<double> w = laneTributaryWeights(m, lane);
    const int loadDof = (lane.dir == "GX") ? 1 : 2;

    EnvelopeValue ev;
    ev.label = "NONLINEAR " + lin.label;
    ev.maxVal = -std::numeric_limits<double>::max();
    ev.minVal = std::numeric_limits<double>::max();
    ev.smax = 0.0;
    ev.smin = 0.0;

    for (std::size_t p = 0; p < lane.nodes.size(); ++p) {
        Model2D mcMax = m;
        Model2D mcMin = m;

        for (std::size_t j = 0; j < lane.nodes.size(); ++j) {
            const double qj = lane.q * w[j];
            const double prod = lane.q * lin.eta[j];
            if (prod > 0.0)
                mcMax.loads.push_back({lane.nodes[j], loadDof, qj});
            else if (prod < 0.0)
                mcMin.loads.push_back({lane.nodes[j], loadDof, qj});
        }

        mcMax.loads.push_back({lane.nodes[p], loadDof, lane.P});
        mcMin.loads.push_back({lane.nodes[p], loadDof, lane.P});

        const ModelResult rMax = solveModel2D(mcMax, type);
        const double vMax = evaluateInfluenceResponse(m, rMax, spec);
        if (vMax > ev.maxVal) { ev.maxVal = vMax; ev.smax = lin.s[p]; }

        const ModelResult rMin = solveModel2D(mcMin, type);
        const double vMin = evaluateInfluenceResponse(m, rMin, spec);
        if (vMin < ev.minVal) { ev.minVal = vMin; ev.smin = lin.s[p]; }
    }
    return ev;
}

inline std::vector<EnvelopeSpec> defaultEnvelopeSpecs(const Model2D& m) {
    std::vector<EnvelopeSpec> specs;
    int id = 0;
    for (const auto& e : m.elements) {
        if (isRigid(e.type)) continue;

        for (int end = 1; end <= 2; ++end) {
            Model2D::InfluenceRequest req;
            req.id = ++id;
            req.kind = "FORCE";
            req.elem = e.id;
            req.end = end;
            req.component = "N";
            specs.push_back({req, "E" + std::to_string(e.id) + "_N" + std::to_string(end)});

            if (e.type == "BEAM") {
                req.id = ++id;
                req.component = "V";
                specs.push_back({req, "E" + std::to_string(e.id) + "_V" + std::to_string(end)});
                req.id = ++id;
                req.component = "M";
                specs.push_back({req, "E" + std::to_string(e.id) + "_M" + std::to_string(end)});
            }
        }

        if (isAxial(e.type)) {
            for (int end = 1; end <= 2; ++end) {
                Model2D::InfluenceRequest req;
                req.id = ++id;
                req.kind = "STRESS";
                req.elem = e.id;
                req.end = end;
                req.component = "AXIAL";
                specs.push_back({req, "E" + std::to_string(e.id) + "_S" + std::to_string(end)});
            }
        } else if (e.type == "BEAM" && e.Wt > 0.0 && e.Wb > 0.0) {
            for (int end = 1; end <= 2; ++end) {
                const std::vector<std::string> comps = {"TOP", "BOT"};
                for (const std::string& comp : comps) {
                    Model2D::InfluenceRequest req;
                    req.id = ++id;
                    req.kind = "STRESS";
                    req.elem = e.id;
                    req.end = end;
                    req.component = comp;
                    specs.push_back({req, "E" + std::to_string(e.id) + "_" + comp + std::to_string(end)});
                }
            }
        }
    }
    return specs;
}

// 把节点位移、单元内力、支座反力写入结果文件。
inline void writeResult(const std::string& path, const Model2D& model,
                        const ModelResult& result,
                        const std::vector<ElementForce>& forces) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法写入结果文件: " + path);

    out << std::setprecision(8);
    out << "# FEM 计算结果\n";
    out << "# 模型: " << (model.title.empty() ? "未命名" : model.title) << "\n\n";

    out << "节点位移 (m, rad)\n";
    out << "节点\tu(x)\t\tv(y)\t\ttheta\n";
    for (std::size_t i = 0; i < result.nodeIds.size(); ++i) {
        const int d = 3 * static_cast<int>(i);
        out << result.nodeIds[i] << '\t'
            << result.U(d + 0) << '\t' << result.U(d + 1) << '\t'
            << result.U(d + 2) << '\n';
    }

    out << "\n单元内力 (N, N*m)\n";
    out << "单元\t类型\tN1(轴力i)\tN2(轴力j)\tV1(剪力)\tM1(弯矩)\tV2(剪力)\tM2(弯矩)\n";
    for (const auto& f : forces) {
        if (isRigid(f.type)) continue;
        out << f.id << '\t' << f.type << '\t' << f.N1 << '\t' << f.N2;
        if (f.type == "BEAM") {
            out << '\t' << f.V1 << '\t' << f.M1 << '\t' << f.V2 << '\t' << f.M2;
        } else {
            out << "\t0\t0\t0\t0";
        }
        out << '\n';
    }

    try {
        const std::vector<ElementStress> stresses = elementStress2D(model, forces);
        out << "\n单元应力 (Pa)\n";
        out << "单元\t类型\tσ1\t\tσ2\t\tσtop1\t\tσbot1\t\tσtop2\t\tσbot2\n";
        for (const auto& s : stresses) {
            if (isRigid(s.type)) continue;
            if (isAxial(s.type)) {
                out << s.id << "\t" << s.type << "\t" << s.s1 << '\t' << s.s2 << "\t-\t-\t-\t-\n";
            } else {
                out << s.id << "\tBEAM\t-\t-\t" << s.stop1 << '\t' << s.sbot1
                    << '\t' << s.stop2 << '\t' << s.sbot2 << '\n';
            }
        }
    } catch (const std::runtime_error& e) {
        out << "\n单元应力: (跳过，原因: " << e.what() << ")\n";
    }

    out << "\n支座反力 (N, N*m)\n";
    out << "节点\tdof\t反力\n";
    double maxAbsR = 0.0;
    for (int d = 0; d < result.R.size(); ++d) maxAbsR = std::max(maxAbsR, std::abs(result.R(d)));
    const double rThresh = 1e-9 * std::max(1.0, maxAbsR); // 相对阈值，避免大模型漏噪声、小模型滤真值
    for (std::size_t i = 0; i < result.nodeIds.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            const double v = result.R(3 * static_cast<int>(i) + d);
            if (std::abs(v) > rThresh)
                out << result.nodeIds[i] << '\t' << (d + 1) << '\t' << v << '\n';
        }
    }
}

} // namespace fem
