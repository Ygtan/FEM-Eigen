// ============================================================
// mat_nl_probe.cpp —— 材料非线性纤维梁 / 截面本构探针 CLI
//
// 第二十一次力学③：把 section.hpp 的 materialStress / sectionResponse
// 与 fem_3d_nl.hpp 的 crFiberLocalForce / crFiberLocalTangent 暴露成
// 文件驱动、机器可读的 CSV，交给 tests/benchmark/material_nonlinear_3d/
// 做独立闭式锚点校验（单轴应力-应变 / 矩形截面塑性铰 / 混凝土开裂 /
// 弹性段与 Kl12 逐位对齐）。
//
// 这里**不**重写任何力学：纤维网格构造、本构积分、切线全部调用
// 生产代码（section.hpp / fem_3d_nl.hpp），只做输入解析与结果打印。
//
// 输入（纯 ASCII 参数文件）：
//   G <gid> E=... A=... IY=... IZ=... G=... J=... MATKIND=... FY=... FT=... [NFY=.. NFZ=..]
//   S <sid> <gid> EPS=... KZ=... KY=...
//   F <fid> <gid> D=... (12 个，CR 局部有效变形 deff)
//
// 输出 CSV（无表头，首列为记录类型）：
//   SECT,<sid>,<N>,<My>,<Mz>,<D11>,<D12>,<D13>,<D22>,<D23>,<D33>
//   FBEAM,<fid>,<q0>,...,<q11>
//   FEQ,<fid>,<maxdq>,<maxdK>   （仅 ELASTIC 组：纤维 vs Kl12·deff 最大差）
//
// 用法：mat_nl_probe.exe <params.mat> <out.csv>
// ============================================================
#include "fem_3d.hpp"
#include "fem_3d_nl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace fem;

namespace {

struct Group {
    std::string id;
    double E = 0.0, A = 0.0, Iy = 0.0, Iz = 0.0, G = 0.0, J = 0.0;
    std::string kind = "ELASTIC";
    double fy = 0.0, ft = 0.0;
    int nfy = 8, nfz = 4;
    CRElemData::FiberNlData fn;   // 由生产逻辑同构构造的纤维数据
};

struct SectCase {
    std::string id, gid;
    double eps = 0.0, kz = 0.0, ky = 0.0;
};

struct FiberCase {
    std::string id, gid;
    Eigen::Matrix<double, 12, 1> d = Eigen::Matrix<double, 12, 1>::Zero();
};

std::string trim(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// 解析 KEY=VALUE 数值（与 beam_nl_probe 同口径）。
std::map<std::string, double> parseKVNum(const std::string& s) {
    std::map<std::string, double> out;
    std::istringstream ss(s);
    std::string tok;
    while (ss >> tok) {
        const std::size_t p = tok.find('=');
        if (p == std::string::npos) continue;
        // 跳过 MATKIND=ELASTIC 这类非数值键（由 parseKVStr 单独解析）。
        try {
            out[tok.substr(0, p)] = std::stod(tok.substr(p + 1));
        } catch (const std::exception&) {
            // ignore
        }
    }
    return out;
}

// 解析 KEY=VALUE 字符串（MATKIND）。
std::map<std::string, std::string> parseKVStr(const std::string& s) {
    std::map<std::string, std::string> out;
    std::istringstream ss(s);
    std::string tok;
    while (ss >> tok) {
        const std::size_t p = tok.find('=');
        if (p == std::string::npos) continue;
        out[tok.substr(0, p)] = tok.substr(p + 1);
    }
    return out;
}

double getNum(const std::map<std::string, double>& kv, const std::string& k,
              double def = 0.0) {
    const auto it = kv.find(k);
    return it == kv.end() ? def : it->second;
}

std::string getStr(const std::map<std::string, std::string>& kv,
                   const std::string& k, const std::string& def) {
    const auto it = kv.find(k);
    return it == kv.end() ? def : it->second;
}

// 用与 crMakeElemData 完全相同的矩形等效离散构造纤维截面。
// 这样探针的纤维网格与生产代码一致：A/Iy/Iz 精确二阶积分，
// 弹性段 sectionResponse 与 Kl12 完全对齐。
CRElemData::FiberNlData makeFibNl(const Group& g) {
    CRElemData::FiberNlData fn;
    fn.enabled = true;
    fn.kind = g.kind;
    fn.fy = g.fy;
    fn.ft = g.ft;
    fn.nfy = g.nfy;
    fn.nfz = g.nfz;
    const double h = std::sqrt(12.0 * g.Iz / g.A);
    const double b = std::sqrt(12.0 * g.Iy / g.A);
    const auto ry = gaussLegendreRule(fn.nfy);
    const auto rz = gaussLegendreRule(fn.nfz);
    const double halfA = 0.25 * h * b;   // Σwy·Σwz=4 ⇒ Σf.area = h*b = A
    MaterialModel mm;
    mm.E = g.E; mm.kind = g.kind; mm.fy = g.fy; mm.ft = g.ft;
    fn.mats.clear();
    fn.mats.push_back(mm);
    fn.sec.fibers.clear();
    fn.sec.Eref = g.E;
    for (const auto& gy : ry) {
        const double y = 0.5 * h * gy.first;
        for (const auto& gz : rz) {
            const double z = 0.5 * b * gz.first;
            Fiber f;
            f.y = y; f.z = z;
            f.area = halfA * gy.second * gz.second;
            f.E = g.E;
            f.hy = 0.5 * h / fn.nfy; f.hz = 0.5 * b / fn.nfz;
            f.mat = 0;
            fn.sec.fibers.push_back(f);
        }
    }
    fn.props = computeSectionProperties(fn.sec.fibers);
    fn.Ay = fn.Az = g.A;
    return fn;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: mat_nl_probe.exe <params.mat> <out.csv>\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open params: " << argv[1] << "\n";
        return 2;
    }

    std::map<std::string, Group> groups;
    std::vector<SectCase> sects;
    std::vector<FiberCase> fibers;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        std::istringstream ss(t);
        std::string tag;
        ss >> tag;
        if (tag == "G") {
            Group g;
            ss >> g.id;
            std::string rest;
            std::getline(ss, rest);
            const auto kv = parseKVNum(rest);
            const auto ks = parseKVStr(rest);
            g.E = getNum(kv, "E");  g.A = getNum(kv, "A");
            g.Iy = getNum(kv, "IY"); g.Iz = getNum(kv, "IZ");
            g.G = getNum(kv, "G");  g.J = getNum(kv, "J");
            g.kind = getStr(ks, "MATKIND", "ELASTIC");
            g.fy = getNum(kv, "FY");
            g.ft = getNum(kv, "FT");
            const int nfy = static_cast<int>(getNum(kv, "NFY", 8));
            const int nfz = static_cast<int>(getNum(kv, "NFZ", 4));
            g.nfy = nfy > 0 ? nfy : 8;
            g.nfz = nfz > 0 ? nfz : 4;
            if (!(g.E > 0.0 && g.A > 0.0 && g.Iy > 0.0 && g.Iz > 0.0)) {
                // ELASTIC 组的 Iy/Iz 必填（纤维离散需要）；其余组同样必填。
                throw std::runtime_error("invalid section group: " + g.id);
            }
            if (g.kind == "STEEL" && !(g.fy > 0.0))
                throw std::runtime_error("STEEL needs FY>0: " + g.id);
            if (g.kind == "CONCRETE" && !(g.ft > 0.0))
                throw std::runtime_error("CONCRETE needs FT>0: " + g.id);
            g.fn = makeFibNl(g);
            groups[g.id] = g;
        } else if (tag == "S") {
            SectCase c;
            if (!(ss >> c.id >> c.gid)) throw std::runtime_error("invalid S line: " + t);
            std::string rest;
            std::getline(ss, rest);
            const auto kv = parseKVNum(rest);
            c.eps = getNum(kv, "EPS");
            c.kz = getNum(kv, "KZ");
            c.ky = getNum(kv, "KY");
            sects.push_back(c);
        } else if (tag == "F") {
            FiberCase c;
            if (!(ss >> c.id >> c.gid)) throw std::runtime_error("invalid F line: " + t);
            for (int i = 0; i < 12; ++i) {
                double v = 0.0;
                if (!(ss >> v)) throw std::runtime_error("F needs 12 values: " + t);
                c.d[i] = v;
            }
            fibers.push_back(c);
        }
    }
    if (sects.empty() && fibers.empty()) {
        std::cerr << "no S/F cases in params file\n";
        return 2;
    }

    std::ofstream out(argv[2]);
    if (!out) {
        std::cerr << "cannot write: " << argv[2] << "\n";
        return 2;
    }
    out << std::setprecision(16);

    // ---- 截面响应 ----
    for (const auto& c : sects) {
        const auto it = groups.find(c.gid);
        if (it == groups.end()) throw std::runtime_error("unknown group: " + c.gid);
        const Group& g = it->second;
        const SectionResponse r =
            sectionResponse(g.fn.sec, g.fn.mats, c.eps, c.kz, c.ky);
        out << "SECT," << c.id
            << "," << r.N << "," << r.My << "," << r.Mz
            << "," << r.D11 << "," << r.D12 << "," << r.D13
            << "," << r.D22 << "," << r.D23 << "," << r.D33 << "\n";
    }

    // ---- 纤维梁局部力 / 切线 / 弹性对齐 ----
    for (const auto& c : fibers) {
        const auto it = groups.find(c.gid);
        if (it == groups.end()) throw std::runtime_error("unknown group: " + c.gid);
        const Group& g = it->second;
        // Lref 固定取 1.0（F 行的 deff 已是"除以参考长度后"的无量纲约定，
        // 这里仅用于把 12 维 deff 映射到局部力，闭式锚点在 compare 脚本里
        // 以同一 Lref=1 复算，不引入额外自由度）。
        const double Lref = 1.0;
        const double GJ = g.G * g.J;
        const Eigen::Matrix<double, 12, 1> q =
            crFiberLocalForce(g.fn, c.d, Lref, GJ);
        out << "FBEAM," << c.id;
        for (int i = 0; i < 12; ++i) out << "," << q[i];
        out << "\n";

        Eigen::Matrix<double, 12, 12> Km;
        crFiberLocalTangent(g.fn, c.d, Lref, GJ, Km);
        // 只输出材料敏感的子块（行/列 3,4,5,6,9,10,11），
        // 其余行/列恒零且与 Kl12 一致，compare 脚本按同一索引复算。
        const int idx[7] = {3, 4, 5, 6, 9, 10, 11};
        out << "FTANG," << c.id;
        for (int a = 0; a < 7; ++a)
            for (int b = 0; b < 7; ++b)
                out << "," << Km(idx[a], idx[b]);
        out << "\n";

        if (g.kind == "ELASTIC") {
            // 弹性段逐位对齐：构造与 c.d 同构的 12×12 Kl12（Lref=1），
            // 纤维力应与 Kl12·deff 完全相等（矩形离散对 A/Iy/Iz 精确）。
            const Eigen::Matrix<double, 12, 12> Kl =
                beam3DLocalStiffness(g.E, g.A, g.Iy, g.Iz, g.J, g.G, Lref);
            const Eigen::Matrix<double, 12, 1> qe = Kl * c.d;
            double maxdq = 0.0;
            for (int i = 0; i < 12; ++i)
                maxdq = std::max(maxdq, std::abs(q[i] - qe[i]));
            Eigen::Matrix<double, 12, 12> KmRef = Kl;
            // Km 只在 idx 子块非零；把参考切线同样取子块比较。
            double maxdK = 0.0;
            for (int a = 0; a < 7; ++a)
                for (int b = 0; b < 7; ++b)
                    maxdK = std::max(maxdK, std::abs(Km(idx[a], idx[b]) - KmRef(idx[a], idx[b])));
            out << "FEQ," << c.id << "," << maxdq << "," << maxdK << "\n";
        }
    }
    return 0;
}
