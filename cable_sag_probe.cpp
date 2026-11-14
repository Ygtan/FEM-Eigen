// ============================================================
// cable_sag_probe.cpp —— 索垂度 / Ernst 等效弹模单元级探针 CLI
//
// 目的：把 fem_3d.hpp / fem_3d_nl.hpp 里的**内核本构闭式**暴露成机器可读
// CSV，交给 tests/benchmark/cable_sag_ansys_3d/ 与 ANSYS 分段 LINK180
// 索链做端到端对标。这里不重新推导任何公式：T 用
// cableTensionFromChord()，Eeq 用 cableErnstModulus()，两者与
// cable_solver.exe 找形 / 施工阶段用的是同一份代码。
//
// 为什么需要分段索链：ANSYS LINK180 是**无垂度**轴向杆，单根单元无法
// 产生 Ernst 软化。把一根索离散成 N 段 LINK180、按水平均布施加节点
// 自重后，NLGEOM 自己会收敛出悬链线，等效轴向柔度即包含垂度效应。
// 因此本 CLI 的 T(L0) 曲线 / Ernst 切线可以被 ANSYS 独立验证。
//
// 输入（纯 ASCII 参数文件，由 BridgeStudio/ansys_bridge_sag.py 同源导出）：
//     E=<Pa>            # 弹性模量
//     A=<m2>            # 截面面积
//     LH=<m>            # 水平跨径（两端等高，弦长 Ln = LH）
//     NSEG=<int>        # ANSYS 分段数（本 CLI 只记录，不参与闭式）
//     C <id> <W> <L0> <amp>
//         W   = 单位水平长度重量 N/m（0 => 退化为直杆）
//         L0  = 无应力（弧）长度 m
//         amp = ANSYS 初始几何正弦扰动幅值 m，仅随 CSV 透传
//
// 输出 CSV（首行表头）：
//     case,W,L0,amp,T,f,Eeq,Eeq_ratio,dTdL0,sag_ratio,s
//       T   = cableTensionFromChord(LH, L0, EA, W, LH)
//       f   = W*LH^2/(8T)                 抛物线垂度
//       Eeq = cableErnstModulus(E,A,W,LH,T,L0)
//       dTdL0 = -(1+T/EA)*Eeq*A/L0      Ln 固定时对 L0 的解析切线
//       s   = W^2*LH^3/(24T^2)         垂度缩短量
//
// 用法：cable_sag_probe.exe <params.sag> <out.csv>
// ============================================================
#include "fem_3d.hpp"
#include "fem_3d_nl.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Case {
    std::string id;
    double W = 0.0;
    double L0 = 0.0;
    double amp = 0.0;
};

std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

bool parseKV(const std::string& line, const std::string& key, double& out) {
    std::string t = trim(line);
    if (t.rfind(key + "=", 0) != 0) return false;
    out = std::stod(t.substr(key.size() + 1));
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: cable_sag_probe.exe <params.sag> <out.csv>\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open params: " << argv[1] << "\n";
        return 2;
    }
    double E = 0.0, A = 0.0, LH = 0.0, NSEG = 0.0;
    std::vector<Case> cases;
    std::string line;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        double v = 0.0;
        if (parseKV(t, "E", v)) { E = v; continue; }
        if (parseKV(t, "A", v)) { A = v; continue; }
        if (parseKV(t, "LH", v)) { LH = v; continue; }
        if (parseKV(t, "NSEG", v)) { NSEG = v; continue; }
        if (t.rfind("C ", 0) == 0 || t.rfind("C\t", 0) == 0) {
            std::istringstream ss(t);
            Case c;
            char tag = 0;
            if (ss >> tag >> c.id >> c.W >> c.L0 >> c.amp) cases.push_back(c);
        }
    }
    if (!(E > 0.0) || !(A > 0.0) || !(LH > 0.0) || cases.empty()) {
        std::cerr << "invalid params: E/A/LH must be positive and >=1 case\n";
        return 2;
    }
    const double EA = E * A;

    std::ofstream out(argv[2]);
    if (!out) {
        std::cerr << "cannot write: " << argv[2] << "\n";
        return 2;
    }
    out << "# cable_sag_probe: kernel closed-form (fem_3d.hpp / fem_3d_nl.hpp)\n";
    out << "# E=" << std::setprecision(15) << E << " A=" << A
        << " LH=" << LH << " NSEG=" << NSEG << "\n";
    out << "case,W,L0,amp,T,f,Eeq,Eeq_ratio,dTdL0,sag_ratio,s\n";
    out << std::setprecision(12);
    for (const auto& c : cases) {
        const double W = c.W, L0 = c.L0;
        if (!(L0 > 0.0)) {
            std::cerr << "case " << c.id << ": L0 must be positive\n";
            return 2;
        }
        // 内核闭式：由弦长反解索力（W=0 时自动退化为直杆）
        const double T = fem::cableTensionFromChord(LH, L0, EA, W, LH);
        const double f = (W > 0.0) ? W * LH * LH / (8.0 * T) : 0.0;
        const double Eeq = fem::cableErnstModulus(E, A, W, LH, T, L0);
        // g(T,L0)=L0(1+T/EA) - cSag/T^2 - Ln = 0，Ln 固定时对 L0 隐式求导：
        //   dT/dL0 = -(1+T/EA) / (L0/EA + 2cSag/T^3)
        //          = -(1+T/EA) * Eeq*A/L0
        const double dTdL0 = -(1.0 + T / EA) * Eeq * A / L0;
        const double s = (W > 0.0) ? W * W * LH * LH * LH / (24.0 * T * T) : 0.0;
        out << c.id << "," << W << "," << L0 << "," << c.amp << ","
            << T << "," << f << "," << Eeq << "," << Eeq / E << ","
            << dTdL0 << "," << f / LH << "," << s << "\n";
    }
    return 0;
}
