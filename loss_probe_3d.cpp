// loss_probe_3d.cpp — 3D 分项预应力损失只读探针（第二十三次 B4）。
//
// 不新增力学：只调用 fem_3d.hpp 的 readModel3D + resolvePrestressLosses3D，
// 把 3D 内核已算出的 σl1~σl6 / σpe / RATIO 与附录 G 反摩擦中间量按 12 位
// 有效数字写成 CSV，供 tests/benchmark/prestress_loss_anti_2d 做
// 「同一组 *PRESTRESS_LOSS 参数 2D/3D 同值」端到端锚点。
//
// 用法：
//   loss_probe_3d.exe <model.dat> [out.csv]
// 输出（无表头）：
//   S,<stage>,<beam>,<KEY>,<value>
//   E,<reason>                      （没有 LOSSREF 束时，非 0 退出）
// 失败（解析 / 内核 throw）写 stderr 并返回 2；不使用任何默认参数兜底。
#include "fem_3d.hpp"

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "用法: loss_probe_3d.exe <model.dat> [out.csv]\n";
        return 2;
    }
    try {
        fem::Model3D m = fem::readModel3D(argv[1]);
        fem::resolvePrestressLosses3D(m);

        std::ofstream f;
        std::ostream* os = &std::cout;
        if (argc >= 3) {
            f.open(argv[2], std::ios::out | std::ios::trunc);
            if (!f) {
                std::cerr << "无法写出: " << argv[2] << "\n";
                return 2;
            }
            os = &f;
        }
        *os << std::setprecision(12);
        int n = 0;
        for (const auto& st : m.stages) {
            for (const auto& ps : st.prestresses) {
                if (!ps.hasLossDetail) continue;
                const fem::PrestressLossResult& d = ps.lossDetail;
                const auto it = m.prestressLosses.find(ps.lossRef);
                const double sc = (it == m.prestressLosses.end())
                    ? 0.0 : it->second.sigmaCon;
                const std::string tag = std::to_string(st.id) + ","
                                      + std::to_string(ps.id) + ",";
                *os << "S," << tag << "SIGMA_CON," << sc << "\n";
                *os << "S," << tag << "SL1," << d.sigmaL1 << "\n";
                *os << "S," << tag << "SL2," << d.sigmaL2 << "\n";
                *os << "S," << tag << "SL3," << d.sigmaL3 << "\n";
                *os << "S," << tag << "SL4," << d.sigmaL4 << "\n";
                *os << "S," << tag << "SL5," << d.sigmaL5 << "\n";
                *os << "S," << tag << "SL6," << d.sigmaL6 << "\n";
                *os << "S," << tag << "SPE," << d.sigmaPe << "\n";
                *os << "S," << tag << "RATIO,"
                    << (sc > 0.0 ? d.sigmaPe / sc : 0.0) << "\n";
                *os << "S," << tag << "LF," << d.lf << "\n";
                *os << "S," << tag << "DSL2D," << d.deltaSigmaD << "\n";
                *os << "S," << tag << "DSL2PEAK," << d.deltaSigmaPeak << "\n";
                *os << "S," << tag << "BRANCH," << d.l2Branch << "\n";
                // 附录 C 收缩徐变口径（仅 cSource>=1 输出，旧路径冻结行逐位不变）。
                if (d.cSource >= 1) {
                    *os << "S," << tag << "C_SOURCE," << d.cSource << "\n";
                    if (d.cSource == 1) {
                        *os << "S," << tag << "ECS_DEV," << d.epsCsDev << "\n";
                        *os << "S," << tag << "PHI_DEV," << d.phiDev << "\n";
                        *os << "S," << tag << "BETAS," << d.betaS << "\n";
                        *os << "S," << tag << "BETAS0," << d.betaS0 << "\n";
                        *os << "S," << tag << "BETAC," << d.betaC << "\n";
                        *os << "S," << tag << "BETAH," << d.betaH << "\n";
                        *os << "S," << tag << "RHREP," << d.cHRepresentative << "\n";
                    }
                }
                ++n;
            }
        }
        if (n == 0) {
            *os << "E,no_lossref_prestress\n";
            return 3;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 2;
    }
}
