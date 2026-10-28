#include "fem_3d.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace fem;

// 模态/屈曲性能基准：生成 nElem 段的三维悬臂梁，计时前 6 阶自振频率与一阶屈曲。
// 用于防止稠密特征求解（O(n^3)）回潮：稀疏求解应近似 O(n·k)，远快于阈值。
int main(int argc, char** argv) {
    setupConsoleUTF8();
    int nElem = 1000;
    if (argc > 1) nElem = std::atoi(argv[1]);
    const double E = 2.06e11, G = 7.9e10, L = 100.0;

    Model3D m;
    for (int i = 0; i <= nElem; ++i)
        m.nodes[i + 1] = Vec3(L * i / nElem, 0.0, 0.0);
    for (int i = 1; i <= nElem; ++i) {
        Element3D e;
        e.id = i; e.type = "BEAM"; e.E = E; e.A = 0.18;
        e.Iy = 0.4; e.Iz = 0.6; e.J = 0.3; e.G = G; e.dens = 7850.0;
        e.n1 = i; e.n2 = i + 1;
        m.elements.push_back(e);
    }
    for (int d = 1; d <= 6; ++d) m.constraints.push_back({1, d, 0.0});   // 悬臂固定端

    const auto t0 = std::chrono::steady_clock::now();
    const ModalResult3D modal = computeNaturalFrequencies3D(m, 6);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    const double f1 = modal.frequencies.empty() ? 0.0 : modal.frequencies[0];
    std::printf("nElem=%d dof=%d 模态耗时=%.1f ms 基频=%.6f Hz\n",
                nElem, 6 * (nElem + 1), ms, f1);

    bool ok = true;
    // 稀疏路径实测约 0.03·nElem ms（模态）/ 0.24·nElem ms（屈曲）；稠密 O(n³) 在
    // nElem=1000 时约需 200 s。阈值取线性预期加足量余量，专门拦稠密回潮。
    const double modalLimitMs  = std::max(2000.0, 10.0 * nElem);
    if (modal.frequencies.empty() || !std::isfinite(f1) || f1 <= 0.0) {
        std::printf("[FAIL] 模态结果无效（f1=%.6g）\n", f1);
        ok = false;
    } else if (ms > modalLimitMs) {
        std::printf("[FAIL] 模态耗时 %.1f ms 超阈值 %.1f ms（疑似稠密回退）\n",
                    ms, modalLimitMs);
        ok = false;
    }

    // 屈曲：固定端 + 自由端轴压，计时一阶屈曲因子
    m.loads.push_back({nElem + 1, 1, -1.0e6});
    const auto b0 = std::chrono::steady_clock::now();
    const std::vector<double> lam = computeBucklingFactors3D(m, 1);
    const auto b1 = std::chrono::steady_clock::now();
    const double bms = std::chrono::duration<double, std::milli>(b1 - b0).count();
    std::printf("nElem=%d dof=%d 屈曲耗时=%.1f ms λ1=%.3f\n",
                nElem, 6 * (nElem + 1), bms, lam.empty() ? 0.0 : lam[0]);

    const double buckleLimitMs = std::max(2000.0, 20.0 * nElem);
    if (lam.empty() || !std::isfinite(lam[0]) || lam[0] <= 0.0) {
        std::printf("[FAIL] 屈曲结果无效\n");
        ok = false;
    } else if (bms > buckleLimitMs) {
        std::printf("[FAIL] 屈曲耗时 %.1f ms 超阈值 %.1f ms（疑似稠密回退）\n",
                    bms, buckleLimitMs);
        ok = false;
    }

    std::printf("%s\n", ok ? "性能基准 PASS" : "性能基准 FAIL");
    return ok ? 0 : 1;
}
