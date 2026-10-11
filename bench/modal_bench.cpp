#include "fem_3d.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>

using namespace fem;

// 模态性能基准：生成 nElem 段的三维悬臂梁，计时前 6 阶自振频率。
// 用于防止稠密特征求解（O(n^3)）回潮。
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
    return 0;
}
