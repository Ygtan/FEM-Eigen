#include "fem_io.hpp"
#include "fem_3d.hpp"

#include <iomanip>
#include <iostream>
#include <string>

using namespace fem;

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "bridge_3d_grillage.dat";
        if (argc > 1) path = argv[1];
        const Model3D model = readModel3D(path);

        const ModalResult3D modal = computeNaturalFrequencies3D(model, 6);
        std::cout << "== 自振频率 ==\n";
        for (std::size_t i = 0; i < modal.frequencies.size(); ++i)
            std::cout << "f" << (i + 1) << " = " << std::setprecision(6)
                      << modal.frequencies[i] << " Hz\n";
        if (!modal.frequencies.empty()) {
            const double mu = impactCoefficientFromFrequency(modal.frequencies[0]);
            std::cout << "基频 f1 = " << modal.frequencies[0]
                      << " Hz，冲击系数 μ = " << mu << "\n";
        }

        const std::vector<double> bf = computeBucklingFactors3D(model, 3);
        if (!bf.empty()) {
            std::cout << "== 屈曲荷载因子（Pcr = λ·P_ref）==\n";
            for (std::size_t i = 0; i < bf.size(); ++i)
                std::cout << "λ" << (i + 1) << " = " << std::setprecision(6) << bf[i] << "\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
