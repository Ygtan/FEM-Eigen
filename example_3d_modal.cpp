#include "fem_2d.hpp"
#include "fem_3d.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace fem;

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "verify_3d_modal.dat";
        if (argc > 1) path = argv[1];
        const Model3D model = readModel3D(path);

        const ModalResult3D modal = computeNaturalFrequencies3D(model, 6);
        if (modal.frequencies.empty()) {
            std::cerr << "[提示] 未求得自振频率（可能无自由自由度或质量为零）\n";
            return 1;
        }
        const double mu = impactCoefficientFromFrequency(modal.frequencies[0]);
        const std::vector<double> bf = computeBucklingFactors3D(model, 3);

        auto emit = [&](std::ostream& os) {
            os << "== 自振频率 ==\n";
            for (std::size_t i = 0; i < modal.frequencies.size(); ++i)
                os << "f" << (i + 1) << " = " << std::setprecision(6)
                   << modal.frequencies[i] << " Hz\n";
            os << "基频 f1 = " << modal.frequencies[0] << " Hz，冲击系数 μ = " << mu << "\n";
            if (!bf.empty()) {
                os << "== 屈曲荷载因子（Pcr = λ·P_ref）==\n";
                for (std::size_t i = 0; i < bf.size(); ++i)
                    os << "λ" << (i + 1) << " = " << std::setprecision(6) << bf[i] << "\n";
            }
        };
        emit(std::cout);
        std::ofstream out(path.substr(0, path.rfind(".dat")) + ".modal.txt");
        emit(out);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
