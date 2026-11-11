#include "fem_2d.hpp"
#include "modal_output_2d.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace fem;

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "verify_modal.dat";
        if (argc > 1) path = argv[1];
        // opt-in：只有显式带 --modal-detail 才在 .modal.txt 末尾追加
        // [MODAL2D] 明细段。不带开关时文件与控制台输出逐字节不变，
        // compare_regression.py 的 verify_modal.modal.txt 黄金基线不刷新。
        bool detail = false;
        for (int a = 2; a < argc; ++a)
            if (std::string(argv[a]) == "--modal-detail") detail = true;
        const Model2D model = readModel2D(path);

        const ModalResult2D modal = computeNaturalFrequencies2D(model, 6);
        if (modal.frequencies.empty()) {
            std::cerr << "[提示] 未求得自振频率（可能无自由自由度或质量为零）\n";
            return 1;
        }
        const double mu = impactCoefficientFromFrequency(modal.frequencies[0]);
        const auto bf = computeBucklingFactors2D(model, 3);

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
        if (detail) out << formatModalSection2D(model, modal);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
