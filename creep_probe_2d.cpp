// creep_probe_2d.cpp — 2D 结构时变徐变/收缩只读探针（第二十四次 B4）。
//
// 不新增力学：只调用 fem_2d.hpp / code_loads.hpp 已公开的函数与字段，
// 把 2D 结构时变链的 JTG 3362 附录 C 连续拟合中间量、
// 逐龄期徐变系数 φ(t) 与收缩应变 εcs(t)、AEMM 折减后的 E_eff
// 按 12 位有效数字写成 CSV，供 tests/benchmark/creep_shrinkage_time_2d
// 做「独立 Python 闭式 oracle 逐键核对」。
//
// 与 σl6 表值链的区别（必须保持）：
//   * 本探针读 code_loads.hpp 的**连续 RH 拟合**（φ0=φRH·βF·βT0，
//     βH 用实际 RH 连续值，εcs0 用 1.55(1−(RH/100)³)）——结构整体时变口径；
//   * σl6 用表值 εcs0/φ0 + RH 代表值 55/80——预应力损失口径，二者不相交。
//
// 用法：
//   creep_probe_2d.exe <model.dat> [out.csv]
// 输出（无表头）：
//   M,PHI0,<value>         名义徐变系数 φ0（连续拟合）
//   M,BETAH,<value>        βH（连续 RH，mm→h/100）
//   M,ECS0,<value>         名义收缩系数 εcs0（连续拟合）
//   S,<age>,PHI,<value>    龄期 age 的徐变系数 φ(age)
//   S,<age>,ECS,<value>    龄期 age 的新发展收缩应变 εcs(age,t0,ts)
//   S,<age>,BETAS,<value>  收缩发展系数 βs(age,ts,h)
//   S,<age>,EEFF,<value>   AEMM 折减 E/(1+χ·φ(age))（χ=0.8）
//   E,<reason>             无徐变模型 / 解析失败时，非 0 退出
// 失败（解析 / throw）写 stderr 并返回 2；不使用任何默认参数兜底。
#include "fem_2d.hpp"
#include "code_loads.hpp"

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "用法: creep_probe_2d.exe <model.dat> [out.csv]\n";
        return 2;
    }
    try {
        fem::Model2D m = fem::readModel2D(argv[1]);
        if (!(m.creepRh > 0.0 && m.creepH > 0.0 && m.creepFcm > 0.0)) {
            std::cerr << "错误: 模型未给全 RH/H/FCM，不进入 JTG 3362 连续拟合路径\n";
            return 3;
        }

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

        const fem::Jtg3362Creep c = fem::jtg3362NominalCreep(
            m.creepRh, m.creepH, m.creepFcm, m.creepT0);
        const double ecs0 = fem::jtg3362ShrinkageNominal(m.creepRh, m.creepFcm);
        *os << "M,PHI0," << c.phi0 << "\n";
        *os << "M,BETAH," << c.betaH << "\n";
        *os << "M,ECS0," << ecs0 << "\n";

        // 取模型的阶段终龄期作为采样点；没有阶段时用 T0+28/90/365 三个默认点。
        std::vector<double> ages;
        for (const auto& st : m.stages)
            if (st.time > 0.0) ages.push_back(st.time);
        if (ages.empty()) ages = {m.creepT0 + 28.0, m.creepT0 + 90.0,
                                  m.creepT0 + 365.0};
        const double ts = m.shrinkageTs > 0.0 ? m.shrinkageTs : m.creepT0;
        // E0 取第一个 BEAM 单元的弹性模量（AEMM 折减的基准模量）。
        double e0 = 0.0;
        for (const auto& e : m.elements)
            if (e.type == "BEAM") { e0 = e.E; break; }
        for (double age : ages) {
            const double phi = fem::creepCoefficientAt(m, age);
            const double ecs = fem::shrinkageStrainAt(m, age);
            const double betas = fem::jtg3362ShrinkageBeta(age, ts, m.creepH);
            const double eeff = e0 / (1.0 + 0.8 * phi);
            *os << "S," << age << ",PHI," << phi << "\n";
            *os << "S," << age << ",ECS," << ecs << "\n";
            *os << "S," << age << ",BETAS," << betas << "\n";
            *os << "S," << age << ",EEFF," << eeff << "\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 2;
    }
}
