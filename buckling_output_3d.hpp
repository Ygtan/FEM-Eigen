#pragma once
// buckling_output_3d.hpp
// ======================
// 把内核**已经算出来**的线性屈曲荷载因子写成文本段，供 BridgeStudio 解析。
//
// 这里不做任何特征值计算：(K + λ·Kg)·φ = 0 由 fem_3d.hpp 的
// computeBucklingResult3D 求解（Spectra 稀疏 shift-invert，失败回退稠密）。
// 本文件只负责格式化同一份 λ 及其诊断，和 modal_output_3d.hpp 是同一种
// "暴露内核已有能力"的做法，不在 Python 侧另组装 K/Kg 重算一遍。
//
// 段格式（每行定长字段，数值 %.10g）：
//   <header>                       "buckling"（summary.txt）/ "[BUCKLING]"（cable.txt）
//   NMODES=<n> NDOF=<nd> BASIS=linear_elastic
//   REF NLOAD=<条数> PREF=<Σ|Fi| N> SUMFZ=<Σ Fz N> NBEAM=<n> NAXIAL=<n>
//       KGNZ=<Kg 非零元数> NFREE=<自由自由度> COMPRESSION=<0|1>
//       CHOLESKY=<0|1> STATUS=<ok|no_compression|no_modes|no_free_dof|cholesky_failed>
//       PREF 是参考荷载的"规模"（各节点荷载绝对值之和，不限方向），
//       λ·PREF 才是临界荷载；SUMFZ 只作诊断，竖向是 dof 3。
//       COMPRESSION 表示 Kg 对角存在受压软化项 —— 没有受压杆件时
//       "算不出 λ" 是物理预期，不是求解失败，STATUS 必须区分开。
//   L <i> <lambda> <Pcr_kN>       i 从 1 起；Pcr = λ·PREF/1000（kN）
//
// 口径（必须写进报告，不能只给一个 λ 数字）：
//   λ 是**参考荷载整体放大倍数**：λ > 1 表示成桥态参考荷载下未失稳，
//   λ ≤ 1 表示参考荷载已达到/超过临界荷载。它不是规范的"稳定安全系数"，
//   规范也没给 λ 的合格限值；线性屈曲不含初始缺陷、材料非线性与几何
//   非线性，是真实稳定承载力的**上限**，不能替代非线性稳定全过程分析。
//
// 索承体系（cable_solver）走共转列式几何非线性成桥态，切线刚度已含成桥态
// 几何刚度，这个基于未变形构型的一次特征值 λ 不适用。那边只写 ERROR 段，
// 下游明确标注"本次未产出"，绝不拿线弹性 λ 冒充非线性稳定结论。

#include "fem_3d.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace fem {

struct BucklingReference3D {
    int nload = 0;
    double prefN = 0.0;     // Σ |Fi|（N）
    double sumFzN = 0.0;    // Σ Fz（N，带符号；dof 3 = 竖向 w）
    int nbeam = 0;
    int naxial = 0;         // TRUSS / CABLE 单元数
};

inline BucklingReference3D bucklingReference3D(const Model3D& m) {
    BucklingReference3D r;
    r.nload = static_cast<int>(m.loads.size());
    for (const auto& [node, dof, val] : m.loads) {
        (void)node;
        r.prefN += std::abs(val);
        if (dof == 3) r.sumFzN += val;
    }
    for (const auto& e : m.elements) {
        if (e.type == "BEAM") ++r.nbeam;
        else if (e.type == "TRUSS" || e.type == "CABLE") ++r.naxial;
    }
    return r;
}

// STATUS 的判定只有这一处，Python 侧只读不算，避免两边各判一遍漂移。
inline std::string bucklingStatus3D(const BucklingResult3D& r) {
    if (!r.factors.empty()) return "ok";
    if (r.nFree == 0) return "no_free_dof";
    if (!r.choleskyOk) return "cholesky_failed";
    if (!r.compression) return "no_compression";
    return "no_modes";
}

inline std::string formatBucklingSection3D(const Model3D& ref,
                                           const DofMap3D& dm,
                                           const BucklingResult3D& result,
                                           const std::string& header) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << header << "\n";
    const BucklingReference3D r = bucklingReference3D(ref);
    ss << "NMODES=" << result.factors.size() << " NDOF=" << dm.nDofs
       << " BASIS=linear_elastic\n";
    ss << "REF NLOAD=" << r.nload << " PREF=" << r.prefN
       << " SUMFZ=" << r.sumFzN
       << " NBEAM=" << r.nbeam << " NAXIAL=" << r.naxial
       << " KGNZ=" << result.kgNonZeros
       << " NFREE=" << result.nFree
       << " COMPRESSION=" << (result.compression ? 1 : 0)
       << " CHOLESKY=" << (result.choleskyOk ? 1 : 0)
       << " STATUS=" << bucklingStatus3D(result) << "\n";
    for (std::size_t i = 0; i < result.factors.size(); ++i)
        ss << "L " << (i + 1) << " " << result.factors[i] << " "
           << result.factors[i] * r.prefN / 1000.0 << "\n";
    return ss.str();
}

// 算不出屈曲时也要写一段，不能让下游"没解析到就当没这回事"——不静默。
inline std::string bucklingErrorSection3D(const std::string& header,
                                          const std::string& code) {
    std::ostringstream ss;
    ss << header << "\n"
       << "NMODES=0 NDOF=0 BASIS=linear_elastic\n"
       << "ERROR=" << code << "\n";
    return ss.str();
}

} // namespace fem
