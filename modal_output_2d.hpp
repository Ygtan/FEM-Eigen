#pragma once
// modal_output_2d.hpp
// ==================
// 把 2D 内核**已经算出**的自振结果写成 opt-in 明细段，供 BridgeStudio 解析。
//
// 这里不做任何特征值计算：K φ = λ M φ 由 fem_2d.cpp 的
// computeNaturalFrequencies2D() 求解，ModalResult2D.modes 已含全阶全局振型；
// 本文件只用 assembleMassMatrix2D(m) 重建同一份 M，再做向量内积：
//     Γ_d = φᵀ M r_d,  eff_d = Γ_d² / (φᵀ M φ),  mtot_d = r_dᵀ M r_d
// 与 3D 的 modal_output_3d.hpp::formatModalSection3D 逐项同构，因此报告里的
// 频率 / 周期 / 振型 / 参与质量必然同源，不会出现「两套频率」。
//
// 方向口径（BASIS=xy，用户 2026-11-11 拍板）：
//   2D 每节点 3 个自由度：1=u(x) 水平、2=v(y) 竖向、3=θz 转角。
//   刚体平动只有 X / Y 两个方向；2D **没有出平面自由度**，故不设 Z 向
//   参与质量（是不存在，不是 0，不写 mz 字段，不补零冒充 3D）。
//   θz 是转角自由度、不是刚体平动，也不设参与质量；S 行的第三个分量只是
//   振型分量的如实展示，不得当成参与质量引用。
//   周期 T = 1/f 由内核写出，Python 侧不得现算。
//
// 段格式（opt-in，追加在既有「== 自振频率 ==」「== 屈曲荷载因子 ==」之后）：
//   [MODAL2D]
//   NMODES=<n> NDOF=<nd> BASIS=xy
//   NOTE=<方向口径与边界说明>
//   MASS X=<kg> Y=<kg>              各方向平动总质量 r_dᵀ M r_d
//   F <i> <f_Hz> <T_s> <mx> <my> <mnorm>
//                                  i 从 1 起；mx/my 为该阶 X/Y 方向有效
//                                  参与质量 Γ_d²/(φᵀ M φ)（kg）；
//                                  mnorm = φᵀ M φ，M-正交归一时 ≈ 1。
//   S <i> <node> <ux> <uy> <theta> 该阶该节点的三个自由度振型分量
//
// 自由度顺序与 fem_2d.cpp 的 assembleLinear2D / assembleMassMatrix2D 一致：
// 3*index[node] + (dof-1)，index 为 m.nodes 按 id 升序的压缩编号。
// 段只在 example_modal.exe --modal-detail 时写出；不加开关时 .modal.txt 逐位
// 不变，compare_regression.py 的 verify_modal.modal.txt 黄金基线无需刷新。

#include "fem_2d.hpp"

#include <algorithm>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fem {

inline std::string formatModalSection2D(const Model2D& m,
                                        const ModalResult2D& modal) {
    std::ostringstream ss;
    ss << std::setprecision(10);
    ss << "[MODAL2D]\n";

    // 与 computeNaturalFrequencies2D / assembleMassMatrix2D 完全相同的索引。
    std::map<int, int> index;
    std::vector<int> ids;
    for (const auto& [id, c] : m.nodes) {
        index[id] = static_cast<int>(ids.size());
        ids.push_back(id);
    }
    const int nd = 3 * static_cast<int>(ids.size());
    const int nm = static_cast<int>(modal.frequencies.size());
    ss << "NMODES=" << nm << " NDOF=" << nd << " BASIS=xy\n";
    ss << "NOTE=2D 无出平面自由度，Z 向参与质量不存在（不是 0）；theta 为转角"
          "自由度，只作振型展示、不设参与质量；T=1/f 由内核写出，Python 侧"
          "不得现算。\n";

    // 不静默：段写了但没阶、振型维度对不上，都显式给 ERROR，不摆空表。
    if (nm == 0) {
        ss << "ERROR=no_modes\n";
        return ss.str();
    }
    if (static_cast<int>(modal.modes.size()) < nm) {
        ss << "ERROR=missing_mode_vectors\n";
        return ss.str();
    }
    for (int i = 0; i < nm; ++i) {
        if (static_cast<int>(modal.modes[i].size()) != nd) {
            ss << "ERROR=mode_size_mismatch\n";
            return ss.str();
        }
    }

    const SparseMat M = assembleMassMatrix2D(m);
    DenseVec r[2], Mr[2];
    double mtot[2] = {0.0, 0.0};
    for (int d = 0; d < 2; ++d) {
        r[d] = DenseVec::Zero(nd);
        for (int nid : ids) r[d][3 * index.at(nid) + d] = 1.0;
        Mr[d] = M * r[d];
        mtot[d] = r[d].dot(Mr[d]);
    }
    ss << "MASS X=" << mtot[0] << " Y=" << mtot[1] << "\n";

    for (int i = 0; i < nm; ++i) {
        const double f = modal.frequencies[i];
        const double period = (f > 0.0) ? 1.0 / f : 0.0;
        const DenseVec& mode = modal.modes[i];
        const DenseVec Mphi = M * mode;
        const double mnorm = mode.dot(Mphi);
        double eff[2] = {0.0, 0.0};
        if (mnorm > 0.0) {
            for (int d = 0; d < 2; ++d) {
                const double gamma = mode.dot(Mr[d]);
                eff[d] = gamma * gamma / mnorm;
            }
        }
        ss << "F " << (i + 1) << " " << f << " " << period << " "
           << eff[0] << " " << eff[1] << " " << mnorm << "\n";
        for (int nid : ids) {
            const int g = 3 * index.at(nid);
            ss << "S " << (i + 1) << " " << nid << " "
               << mode[g] << " " << mode[g + 1] << " " << mode[g + 2] << "\n";
        }
    }
    return ss.str();
}

} // namespace fem
