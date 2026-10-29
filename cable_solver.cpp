// cable_solver.cpp —— 斜拉桥成桥态求解器的命令行接缝
// =====================================================
// 为什么单独做一个 exe，而不是塞进 BridgeStudio 的 Python：
//   成桥态的三种工况（静力 / 定索力找形 / 成桥态闭环）都要用到 P3 内核的
//   共转列式非线性、索垂度、影响矩阵、最小弯曲能量法。Python 侧只该做两件
//   事——**生成 .dat** 和**解析结果**，任何力学都在这里。
//
// 输入：FEM-Eigen 3D new format .dat（readModel3D 可读）＋ 两个扩展段：
//
//     *CABLE_BRIDGE
//     MODE=COMPLETION          ; COMPLETION | TENSION | STATIC
//     NLSTEPS=3                ; 非线性荷载步
//     NEWTON=30                ; 每个荷载步的 Newton 上限
//     MAXITER=20               ; 外层迭代（索长/几何闭合）上限
//     TOL=1e-8                 ; 几何闭合容差 m
//     TOLT=1e-7                ; 索力相对容差
//     METHOD=MINENERGY|RIGID   ; 自动定索力方法（仅 COMPLETION 且未全给索力时用）
//     MINT=0.02                ; 保力下限 = 该比例 × 最大索力
//     VERBOSE=0
//
//     *CABLE
//     <elemId> <梁侧锚点nodeId> [目标索力N]     ; 每行一根；索力缺省/为 0 = 自动
//
//   若省略 *CABLE 段，则自动取模型里全部 CABLE 单元，锚点取 z 较低的那端。
//
// 输出：与 .dat 同目录 <stem>_cable.txt，分段纯文本（Python 侧顺序解析）：
//     mode=... nodes=... elements=... cables=... dofs=...
//     converged=... iterations=... elapsed_ms=...
//     [CABLE]  elem anchor T L0 chord Ltarget dL dz
//     [CAMBER] node dx dy dz        制造线形相对设计线形的偏移（m）
//     [ERECTED] node x y z          制造线形绝对坐标
//     [U]      node dof value       成桥位移（相对制造线形）
//     [END]
//
// 编译：随 CMake 一起（目标 cable_solver）。

#include "fem_3d_nl.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace fem;

namespace {

double nowMs() {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

struct CableReq {
    int elem = 0;
    int anchor = 0;        // 梁侧锚点；0 = 由几何自动判定
    double tension = 0.0;  // <=0 = 自动定索力
};

struct Cfg {
    std::string mode = "COMPLETION";
    int nlSteps = 3;
    int newton = 30;
    int maxIter = 20;
    int maxLenIter = 15;
    double tol = 1e-8;
    double tolTension = 1e-7;
    bool minEnergy = true;          // false = 刚性支承连续梁法
    double minTension = 0.0;
    double minTensionFraction = 0.02;
    double deckLoad = 0.0;          // 桥面恒载合计 N（仅用于校核报告，0=未知）
    bool verbose = false;
    std::vector<CableReq> cables;
};

// 只扫 *CABLE_BRIDGE / *CABLE 两个扩展段；其余交给 readModel3D。
// 单独扫一遍比往 readModel3D 里塞桥型专用逻辑干净——内核的 .dat 格式
// 应该是"几何 + 载荷"，而"要哪根索拉多少力"是求解工况，不是模型。
Cfg readCfg(const std::string& path) {
    Cfg c;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开输入文件: " + path);
    std::string line, sec;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        if (line[0] == '*') {
            std::string kw = line.substr(1);
            const std::size_t sp = kw.find_first_of(" \t,");
            if (sp != std::string::npos) kw = kw.substr(0, sp);
            sec = upper(kw);
            continue;
        }
        if (sec == "CABLE_BRIDGE") {
            std::istringstream ss(line);
            std::string kv;
            while (ss >> kv) {
                const std::size_t p = kv.find('=');
                if (p == std::string::npos) continue;
                const std::string key = upper(kv.substr(0, p));
                const std::string val = kv.substr(p + 1);
                if (key == "MODE") c.mode = upper(val);
                else if (key == "NLSTEPS") c.nlSteps = std::max(1, std::stoi(val));
                else if (key == "NEWTON") c.newton = std::max(1, std::stoi(val));
                else if (key == "MAXITER") c.maxIter = std::max(1, std::stoi(val));
                else if (key == "MAXLENITER") c.maxLenIter = std::max(1, std::stoi(val));
                else if (key == "TOL") c.tol = std::stod(val);
                else if (key == "TOLT") c.tolTension = std::stod(val);
                else if (key == "METHOD") {
                    c.minEnergy = (upper(val) != "RIGID");
                } else if (key == "MINT") c.minTensionFraction = std::stod(val);
                else if (key == "MINABS") c.minTension = std::stod(val);
                else if (key == "DECKLOAD") c.deckLoad = std::stod(val);
                else if (key == "VERBOSE") c.verbose = (std::stod(val) != 0.0);
                else std::fprintf(stderr, "[警告] *CABLE_BRIDGE 未知关键字: %s\n", kv.c_str());
            }
        } else if (sec == "CABLE") {
            std::istringstream ss(line);
            CableReq r;
            if (ss >> r.elem >> r.anchor) {
                ss >> r.tension;          // 可缺省
                c.cables.push_back(r);
            }
        }
    }
    return c;
}

const Element3D* findElem(const Model3D& m, int id) {
    for (const auto& e : m.elements) if (e.id == id) return &e;
    return nullptr;
}

// 梁侧锚点 = 索较低的那端（斜拉索从塔顶下到梁面）。塔顶始终高于梁面，
// 这条规则对竖索/斜索都成立；真要反过来用 *CABLE 段显式给锚点即可。
int autoAnchor(const Model3D& m, const Element3D& e) {
    const double z1 = m.nodes.at(e.n1).z(), z2 = m.nodes.at(e.n2).z();
    return (z1 <= z2) ? e.n1 : e.n2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("用法: %s <model.dat>\n"
                    "  读 *CABLE_BRIDGE / *CABLE，输出 <stem>_cable.txt\n", argv[0]);
        return 1;
    }
    const std::string dat = argv[1];
    try {
        Cfg cfg = readCfg(dat);
        Model3D m = readModel3D(dat);          // 参考构形（只读副本，后续传值）

        if (cfg.cables.empty()) {
            for (const auto& e : m.elements) {
                if (e.type != "CABLE") continue;
                cfg.cables.push_back({e.id, autoAnchor(m, e), 0.0});
            }
        }
        if (cfg.cables.empty())
            throw std::runtime_error("模型中没有 CABLE 单元，也没有 *CABLE 段");

        for (auto& r : cfg.cables) {
            const Element3D* e = findElem(m, r.elem);
            if (!e) throw std::runtime_error("索单元不存在: " + std::to_string(r.elem));
            if (e->type != "CABLE")
                throw std::runtime_error("单元 " + std::to_string(r.elem) + " 不是 CABLE");
            if (r.anchor <= 0) r.anchor = autoAnchor(m, *e);
            if (r.anchor != e->n1 && r.anchor != e->n2)
                throw std::runtime_error("锚点 " + std::to_string(r.anchor)
                                         + " 不是索单元 " + std::to_string(e->id) + " 的端点");
        }

        const DofMap3D dm = buildDofMap3D(m);
        const std::size_t nc = cfg.cables.size();

        // 总恒载（竖向向下为正），用于全局平衡校核
        double Wtot = 0.0;
        for (const auto& [node, dof, val] : m.loads) if (dof == 3) Wtot -= val;

        NonlinearOptions3D nl;
        nl.loadSteps = cfg.nlSteps;
        nl.maxIter = cfg.newton;
        nl.tolForce = 1e-8;
        nl.verbose = cfg.verbose;

        bool converged = false;
        int iterations = 0;
        double closureErr = 0.0, tensionErr = 0.0, elapsed = 0.0;
        std::vector<double> T(nc, 0.0), L0(nc, 0.0), chord(nc, 0.0), Ltar(nc, 0.0);
        std::vector<double> u;
        std::map<int, Vec3> erected;
        std::vector<double> camber;
        std::string method = "given";

        const double t0 = nowMs();
        if (cfg.mode == "STATIC") {
            // 只算：索的 L0 来自 .dat（每根索的 *R W=/L0= 或 *EN 覆盖）。
            const NonlinearResult3D r = solveModel3DNonlinear(m, nl);
            converged = r.converged;
            iterations = r.totalIterations;
            u = r.u;
            for (std::size_t k = 0; k < nc; ++k) {
                std::size_t idx = 0;
                while (idx < r.elemIds.size() && r.elemIds[idx] != cfg.cables[k].elem) ++idx;
                if (idx >= r.elemIds.size())
                    throw std::runtime_error("结果缺少索单元 " +
                                             std::to_string(cfg.cables[k].elem));
                T[k] = r.axial[idx];
                chord[k] = r.Ln[idx];
                const Element3D* e = findElem(m, cfg.cables[k].elem);
                L0[k] = (e->L0 > 0.0) ? e->L0
                                      : (m.nodes.at(e->n2) - m.nodes.at(e->n1)).norm();
                Ltar[k] = (m.nodes.at(e->n2) - m.nodes.at(e->n1)).norm();
            }
        } else if (cfg.mode == "TENSION") {
            std::vector<CableTarget3D> tg;
            for (const auto& c : cfg.cables) {
                if (c.tension <= 0.0)
                    throw std::runtime_error(
                        "MODE=TENSION 要求 *CABLE 段给全目标索力（单元 "
                        + std::to_string(c.elem) + " 未给）");
                tg.push_back({c.elem, c.tension});
            }
            CableLengthOptions3D co;
            co.maxIter = cfg.maxLenIter;
            co.tolTension = cfg.tolTension;
            co.nl = nl;
            const CableLengthResult3D r = iterateCableUnstrainedLengths(m, tg, co);
            converged = r.converged;
            iterations = r.iterations;
            T = r.tension;
            L0 = r.L0;
            chord = r.chord;
            Ltar = r.Ltarget;
            u = r.u;
            tensionErr = r.maxTensionErr;
            closureErr = r.maxGeomErr;
            method = "target";
        } else if (cfg.mode == "COMPLETION") {
            std::vector<CompletionCable3D> cc;
            for (const auto& c : cfg.cables) cc.push_back({c.elem, c.anchor, c.tension});
            CompletionStateOptions3D co;
            co.maxIter = cfg.maxIter;
            co.tolClosure = cfg.tol;
            co.verbose = cfg.verbose;
            co.cableForce = cfg.minEnergy
                                ? CompletionStateOptions3D::CableForceMethod::MinBendingEnergy
                                : CompletionStateOptions3D::CableForceMethod::RigidSupport;
            co.minEnergy.minTension = cfg.minTension;
            co.minEnergy.minTensionFraction = cfg.minTensionFraction;
            co.minEnergy.verbose = cfg.verbose;
            co.nl = nl;
            bool anyAuto = false, allAuto = true;
            for (const auto& c : cfg.cables) {
                if (c.tension > 0.0) allAuto = false; else anyAuto = true;
            }
            const CompletionStateResult3D r = solveCompletionState3D(m, cc, co);
            converged = r.converged;
            iterations = r.iterations;
            T = r.tension;
            L0 = r.L0;
            chord = r.chord;
            Ltar = r.Ltarget;
            u = r.u;
            erected = r.erected;
            camber = r.camber;
            closureErr = r.maxClosureErr;
            tensionErr = r.maxTensionErr;
            method = allAuto ? (cfg.minEnergy ? "MinBendingEnergy" : "RigidSupport")
                             : (anyAuto ? "Mixed" : "given");
        } else {
            throw std::runtime_error("未知 MODE: " + cfg.mode
                                     + "（应为 COMPLETION / TENSION / STATIC）");
        }
        elapsed = nowMs() - t0;

        // 全局平衡：索力竖向分量之和 vs 总恒载
        double sumT = 0.0, sumV = 0.0, tmin = 1e300, tmax = -1e300;
        for (std::size_t k = 0; k < nc; ++k) {
            const Element3D* e = findElem(m, cfg.cables[k].elem);
            const int A = cfg.cables[k].anchor;
            const int B = (e->n1 == A) ? e->n2 : e->n1;
            const Vec3 ev = (m.nodes.at(B) - m.nodes.at(A)).normalized();
            sumT += T[k];
            sumV += T[k] * ev.z();
            tmin = std::min(tmin, T[k]);
            tmax = std::max(tmax, T[k]);
        }

        // ---- 输出 ----
        std::string outPath = dat;
        const std::size_t dot = outPath.find_last_of('.');
        const std::size_t sep = outPath.find_last_of("/\\");
        if (dot != std::string::npos && (sep == std::string::npos || dot > sep))
            outPath = outPath.substr(0, dot);
        outPath += "_cable.txt";

        std::FILE* f = std::fopen(outPath.c_str(), "w");
        if (!f) throw std::runtime_error("无法写入 " + outPath);
        std::fprintf(f, "cable_solver v1\n");
        std::fprintf(f, "mode=%s\n", cfg.mode.c_str());
        std::fprintf(f, "tension_method=%s\n", method.c_str());
        std::fprintf(f, "model nodes=%d elements=%d cables=%d dofs=%d\n",
                     static_cast<int>(m.nodes.size()),
                     static_cast<int>(m.elements.size()),
                     static_cast<int>(nc), dm.nDofs);
        std::fprintf(f, "converged=%d iterations=%d elapsed_ms=%.1f\n",
                     converged ? 1 : 0, iterations, elapsed);
        std::fprintf(f, "closure_err=%.6e tension_err=%.6e\n", closureErr, tensionErr);
        std::fprintf(f, "tension_min_kN=%.6g tension_max_kN=%.6g tension_sum_kN=%.6g\n",
                     tmin / 1e3, tmax / 1e3, sumT / 1e3);
        std::fprintf(f, "vertical_sum_kN=%.6g total_load_kN=%.6g balance_pct=%.4f\n",
                     sumV / 1e3, Wtot / 1e3,
                     Wtot != 0.0 ? 100.0 * (sumV - Wtot) / std::abs(Wtot) : 0.0);
        // 桥面恒载校核：索的竖向分量应当 ≤ 桥面恒载，差额由支座承担。
        // 用 total_load（含塔自重）算比例会得到"索只承担一半"的假象——塔自重
        // 由塔底直接传走，永远不经过索。
        if (cfg.deckLoad > 0.0)
            std::fprintf(f, "deck_load_kN=%.6g deck_share_pct=%.2f support_share_kN=%.6g\n",
                         cfg.deckLoad / 1e3, 100.0 * sumV / cfg.deckLoad,
                         (cfg.deckLoad - sumV) / 1e3);

        std::fprintf(f, "[CABLE]\n");
        for (std::size_t k = 0; k < nc; ++k)
            std::fprintf(f, "%d %d %.10g %.10g %.10g %.10g %.10g\n",
                         cfg.cables[k].elem, cfg.cables[k].anchor,
                         T[k], L0[k], chord[k], Ltar[k], chord[k] - Ltar[k]);

        // 预抬高只对节点平动有意义：camber/u 是按自由度排的向量，靠 DofMap 取 dz
        auto dofOf = [&](int node, int a) {
            auto it = dm.index.find(node);
            if (it == dm.index.end()) return -1;
            const int cnt = dm.dofCount[it->second];
            if (a >= cnt) return -1;
            return dm.dofStart[it->second] + a;
        };
        if (!camber.empty()) {
            std::fprintf(f, "[CAMBER]\n");
            for (const auto& [id, x] : m.nodes) {
                (void)x;
                const int dx = dofOf(id, 0), dy = dofOf(id, 1), dz = dofOf(id, 2);
                std::fprintf(f, "%d %.10g %.10g %.10g\n", id,
                             dx >= 0 ? camber[dx] : 0.0,
                             dy >= 0 ? camber[dy] : 0.0,
                             dz >= 0 ? camber[dz] : 0.0);
            }
        }
        if (!erected.empty()) {
            std::fprintf(f, "[ERECTED]\n");
            for (const auto& [id, p] : erected)
                std::fprintf(f, "%d %.10g %.10g %.10g\n", id, p.x(), p.y(), p.z());
        }
        if (!u.empty()) {
            std::fprintf(f, "[U]\n");
            for (const auto& [id, x] : m.nodes) {
                (void)x;
                for (int a = 0; a < 6; ++a) {
                    const int d = dofOf(id, a);
                    if (d >= 0) std::fprintf(f, "%d %d %.10g\n", id, a + 1, u[d]);
                }
            }
        }
        std::fprintf(f, "[END]\n");
        std::fclose(f);

        std::printf("cable_solver: mode=%s cables=%d converged=%d iterations=%d "
                    "%.1f ms\n", cfg.mode.c_str(), static_cast<int>(nc),
                    converged ? 1 : 0, iterations, elapsed);
        std::printf("  索力 %.4g ~ %.4g kN（Σ=%.4g kN，ΣV=%.4g kN vs 恒载 %.4g kN，偏差 %.1f%%）\n",
                    tmin / 1e3, tmax / 1e3, sumT / 1e3, sumV / 1e3, Wtot / 1e3,
                    Wtot != 0.0 ? 100.0 * (sumV - Wtot) / std::abs(Wtot) : 0.0);
        std::printf("  结果: %s\n", outPath.c_str());
        return converged ? 0 : 3;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "cable_solver 失败: %s\n", e.what());
        return 2;
    }
}
