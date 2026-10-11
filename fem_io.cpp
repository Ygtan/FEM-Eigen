#include "fem_io.hpp"

namespace fem {

bool isAxial(const std::string& type) {
    return type == "TRUSS" || type == "CABLE" || type == "STRUT";
}

bool isRigid(const std::string& type) {
    return type == "RIGID" || type == "RIGID_PIN";
}

DenseMat rigidPinStiffness(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                  double penalty) {
    const Eigen::Vector2d r = p2 - p1; // 梁节点 -> 筋节点

    // DOF: [u_beam, v_beam, θ_beam, u_tendon, v_tendon, θ_tendon]
    // 前两个约束为：u_t - u_b + r_y·θ_b = 0
    //               v_t - v_b - r_x·θ_b = 0
    Eigen::Matrix<double, 2, 6> C = Eigen::Matrix<double, 2, 6>::Zero();
    C(0, 0) = -1.0; C(0, 2) =  r(1); C(0, 3) = 1.0;
    C(1, 1) = -1.0; C(1, 2) = -r(0); C(1, 4) = 1.0;

    return penalty * (C.transpose() * C);
}

double axialForce(const Model2D::Element& e, const DenseVec& U,
                         int i1, int i2,
                         const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const Eigen::Vector2d c = d / L;
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const double dL = (u2 - u1).dot(c);
    return e.E * e.A / L * dL + e.E * e.A * e.istrn;
}

std::vector<int> elementGdofs2D(const std::string& type, int i1, int i2) {
    if (isAxial(type))
        return {3 * i1 + 0, 3 * i1 + 1, 3 * i2 + 0, 3 * i2 + 1};
    return {3 * i1 + 0, 3 * i1 + 1, 3 * i1 + 2,
            3 * i2 + 0, 3 * i2 + 1, 3 * i2 + 2};
}

Eigen::Vector2d resolveDir(const std::string& dir,
                                  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const std::string t = detail::toUpper(detail::trim(dir));
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const double c = d(0) / L, s = d(1) / L;
    if (t == "GX") return {1.0, 0.0};
    if (t == "GY") return {0.0, 1.0};
    if (t == "LX") return {c, s};
    if (t == "LY") return {-s, c};
    if (t.size() >= 2 && (t[0] == 'G' || t[0] == 'L')) {
        const double th = std::stod(t.substr(1)) * kPi / 180.0;
        const double co = std::cos(th), si = std::sin(th);
        if (t[0] == 'G') return {co, si};
        return {c * co - s * si, s * co + c * si}; // 局部 -> 全局
    }
    throw std::runtime_error("未知荷载方向: " + dir);
}

DenseVec distributedLoadEquivalent(const std::string& type,
                                          const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                          const Eigen::Vector2d& gdir,
                                          double q1, double q2) {
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const double c = d(0) / L, s = d(1) / L;

    if (isAxial(type)) {
        DenseVec fe(4);
        fe << L / 6.0 * (2.0 * q1 + q2) * gdir(0),
              L / 6.0 * (2.0 * q1 + q2) * gdir(1),
              L / 6.0 * (q1 + 2.0 * q2) * gdir(0),
              L / 6.0 * (q1 + 2.0 * q2) * gdir(1);
        return fe;
    }

    // BEAM：全局方向分解到局部轴向/横向。
    const double lx =  c * gdir(0) + s * gdir(1);
    const double ly = -s * gdir(0) + c * gdir(1);
    const double qx1 = q1 * lx, qx2 = q2 * lx;
    const double qy1 = q1 * ly, qy2 = q2 * ly;

    DenseVec fl(6);
    fl << L / 6.0 * (2.0 * qx1 + qx2),
          L / 20.0 * (7.0 * qy1 + 3.0 * qy2),
          L * L / 60.0 * (3.0 * qy1 + 2.0 * qy2),
          L / 6.0 * (qx1 + 2.0 * qx2),
          L / 20.0 * (3.0 * qy1 + 7.0 * qy2),
          -L * L / 60.0 * (2.0 * qy1 + 3.0 * qy2);

    EulerBeam2D be{1.0, 1.0, 1.0, p1, p2}; // 只需其坐标变换
    return be.transformation().transpose() * fl;
}

double localTransverseUDL(const Model2D& m, const Model2D::Element& e,
                          const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    if (L <= 0.0) return 0.0;
    const double c = d(0) / L, s = d(1) / L;

    double qy = 0.0;   // 局部横向（+局部 y，向下为负）
    auto add = [&](const Eigen::Vector2d& gdir, double q) {
        const double ly = -s * gdir(0) + c * gdir(1);
        qy += q * ly;
    };

    // 单元均布/梯形荷载（取两端平均，近似为均布）
    for (const auto& el : m.elemLoads) {
        if (el.elem != e.id) continue;
        add(resolveDir(el.dir, p1, p2), 0.5 * (el.q1 + el.q2));
    }

    // 自重（自动 + 显式）
    const double gmag = std::sqrt(m.gx * m.gx + m.gy * m.gy);
    if (gmag > 0.0) {
        const Eigen::Vector2d gdir(m.gx / gmag, m.gy / gmag);
        const double wAuto = (e.dens * e.A + e.addmass) * gmag;
        if (wAuto > 0.0) add(gdir, wAuto);
        for (const auto& sw : m.selfWeights) {
            if (!sw.all && sw.elem != e.id) continue;
            double w = 0.0;
            if (sw.rule == "RHO") w = sw.value * gmag * e.A;
            else if (sw.rule == "MASS") w = sw.value * gmag;
            add(gdir, w);
        }
    }

    return -qy;   // 向下为正（产生下缘受拉的正弯矩）
}

double beamShearParameter(const Model2D::Element& e, double L) {
    if (e.nu <= 0.0 || e.K <= 0.0 || e.A <= 0.0) return 0.0;
    return 24.0 * (1.0 + e.nu) * e.I / (e.K * e.A * L * L);
}

BeamRelease2D releaseBeam2D(const DenseMat& Kl, const DenseVec& Fq,
                            const std::vector<int>& releasedLocal) {
    BeamRelease2D out;
    out.relLocal = releasedLocal;
    for (int i = 0; i < 6; ++i)
        if (std::find(releasedLocal.begin(), releasedLocal.end(), i) == releasedLocal.end())
            out.freeLocal.push_back(i);
    if (releasedLocal.empty()) {
        out.Kl = Kl;
        out.Fq = Fq;
        return out;
    }

    const int nf = static_cast<int>(out.freeLocal.size());
    const int nr = static_cast<int>(releasedLocal.size());
    DenseMat Kff(nf, nf), Kfr(nf, nr), Krr(nr, nr);
    DenseVec Ff(nf), Fr(nr);
    for (int a = 0; a < nf; ++a) {
        Ff(a) = Fq(out.freeLocal[a]);
        for (int b = 0; b < nf; ++b) Kff(a, b) = Kl(out.freeLocal[a], out.freeLocal[b]);
        for (int b = 0; b < nr; ++b) Kfr(a, b) = Kl(out.freeLocal[a], releasedLocal[b]);
    }
    for (int a = 0; a < nr; ++a) {
        Fr(a) = Fq(releasedLocal[a]);
        for (int b = 0; b < nr; ++b) Krr(a, b) = Kl(releasedLocal[a], releasedLocal[b]);
    }

    out.KrrInv = Krr.inverse();
    const DenseMat Krf = Kfr.transpose();
    out.KrrInvKrf = out.KrrInv * Krf;           // K_rr^{-1} K_rf
    out.Kl = Kff - Kfr * out.KrrInvKrf;         // 凝聚刚度
    out.Fq = Ff - Kfr * (out.KrrInv * Fr);      // 凝聚等效荷载
    return out;
}

namespace {
// 带端部释放的梁：凝聚后的全局刚度 / 全局自由度 / 全局等效荷载。
struct CondensedBeamGlobal {
    DenseMat Ke;
    std::vector<int> gdofs;
    DenseVec feq;
};
CondensedBeamGlobal condenseBeamGlobal(const Model2D::Element& e,
                                       const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                                       int i1, int i2, const DenseVec& fe) {
    EulerBeam2D be{e.E, e.A, e.I, p1, p2};
    const DenseMat T = be.transformation();
    const DenseVec Fq = T * fe;   // 全局 -> 局部
    const double Phi = beamShearParameter(e, (p2 - p1).norm());
    const BeamRelease2D rel = releaseBeam2D(be.localStiffness(Phi), Fq, e.releases);
    const int nf = static_cast<int>(rel.freeLocal.size());
    DenseMat Tc(nf, 6);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < 6; ++b)
            Tc(a, b) = T(rel.freeLocal[a], b);
    CondensedBeamGlobal out;
    out.Ke = Tc.transpose() * rel.Kl * Tc;
    out.feq = Tc.transpose() * rel.Fq;
    for (int d : rel.freeLocal)
        out.gdofs.push_back(d < 3 ? 3 * i1 + d : 3 * i2 + (d - 3));
    return out;
}

// 带端部释放的梁：反算 6 个局部杆端力（释放端内力强制为 0）。
DenseVec releasedBeamLocalForces(const Model2D& m, const Model2D::Element& e,
                                 const DenseVec& U, int i1, int i2,
                                 const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    EulerBeam2D be{e.E, e.A, e.I, p1, p2};
    const DenseMat T = be.transformation();
    const DenseMat Kl = be.localStiffness(beamShearParameter(e, (p2 - p1).norm()));
    const DenseVec Fq = T * elementEquivalentLoad(m, e, p1, p2);
    const BeamRelease2D rel = releaseBeam2D(Kl, Fq, e.releases);

    Eigen::VectorXd ug(6);
    for (int b = 0; b < 6; ++b)
        ug(b) = (b < 3) ? U(3 * i1 + b) : U(3 * i2 + (b - 3));

    const int nf = static_cast<int>(rel.freeLocal.size());
    const int nr = static_cast<int>(rel.relLocal.size());
    Eigen::VectorXd uf(nf);
    for (int a = 0; a < nf; ++a) {
        const int d = rel.freeLocal[a];
        double s = 0.0;
        for (int b = 0; b < 6; ++b) s += T(d, b) * ug(b);
        uf(a) = s;
    }
    Eigen::VectorXd Fr(nr), KrfUf(nr);
    for (int a = 0; a < nr; ++a) {
        Fr(a) = Fq(rel.relLocal[a]);
        double s = 0.0;
        for (int b = 0; b < nf; ++b)
            s += Kl(rel.freeLocal[b], rel.relLocal[a]) * uf(b);
        KrfUf(a) = s;
    }
    const Eigen::VectorXd ur = rel.KrrInv * (Fr - KrfUf);

    Eigen::VectorXd ul(6);
    for (int a = 0; a < nf; ++a) ul(rel.freeLocal[a]) = uf(a);
    for (int a = 0; a < nr; ++a) ul(rel.relLocal[a]) = ur(a);

    DenseVec fl = Kl * ul - Fq;
    for (int r : rel.relLocal) fl(r) = 0.0;
    return fl;
}
} // namespace

DenseVec uniformTempEquivalent(const std::string& type, double E, double A, double alpha,
                                      double dT, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const double NT = E * A * alpha * dT;
    const Eigen::Vector2d d = p2 - p1;
    const double L = d.norm();
    const Eigen::Vector2d c = d / L;
    if (isAxial(type)) {
        DenseVec fe(4);
        fe << -NT * c(0), -NT * c(1), NT * c(0), NT * c(1);
        return fe;
    }
    DenseVec fl = DenseVec::Zero(6);
    fl(0) = -NT;
    fl(3) = NT;
    EulerBeam2D be{1.0, 1.0, 1.0, p1, p2};
    return be.transformation().transpose() * fl;
}

DenseVec gradientTempEquivalent(double E, double I, double alpha, double h, double dT,
                                       const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    if (h <= 0.0)
        throw std::runtime_error("梁温度梯度需要 *ELEMENT 里的 H");
    const double MT = E * I * alpha * dT / h;
    DenseVec fl = DenseVec::Zero(6);
    fl(2) = -MT;
    fl(5) = MT;
    EulerBeam2D be{1.0, 1.0, 1.0, p1, p2};
    return be.transformation().transpose() * fl;
}

std::pair<double, double> tempProfileResultants(double E, double A, double H,
                                                       double alpha,
                                                       double T1, double T2,
                                                       double A0, double B) {
    if (H <= 0.0) throw std::runtime_error("TEMP_PROFILE 需要梁高 H");
    if (A0 <= 0.0) throw std::runtime_error("TEMP_PROFILE A 必须大于 0");
    const double b = A / H;
    const double B0 = (B <= 0.0 || B > H) ? H : B;
    if (B0 < A0) throw std::runtime_error("TEMP_PROFILE B 不能小于 A");

    auto T = [&](double z) -> double {
        if (z <= 0.0) return T1;
        if (z <= A0)
            return T1 + (T2 - T1) * (z / A0);
        if (z <= B0)
            return T2 + (0.0 - T2) * ((z - A0) / (B0 - A0));
        return 0.0;
    };

    const int n = 400;
    double N = 0.0, M = 0.0;
    for (int i = 0; i < n; ++i) {
        const double z0 = i * H / n;
        const double z1 = (i + 1) * H / n;
        const double zm = 0.5 * (z0 + z1);
        const double y = H * 0.5 - zm; // 向上为正
        const double dz = z1 - z0;
        N += b * T(zm) * dz;
        M += b * T(zm) * y * dz;
    }
    return {E * alpha * N, E * alpha * M};
}

DenseVec tempProfileEquivalent(const Model2D::Element& e,
                                      const Eigen::Vector2d& p1,
                                      const Eigen::Vector2d& p2,
                                      double N, double M) {
    DenseVec fl = DenseVec::Zero(6);
    fl(0) = -N;
    fl(3) = N;
    fl(2) = -M;
    fl(5) = M;
    EulerBeam2D be{e.E, e.A, e.I, p1, p2};
    return be.transformation().transpose() * fl;
}

void addElementLoad(const std::vector<int>& gdofs, const DenseVec& fe, DenseVec& F) {
    for (std::size_t i = 0; i < gdofs.size(); ++i)
        F(gdofs[i]) += fe(i);
}

DenseVec elementEquivalentLoad(const Model2D& m, const Model2D::Element& e,
                                      const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    if (isRigid(e.type))
        return DenseVec::Zero(6);
    const int nd = isAxial(e.type) ? 4 : 6;
    DenseVec fe = DenseVec::Zero(nd);

    for (const auto& el : m.elemLoads) {
        if (el.elem != e.id) continue;
        const Eigen::Vector2d dir = resolveDir(el.dir, p1, p2);
        fe += distributedLoadEquivalent(e.type, p1, p2, dir, el.q1, el.q2);
    }

    const double gmag = std::sqrt(m.gx * m.gx + m.gy * m.gy);
    Eigen::Vector2d gdir(0.0, 0.0);
    if (gmag > 0.0) gdir = Eigen::Vector2d(m.gx / gmag, m.gy / gmag);

    // 自动自重：材料密度×截面积 + 实常数附加质量（线密度 kg/m）。
    // e.addmass 既参与静力自重，也保留在 Element 里供后续动力分析使用。
    const double autoLineMass = e.dens * e.A + e.addmass;
    if (autoLineMass > 0.0) {
        const double w = autoLineMass * gmag;
        fe += distributedLoadEquivalent(e.type, p1, p2, gdir, w, w);
    }

    // 兼容旧格式的显式 *SELF_WEIGHT。
    for (const auto& sw : m.selfWeights) {
        if (!sw.all && sw.elem != e.id) continue;
        double w = 0.0;
        if (sw.rule == "RHO") w = sw.value * gmag * e.A;
        else if (sw.rule == "MASS") w = sw.value * gmag;
        else throw std::runtime_error("未知自重规则: " + sw.rule);
        fe += distributedLoadEquivalent(e.type, p1, p2, gdir, w, w);
    }

    if (m.uniformTemp != 0.0)
        fe += uniformTempEquivalent(e.type, e.E, e.A, e.alpha, m.uniformTemp, p1, p2);

    for (const auto& et : m.elemTemps) {
        if (et.elem != e.id) continue;
        if (et.kind == "UNIFORM")
            fe += uniformTempEquivalent(e.type, e.E, e.A, e.alpha, et.dT, p1, p2);
        else if (et.kind == "GRADIENT")
            fe += gradientTempEquivalent(e.E, e.I, e.alpha, e.h, et.dT, p1, p2);
    }

    // 收缩应变：正值为缩短，按负温度应变等效为节点荷载。
    for (const auto& sh : m.shrinkageSpecs) {
        if (!sh.all && std::find(sh.elems.begin(), sh.elems.end(), e.id) == sh.elems.end())
            continue;
        fe += uniformTempEquivalent(e.type, e.E, e.A, 1.0, -sh.ecs, p1, p2);
    }

    // 规范竖向温度梯度：仅用于梁单元。
    for (const auto& tp : m.tempProfiles) {
        if (e.type != "BEAM") continue;
        if (!tp.all && std::find(tp.elems.begin(), tp.elems.end(), e.id) == tp.elems.end())
            continue;
        const double h = (tp.H > 0.0) ? tp.H : e.h;
        const auto [N, M] = tempProfileResultants(e.E, e.A, h, e.alpha,
                                                  tp.T1, tp.T2, tp.A, tp.B);
        fe += tempProfileEquivalent(e, p1, p2, N, M);
    }

    // 初应变 ISTRN（cable/strut/truss 预应力）：等效节点力 ±EA·ISTRN 沿杆轴。
    if (e.istrn != 0.0 && isAxial(e.type)) {
        const Eigen::Vector2d d = p2 - p1;
        const double L = d.norm();
        const Eigen::Vector2d c = d / L;
        const double P = e.E * e.A * e.istrn;
        fe(0) +=  P * c(0); fe(1) +=  P * c(1);   // 节点 i 沿 +c 方向
        fe(2) += -P * c(0); fe(3) += -P * c(1);   // 节点 j 沿 -c 方向
    }
    return fe;
}

double axialForceNonlinear(const Model2D::Element& e, const DenseVec& U,
                                  int i1, int i2,
                                  const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d0 = p2 - p1;
    const double L0 = d0.norm();
    if (L0 <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const Eigen::Vector2d d = d0 + u2 - u1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
    return e.E * e.A * ((L - L0) / L0 + e.istrn);
}

DenseVec axialInternalForceNonlinear(const Model2D::Element& e, const DenseVec& U,
                                             int i1, int i2,
                                             const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d0 = p2 - p1;
    const double L0 = d0.norm();
    if (L0 <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const Eigen::Vector2d d = d0 + u2 - u1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
    const Eigen::Vector2d c = d / L;
    const double N = e.E * e.A * ((L - L0) / L0 + e.istrn);

    DenseVec fe(4);
    fe << -N * c(0), -N * c(1), N * c(0), N * c(1);
    return fe;
}

DenseMat axialTangentNonlinear(const Model2D::Element& e, const DenseVec& U,
                                      int i1, int i2,
                                      const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    const Eigen::Vector2d d0 = p2 - p1;
    const double L0 = d0.norm();
    if (L0 <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
    const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
    const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
    const Eigen::Vector2d d = d0 + u2 - u1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
    const Eigen::Vector2d c = d / L;
    const double N = e.E * e.A * ((L - L0) / L0 + e.istrn);

    const double kMat = e.E * e.A / L0;
    const double kGeo = N / L;
    const Eigen::Matrix2d A =
        (kMat - kGeo) * (c * c.transpose()) + kGeo * Eigen::Matrix2d::Identity();

    DenseMat Kt = DenseMat::Zero(4, 4);
    Kt.block<2, 2>(0, 0) = A;
    Kt.block<2, 2>(0, 2) = -A;
    Kt.block<2, 2>(2, 0) = -A;
    Kt.block<2, 2>(2, 2) = A;
    return Kt;
}

DenseVec elementExternalEquivalentLoad(const Model2D& m, const Model2D::Element& e,
                                              const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    DenseVec fe = elementEquivalentLoad(m, e, p1, p2);
    if (e.istrn != 0.0 && isAxial(e.type)) {
        const Eigen::Vector2d d = p2 - p1;
        const double L = d.norm();
        if (L <= 0.0) throw std::runtime_error("轴向单元初始长度为零");
        const Eigen::Vector2d c = d / L;
        const double P = e.E * e.A * e.istrn;
        fe(0) -= P * c(0);
        fe(1) -= P * c(1);
        fe(2) += P * c(0);
        fe(3) += P * c(1);
    }
    return fe;
}

void processPrestress(Model2D& m) {
    if (m.prestresses.empty()) return;

    // 预建单元索引，避免 getElem 的 O(n) 线性扫描。
    std::unordered_map<int, Model2D::Element*> elemById;
    elemById.reserve(m.elements.size() * 2);
    for (auto& e : m.elements) elemById[e.id] = &e;

    // 预建约束索引，避免 isConstrained 的 O(n_bc) 线性扫描。
    std::unordered_set<std::uint64_t> constrainedSet;
    constrainedSet.reserve(m.constraints.size() * 2);
    for (const auto& c : m.constraints)
        constrainedSet.insert(static_cast<std::uint64_t>(c.node) * 10ull + static_cast<std::uint64_t>(c.dof));

    auto getElem = [&](int id) -> Model2D::Element* {
        auto it = elemById.find(id);
        return (it != elemById.end()) ? it->second : nullptr;
    };

    auto isConstrained = [&](int node, int dof) -> bool {
        return constrainedSet.find(static_cast<std::uint64_t>(node) * 10ull + static_cast<std::uint64_t>(dof))
               != constrainedSet.end();
    };

    int nextNode = 0, nextElem = 0, nextType = 0, nextMat = 0, nextReal = 0;
    for (const auto& [id, c] : m.nodes) nextNode = std::max(nextNode, id);
    for (const auto& e : m.elements) nextElem = std::max(nextElem, e.id);
    for (const auto& [id, v] : m.elementTypes) nextType = std::max(nextType, id);
    for (const auto& [id, v] : m.materials) nextMat = std::max(nextMat, id);
    for (const auto& [id, v] : m.realConsts) nextReal = std::max(nextReal, id);

    for (const auto& ps : m.prestresses) {
        if (ps.elems.empty())
            throw std::runtime_error("*PRESTRESS 的 ELEM 列表为空");
        double loss = ps.loss;
        if (ps.lossRef > 0) {
            const auto lit = m.prestressLosses.find(ps.lossRef);
            if (lit == m.prestressLosses.end())
                throw std::runtime_error("*PRESTRESS LOSSREF 引用了不存在的 *PRESTRESS_LOSS: "
                                         + std::to_string(ps.lossRef));
            loss = effectivePrestressRatio(lit->second);
        }
        const double Pe = ps.P * loss;

        // 按输入顺序建立梁节点链。
        std::vector<int> beamNodes;
        int prevEnd = 0;
        for (std::size_t i = 0; i < ps.elems.size(); ++i) {
            Model2D::Element* e = getElem(ps.elems[i]);
            if (!e) throw std::runtime_error("*PRESTRESS 引用了不存在的单元: " + std::to_string(ps.elems[i]));
            if (e->type != "BEAM")
                throw std::runtime_error("*PRESTRESS 目前只支持施加到 BEAM 单元");
            if (i == 0) {
                beamNodes.push_back(e->conn[0]);
                beamNodes.push_back(e->conn[1]);
            } else {
                if (e->conn[0] != prevEnd)
                    throw std::runtime_error("*PRESTRESS 的 ELEM 必须按节点顺序连续");
                beamNodes.push_back(e->conn[1]);
            }
            prevEnd = e->conn[1];
        }

        const int startNode = beamNodes.front();
        const int endNode = beamNodes.back();
        const Eigen::Vector2d pstart = m.nodes.at(startNode);
        const Eigen::Vector2d pend = m.nodes.at(endNode);
        const Eigen::Vector2d dspan = pend - pstart;
        const double L = dspan.norm();
        if (L <= 0.0) throw std::runtime_error("*PRESTRESS 梁组长度为零");
        const Eigen::Vector2d t = dspan / L;
        const Eigen::Vector2d n(-t(1), t(0)); // 局部 +y

        std::string method = ps.method;
        if (method == "AUTO")
            method = (m.analysisType == "NONLINEAR") ? "TENDON" : "EQUIV";

        if (method == "EQUIV") {
            if (ps.F != 0.0) {
                const double w = 8.0 * Pe * ps.F / (L * L);
                for (int elem : ps.elems)
                    m.elemLoads.push_back({elem, "LY", w, w});
            }

            if (ps.E1 != 0.0)
                m.loads.push_back({startNode, 3, Pe * ps.E1});
            if (ps.E2 != 0.0)
                m.loads.push_back({endNode, 3, -Pe * ps.E2});

            const bool startAxialFixed = isConstrained(startNode, 1);
            const bool endAxialFixed = isConstrained(endNode, 1);
            if (startAxialFixed && !endAxialFixed)
                m.loads.push_back({endNode, 1, -Pe});
            else if (!startAxialFixed && endAxialFixed)
                m.loads.push_back({startNode, 1, Pe});
            else if (startAxialFixed && endAxialFixed)
                throw std::runtime_error("*PRESTRESS EQUIV：两端 dof1 均固定时无法用外荷载模拟弹性缩短，请使用 METHOD=TENDON");
            else
                throw std::runtime_error("*PRESTRESS EQUIV：梁组至少需要一端约束 dof1");
            continue;
        }

        if (method != "TENDON")
            throw std::runtime_error("*PRESTRESS 未知 METHOD: " + ps.method);
        if (ps.Ap <= 0.0)
            throw std::runtime_error("*PRESTRESS TENDON 需要 APS/AP 预应力筋面积");

        const double istrn = Pe / (ps.Ep * ps.Ap);

        const int tendonType = ++nextType;
        const int linkType = ++nextType;
        const int tendonMat = ++nextMat;
        const int linkMat = ++nextMat;
        const int tendonReal = ++nextReal;
        const int linkReal = ++nextReal;

        m.elementTypes[tendonType] = "CABLE";
        m.elementTypes[linkType] = "RIGID_PIN";

        Model2D::Material tm;
        tm.E = ps.Ep;
        tm.alpha = 0.0;
        tm.dens = 0.0;
        m.materials[tendonMat] = tm;

        Model2D::Material lm;
        lm.E = 1e11;
        lm.alpha = 0.0;
        lm.dens = 0.0;
        m.materials[linkMat] = lm;

        Model2D::RealConst tr;
        tr.A = ps.Ap;
        tr.istrn = istrn;
        m.realConsts[tendonReal] = tr;

        Model2D::RealConst lr;
        lr.A = 1e-2;
        lr.I = 1e-5;
        m.realConsts[linkReal] = lr;

        std::vector<int> tendonNodes;
        tendonNodes.reserve(beamNodes.size());
        for (std::size_t k = 0; k < beamNodes.size(); ++k) {
            const int bn = beamNodes[k];
            const Eigen::Vector2d pb = m.nodes.at(bn);
            const double s = std::max(0.0, std::min(L, (pb - pstart).dot(t)));
            const double xi = L > 0.0 ? s / L : 0.0;
            const double e = (1.0 - xi) * ps.E1 + xi * ps.E2 + 4.0 * ps.F * xi * (1.0 - xi);
            const Eigen::Vector2d pt = pb - e * n;

            if ((pt - pb).norm() < 1e-8) {
                tendonNodes.push_back(bn); // 锚固点或偏心为 0 时直接共用梁节点
                continue;
            }

            const int tn = ++nextNode;
            m.nodes[tn] = pt;
            tendonNodes.push_back(tn);
            // 筋节点只有 u/v 参与受力；θ 为无刚度自由度，固定为 0。
            m.constraints.push_back({tn, 3, 0.0});

            Model2D::Element link;
            link.id = ++nextElem;
            link.type = "RIGID_PIN";
            link.E = lm.E;
            link.A = lr.A;
            link.I = lr.I;
            link.conn = {bn, tn};
            m.elements.push_back(link);
        }

        for (std::size_t k = 0; k + 1 < tendonNodes.size(); ++k) {
            Model2D::Element cable;
            cable.id = ++nextElem;
            cable.type = "CABLE";
            cable.E = ps.Ep;
            cable.A = ps.Ap;
            cable.istrn = istrn;
            cable.conn = {tendonNodes[k], tendonNodes[k + 1]};
            m.elements.push_back(cable);
        }
    }
}

void snapshotBaseModulus(Model2D& m) {
    for (const auto& e : m.elements)
        m.baseModulus[e.id] = e.E;
}

double baseModulusOf(const Model2D& m, const Model2D::Element& e) {
    auto it = m.baseModulus.find(e.id);
    return (it != m.baseModulus.end()) ? it->second : e.E;
}

double creepReductionFactor(double chi, double dphi) {
    const double r = chi * std::max(0.0, dphi);
    return 1.0 / (1.0 + r);
}

double interpolateTimeCurve(const Model2D::TimeCurve& curve, double time) {
    if (curve.points.empty())
        throw std::runtime_error("时变曲线没有数据点");
    std::vector<std::pair<double, double>> pts = curve.points;
    std::sort(pts.begin(), pts.end());
    if (time <= pts.front().first) return pts.front().second;
    if (time >= pts.back().first) return pts.back().second;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (time <= pts[i].first) {
            const double t0 = pts[i - 1].first, t1 = pts[i].first;
            const double v0 = pts[i - 1].second, v1 = pts[i].second;
            return v0 + (v1 - v0) * (time - t0) / (t1 - t0);
        }
    }
    return pts.back().second;
}

double creepCoefficientAt(const Model2D& m, double time) {
    if (m.creepModel == "CURVE" && !m.creepCurve.points.empty())
        return interpolateTimeCurve(m.creepCurve, time);
    const double dt = std::max(0.0, time - m.creepT0);
    return m.creepPhi0 * std::pow(dt / (m.creepBetaH + dt), m.creepPower);
}

double shrinkageStrainAt(const Model2D& m, double time) {
    if (m.shrinkageModel == "CURVE" && !m.shrinkageCurve.points.empty())
        return interpolateTimeCurve(m.shrinkageCurve, time);
    const double dt = std::max(0.0, time - m.shrinkageTs);
    return m.shrinkageEcsInf * std::pow(dt / (m.shrinkageAlphaDs + dt), m.shrinkagePower);
}

void processTimeEffects(Model2D& m,
                               const std::map<int, double>* loadAge,
                               double now) {
    if (m.creepSpecs.empty()) return;
    const bool useAge = (now > 0.0);

    for (auto& e : m.elements) {
        if (isRigid(e.type)) continue;

        double bestChi = 0.0, bestDphi = 0.0;
        bool hit = false;
        for (const auto& cs : m.creepSpecs) {
            if (!cs.all && std::find(cs.elems.begin(), cs.elems.end(), e.id) == cs.elems.end())
                continue;
            double dphi = cs.phi;
            if (useAge) {
                double t0 = 0.0;
                if (loadAge) {
                    auto it = loadAge->find(e.id);
                    if (it != loadAge->end()) t0 = it->second;
                }
                // AEMM：只计加载龄期之后发展的那部分徐变
                dphi = creepCoefficientAt(m, now) - creepCoefficientAt(m, t0);
            }
            const double red = cs.chi * std::max(0.0, dphi);
            if (!hit || red > bestChi * bestDphi) {
                bestChi = cs.chi;
                bestDphi = std::max(0.0, dphi);
                hit = true;
            }
        }
        if (!hit) continue;

        const double E0 = baseModulusOf(m, e);
        e.E = E0 * creepReductionFactor(bestChi, bestDphi);
    }
}

Model2D readModel2D(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开文件: " + path);

    const std::vector<std::string> lines = detail::preprocess(in);

    Model2D m;
    std::string section;
    Model2D::Element legacyTemplate; // 旧式 *ELEMENT 段的内联模板
    int activeType = 0, activeMat = 0, activeReal = 0;

    for (const std::string& line : lines) {
        const std::string t = detail::trim(line);
        if (t.empty() || t[0] == '#' || t[0] == '!') continue;

        if (t[0] == '*') {
            // 关键字 = 第一个空白或逗号之前的部分；其余为参数。
            const std::size_t sep = t.find_first_of(", \t");
            const std::string kw = detail::toUpper(sep == std::string::npos ? t : t.substr(0, sep));
            const std::string rest = sep == std::string::npos ? std::string() : detail::trim(t.substr(sep + 1));

            // ---- 命令：*TYPE / *MAT / *REAL / *EN ----
            if (kw == "*TYPE" || kw == "*MAT" || kw == "*REAL" || kw == "*EN") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                if (kw == "*TYPE") {
                    ss >> activeType;
                } else if (kw == "*MAT") {
                    ss >> activeMat;
                } else if (kw == "*REAL") {
                    ss >> activeReal;
                } else { // *EN,eid,nodei,nodej
                    int eid = 0, ni = 0, nj = 0;
                    if (!(ss >> eid >> ni >> nj))
                        throw std::runtime_error("*EN 行格式错误: " + t);
                    if (activeType <= 0 || activeMat <= 0 || activeReal <= 0)
                        throw std::runtime_error("*EN 前需先激活 *TYPE/*MAT/*REAL");
                    auto tit = m.elementTypes.find(activeType);
                    auto mit = m.materials.find(activeMat);
                    auto rit = m.realConsts.find(activeReal);
                    if (tit == m.elementTypes.end() || mit == m.materials.end() || rit == m.realConsts.end())
                        throw std::runtime_error("*EN 引用的 TYPE/MAT/REAL 未定义");
                    Model2D::Element e;
                    e.id = eid;
                    e.type = tit->second;
                    e.E = mit->second.E;
                    e.alpha = mit->second.alpha;
                    e.dens = mit->second.dens;
                    e.nu = mit->second.nu;
                    e.A = rit->second.A;
                    e.I = rit->second.I;
                    e.Wt = rit->second.Wt;
                    e.Wb = rit->second.Wb;
                    e.h = rit->second.h;
                    e.addmass = rit->second.addmass;
                    e.istrn = rit->second.istrn;
                    e.K = rit->second.K;
                    e.conn = {ni, nj};
                    if (e.type == "BEAM" && e.I <= 0.0)
                        throw std::runtime_error("梁单元缺少实常数 I");
                    if (e.type == "BEAM" && (e.Wt <= 0.0 || e.Wb <= 0.0))
                        throw std::runtime_error("梁单元缺少 WT/WB（截面抵抗矩，应力计算必需）");
                    if (isAxial(e.type) && e.A <= 0.0)
                        throw std::runtime_error("桁架/索/压杆单元缺少实常数 A");
                    m.elements.push_back(e);
                }
                continue;
            }

            // *TEMPERATURE 既支持命令式（同行带值），也支持段式（下一行给值）。
            if (kw == "*TEMPERATURE" && !rest.empty()) {
                m.uniformTemp = std::stod(rest);
                continue;
            }

            // *ANALYSIS：分析类型。LINEAR=方案A（默认），NONLINEAR=方案C。
            // 可选参数：TOL=<相对容差> MAXITER=<最大迭代次数>。
            if (kw == "*ANALYSIS") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                auto p = detail::parseKeyValues(args);
                std::string v;
                std::istringstream probe(args);
                std::string first;
                probe >> first;
                const std::string firstUp = detail::toUpper(detail::trim(first));
                if (firstUp == "LINEAR" || firstUp == "NONLINEAR") {
                    v = firstUp;
                } else {
                    auto it = p.find("TYPE");
                    if (it != p.end())
                        v = detail::toUpper(detail::trim(it->second));
                }

                if (v == "LINEAR" || v == "NONLINEAR")
                    m.analysisType = v;
                else
                    throw std::runtime_error("*ANALYSIS 未知类型: " + rest);

                auto it = p.find("TOL");
                if (it != p.end()) m.nonlinearTol = std::stod(it->second);
                it = p.find("MAXITER");
                if (it != p.end()) m.nonlinearMaxIter = std::stoi(it->second);
                // STEPS：荷载分级步数（默认 1，即不分级）。
                it = p.find("STEPS");
                if (it != p.end()) m.nonlinearLoadSteps = std::stoi(it->second);
                // LINESEARCH：ON/OFF，ΔU 阻尼回溯（默认 ON）。
                it = p.find("LINESEARCH");
                if (it != p.end()) {
                    const std::string s = detail::toUpper(detail::trim(it->second));
                    if (s == "ON" || s == "1" || s == "TRUE" || s == "YES")
                        m.nonlinearLineSearch = true;
                    else if (s == "OFF" || s == "0" || s == "FALSE" || s == "NO")
                        m.nonlinearLineSearch = false;
                    else
                        throw std::runtime_error("*ANALYSIS LINESEARCH 只能为 ON 或 OFF");
                }
                // THROW：不收敛时是否抛异常（默认 ON）。OFF 则返回部分成功结果。
                it = p.find("THROW");
                if (it != p.end()) {
                    const std::string s = detail::toUpper(detail::trim(it->second));
                    if (s == "ON" || s == "1" || s == "TRUE" || s == "YES")
                        m.nonlinearThrowOnFailure = true;
                    else if (s == "OFF" || s == "0" || s == "FALSE" || s == "NO")
                        m.nonlinearThrowOnFailure = false;
                    else
                        throw std::runtime_error("*ANALYSIS THROW 只能为 ON 或 OFF");
                }
                if (m.nonlinearTol <= 0.0)
                    throw std::runtime_error("*ANALYSIS TOL 必须大于 0");
                if (m.nonlinearMaxIter < 1)
                    throw std::runtime_error("*ANALYSIS MAXITER 必须至少为 1");
                if (m.nonlinearLoadSteps < 1)
                    throw std::runtime_error("*ANALYSIS STEPS 必须至少为 1");
                continue;
            }

            if (kw == "*PRESTRESS" || kw == "*PRESTRESS_LOSS") {
                section = kw;
                continue;
            }

            if (kw == "*INFLUENCE" || kw == "*LANE_LOAD") {
                section = kw;
                continue;
            }

            if (kw == "*CE_PARAMETER") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                auto p = detail::parseKeyValues(args);
                auto it = p.find("PENALTY");
                if (it != p.end()) m.cePenalty = std::stod(it->second);
                it = p.find("ALPHA");
                if (it != p.end()) m.cePenalty = std::stod(it->second);
                if (m.cePenalty <= 0.0)
                    throw std::runtime_error("*CE_PARAMETER PENALTY 必须大于 0");
                continue;
            }

            if (kw == "*CE") {
                section = kw;
                continue;
            }

            if (kw == "*CREEP_MODEL" || kw == "*SHRINKAGE_MODEL") {
                std::string args = rest;
                std::replace(args.begin(), args.end(), ',', ' ');
                auto p = detail::parseKeyValues(args);
                if (kw == "*CREEP_MODEL") {
                    auto it = p.find("MODEL"); if (it != p.end()) m.creepModel = detail::toUpper(it->second);
                    it = p.find("PHI0"); if (it != p.end()) m.creepPhi0 = std::stod(it->second);
                    it = p.find("T0"); if (it != p.end()) m.creepT0 = std::stod(it->second);
                    it = p.find("BETA_H"); if (it != p.end()) m.creepBetaH = std::stod(it->second);
                    it = p.find("POWER"); if (it != p.end()) m.creepPower = std::stod(it->second);
                } else {
                    auto it = p.find("MODEL"); if (it != p.end()) m.shrinkageModel = detail::toUpper(it->second);
                    it = p.find("ECS_INF"); if (it != p.end()) m.shrinkageEcsInf = std::stod(it->second);
                    it = p.find("TS"); if (it != p.end()) m.shrinkageTs = std::stod(it->second);
                    it = p.find("ALPHA_DS"); if (it != p.end()) m.shrinkageAlphaDs = std::stod(it->second);
                    it = p.find("POWER"); if (it != p.end()) m.shrinkagePower = std::stod(it->second);
                }
                continue;
            }

            if (kw == "*CREEP" || kw == "*SHRINKAGE") {
                section = kw;
                continue;
            }

            if (kw == "*TEMP_PROFILE") {
                section = kw;
                continue;
            }

            if (kw == "*STAGE" || kw == "*STAGE_LOAD" || kw == "*STAGE_ELEMENT_LOAD" ||
                kw == "*STAGE_PRESTRESS" || kw == "*STAGE_BOUNDARY" || kw == "*STAGE_UNBOUNDARY" ||
                kw == "*STAGE_CREEP" || kw == "*STAGE_SHRINKAGE" || kw == "*STAGE_UNLOAD" ||
                kw == "*STAGE_TARGET" || kw == "*STAGE_CONTROL" || kw == "*STAGE_TEMP_PROFILE" ||
                kw == "*CREEP_CURVE" || kw == "*SHRINKAGE_CURVE" || kw == "*STAGE_CONSTRAINT" ||
                kw == "*CAMBER" || kw == "*FEEDBACK" || kw == "*CREEP_MODEL" || kw == "*SHRINKAGE_MODEL") {
                section = kw;
                continue;
            }

            section = kw;

            if (kw == "*ELEMENT") {
                // 旧式：*ELEMENT, TYPE=..., E=..., A=..., I=..., ALPHA=..., H=..., WT=..., WB=...
                auto p = detail::parseParams(rest);
                legacyTemplate = Model2D::Element();
                legacyTemplate.type = detail::toUpper(p["TYPE"]);
                legacyTemplate.E = detail::toDouble(p, "E");
                legacyTemplate.A = detail::toDouble(p, "A");
                legacyTemplate.I = 0.0;
                legacyTemplate.alpha = 0.0;
                legacyTemplate.h = 0.0;
                legacyTemplate.Wt = 0.0;
                legacyTemplate.Wb = 0.0;
                legacyTemplate.conn.clear();
                auto it = p.find("I");
                if (it != p.end()) legacyTemplate.I = std::stod(it->second);
                it = p.find("ALPHA");
                if (it != p.end()) legacyTemplate.alpha = std::stod(it->second);
                it = p.find("H");
                if (it != p.end()) legacyTemplate.h = std::stod(it->second);
                it = p.find("WT");
                if (it != p.end()) legacyTemplate.Wt = std::stod(it->second);
                it = p.find("WB");
                if (it != p.end()) legacyTemplate.Wb = std::stod(it->second);
            }
            continue;
        }

        // ---- 数据行，按当前 section 解析 ----
        std::istringstream ss(t);
        if (section == "*TITLE") {
            m.title = t;
        } else if (section == "*NODE") {
            int id = 0;
            double x = 0.0, y = 0.0;
            if (!(ss >> id >> x >> y))
                throw std::runtime_error("NODE 行格式错误: " + t);
            m.nodes[id] = {x, y};
        } else if (section == "*ELEMENT") {
            int eid = 0, a = 0, b = 0;
            if (!(ss >> eid >> a >> b))
                throw std::runtime_error("ELEMENT 行格式错误: " + t);
            Model2D::Element e = legacyTemplate;
            e.id = eid;
            e.conn = {a, b};
            if (e.type == "BEAM" && e.I <= 0.0)
                throw std::runtime_error("旧式 *ELEMENT：梁单元缺少 I");
            if (e.type == "BEAM" && (e.Wt <= 0.0 || e.Wb <= 0.0))
                throw std::runtime_error("旧式 *ELEMENT：梁单元缺少 WT/WB（截面抵抗矩，应力计算必需）");
            if (isAxial(e.type) && e.A <= 0.0)
                throw std::runtime_error("旧式 *ELEMENT：桁架/索/压杆单元缺少 A");
            m.elements.push_back(e);
        } else if (section == "*ET") {
            int id = 0;
            std::string kind;
            if (!(ss >> id >> kind))
                throw std::runtime_error("ET 行格式错误: " + t);
            m.elementTypes[id] = detail::toUpper(kind);
        } else if (section == "*MP") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("MP 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Material mat;
            auto it = p.find("E"); if (it != p.end()) mat.E = std::stod(it->second);
            it = p.find("ALPX"); if (it != p.end()) mat.alpha = std::stod(it->second);
            it = p.find("DENS"); if (it == p.end()) it = p.find("DENSITY");   // DENSITY 为 DENS 别名
            if (it != p.end()) mat.dens = std::stod(it->second);
            it = p.find("NU"); if (it != p.end()) mat.nu = std::stod(it->second);
            for (const auto& [k, v] : p)
                if (k != "E" && k != "ALPX" && k != "DENS" && k != "DENSITY" && k != "NU")
                    std::cerr << "[警告] *MP 未知关键字: " << k << "=" << v << "\n";
            m.materials[id] = mat;
        } else if (section == "*R") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("R 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::RealConst rc;
            auto it = p.find("A"); if (it != p.end()) rc.A = std::stod(it->second);
            it = p.find("I"); if (it != p.end()) rc.I = std::stod(it->second);
            it = p.find("WT"); if (it != p.end()) rc.Wt = std::stod(it->second);
            it = p.find("WB"); if (it != p.end()) rc.Wb = std::stod(it->second);
            it = p.find("H"); if (it != p.end()) rc.h = std::stod(it->second);
            it = p.find("ADDMASS"); if (it != p.end()) rc.addmass = std::stod(it->second);
            it = p.find("ISTRN"); if (it != p.end()) rc.istrn = std::stod(it->second);
            it = p.find("K"); if (it != p.end()) rc.K = std::stod(it->second);
            for (const auto& [k, v] : p)
                if (k != "A" && k != "I" && k != "WT" && k != "WB" && k != "H" &&
                    k != "ADDMASS" && k != "ISTRN" && k != "K")
                    std::cerr << "[警告] *R 未知关键字: " << k << "=" << v << "\n";
            m.realConsts[id] = rc;
        } else if (section == "*BOUNDARY") {
            std::string ns, ds, vs;
            if (!(ss >> ns >> ds >> vs))
                throw std::runtime_error("BOUNDARY 行格式错误: " + t);
            const double value = std::stod(vs);
            double k = 0.0;
            std::string ks;
            if (ss >> ks) {
                if (ks.rfind("K=", 0) == 0) ks = ks.substr(2);
                k = std::stod(ks);
            }
            for (int node : detail::parseRange(ns))
                for (int dof : detail::parseRange(ds))
                    m.constraints.push_back({node, dof, value, k});
        } else if (section == "*RELEASE") {
            // 端部释放：<eid> <END=1/I|2/J> <TYPE=N/V/M> [<TYPE>...]
            std::string es, ends;
            if (!(ss >> es >> ends))
                throw std::runtime_error("RELEASE 行格式错误: " + t);
            const std::string eup = detail::toUpper(detail::trim(ends));
            const int end = (eup == "1" || eup == "I") ? 0 : 1;
            if (!(eup == "1" || eup == "I" || eup == "2" || eup == "J"))
                throw std::runtime_error("RELEASE 的 END 必须是 1/I 或 2/J: " + ends);
            const int eid = std::stoi(es);
            auto it = std::find_if(m.elements.begin(), m.elements.end(),
                                   [&](const Model2D::Element& e) { return e.id == eid; });
            if (it == m.elements.end())
                throw std::runtime_error("RELEASE 引用了不存在的单元: " + es);
            if (it->type != "BEAM")
                throw std::runtime_error("端部释放仅支持 BEAM 单元（单元 " + es + "）");
            std::string typ;
            while (ss >> typ) {
                const std::string tt = detail::toUpper(detail::trim(typ));
                int dof = -1;
                if (tt == "N" || tt == "AXIAL" || tt == "U") dof = 0;
                else if (tt == "V" || tt == "SHEAR") dof = 1;
                else if (tt == "M" || tt == "MOMENT" || tt == "ROT") dof = 2;
                else throw std::runtime_error("RELEASE 未知释放类型: " + typ + "（应为 N/V/M）");
                const int ld = end * 3 + dof;
                if (std::find(it->releases.begin(), it->releases.end(), ld) == it->releases.end())
                    it->releases.push_back(ld);
            }
        } else if (section == "*CLOAD") {
            std::string ns, ds, vs;
            if (!(ss >> ns >> ds >> vs))
                throw std::runtime_error("CLOAD 行格式错误: " + t);
            const double value = std::stod(vs);
            for (int node : detail::parseRange(ns))
                for (int dof : detail::parseRange(ds))
                    m.loads.push_back({node, dof, value});
        } else if (section == "*ELEMENT_LOAD") {
            std::string es, dir;
            double q1 = 0.0, q2 = 0.0;
            if (!(ss >> es >> dir >> q1 >> q2))
                throw std::runtime_error("ELEMENT_LOAD 行格式错误: " + t);
            for (int elem : detail::parseRange(es))
                m.elemLoads.push_back({elem, dir, q1, q2});
        } else if (section == "*GRAVITY") {
            double gx = 0.0, gy = 0.0;
            if (!(ss >> gx >> gy))
                throw std::runtime_error("GRAVITY 行格式错误: " + t);
            m.gx = gx;
            m.gy = gy;
        } else if (section == "*SELF_WEIGHT") {
            std::string who, rule;
            double value = 0.0;
            if (!(ss >> who >> rule >> value))
                throw std::runtime_error("SELF_WEIGHT 行格式错误: " + t);
            Model2D::SelfWeight sw;
            sw.rule = detail::toUpper(rule);
            sw.value = value;
            if (detail::toUpper(who) == "ALL") sw.all = true;
            else sw.elem = std::stoi(who);
            m.selfWeights.push_back(sw);
        } else if (section == "*TEMPERATURE") {
            double dT = 0.0;
            if (!(ss >> dT))
                throw std::runtime_error("TEMPERATURE 行格式错误: " + t);
            m.uniformTemp = dT;
        } else if (section == "*ELEMENT_TEMPERATURE") {
            std::string es, kind;
            double dT = 0.0;
            if (!(ss >> es >> kind >> dT))
                throw std::runtime_error("ELEMENT_TEMPERATURE 行格式错误: " + t);
            for (int elem : detail::parseRange(es))
                m.elemTemps.push_back({elem, detail::toUpper(kind), dT});
        } else if (section == "*PRESTRESS_LOSS") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("PRESTRESS_LOSS 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            PrestressLossInput li;
            auto g = [&](const char* key, double& dst) {
                auto it = p.find(key);
                if (it != p.end()) dst = std::stod(it->second);
            };
            g("SIGMA_CON", li.sigmaCon);
            g("EP", li.Ep);
            g("EC", li.Ec);
            g("FPK", li.fpk);
            g("MU", li.mu);
            g("K", li.kappa);
            g("THETA", li.theta);
            g("X", li.x);
            g("DELTA_L", li.sumDeltaL);
            g("L", li.tendonLength);
            g("T1", li.t1);
            g("T2", li.t2);
            g("SIGMA_PC4", li.sigmaPc4);
            g("SIGMA_PC6", li.sigmaPc6);
            g("PSI", li.psi);
            g("ZETA", li.zeta);
            g("ECS", li.epsCs);
            g("PHI", li.phi);
            g("RHO", li.rho);
            g("RHOPS", li.rhoPs);
            auto it = p.find("POST");
            if (it != p.end()) li.postTensioned = (std::stod(it->second) != 0.0);
            auto itRebar = p.find("REBAR");
            if (itRebar != p.end()) li.threadedBar = (std::stod(itRebar->second) != 0.0);
            m.prestressLosses[id] = li;
        } else if (section == "*PRESTRESS") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("PRESTRESS 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Prestress ps;
            ps.id = id;
            auto it = p.find("ELEM");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS 缺少 ELEM");
            ps.elems = detail::parseIntList(it->second);
            it = p.find("P");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS 缺少 P");
            ps.P = std::stod(it->second);
            it = p.find("F"); if (it != p.end()) ps.F = std::stod(it->second);
            it = p.find("E1"); if (it != p.end()) ps.E1 = std::stod(it->second);
            it = p.find("E2"); if (it != p.end()) ps.E2 = std::stod(it->second);
            const bool hasLoss = (p.find("LOSS") != p.end());
            const bool hasLossRef = (p.find("LOSSREF") != p.end());
            if (hasLoss == hasLossRef)
                throw std::runtime_error("*PRESTRESS 必须且只能指定 LOSS 或 LOSSREF 之一");
            if (hasLoss) ps.loss = std::stod(p.find("LOSS")->second);
            if (hasLossRef) ps.lossRef = std::stoi(p.find("LOSSREF")->second);
            it = p.find("EPS"); if (it != p.end()) ps.Ep = std::stod(it->second);
            it = p.find("AP"); if (it != p.end()) ps.Ap = std::stod(it->second);
            it = p.find("APS"); if (it != p.end()) ps.Ap = std::stod(it->second);
            it = p.find("METHOD"); if (it != p.end()) ps.method = detail::toUpper(it->second);
            if (ps.P <= 0.0 || ps.Ep <= 0.0)
                throw std::runtime_error("*PRESTRESS 的 P/EPS 必须大于 0");
            if (ps.method != "AUTO" && ps.method != "EQUIV" && ps.method != "TENDON")
                throw std::runtime_error("*PRESTRESS 未知 METHOD: " + ps.method);
            m.prestresses.push_back(ps);
        } else if (section == "*INFLUENCE") {
            int id = 0;
            std::string kind;
            if (!(ss >> id >> kind))
                throw std::runtime_error("INFLUENCE 行格式错误: " + t);
            Model2D::InfluenceRequest ir;
            ir.id = id;
            ir.kind = detail::toUpper(kind);
            if (ir.kind == "DISP" || ir.kind == "REACTION") {
                if (!(ss >> ir.node >> ir.dof))
                    throw std::runtime_error("INFLUENCE DISP/REACTION 需要 node dof");
            } else if (ir.kind == "FORCE" || ir.kind == "STRESS") {
                std::string comp;
                if (!(ss >> ir.elem >> comp >> ir.end))
                    throw std::runtime_error("INFLUENCE FORCE/STRESS 需要 elem component end");
                ir.component = detail::toUpper(comp);
            } else {
                throw std::runtime_error("INFLUENCE 未知类型: " + kind);
            }
            if (ir.end != 1 && ir.end != 2)
                throw std::runtime_error("INFLUENCE end 只能为 1 或 2");
            m.influences.push_back(ir);
        } else if (section == "*LANE_LOAD") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("LANE_LOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::LaneLoad ll;
            ll.id = id;
            auto it = p.find("NODES");
            if (it == p.end()) throw std::runtime_error("*LANE_LOAD 缺少 NODES");
            ll.nodes = detail::parseIntList(it->second);
            it = p.find("Q"); if (it != p.end()) ll.q = std::stod(it->second);
            it = p.find("P"); if (it != p.end()) ll.P = std::stod(it->second);
            it = p.find("DIR"); if (it != p.end()) ll.dir = detail::toUpper(it->second);
            it = p.find("METHOD"); if (it != p.end()) ll.method = detail::toUpper(it->second);
            it = p.find("CODE"); if (it != p.end()) ll.autoCode = (detail::toUpper(it->second) == "JTG_D60");
            if (ll.dir != "GX" && ll.dir != "GY")
                throw std::runtime_error("*LANE_LOAD DIR 只能为 GX 或 GY");
            if (ll.method != "LINEAR" && ll.method != "NONLINEAR")
                throw std::runtime_error("*LANE_LOAD METHOD 只能为 LINEAR 或 NONLINEAR");
            m.laneLoads.push_back(ll);
        } else if (section == "*CE") {
            Model2D::ConstraintEquation ce;
            if (!(ss >> ce.id >> ce.constant))
                throw std::runtime_error("CE 行格式错误: " + t);
            int node = 0, dof = 0;
            double coef = 0.0;
            while (ss >> node >> dof >> coef) {
                if (dof < 1 || dof > 3)
                    throw std::runtime_error("CE 自由度必须在 1~3 之间");
                ce.terms.push_back({node, dof, coef});
            }
            if (ce.terms.empty())
                throw std::runtime_error("CE 至少需要一个自由度项");
            m.ces.push_back(ce);
        } else if (section == "*CREEP") {
            std::string who;
            if (!(ss >> who))
                throw std::runtime_error("CREEP 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::CreepSpec cs;
            if (detail::toUpper(who) == "ALL") cs.all = true;
            else cs.elems = detail::parseIntList(who);
            auto it = p.find("PHI");
            if (it == p.end()) throw std::runtime_error("*CREEP 缺少 PHI");
            cs.phi = std::stod(it->second);
            it = p.find("CHI"); if (it != p.end()) cs.chi = std::stod(it->second);
            if (cs.phi < 0.0 || cs.chi <= 0.0)
                throw std::runtime_error("*CREEP PHI/CHI 参数错误");
            m.creepSpecs.push_back(cs);
        } else if (section == "*SHRINKAGE") {
            std::string who;
            if (!(ss >> who))
                throw std::runtime_error("SHRINKAGE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::ShrinkageSpec ss_;
            if (detail::toUpper(who) == "ALL") ss_.all = true;
            else ss_.elems = detail::parseIntList(who);
            auto it = p.find("ECS");
            if (it == p.end()) throw std::runtime_error("*SHRINKAGE 缺少 ECS");
            ss_.ecs = std::stod(it->second);
            m.shrinkageSpecs.push_back(ss_);
        } else if (section == "*TEMP_PROFILE") {
            std::string who;
            if (!(ss >> who))
                throw std::runtime_error("TEMP_PROFILE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::TempProfileSpec tp;
            if (detail::toUpper(who) == "ALL") tp.all = true;
            else tp.elems = detail::parseIntList(who);
            auto it = p.find("T1");
            if (it != p.end()) tp.T1 = std::stod(it->second);
            it = p.find("T2");
            if (it != p.end()) tp.T2 = std::stod(it->second);
            it = p.find("A"); if (it != p.end()) tp.A = std::stod(it->second);
            it = p.find("B"); if (it != p.end()) tp.B = std::stod(it->second);
            it = p.find("H"); if (it != p.end()) tp.H = std::stod(it->second);
            if (tp.A <= 0.0)
                throw std::runtime_error("*TEMP_PROFILE A 必须大于 0");
            m.tempProfiles.push_back(tp);
        } else if (section == "*STAGE") {
            int id = 0;
            if (!(ss >> id))
                throw std::runtime_error("STAGE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Stage st;
            st.id = id;
            auto it = p.find("ACTIVE");
            if (it != p.end()) st.activate = detail::parseIntList(it->second);
            it = p.find("DEACTIVE");
            if (it != p.end()) st.deactivate = detail::parseIntList(it->second);
            it = p.find("TIME");
            if (it != p.end()) st.time = std::stod(it->second);
            m.stages.push_back(st);
        } else if (section == "*STAGE_LOAD") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_LOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            auto itV = p.find("VALUE");
            if (itN == p.end() || itD == p.end() || itV == p.end())
                throw std::runtime_error("*STAGE_LOAD 需要 NODE/DOF/VALUE");
            for (int node : detail::parseIntList(itN->second))
                for (int dof : detail::parseIntList(itD->second))
                    m.stages[stageId - 1].loads.push_back({node, dof, std::stod(itV->second)});
        } else if (section == "*STAGE_ELEMENT_LOAD") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_ELEMENT_LOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itE = p.find("ELEM");
            auto itDir = p.find("DIR");
            auto itQ1 = p.find("Q1");
            auto itQ2 = p.find("Q2");
            if (itE == p.end() || itDir == p.end() || itQ1 == p.end() || itQ2 == p.end())
                throw std::runtime_error("*STAGE_ELEMENT_LOAD 需要 ELEM/DIR/Q1/Q2");
            for (int elem : detail::parseIntList(itE->second))
                m.stages[stageId - 1].elemLoads.push_back({elem, itDir->second, std::stod(itQ1->second), std::stod(itQ2->second)});
        } else if (section == "*STAGE_PRESTRESS") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_PRESTRESS 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::Prestress ps;
            auto it = p.find("ELEM");
            if (it == p.end()) throw std::runtime_error("*STAGE_PRESTRESS 缺少 ELEM");
            ps.elems = detail::parseIntList(it->second);
            it = p.find("P");
            if (it == p.end()) throw std::runtime_error("*STAGE_PRESTRESS 缺少 P");
            ps.P = std::stod(it->second);
            it = p.find("F"); if (it != p.end()) ps.F = std::stod(it->second);
            it = p.find("E1"); if (it != p.end()) ps.E1 = std::stod(it->second);
            it = p.find("E2"); if (it != p.end()) ps.E2 = std::stod(it->second);
            it = p.find("LOSS");
            if (it == p.end()) throw std::runtime_error("*PRESTRESS / *STAGE_PRESTRESS 缺少 LOSS（损失系数，必须显式指定）");
            ps.loss = std::stod(it->second);
            it = p.find("EPS"); if (it != p.end()) ps.Ep = std::stod(it->second);
            it = p.find("APS"); if (it != p.end()) ps.Ap = std::stod(it->second);
            ps.method = "EQUIV"; // 施工阶段第一版只使用等效荷载法
            m.stages[stageId - 1].prestresses.push_back(ps);
        } else if (section == "*STAGE_BOUNDARY" || section == "*STAGE_UNBOUNDARY") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_BOUNDARY/UNBOUNDARY 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            if (itN == p.end() || itD == p.end())
                throw std::runtime_error("*STAGE_BOUNDARY/UNBOUNDARY 需要 NODE/DOF");
            const auto nodes = detail::parseIntList(itN->second);
            const auto dofs = detail::parseIntList(itD->second);
            for (int node : nodes) {
                for (int dof : dofs) {
                    if (section == "*STAGE_BOUNDARY") {
                        double value = 0.0;
                        auto itV = p.find("VALUE");
                        if (itV != p.end()) value = std::stod(itV->second);
                        m.stages[stageId - 1].addedConstraints.push_back({node, dof, value});
                    } else {
                        m.stages[stageId - 1].removedConstraints.push_back({node, dof});
                    }
                }
            }
        } else if (section == "*STAGE_CREEP" || section == "*STAGE_SHRINKAGE") {
            int stageId = 0;
            std::string who;
            if (!(ss >> stageId >> who))
                throw std::runtime_error("STAGE_CREEP/SHRINKAGE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            if (section == "*STAGE_CREEP") {
                Model2D::CreepSpec cs;
                if (detail::toUpper(who) == "ALL") cs.all = true;
                else cs.elems = detail::parseIntList(who);
                auto it = p.find("PHI");
                if (it == p.end()) throw std::runtime_error("*STAGE_CREEP 缺少 PHI");
                cs.phi = std::stod(it->second);
                it = p.find("CHI"); if (it != p.end()) cs.chi = std::stod(it->second);
                m.stages[stageId - 1].creepSpecs.push_back(cs);
            } else {
                Model2D::ShrinkageSpec sh;
                if (detail::toUpper(who) == "ALL") sh.all = true;
                else sh.elems = detail::parseIntList(who);
                auto it = p.find("ECS");
                if (it == p.end()) throw std::runtime_error("*STAGE_SHRINKAGE 缺少 ECS");
                sh.ecs = std::stod(it->second);
                m.stages[stageId - 1].shrinkageSpecs.push_back(sh);
            }
        } else if (section == "*STAGE_TEMP_PROFILE") {
            int stageId = 0;
            std::string who;
            if (!(ss >> stageId >> who))
                throw std::runtime_error("STAGE_TEMP_PROFILE 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            Model2D::TempProfileSpec tp;
            if (detail::toUpper(who) == "ALL") tp.all = true;
            else tp.elems = detail::parseIntList(who);
            auto it = p.find("T1"); if (it != p.end()) tp.T1 = std::stod(it->second);
            it = p.find("T2"); if (it != p.end()) tp.T2 = std::stod(it->second);
            it = p.find("A"); if (it != p.end()) tp.A = std::stod(it->second);
            it = p.find("B"); if (it != p.end()) tp.B = std::stod(it->second);
            it = p.find("H"); if (it != p.end()) tp.H = std::stod(it->second);
            if (tp.A <= 0.0)
                throw std::runtime_error("*STAGE_TEMP_PROFILE A 必须大于 0");
            m.stages[stageId - 1].tempProfiles.push_back(tp);
        } else if (section == "*CREEP_CURVE" || section == "*SHRINKAGE_CURVE") {
            double time = 0.0, value = 0.0;
            if (!(ss >> time >> value))
                throw std::runtime_error("CREEP_CURVE/SHRINKAGE_CURVE 行格式错误: " + t);
            if (section == "*CREEP_CURVE") m.creepCurve.points.push_back({time, value});
            else m.shrinkageCurve.points.push_back({time, value});
        } else if (section == "*CREEP_MODEL" || section == "*SHRINKAGE_MODEL") {
            auto p = detail::parseKeyValues(t);
            if (section == "*CREEP_MODEL") {
                auto it = p.find("MODEL"); if (it != p.end()) m.creepModel = detail::toUpper(it->second);
                it = p.find("PHI0"); if (it != p.end()) m.creepPhi0 = std::stod(it->second);
                it = p.find("T0"); if (it != p.end()) m.creepT0 = std::stod(it->second);
                it = p.find("BETA_H"); if (it != p.end()) m.creepBetaH = std::stod(it->second);
                it = p.find("POWER"); if (it != p.end()) m.creepPower = std::stod(it->second);
            } else {
                auto it = p.find("MODEL"); if (it != p.end()) m.shrinkageModel = detail::toUpper(it->second);
                it = p.find("ECS_INF"); if (it != p.end()) m.shrinkageEcsInf = std::stod(it->second);
                it = p.find("TS"); if (it != p.end()) m.shrinkageTs = std::stod(it->second);
                it = p.find("ALPHA_DS"); if (it != p.end()) m.shrinkageAlphaDs = std::stod(it->second);
                it = p.find("POWER"); if (it != p.end()) m.shrinkagePower = std::stod(it->second);
            }
        } else if (section == "*STAGE_UNLOAD") {
            int stageId = 0;
            if (!(ss >> stageId))
                throw std::runtime_error("STAGE_UNLOAD 行格式错误: " + t);
            std::string rest2;
            std::getline(ss, rest2);
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            auto itE = p.find("ELEM");
            if (itN != p.end() && itD != p.end()) {
                for (int node : detail::parseIntList(itN->second))
                    for (int dof : detail::parseIntList(itD->second))
                        m.stages[stageId - 1].unloadNodeDofs.push_back({node, dof});
            } else if (itE != p.end()) {
                const auto elems = detail::parseIntList(itE->second);
                m.stages[stageId - 1].unloadElems.insert(m.stages[stageId - 1].unloadElems.end(),
                                                        elems.begin(), elems.end());
            } else {
                throw std::runtime_error("*STAGE_UNLOAD 需要 NODE/DOF 或 ELEM");
            }
        } else if (section == "*STAGE_TARGET") {
            std::string rest2 = t;
            auto p = detail::parseKeyValues(rest2);
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            auto itV = p.find("VALUE");
            Model2D::StageTarget tg;
            if (itV == p.end())
                throw std::runtime_error("*STAGE_TARGET 需要 VALUE");
            tg.value = std::stod(itV->second);
            auto itW = p.find("WEIGHT");
            if (itW != p.end()) tg.weight = std::stod(itW->second);
            auto itK = p.find("TYPE");
            tg.kind = itK != p.end() ? detail::toUpper(itK->second) : "DISP";
            if (tg.kind == "DISP" || tg.kind == "REACTION") {
                if (itN == p.end() || itD == p.end())
                    throw std::runtime_error("*STAGE_TARGET DISP/REACTION 需要 NODE/DOF");
                tg.node = std::stoi(itN->second);
                tg.dof = std::stoi(itD->second);
            } else if (tg.kind == "FORCE" || tg.kind == "STRESS") {
                auto itE = p.find("ELEM");
                auto itC = p.find("COMPONENT");
                auto itEnd = p.find("END");
                if (itE == p.end() || itC == p.end())
                    throw std::runtime_error("*STAGE_TARGET FORCE/STRESS 需要 ELEM/COMPONENT");
                tg.elem = std::stoi(itE->second);
                tg.component = detail::toUpper(itC->second);
                if (itEnd != p.end()) tg.end = std::stoi(itEnd->second);
            } else {
                throw std::runtime_error("*STAGE_TARGET 未知 TYPE: " + tg.kind);
            }
            m.stageTargets.push_back(tg);
        } else if (section == "*STAGE_CONTROL") {
            auto p = detail::parseKeyValues(t);
            auto itS = p.find("STAGE");
            auto itP = p.find("PRESTRESS");
            if (itS == p.end() || itP == p.end())
                throw std::runtime_error("*STAGE_CONTROL 需要 STAGE/PRESTRESS");
            Model2D::StageControl ctrl;
            ctrl.stage = std::stoi(itS->second);
            ctrl.prestressIndex = std::stoi(itP->second);
            auto itL = p.find("LOWER");
            if (itL != p.end()) ctrl.lower = std::stod(itL->second);
            auto itU = p.find("UPPER");
            if (itU != p.end()) ctrl.upper = std::stod(itU->second);
            if (ctrl.lower > ctrl.upper)
                throw std::runtime_error("*STAGE_CONTROL LOWER 不能大于 UPPER");
            m.stageControls.push_back(ctrl);
        } else if (section == "*STAGE_CONSTRAINT") {
            auto p = detail::parseKeyValues(t);
            Model2D::StageConstraint sc;
            auto itV = p.find("VALUE");
            if (itV != p.end()) sc.target.value = std::stod(itV->second);
            auto itW = p.find("WEIGHT");
            if (itW != p.end()) sc.target.weight = std::stod(itW->second);
            auto itK = p.find("TYPE");
            sc.target.kind = itK != p.end() ? detail::toUpper(itK->second) : "DISP";
            auto itN = p.find("NODE");
            auto itD = p.find("DOF");
            if (sc.target.kind == "DISP" || sc.target.kind == "REACTION") {
                if (itN == p.end() || itD == p.end())
                    throw std::runtime_error("*STAGE_CONSTRAINT DISP/REACTION 需要 NODE/DOF");
                sc.target.node = std::stoi(itN->second);
                sc.target.dof = std::stoi(itD->second);
            } else if (sc.target.kind == "FORCE" || sc.target.kind == "STRESS") {
                auto itE = p.find("ELEM");
                auto itC = p.find("COMPONENT");
                auto itEnd = p.find("END");
                if (itE == p.end() || itC == p.end())
                    throw std::runtime_error("*STAGE_CONSTRAINT FORCE/STRESS 需要 ELEM/COMPONENT");
                sc.target.elem = std::stoi(itE->second);
                sc.target.component = detail::toUpper(itC->second);
                if (itEnd != p.end()) sc.target.end = std::stoi(itEnd->second);
            } else {
                throw std::runtime_error("*STAGE_CONSTRAINT 未知 TYPE");
            }
            auto itL = p.find("LOWER");
            if (itL != p.end()) sc.lower = std::stod(itL->second);
            auto itU = p.find("UPPER");
            if (itU != p.end()) sc.upper = std::stod(itU->second);
            m.stageConstraints.push_back(sc);
        } else if (section == "*CAMBER") {
            int node = 0, dof = 0;
            double value = 0.0;
            if (!(ss >> node >> dof))
                throw std::runtime_error("CAMBER 行格式错误: " + t);
            ss >> value;
            if (dof < 1 || dof > 3)
                throw std::runtime_error("CAMBER dof 必须在 1~3 之间");
            m.cambers.push_back({node, dof, value});
        } else if (section == "*FEEDBACK") {
            int node = 0, dof = 0;
            double measured = 0.0;
            if (!(ss >> node >> dof >> measured))
                throw std::runtime_error("FEEDBACK 行格式错误: " + t);
            m.feedbacks.push_back({node, dof, measured});
        } else if (section == "*END") {
            break;
        }
    }

    if (m.nodes.empty())
        throw std::runtime_error("输入文件没有定义任何节点");
    processPrestress(m);
    // 先记录未折减的基准模量，再做徐变折减。
    // 之后每一次徐变处理都从这个基准重算，保证幂等。
    snapshotBaseModulus(m);
    processTimeEffects(m);

    // 危险默认值检查：若存在重力/自重荷载但单元密度为 0，给出警告。
    const bool hasGravity = (m.gx != 0.0 || m.gy != 0.0 || !m.selfWeights.empty());
    if (hasGravity) {
        for (const auto& e : m.elements) {
            if (e.dens <= 0.0 && e.addmass <= 0.0 && !isRigid(e.type)) {
                std::cerr << "[警告] 单元 " << e.id << " 的密度 dens=0 且 addmass=0，"
                          << "但模型存在 *GRAVITY 或 *SELF_WEIGHT，该单元将不产生自重。"
                          << "若需自重，请在 *MAT 或 *ELEMENT 中指定 DENS。\n";
            }
        }
    }

    applyJtgD60LaneLoads2D(m);
    return m;
}

void appendConstraintPenaltyTriplets(const Model2D& m,
                                            const std::map<int, int>& index,
                                            int nDofs,
                                            double penalty,
                                            std::vector<Eigen::Triplet<double>>& triplets) {
    for (const auto& ce : m.ces) {
        std::map<int, double> c;
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            c[d] += term.coef;
        }
        for (const auto& [i, ci] : c) {
            for (const auto& [j, cj] : c) {
                const double v = penalty * ci * cj;
                if (v != 0.0)
                    triplets.emplace_back(i, j, v);
            }
        }
    }
}

double effectiveCePenalty(const std::vector<Eigen::Triplet<double>>& triplets,
                                 double userPenalty) {
    if (userPenalty > 0.0) return userPenalty;
    double maxDiag = 0.0;
    for (const auto& t : triplets) {
        if (t.row() == t.col())
            maxDiag = std::max(maxDiag, std::abs(t.value()));
    }
    if (maxDiag <= 0.0) maxDiag = 1.0;
    return 1e8 * maxDiag;
}

DenseVec constraintPenaltyConstantForce(const Model2D& m,
                                               const std::map<int, int>& index,
                                               int nDofs,
                                               double penalty) {
    DenseVec Fce = DenseVec::Zero(nDofs);
    for (const auto& ce : m.ces) {
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            Fce(d) += penalty * ce.constant * term.coef;
        }
    }
    return Fce;
}

DenseVec constraintPenaltyResidualForce(const Model2D& m,
                                               const std::map<int, int>& index,
                                               int nDofs,
                                               double penalty,
                                               const DenseVec& U) {
    DenseVec Fce = DenseVec::Zero(nDofs);
    for (const auto& ce : m.ces) {
        double g = -ce.constant;
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            g += term.coef * U(d);
        }
        const double a = penalty * g;
        for (const auto& term : ce.terms) {
            const int d = 3 * index.at(term.node) + (term.dof - 1);
            Fce(d) += a * term.coef;
        }
    }
    return Fce;
}

void autoConstrainFreeRotations(const SparseMat& K, int nDofs,
                                       const std::vector<int>& baseCdofs,
                                       const std::vector<double>& baseCvals,
                                       std::vector<int>& cdofs,
                                       std::vector<double>& cvals) {
    cdofs = baseCdofs;
    cvals = baseCvals;
    std::vector<char> isC(nDofs, 0);
    for (int d : cdofs) isC[d] = 1;
    for (int d = 2; d < nDofs; d += 3) {   // 每节点的第 3 个自由度 = 转角 θ
        if (isC[d]) continue;
        if (K.coeff(d, d) == 0.0) {
            cdofs.push_back(d);
            cvals.push_back(0.0);
            isC[d] = 1;
        }
    }
}

ModelResult solveModel2DNonlinear(const Model2D& m,
                                         SolverType type,
                                         const DenseVec* U0,
                                         const std::vector<char>* active0,
                                         NonlinearDiag* diagOut,
                                         ReusableSolver* reuse) {
    for (const auto& e : m.elements)
        if (!e.releases.empty())
            throw std::runtime_error("端部释放暂不支持 NONLINEAR 分析（单元 " +
                                     std::to_string(e.id) + "）");
    std::vector<int> ids;
    std::map<int, int> index;
    for (const auto& [id, c] : m.nodes) {
        index[id] = static_cast<int>(ids.size());
        ids.push_back(id);
    }
    const int nNodes = static_cast<int>(ids.size());
    const int nDofs = 3 * nNodes;

    // 约束自由度与指定位移。
    std::vector<int> cdofs;
    std::vector<double> cvals;
    std::vector<char> isConstrained(nDofs, 0);
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        if (c.k > 0.0) continue;   // 弹性支承不作硬约束
        const int d = 3 * index.at(c.node) + (c.dof - 1);
        cdofs.push_back(d);
        cvals.push_back(c.value);
        isConstrained[d] = 1;
    }

    DenseVec U = DenseVec::Zero(nDofs);
    if (U0) {
        if (U0->size() != nDofs)
            throw std::runtime_error("施工阶段初始位移向量长度不匹配");
        U = *U0;
    }
    for (std::size_t i = 0; i < cdofs.size(); ++i)
        U(cdofs[i]) = cvals[i];

    // cable/strut 初始全部激活，后续根据当前构形切换。
    std::vector<char> active(m.elements.size(), 1);
    if (active0 && active0->size() == active.size())
        active = *active0;

    const double tol = m.nonlinearTol;
    const int maxIter = m.nonlinearMaxIter;

    auto axialActiveTol = [](const Model2D::Element& e) {
        return 1e-12 * std::max(1.0, std::abs(e.E * e.A));
    };

    // 荷载分级比例：assembleState 只把**外荷载**乘以它。
    // 1.0 = 全部荷载（默认，与旧行为一致）。
    double loadScale = 1.0;

    // 把指定位移"钉"回向量里。荷载分级/线搜索期间每一步都要保证约束成立。
    auto applyPrescribed = [&](DenseVec& V) {
        for (std::size_t i = 0; i < cdofs.size(); ++i)
            V(cdofs[i]) = cvals[i];
    };

    auto assembleState = [&](const std::vector<char>& act, const DenseVec& Ucur,
                             DenseVec& Fint, DenseVec& Fext,
                             std::vector<Eigen::Triplet<double>>* triplets) {
        Fint = DenseVec::Zero(nDofs);
        Fext = DenseVec::Zero(nDofs);
        if (triplets) triplets->clear();

        // 节点集中荷载始终施加（按当前荷载分级比例缩放）。
        for (const auto& l : m.loads) {
            if (l.dof < 1 || l.dof > 3)
                throw std::runtime_error("载荷自由度必须在 1~3 之间");
            Fext(3 * index.at(l.node) + (l.dof - 1)) += loadScale * l.value;
        }

        for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
            if (!act[ei]) continue;
            const auto& e = m.elements[ei];
            if (e.conn.size() != 2)
                throw std::runtime_error("二维单元必须有两个节点");
            const int i1 = index.at(e.conn[0]);
            const int i2 = index.at(e.conn[1]);
            const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
            const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

            if (isAxial(e.type)) {
                const std::vector<int> g = elementGdofs2D(e.type, i1, i2);
                DenseVec fext = elementExternalEquivalentLoad(m, e, p1, p2);
                fext *= loadScale;
                const DenseVec fint = axialInternalForceNonlinear(e, Ucur, i1, i2, p1, p2);
                addElementLoad(g, fext, Fext);
                addElementLoad(g, fint, Fint);
                if (triplets) {
                    const DenseMat Ke = axialTangentNonlinear(e, Ucur, i1, i2, p1, p2);
                    scatter(Ke, g, *triplets);
                }
            } else if (e.type == "RIGID_PIN") {
                const std::vector<int> g = elementGdofs2D(e.type, i1, i2);
                Eigen::VectorXd ug(6);
                ug << Ucur(3 * i1 + 0), Ucur(3 * i1 + 1), Ucur(3 * i1 + 2),
                      Ucur(3 * i2 + 0), Ucur(3 * i2 + 1), Ucur(3 * i2 + 2);
                const DenseMat Ke = rigidPinStiffness(p1, p2);
                addElementLoad(g, Ke * ug, Fint);
                addElementLoad(g, DenseVec::Zero(6), Fext);
                if (triplets) scatter(Ke, g, *triplets);
            } else if (e.type == "BEAM" || e.type == "RIGID") {
                EulerBeam2D be{e.E, e.A, e.I, p1, p2};
                const std::vector<int> g = elementGdofs2D(e.type, i1, i2);
                Eigen::VectorXd ug(6);
                ug << Ucur(3 * i1 + 0), Ucur(3 * i1 + 1), Ucur(3 * i1 + 2),
                      Ucur(3 * i2 + 0), Ucur(3 * i2 + 1), Ucur(3 * i2 + 2);
                const double Phi = beamShearParameter(e, (p2 - p1).norm());
                const DenseMat Ke = be.stiffness(Phi);
                DenseVec fext = elementExternalEquivalentLoad(m, e, p1, p2);
                fext *= loadScale;
                DenseMat Kt = Ke;
                DenseVec fint = Ke * ug;
                if (e.type == "BEAM") {
                    // 一阶 P-Δ / 几何刚度：轴力 N 会软化（受压）或增强（受拉）弯曲刚度。
                    // N 取当前构形下的物理轴力 = 机械轴力 + 局部等效轴力（含温度/自重轴向分量）。
                    const DenseMat T = be.transformation();
                    const Eigen::VectorXd ul = T * ug;
                    const double L = (p2 - p1).norm();
                    const Eigen::VectorXd fel = T * fext;   // 已按 loadScale 缩放
                    const double N = e.E * e.A / L * (ul(3) - ul(0)) + fel(0);
                    const DenseMat Kgg = T.transpose() * beamGeometricStiffness2D(N, L) * T;
                    Kt += Kgg;
                    fint = Kt * ug;
                }
                addElementLoad(g, fext, Fext);
                addElementLoad(g, fint, Fint);
                if (triplets) scatter(Kt, g, *triplets);
            } else {
                throw std::runtime_error("未知单元类型: " + e.type);
            }
        }

        // 弹性支承（k>0）：Fint += k·u，Fext += k·value，切线 += k。
        for (const auto& c : m.constraints) {
            if (c.k <= 0.0) continue;
            const int d = 3 * index.at(c.node) + (c.dof - 1);
            Fint(d) += c.k * Ucur(d);
            Fext(d) += c.k * c.value;
            if (triplets) triplets->emplace_back(d, d, c.k);
        }
    };

    DenseVec Fint, Fext;
    NonlinearDiag diag;
    std::vector<double> zeroPrescribed(cdofs.size(), 0.0);
    bool haveInitialResidual = false;

    // 给定构形与激活状态，返回自由自由度的残差范数（只组装内力/外力，不算切线）。
    // 线搜索每一步都用它判断"这一步有没有让残差下降"。
    auto residualNorm = [&](const std::vector<char>& act, const DenseVec& Ucur) -> double {
        DenseVec fi, fe;
        assembleState(act, Ucur, fi, fe, nullptr);
        DenseVec rr = fe - fi;
        if (!m.ces.empty())
            rr -= constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, Ucur);
        double s = 0.0;
        for (int d = 0; d < nDofs; ++d)
            if (!isConstrained[d]) s += rr(d) * rr(d);
        return std::sqrt(s);
    };

    // 依据当前构形更新 cable/strut 的单向激活状态。
    auto updateActive = [&](const DenseVec& Ucur) -> std::vector<char> {
        std::vector<char> na(m.elements.size(), 1);
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
            const auto& e = m.elements[ei];
            if (e.type != "CABLE" && e.type != "STRUT") continue;
            const int i1 = index.at(e.conn[0]);
            const int i2 = index.at(e.conn[1]);
            const double N = axialForceNonlinear(e, Ucur, i1, i2,
                                                 m.nodes.at(e.conn[0]),
                                                 m.nodes.at(e.conn[1]));
            const double atol = axialActiveTol(e);
            na[ei] = (e.type == "CABLE") ? (N > atol) : (N < -atol);
        }
        return na;
    };

    // 线性路径已经会为"只连桁架的节点"自动约束多余转角自由度（P1-3 的二维版）。
    // 非线性路径原先没有这一步，纯桁架/索模型必须先做一次装配才知道哪些 θ 没有刚度，
    // 因此这里在初始构形上预装配一次，补上这些约束后再进入迭代。
    {
        std::vector<Eigen::Triplet<double>> trip0;
        DenseVec fi0, fe0;
        assembleState(active, U, fi0, fe0, &trip0);
        if (!m.ces.empty()) {
            const double penalty = effectiveCePenalty(trip0, m.cePenalty);
            appendConstraintPenaltyTriplets(m, index, nDofs, penalty, trip0);
        }
        SparseMat K0(nDofs, nDofs);
        K0.setFromTriplets(trip0.begin(), trip0.end());
        K0.makeCompressed();

        std::vector<int> cd2;
        std::vector<double> cv2;
        autoConstrainFreeRotations(K0, nDofs, cdofs, cvals, cd2, cv2);
        if (cd2.size() != cdofs.size()) {
            cdofs.swap(cd2);
            cvals.swap(cv2);
            isConstrained.assign(nDofs, 0);
            for (int d : cdofs) isConstrained[d] = 1;
            zeroPrescribed.assign(cdofs.size(), 0.0);
            applyPrescribed(U);
        }
    }

    // ---- 外层：荷载分级；内层：Newton-Raphson ----
    const int nSteps = std::max(1, m.nonlinearLoadSteps);
    double lam = 0.0;                           // 已成功施加的荷载比例
    double dlam = 1.0 / static_cast<double>(nSteps);
    const double minLam = 1.0 / 4096.0;         // 再细分下去没有意义
    const int maxLineSearch = 12;               // 步长最多折半 12 次

    int finalIter = 0;
    double finalResidual = 0.0;
    bool converged = false;

    while (lam < 1.0 - 1e-12) {
        const double lamTarget = std::min(1.0, lam + dlam);
        loadScale = lamTarget;

        // 本步起点：不收敛时要回滚到它，再用更小的步长重来。
        const DenseVec Usaved = U;
        const std::vector<char> actSaved = active;

        bool stepOk = false;
        int stepIter = 0;

        for (int iter = 0; iter <= maxIter; ++iter) {
            std::vector<Eigen::Triplet<double>> triplets;
            assembleState(active, U, Fint, Fext, &triplets);
            if (!m.ces.empty()) {
                const double penalty = effectiveCePenalty(triplets, m.cePenalty);
                appendConstraintPenaltyTriplets(m, index, nDofs, penalty, triplets);
            }

            SparseMat Kt(nDofs, nDofs);
            Kt.setFromTriplets(triplets.begin(), triplets.end());
            Kt.makeCompressed();

            // 自由自由度残差范数，用于收敛判断。
            DenseVec r = Fext - Fint;
            if (!m.ces.empty())
                r -= constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, U);
            double rFree2 = 0.0, fFree2 = 0.0, rMaxAbs = 0.0;
            int rMaxDof = -1;
            double kDiagMax = 0.0;
            double kDiagMin = std::numeric_limits<double>::infinity();
            for (int d = 0; d < nDofs; ++d) {
                if (isConstrained[d]) continue;
                rFree2 += r(d) * r(d);
                fFree2 += Fext(d) * Fext(d);
                if (std::abs(r(d)) > rMaxAbs) { rMaxAbs = std::abs(r(d)); rMaxDof = d; }
                const double kd = std::abs(Kt.coeff(d, d));
                if (kd > kDiagMax) kDiagMax = kd;
                if (kd < kDiagMin) kDiagMin = kd;
            }
            const double rFree = std::sqrt(rFree2);
            const double fFree = std::sqrt(fFree2);
            const double resRef = tol * std::max(1.0, fFree);
            stepIter = iter + 1;

            diag.residualHistory.push_back(rFree);
            diag.lastResidual = rFree;
            diag.residualRef = resRef;
            if (!haveInitialResidual) { diag.initialResidual = rFree; haveInitialResidual = true; }
            if (kDiagMin > 0.0) diag.diagCondEst = kDiagMax / kDiagMin;
            if (rMaxDof >= 0) {
                diag.maxResidualNode = ids[rMaxDof / 3];
                diag.maxResidualDof = (rMaxDof % 3) + 1;
                diag.maxResidualValue = rMaxAbs;
            }

            // 自由 DOF 若没有刚度，直接给出可诊断错误。
            for (int d = 0; d < nDofs; ++d) {
                if (!isConstrained[d] && Kt.coeff(d, d) == 0.0) {
                    const int ni = d / 3;
                    const int ld = d % 3;
                    throw std::runtime_error(
                        "方案C：存在零刚度自由度：节点 " + std::to_string(ids[ni]) +
                        " dof " + std::to_string(ld + 1) +
                        "（未约束且未被任何激活单元连接）");
                }
            }

            // 残差判据：当前构形已平衡，无需再解增量。
            if (rFree <= resRef) { stepOk = true; break; }

            // 求解增量：Kt·ΔU = Fext - Fint。
            // Newton 迭代中切线矩阵的**非零结构恒定**，只有数值在变，
            // 因此符号分解只需做一次（reuse 由调用方在施工阶段/扫描间传递）。
            DenseVec rhs = r;
            applyDirichlet(Kt, rhs, cdofs, zeroPrescribed);
            DenseVec dU;
            if (reuse && reuse->reusable()) {
                reuse->factorizeReuse(Kt);
                dU = reuse->solve(rhs);
            } else {
                dU = solve(Kt, rhs, type);
            }
            if (!dU.allFinite()) {
                diag.failureReason = "切线方程解出非有限增量 ΔU（矩阵奇异或严重病态）";
                break;
            }

            // ---- 线搜索（backtracking）：α ∈ {1, 1/2, 1/4, ...} ----
            // 索松弛、跳跃失稳时全步长 ΔU 会让残差不降反升甚至发散。
            // 正常情况下 α=1 立刻被接受，结果与不启用线搜索时一致。
            double alpha = 1.0;
            DenseVec Ubest = U + dU;
            applyPrescribed(Ubest);
            double rBest = residualNorm(active, Ubest);
            if (m.nonlinearLineSearch && !(rBest < rFree)) {
                for (int ls = 0; ls < maxLineSearch; ++ls) {
                    alpha *= 0.5;
                    ++diag.lineSearchCuts;
                    Ubest = U + alpha * dU;
                    applyPrescribed(Ubest);
                    rBest = residualNorm(active, Ubest);
                    if (rBest < rFree) break;
                }
            }

            U = Ubest;

            // 依据新构形更新 cable/strut 激活状态。
            std::vector<char> newActive = updateActive(U);
            const bool statusChanged = (newActive != active);
            active.swap(newActive);

            const double duOk = (alpha * dU.norm()) <= tol * std::max(1.0, U.norm());
            const bool resOk = rBest <= resRef;

            // 发散检测：残差涨到初始值的 1e8 倍以上，再迭代没有意义。
            if (!U.allFinite() ||
                (diag.initialResidual > 0.0 && rBest > 1e8 * diag.initialResidual)) {
                diag.failureReason = "残差发散（超过初始残差的 1e8 倍）";
                break;
            }

            if (!statusChanged && duOk && resOk) { stepOk = true; break; }
        }

        finalIter += stepIter;

        if (stepOk) {
            lam = lamTarget;
            ++diag.loadSteps;
            continue;
        }

        // 本荷载步失败：回滚，步长减半重来（自适应荷载分级）。
        U = Usaved;
        active = actSaved;
        if (diag.failureReason.empty())
            diag.failureReason = "单个荷载步内达到最大迭代次数 " + std::to_string(maxIter);
        dlam *= 0.5;
        ++diag.stepRefinements;
        if (dlam < minLam) break;
    }

    converged = (lam >= 1.0 - 1e-12);
    diag.iterations = finalIter;
    diag.loadRatio = lam;
    diag.converged = converged;

    // 用最终位移和最终激活状态重新组装一次，得到精确内力/反力/残差。
    loadScale = 1.0;
    assembleState(active, U, Fint, Fext, nullptr);
    DenseVec rFinal = Fext - Fint;
    if (!m.ces.empty())
        rFinal -= constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, U);
    double rf2 = 0.0;
    for (int d = 0; d < nDofs; ++d)
        if (!isConstrained[d]) rf2 += rFinal(d) * rFinal(d);
    finalResidual = std::sqrt(rf2);
    diag.lastResidual = finalResidual;

    DenseVec R = Fint - Fext; // 约束反力 = 内部节点力 - 外部节点力
    if (!m.ces.empty())
        R += constraintPenaltyResidualForce(m, index, nDofs, m.cePenalty, U);

    if (!converged && m.nonlinearThrowOnFailure)
        throw std::runtime_error(diag.describe());

    ModelResult result;
    result.U = U;
    result.R = R;
    result.nodeIds = ids;
    result.coords.reserve(nNodes);
    for (int id : ids) result.coords.push_back(m.nodes.at(id));
    result.nonlinearIterations = finalIter;
    result.nonlinearResidual = finalResidual;
    result.nonlinearConverged = converged;
    if (diagOut) *diagOut = diag;
    return result;
}

void assembleLinear2D(const Model2D& m,
                             const std::vector<char>& active,
                             const std::map<int, int>& index,
                             int nDofs,
                             SparseMat& K, DenseVec& F) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
        if (!active[ei]) continue;
        const auto& e = m.elements[ei];
        if (e.conn.size() != 2)
            throw std::runtime_error("二维单元必须有两个节点");
        const int i1 = index.at(e.conn[0]);
        const int i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

        if (isAxial(e.type)) {
            TrussElement te{e.E, e.A, 2,
                            Eigen::Vector3d(p1(0), p1(1), 0.0),
                            Eigen::Vector3d(p2(0), p2(1), 0.0)};
            const DenseMat Ke = te.stiffness(); // 4x4
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i2, 3 * i2 + 1};
            scatter(Ke, g, triplets);
        } else if (e.type == "RIGID_PIN") {
            const DenseMat Ke = rigidPinStiffness(p1, p2);
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i1 + 2,
                                        3 * i2, 3 * i2 + 1, 3 * i2 + 2};
            scatter(Ke, g, triplets);
        } else if (e.type == "BEAM" || e.type == "RIGID") {
            if (e.type == "BEAM" && !e.releases.empty()) {
                const CondensedBeamGlobal cb = condenseBeamGlobal(e, p1, p2, i1, i2,
                                                                  DenseVec::Zero(6));
                scatter(cb.Ke, cb.gdofs, triplets);
            } else {
                EulerBeam2D be{e.E, e.A, e.I, p1, p2};
                const double Phi = beamShearParameter(e, (p2 - p1).norm());
                const DenseMat Ke = be.stiffness(Phi); // 6x6
                const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i1 + 2,
                                            3 * i2, 3 * i2 + 1, 3 * i2 + 2};
                scatter(Ke, g, triplets);
            }
        } else {
            throw std::runtime_error("未知单元类型: " + e.type);
        }
    }
    if (!m.ces.empty()) {
        const double penalty = effectiveCePenalty(triplets, m.cePenalty);
        appendConstraintPenaltyTriplets(m, index, nDofs, penalty, triplets);
    }
    K.resize(nDofs, nDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();

    F = DenseVec::Zero(nDofs);
    for (const auto& l : m.loads) {
        if (l.dof < 1 || l.dof > 3)
            throw std::runtime_error("载荷自由度必须在 1~3 之间");
        F(3 * index.at(l.node) + (l.dof - 1)) += l.value;
    }
    for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
        if (!active[ei]) continue;
        const auto& e = m.elements[ei];
        const int i1 = index.at(e.conn[0]), i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);
        const DenseVec fe = elementEquivalentLoad(m, e, p1, p2);
        if (e.type == "BEAM" && !e.releases.empty()) {
            const CondensedBeamGlobal cb = condenseBeamGlobal(e, p1, p2, i1, i2, fe);
            addElementLoad(cb.gdofs, cb.feq, F);
        } else {
            addElementLoad(elementGdofs2D(e.type, i1, i2), fe, F);
        }
    }
    if (!m.ces.empty())
        F += constraintPenaltyConstantForce(m, index, nDofs, m.cePenalty);

    // 弹性支承（k>0）：刚度对角加 k，等效荷载加 k·value；不作为硬约束。
    for (const auto& c : m.constraints) {
        if (c.k <= 0.0) continue;
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        const int d = 3 * index.at(c.node) + (c.dof - 1);
        K.coeffRef(d, d) += c.k;
        F(d) += c.k * c.value;
    }
}

Assembled2D assembleModel2DLinear(const Model2D& m) {
    Assembled2D a;
    for (const auto& [id, c] : m.nodes) {
        a.index[id] = static_cast<int>(a.ids.size());
        a.ids.push_back(id);
    }
    a.nDofs = 3 * static_cast<int>(a.ids.size());

    const std::vector<char> active(m.elements.size(), 1);
    SparseMat K;
    DenseVec F;
    assembleLinear2D(m, active, a.index, a.nDofs, K, F);
    a.K0 = K;
    a.F0 = F;

    a.isConstrained.assign(a.nDofs, 0);
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        if (c.k > 0.0) continue;   // 弹性支承已在装配中加弹簧，不作硬约束
        const int d = 3 * a.index.at(c.node) + (c.dof - 1);
        a.cdofs.push_back(d);
        a.cvals.push_back(c.value);
        a.isConstrained[d] = 1;
    }
    std::vector<int> cd;
    std::vector<double> cv;
    autoConstrainFreeRotations(K, a.nDofs, a.cdofs, a.cvals, cd, cv);
    a.cdofs = cd;
    a.cvals = cv;
    a.isConstrained.assign(a.nDofs, 0);
    for (int d : a.cdofs) a.isConstrained[d] = 1;

    a.K = K;
    a.F = F;
    applyDirichlet(a.K, a.F, a.cdofs, a.cvals);
    return a;
}

SparseMat assembleMassMatrix2D(const Model2D& m) {
    std::map<int, int> index;
    for (const auto& [id, c] : m.nodes) index[id] = static_cast<int>(index.size());
    const int nDofs = 3 * static_cast<int>(index.size());

    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        if (e.conn.size() != 2) continue;
        const double rhoA = e.dens * e.A + e.addmass;
        if (rhoA <= 0.0) continue;
        const int i1 = index.at(e.conn[0]), i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

        if (isAxial(e.type)) {
            const double L = (p2 - p1).norm();
            const DenseMat Me = truss2DConsistentMassMatrix(rhoA, L);
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i2, 3 * i2 + 1};
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    if (Me(a, b) != 0.0) trips.emplace_back(g[a], g[b], Me(a, b));
        } else if (e.type == "BEAM") {
            EulerBeam2D be{e.E, e.A, e.I, p1, p2};
            const DenseMat Me = be.consistentMassMatrixGlobal(rhoA);
            const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i1 + 2,
                                        3 * i2, 3 * i2 + 1, 3 * i2 + 2};
            for (int a = 0; a < 6; ++a)
                for (int b = 0; b < 6; ++b)
                    if (Me(a, b) != 0.0) trips.emplace_back(g[a], g[b], Me(a, b));
        }
        // RIGID / RIGID_PIN 为建模约束，不贡献质量。
    }
    SparseMat M(nDofs, nDofs);
    M.setFromTriplets(trips.begin(), trips.end());
    M.makeCompressed();
    return M;
}

ModalResult2D computeNaturalFrequencies2D(const Model2D& m, int nModes) {
    std::map<int, int> index;
    std::vector<int> ids;
    for (const auto& [id, c] : m.nodes) { index[id] = static_cast<int>(ids.size()); ids.push_back(id); }
    const int nDofs = 3 * static_cast<int>(ids.size());

    // K（未施加 Dirichlet）与 M
    std::vector<char> active(m.elements.size(), 1);
    SparseMat K;
    DenseVec F;
    assembleLinear2D(m, active, index, nDofs, K, F);
    const SparseMat M = assembleMassMatrix2D(m);
    if (M.nonZeros() == 0)
        throw std::runtime_error("模态分析：全部单元质量为零，请在 *MP 给 DENS（或 *R 给 ADDMASS）");

    // 约束自由度 + 自动约束零刚度转角
    std::vector<int> cdofs;
    std::vector<double> cvals;
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3) continue;
        if (c.k > 0.0) continue;   // 弹性支承不作硬约束（弹簧已在 K 中）
        cdofs.push_back(3 * index.at(c.node) + (c.dof - 1));
        cvals.push_back(c.value);
    }
    std::vector<int> cd2;
    std::vector<double> cv2;
    autoConstrainFreeRotations(K, nDofs, cdofs, cvals, cd2, cv2);
    std::vector<char> isC(nDofs, 0);
    for (int d : cd2) isC[d] = 1;

    std::vector<int> free;
    for (int d = 0; d < nDofs; ++d) if (!isC[d]) free.push_back(d);
    const int nf = static_cast<int>(free.size());

    ModalResult2D out;
    if (nf == 0) return out;
    const int n = std::min(nModes, nf);

#ifdef FEM_HAS_SPECTRA
    const SparseModalResult sm = solveSparseModal(K, M, free, nDofs, n);
    if (!sm.eigenvalues.empty()) {
        for (int i = 0; i < static_cast<int>(sm.eigenvalues.size()); ++i) {
            const double lam = sm.eigenvalues[i];
            const double f = (lam > 0.0) ? std::sqrt(lam) / (2.0 * kPi) : 0.0;
            out.frequencies.push_back(f);
            DenseVec mode = DenseVec::Zero(nDofs);
            for (int a = 0; a < nf; ++a) mode(free[a]) = sm.modes[i][free[a]];
            out.modes.push_back(mode);
        }
        return out;
    }
    // 回退稠密
#endif
    Eigen::MatrixXd Kd(nf, nf), Md(nf, nf);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            Kd(a, b) = K.coeff(free[a], free[b]);
            Md(a, b) = M.coeff(free[a], free[b]);
        }
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Kd, Md);
    for (int i = 0; i < n; ++i) {
        const double lam = ges.eigenvalues()(i);
        const double f = (lam > 0.0) ? std::sqrt(lam) / (2.0 * kPi) : 0.0;
        out.frequencies.push_back(f);
        DenseVec mode = DenseVec::Zero(nDofs);
        for (int a = 0; a < nf; ++a) mode(free[a]) = ges.eigenvectors()(a, i);
        out.modes.push_back(mode);
    }
    return out;
}

std::vector<double> computeBucklingFactors2D(const Model2D& m, int nModes) {
    const ModelResult r = solveModel2D(m);   // 参考荷载下的线性位移

    std::map<int, int> index;
    std::vector<int> ids;
    for (const auto& [id, c] : m.nodes) { index[id] = static_cast<int>(ids.size()); ids.push_back(id); }
    const int nDofs = 3 * static_cast<int>(ids.size());

    std::vector<char> active(m.elements.size(), 1);
    SparseMat K;
    DenseVec F;
    assembleLinear2D(m, active, index, nDofs, K, F);

    // 几何刚度：由各梁的轴力（受拉正）组装
    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        if (e.type != "BEAM") continue;
        const int i1 = index.at(e.conn[0]), i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);
        const double L = (p2 - p1).norm();
        EulerBeam2D be{e.E, e.A, e.I, p1, p2};
        const DenseMat T = be.transformation();
        Eigen::VectorXd ug(6);
        ug << r.U(3 * i1), r.U(3 * i1 + 1), r.U(3 * i1 + 2),
              r.U(3 * i2), r.U(3 * i2 + 1), r.U(3 * i2 + 2);
        const Eigen::VectorXd ul = T * ug;
        const double N = e.E * e.A / L * (ul(3) - ul(0));
        const DenseMat Kg = T.transpose() * beamGeometricStiffness2D(N, L) * T;
        const std::vector<int> g = {3 * i1, 3 * i1 + 1, 3 * i1 + 2,
                                    3 * i2, 3 * i2 + 1, 3 * i2 + 2};
        for (int a = 0; a < 6; ++a)
            for (int b = 0; b < 6; ++b)
                if (Kg(a, b) != 0.0) trips.emplace_back(g[a], g[b], Kg(a, b));
    }
    SparseMat Kg(nDofs, nDofs);
    Kg.setFromTriplets(trips.begin(), trips.end());
    Kg.makeCompressed();

    std::vector<int> cdofs;
    std::vector<double> cvals;
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3) continue;
        if (c.k > 0.0) continue;   // 弹性支承不作硬约束（弹簧已在 K 中）
        cdofs.push_back(3 * index.at(c.node) + (c.dof - 1));
        cvals.push_back(c.value);
    }
    std::vector<int> cd2;
    std::vector<double> cv2;
    autoConstrainFreeRotations(K, nDofs, cdofs, cvals, cd2, cv2);
    std::vector<char> isC(nDofs, 0);
    for (int d : cd2) isC[d] = 1;

    std::vector<int> free;
    for (int d = 0; d < nDofs; ++d) if (!isC[d]) free.push_back(d);
    const int nf = static_cast<int>(free.size());
    if (nf > 3000)
        std::cerr << "[警告] 自由度数 " << nf
                  << " 较大，屈曲分析当前为稠密特征求解，可能很慢或内存不足；建议改用稀疏子空间法\n";

    std::vector<double> out;
    if (nf == 0) return out;
    Eigen::MatrixXd Kd(nf, nf), Kgd(nf, nf);
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            Kd(a, b) = K.coeff(free[a], free[b]);
            Kgd(a, b) = Kg.coeff(free[a], free[b]);
        }

    // (K + λ·Kg)·φ = 0。用 K = L·L^T 的 Cholesky 变换化为标准对称特征值：
    //   A·ψ = μ·ψ，A = L^{-1}·Kg·L^{-T}，λ = -1/μ（μ<0 为屈曲模态）。
    // （Kg 轴向零刚度 → 半正定，不能直接喂给广义自伴求解器。）
    Eigen::LLT<Eigen::MatrixXd> llt(Kd);
    if (llt.info() != Eigen::Success) return {};
    const Eigen::MatrixXd L = llt.matrixL();
    const Eigen::MatrixXd Linv = L.inverse();
    const Eigen::MatrixXd A = Linv * Kgd * Linv.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(A);
    const int n = std::min(nModes, nf);
    std::vector<double> mu;
    for (int i = 0; i < ses.eigenvalues().size() && static_cast<int>(mu.size()) < n; ++i) {
        const double m = ses.eigenvalues()(i);   // 升序，最负的屈曲模态在前
        if (m < -1e-14 && std::isfinite(m)) mu.push_back(m);
    }
    for (double m : mu) out.push_back(-1.0 / m);
    return out;
}

double impactCoefficientFromFrequency(double fHz) {
    if (fHz < 1.5) return 0.05;
    if (fHz <= 14.0) return 0.1767 * std::log(fHz) - 0.0157;
    return 0.45;
}

ModelResult makeResult2D(const Model2D& m, const Assembled2D& a,
                                const DenseVec& U, const DenseVec& F0) {
    ModelResult r;
    r.U = U;
    r.R = a.K0 * U - F0;
    r.nodeIds = a.ids;
    r.coords.reserve(a.ids.size());
    for (int id : a.ids) r.coords.push_back(m.nodes.at(id));
    return r;
}

ModelResult solveModel2D(const Model2D& m,
                                SolverType type,
                                const DenseVec* U0,
                                const std::vector<char>* active0,
                                NonlinearDiag* diagOut,
                                ReusableSolver* reuse) {
    if (m.analysisType == "NONLINEAR")
        return solveModel2DNonlinear(m, type, U0, active0, diagOut, reuse);

    // 节点号 -> 0-based 索引（std::map 按键有序）。
    std::vector<int> ids;
    std::map<int, int> index;
    for (const auto& [id, c] : m.nodes) {
        index[id] = static_cast<int>(ids.size());
        ids.push_back(id);
    }
    const int nNodes = static_cast<int>(ids.size());
    const int nDofs = 3 * nNodes;

    // cable/strut 激活状态（初值全激活）。
    std::vector<char> active(m.elements.size(), 1);

    // 约束自由度。
    std::vector<int> cdofs;
    std::vector<double> cvals;
    for (const auto& c : m.constraints) {
        if (c.dof < 1 || c.dof > 3)
            throw std::runtime_error("约束自由度必须在 1~3 之间");
        if (c.k > 0.0) continue;   // 弹性支承不作硬约束（弹簧已在 K 中）
        cdofs.push_back(3 * index.at(c.node) + (c.dof - 1));
        cvals.push_back(c.value);
    }

    SparseMat K0;
    DenseVec F0, U;
    const int maxIter = 100;

    for (int iter = 0; ; ++iter) {
        // 组装 K / F（只含激活单元）。
        SparseMat K;
        DenseVec F;
        assembleLinear2D(m, active, index, nDofs, K, F);

        K0 = K;   // 施加约束前，用于反力计算
        F0 = F;
        std::vector<int> cd;
        std::vector<double> cv;
        autoConstrainFreeRotations(K, nDofs, cdofs, cvals, cd, cv);
        applyDirichlet(K, F, cd, cv);

        // 诊断：零对角线自由度。
        for (int d = 0; d < nDofs; ++d) {
            if (K.coeff(d, d) == 0.0) {
                const int ni = d / 3;
                const int ld = d % 3;
                throw std::runtime_error(
                    "存在零刚度自由度：节点 " + std::to_string(ids[ni]) +
                    " dof " + std::to_string(ld + 1) +
                    "（未约束且未被任何单元连接）");
            }
        }

        // 传入 reuse 时复用符号分解（结构不变则只做数值分解）。
        // 施工阶段推进、影响线扫描这类"结构不变、只有数值变"的场景收益最直接。
        if (reuse && reuse->reusable()) {
            reuse->factorizeReuse(K);
            U = reuse->solve(F);
        } else {
            U = solve(K, F, type);
        }
        validateSolution(K, F, U, "二维求解");

        // 检查 cable/strut 的开关状态。
        bool changed = false;
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei) {
            if (!active[ei]) continue;
            const auto& e = m.elements[ei];
            if (e.type != "CABLE" && e.type != "STRUT") continue;
            const double N = axialForce(e, U, index.at(e.conn[0]), index.at(e.conn[1]),
                                        m.nodes.at(e.conn[0]), m.nodes.at(e.conn[1]));
            if (e.type == "CABLE" && N <= 0.0) { active[ei] = 0; changed = true; }
            else if (e.type == "STRUT" && N >= 0.0) { active[ei] = 0; changed = true; }
        }
        if (!changed) break;
        if (iter + 1 >= maxIter)
            throw std::runtime_error("cable/strut 开关迭代未收敛");
    }

    const DenseVec R = K0 * U - F0;

    ModelResult result;
    result.U = U;
    result.R = R;
    result.nodeIds = ids;
    result.coords.reserve(nNodes);
    for (int id : ids) result.coords.push_back(m.nodes.at(id));
    return result;
}

std::vector<ElementForce> elementForces2D(const Model2D& m, const DenseVec& U) {
    std::map<int, int> index;
    int n = 0;
    for (const auto& [id, c] : m.nodes) index[id] = n++;

    std::vector<ElementForce> forces;
    forces.reserve(m.elements.size());

    for (const auto& e : m.elements) {
        const int i1 = index.at(e.conn[0]);
        const int i2 = index.at(e.conn[1]);
        const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
        const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);

        ElementForce f;
        f.id = e.id;
        f.type = e.type;

        if (isRigid(e.type)) {
            forces.push_back(f);
            continue;
        }

        if (isAxial(e.type)) {
            if (m.analysisType == "NONLINEAR") {
                const Eigen::Vector2d d0 = p2 - p1;
                const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
                const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
                const Eigen::Vector2d d = d0 + u2 - u1;
                const double L = d.norm();
                if (L <= 0.0) throw std::runtime_error("轴向单元当前长度为零");
                const Eigen::Vector2d c = d / L;
                const DenseVec fint = axialInternalForceNonlinear(e, U, i1, i2, p1, p2);
                const DenseVec fext = elementExternalEquivalentLoad(m, e, p1, p2);
                const double fi = fint(0) * c(0) + fint(1) * c(1);
                const double fj = fint(2) * c(0) + fint(3) * c(1);
                const double ei = fext(0) * c(0) + fext(1) * c(1);
                const double ej = fext(2) * c(0) + fext(3) * c(1);
                f.N1 = -fi + ei;
                f.N2 = fj - ej;
            } else {
                const Eigen::Vector2d d = p2 - p1;
                const double L = d.norm();
                const Eigen::Vector2d c = d / L;
                const Eigen::Vector2d u1(U(3 * i1 + 0), U(3 * i1 + 1));
                const Eigen::Vector2d u2(U(3 * i2 + 0), U(3 * i2 + 1));
                const double Nrestore = (e.E * e.A / L) * (u2 - u1).dot(c); // 恢复力（受拉）
                const DenseVec fe = elementEquivalentLoad(m, e, p1, p2); // [i,j] 全局
                const double fe_i = fe(0) * c(0) + fe(1) * c(1);
                const double fe_j = fe(2) * c(0) + fe(3) * c(1);
                f.N1 = Nrestore + fe_i;  // 节点1端轴力（受拉为正）
                f.N2 = Nrestore - fe_j;  // 节点2端轴力（受拉为正）
            }
            // 索/压杆松弛时内力按 0 输出。
            if (e.type == "CABLE" && f.N1 <= 0.0) { f.N1 = 0.0; f.N2 = 0.0; }
            else if (e.type == "STRUT" && f.N1 >= 0.0) { f.N1 = 0.0; f.N2 = 0.0; }
        } else if (e.type == "BEAM") {
            if (!e.releases.empty()) {
                const DenseVec fphys = releasedBeamLocalForces(m, e, U, i1, i2, p1, p2);
                f.N1 = -fphys(0); f.N2 = fphys(3);
                f.V1 = fphys(1);  f.V2 = -fphys(4);
                f.M1 = -fphys(2); f.M2 = fphys(5);
                forces.push_back(f);
                continue;
            }
            EulerBeam2D be{e.E, e.A, e.I, p1, p2};
            const DenseMat T = be.transformation();
            const DenseMat Kl = be.localStiffness(beamShearParameter(e, (p2 - p1).norm()));
            Eigen::VectorXd ug(6);
            ug << U(3 * i1 + 0), U(3 * i1 + 1), U(3 * i1 + 2),
                  U(3 * i2 + 0), U(3 * i2 + 1), U(3 * i2 + 2);
            const Eigen::VectorXd ul = T * ug;   // 局部位移
            const Eigen::VectorXd fl = Kl * ul;  // 恢复杆端力（局部）
            const DenseVec fe = elementEquivalentLoad(m, e, p1, p2); // 全局等效荷载
            const Eigen::VectorXd fel = T * fe;   // 全局 -> 局部
            const Eigen::VectorXd fphys = fl - fel; // 物理杆端力
            f.N1 = -fphys(0); // 节点1端轴力（受拉）
            f.N2 = fphys(3);  // 节点2端轴力（受拉）
            f.V1 = fphys(1);   // 截面剪力（左端）
            f.V2 = -fphys(4);  // 截面剪力（右端）
            f.M1 = -fphys(2);  // 截面弯矩（左端，受拉为正/下侧受拉）
            f.M2 = fphys(5);   // 截面弯矩（右端）
        }
        forces.push_back(f);
    }
    return forces;
}

std::vector<ElementStress> elementStress2D(const Model2D& m,
                                                  const std::vector<ElementForce>& forces) {
    std::map<int, const Model2D::Element*> byId;
    for (const auto& e : m.elements) byId[e.id] = &e;

    std::vector<ElementStress> out;
    out.reserve(forces.size());
    for (const auto& f : forces) {
        const Model2D::Element* e = byId.at(f.id);
        ElementStress s;
        s.id = f.id;
        s.type = f.type;
        if (isRigid(f.type)) {
            // 刚性连接杆不参与应力输出。
        } else if (isAxial(f.type)) {
            s.s1 = f.N1 / e->A;
            s.s2 = f.N2 / e->A;
        } else { // BEAM
            if (e->Wt <= 0.0 || e->Wb <= 0.0)
                throw std::runtime_error("梁单元应力需要 *ELEMENT 里的 WT 和 WB");
            s.stop1 = f.N1 / e->A - f.M1 / e->Wt;
            s.sbot1 = f.N1 / e->A + f.M1 / e->Wb;
            s.stop2 = f.N2 / e->A - f.M2 / e->Wt;
            s.sbot2 = f.N2 / e->A + f.M2 / e->Wb;
        }
        out.push_back(s);
    }
    return out;
}

std::vector<StageResult> solveModel2DStages(const Model2D& m,
                                                   SolverType type) {
    if (m.stages.empty())
        throw std::runtime_error("模型没有定义施工阶段 *STAGE");

    std::vector<Model2D::Stage> stages = m.stages;
    std::sort(stages.begin(), stages.end(),
              [](const Model2D::Stage& a, const Model2D::Stage& b) { return a.id < b.id; });

    std::vector<StageResult> results;
    std::vector<char> active(m.elements.size(), 0);
    std::vector<Model2D::Element> workingElements = m.elements;

    // 徐变折减基准：readModel2D 已填好；但程序化建模（不经文件解析）时是空的，
    // 这里兜底补一次——否则 baseModulusOf 会退化成 e.E，
    // 而 e.E 在上一阶段已被折减过，折减就会被逐级累乘（正是 P1-5 要修的 bug）。
    std::map<int, double> baseE = m.baseModulus;
    if (baseE.empty())
        for (const auto& e : m.elements) baseE[e.id] = e.E;

    // 各阶段的总体矩阵结构只在"单元激活/钝化"时才变，徐变折减与荷载变化都不改结构。
    // 用可复用句柄把符号分解跨阶段保留下来（P1-7）。
    ReusableSolver reuse(type);
    std::vector<Model2D::Load> cumulativeLoads;
    std::vector<Model2D::ElemLoad> cumulativeElemLoads;
    DenseVec prevU;
    bool hasPrevU = false;
    double prevStageTime = 0.0;
    std::map<int, double> prevAxialForce;
    std::map<int, double> prevM1, prevM2;
    std::vector<char> prevActive;
    bool hasPrevActive = false;
    // 单元加载龄期（天）：单元被激活的那个阶段的龄期。
    // 逐阶段 AEMM 用它算 Δφ = φ(t_now) - φ(t_load)。
    // 未记录（例如阶段未给 TIME）的单元按 0 处理，即"从龄期 0 开始持荷"。
    std::map<int, double> loadAge;

    for (const auto& st : stages) {
        for (int id : st.activate) {
            bool found = false;
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) { active[ei] = 1; found = true; break; }
            if (!found) throw std::runtime_error("*STAGE 激活了不存在的单元: " + std::to_string(id));
            // *STAGE 的 TIME 是"该阶段末的龄期"，单元在阶段开始时就安装并受力，
            // 因此加载龄期取**上一阶段末**的龄期（首个阶段为 0）。
            // 若误取 st.time，则"第 1 阶段 TIME=100"会得到 Δφ=φ(100)-φ(100)=0，
            // 徐变被整个抹掉。
            loadAge[id] = prevStageTime;
        }
        for (int id : st.deactivate) {
            bool found = false;
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) { active[ei] = 0; found = true; break; }
            if (!found) throw std::runtime_error("*STAGE 钝化了不存在的单元: " + std::to_string(id));
            // 钝化再激活视为重新加载，龄期重新起算。
            loadAge.erase(id);
        }

        double creepDPhi = 0.0;

        Model2D mc = m;
        mc.baseModulus = baseE;
        mc.elements.clear();
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
            if (active[ei]) mc.elements.push_back(workingElements[ei]);

        // 本阶段先移除指定施工临时荷载。
        cumulativeLoads.erase(
            std::remove_if(cumulativeLoads.begin(), cumulativeLoads.end(),
                           [&](const Model2D::Load& l) {
                               for (const auto& [node, dof] : st.unloadNodeDofs)
                                   if (l.node == node && l.dof == dof) return true;
                               return false;
                           }),
            cumulativeLoads.end());
        cumulativeElemLoads.erase(
            std::remove_if(cumulativeElemLoads.begin(), cumulativeElemLoads.end(),
                           [&](const Model2D::ElemLoad& el) {
                               return std::find(st.unloadElems.begin(), st.unloadElems.end(), el.elem) != st.unloadElems.end();
                           }),
            cumulativeElemLoads.end());

        cumulativeLoads.insert(cumulativeLoads.end(), st.loads.begin(), st.loads.end());
        cumulativeElemLoads.insert(cumulativeElemLoads.end(), st.elemLoads.begin(), st.elemLoads.end());
        mc.loads = cumulativeLoads;
        mc.elemLoads = cumulativeElemLoads;
        mc.prestresses = st.prestresses;
        // 本阶段是否走"按时变曲线逐阶段 AEMM"路径。
        bool creepByAge = false;
        if (st.time > 0.0 && (m.creepModel != "CURVE" || !m.creepCurve.points.empty())) {
            Model2D::CreepSpec cs;
            cs.all = true;
            const double phiNow = creepCoefficientAt(m, st.time);
            const double phiPrev = prevStageTime > 0.0
                                       ? creepCoefficientAt(m, prevStageTime)
                                       : 0.0;
            creepDPhi = (prevStageTime > 0.0) ? (phiNow - phiPrev) : phiNow;
            cs.phi = phiNow;   // 仅作为无龄期信息时的兜底值
            cs.chi = 0.8;
            mc.creepSpecs = {cs};
            creepByAge = true;
        } else {
            mc.creepSpecs = st.creepSpecs;
        }
        if (st.time > 0.0 && (m.shrinkageModel != "CURVE" || !m.shrinkageCurve.points.empty())) {
            Model2D::ShrinkageSpec sh;
            sh.all = true;
            const double ecsNow = shrinkageStrainAt(m, st.time);
            const double ecsPrev = prevStageTime > 0.0
                                       ? shrinkageStrainAt(m, prevStageTime)
                                       : 0.0;
            sh.ecs = (prevStageTime > 0.0) ? (ecsNow - ecsPrev) : ecsNow;
            mc.shrinkageSpecs = {sh};
        } else {
            mc.shrinkageSpecs = st.shrinkageSpecs;
        }
        mc.tempProfiles = st.tempProfiles;
        mc.uniformTemp = 0.0;
        mc.elemTemps.clear();
        mc.selfWeights.clear();

        // 施工阶段支座变化：先移除，再添加。
        auto constraints = mc.constraints;
        constraints.erase(
            std::remove_if(constraints.begin(), constraints.end(),
                           [&](const Model2D::Constraint& c) {
                               for (const auto& [node, dof] : st.removedConstraints)
                                   if (c.node == node && c.dof == dof) return true;
                               return false;
                           }),
            constraints.end());
        constraints.insert(constraints.end(), st.addedConstraints.begin(), st.addedConstraints.end());
        mc.constraints = constraints;

        // 徐变折减：逐阶段 AEMM。
        // 关键：模量始终从 baseModulus 重算，而不是在上一阶段折减后的
        // 基础上再乘一次 —— 旧实现那样做会让折减被逐级累乘
        // （E0/((1+χφ1)(1+χφ2)…)，而不是 AEMM 的 E0/(1+χ(φ(t)-φ(t0)))）。
        if (!mc.creepSpecs.empty())
            processTimeEffects(mc, creepByAge ? &loadAge : nullptr,
                               creepByAge ? st.time : 0.0);

        if (creepDPhi != 0.0 && !prevAxialForce.empty()) {
            for (const auto& e : mc.elements) {
                auto it = prevAxialForce.find(e.id);
                if (it == prevAxialForce.end() || it->second == 0.0) continue;
                const double P = creepDPhi * it->second;
                const Eigen::Vector2d p1 = m.nodes.at(e.conn[0]);
                const Eigen::Vector2d p2 = m.nodes.at(e.conn[1]);
                const Eigen::Vector2d d = p2 - p1;
                const double L = d.norm();
                if (L <= 0.0) continue;
                const Eigen::Vector2d c = d / L;
                if (isAxial(e.type)) {
                    mc.loads.push_back({e.conn[0], 1, -P * c(0)});
                    mc.loads.push_back({e.conn[0], 2, -P * c(1)});
                    mc.loads.push_back({e.conn[1], 1, P * c(0)});
                    mc.loads.push_back({e.conn[1], 2, P * c(1)});
                } else if (e.type == "BEAM") {
                    DenseVec fl = DenseVec::Zero(6);
                    fl(0) = -P;
                    fl(3) = P;
                    if (prevM1.count(e.id)) fl(2) = -creepDPhi * prevM1[e.id];
                    if (prevM2.count(e.id)) fl(5) = creepDPhi * prevM2[e.id];
                    EulerBeam2D be{e.E, e.A, e.I, p1, p2};
                    const DenseVec fe = be.transformation().transpose() * fl;
                    for (int d = 0; d < 3; ++d) {
                        mc.loads.push_back({e.conn[0], d + 1, fe(d)});
                        mc.loads.push_back({e.conn[1], d + 1, fe(d + 3)});
                    }
                }
            }
        }

        if (!mc.prestresses.empty())
            processPrestress(mc);

        if (st.time > 0.0)
            prevStageTime = st.time;

        cumulativeLoads = mc.loads;
        cumulativeElemLoads = mc.elemLoads;
        for (const auto& me : mc.elements) {
            for (auto& we : workingElements)
                if (we.id == me.id) { we.E = me.E; break; }
        }

        // 未激活节点暂按固定处理，避免出现零刚度自由度。
        std::map<int, bool> connected;
        for (const auto& e : mc.elements) {
            connected[e.conn[0]] = true;
            connected[e.conn[1]] = true;
        }
        for (const auto& [id, c] : m.nodes) {
            if (connected[id]) continue;
            for (int dof = 1; dof <= 3; ++dof)
                mc.constraints.push_back({id, dof, 0.0});
        }

        const DenseVec* U0 = (m.analysisType == "NONLINEAR" && hasPrevU) ? &prevU : nullptr;
        const std::vector<char>* A0 = (m.analysisType == "NONLINEAR" && hasPrevActive) ? &prevActive : nullptr;
        const ModelResult r = solveModel2D(mc, type, U0, A0, nullptr, &reuse);
        prevU = r.U;
        hasPrevU = true;
        StageResult sr;
        sr.stage = st.id;
        sr.U = r.U;
        sr.R = r.R;
        sr.forces = elementForces2D(mc, r.U);
        if (m.analysisType == "NONLINEAR") {
            prevActive.assign(mc.elements.size(), 1);
            for (std::size_t ei = 0; ei < mc.elements.size(); ++ei) {
                const auto& e = mc.elements[ei];
                if (e.type == "CABLE") prevActive[ei] = (sr.forces[ei].N1 > 0.0) ? 1 : 0;
                else if (e.type == "STRUT") prevActive[ei] = (sr.forces[ei].N1 < 0.0) ? 1 : 0;
                else prevActive[ei] = 1;
            }
            hasPrevActive = true;
        }
        for (const auto& f : sr.forces) {
            prevAxialForce[f.id] = f.N1;
            if (f.type == "BEAM") {
                prevM1[f.id] = f.M1;
                prevM2[f.id] = f.M2;
            }
        }
        try {
            sr.stresses = elementStress2D(mc, sr.forces);
        } catch (const std::exception&) {
            sr.stresses.clear();
        }
        results.push_back(sr);
    }
    return results;
}

std::vector<CamberResult> computeCamber(const Model2D& m,
                                               SolverType type) {
    const auto stages = solveModel2DStages(m, type);
    if (stages.empty()) return {};
    const DenseVec U = stages.back().U;
    std::map<int, int> index;
    int n = 0;
    for (const auto& [id, c] : m.nodes) index[id] = n++;

    std::vector<CamberResult> out;
    out.reserve(m.cambers.size());
    for (const auto& c : m.cambers) {
        const int d = 3 * index.at(c.node) + (c.dof - 1);
        const double disp = U(d);
        const Eigen::Vector2d coord = m.nodes.at(c.node);
        const double orig = (c.dof == 1) ? coord(0) : (c.dof == 2) ? coord(1) : 0.0;
        const double camber = c.value - disp;
        out.push_back({c.node, c.dof, camber, disp, orig, orig + camber});
    }
    return out;
}

CamberIterationResult iterateCamber(const Model2D& m,
                                           SolverType type,
                                           double tol,
                                           int maxIter) {
    if (m.cambers.empty())
        throw std::runtime_error("非线性预拱度迭代需要 *CAMBER");
    const int n = static_cast<int>(m.cambers.size());
    DenseVec x = DenseVec::Zero(n);

    std::map<int, int> nid;
    int nn = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = nn++;

    auto eval = [&](const DenseVec& s) -> DenseVec {
        Model2D mc = m;
        for (int k = 0; k < n; ++k) {
            const auto& cp = m.cambers[k];
            auto& coord = mc.nodes.at(cp.node);
            if (cp.dof == 1) coord(0) += s(k);
            else if (cp.dof == 2) coord(1) += s(k);
            else throw std::runtime_error("非线性预拱度目前只支持 UX/UY");
        }
        const auto stages = solveModel2DStages(mc, type);
        if (stages.empty()) throw std::runtime_error("施工阶段结果为空");
        const DenseVec U = stages.back().U;
        DenseVec y(n);
        for (int k = 0; k < n; ++k) {
            const int d = 3 * nid.at(m.cambers[k].node) + (m.cambers[k].dof - 1);
            y(k) = s(k) + U(d);
        }
        return y;
    };

    for (int iter = 0; iter < maxIter; ++iter) {
        const DenseVec y = eval(x);
        DenseVec r(n);
        for (int k = 0; k < n; ++k) r(k) = y(k) - m.cambers[k].value;
        if (r.norm() < tol) {
            CamberIterationResult res;
            res.cambers = x;
            res.finalDisplacements = y;
            res.iterations = iter;
            Model2D mc = m;
            for (int k = 0; k < n; ++k) {
                auto& coord = mc.nodes.at(m.cambers[k].node);
                if (m.cambers[k].dof == 1) coord(0) += x(k);
                else coord(1) += x(k);
            }
            res.stages = solveModel2DStages(mc, type);
            return res;
        }

        Eigen::MatrixXd J(n, n);
        for (int j = 0; j < n; ++j) {
            const double h = 1e-6 * std::max(1.0, std::abs(x(j)));
            DenseVec xp = x, xm = x;
            xp(j) += h;
            xm(j) -= h;
            J.col(j) = (eval(xp) - eval(xm)) / (2.0 * h);
        }
        const double lambda = 1e-8;
        Eigen::MatrixXd A = J.transpose() * J;
        for (int j = 0; j < n; ++j) A(j, j) += lambda * std::max(1.0, A(j, j));
        const Eigen::VectorXd dx = A.fullPivLu().solve(-(J.transpose() * r));
        if (!dx.allFinite())
            throw std::runtime_error("非线性预拱度迭代线性求解失败");
        x += dx;
    }
    throw std::runtime_error("非线性预拱度迭代未收敛");
}

FeedbackCorrection computeFeedbackCorrection(const Model2D& m,
                                                     SolverType type) {
    if (m.feedbacks.empty())
        throw std::runtime_error("反馈修正需要 *FEEDBACK");
    if (m.stageControls.empty())
        throw std::runtime_error("反馈修正需要 *STAGE_CONTROL");
    const auto& fb = m.feedbacks.front();
    std::map<int, int> nid;
    int n = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = n++;

    auto response = [&](double scale) {
        Model2D mc = m;
        const auto& ctrl = m.stageControls.front();
        for (auto& st : mc.stages) {
            if (st.id == ctrl.stage) {
                st.prestresses[ctrl.prestressIndex - 1].P *= scale;
                break;
            }
        }
        const auto stages = solveModel2DStages(mc, type);
        const int d = 3 * nid.at(fb.node) + (fb.dof - 1);
        return stages.back().U(d);
    };

    const double s0 = 1.0;
    const double s1 = 1.01;
    const double y0 = response(s0);
    const double y1 = response(s1);
    const double deriv = (y1 - y0) / (s1 - s0);
    FeedbackCorrection out;
    out.current = y0;
    out.derivative = deriv;
    out.correction = (std::abs(deriv) > 1e-20) ? (fb.measured - y0) / deriv : 0.0;
    return out;
}

std::vector<StageResult> solveModel2DStagesBackward(const Model2D& m,
                                                           SolverType type) {
    if (m.stages.empty())
        throw std::runtime_error("模型没有定义施工阶段 *STAGE");

    std::vector<Model2D::Stage> stages = m.stages;
    std::sort(stages.begin(), stages.end(),
              [](const Model2D::Stage& a, const Model2D::Stage& b) { return a.id < b.id; });

    std::vector<char> active(m.elements.size(), 0);
    std::vector<Model2D::Load> cumLoads;
    std::vector<Model2D::ElemLoad> cumElemLoads;
    for (const auto& st : stages) {
        for (int id : st.activate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 1;
        }
        for (int id : st.deactivate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 0;
        }
        cumLoads.insert(cumLoads.end(), st.loads.begin(), st.loads.end());
        cumElemLoads.insert(cumElemLoads.end(), st.elemLoads.begin(), st.elemLoads.end());
    }

    auto removeLoads = [&](std::vector<Model2D::Load>& vec, const std::vector<Model2D::Load>& rem) {
        vec.erase(std::remove_if(vec.begin(), vec.end(),
                                 [&](const Model2D::Load& l) {
                                     for (const auto& r : rem)
                                         if (l.node == r.node && l.dof == r.dof && l.value == r.value)
                                             return true;
                                     return false;
                                 }),
                  vec.end());
    };
    auto removeElemLoads = [&](std::vector<Model2D::ElemLoad>& vec,
                               const std::vector<Model2D::ElemLoad>& rem) {
        vec.erase(std::remove_if(vec.begin(), vec.end(),
                                 [&](const Model2D::ElemLoad& el) {
                                     for (const auto& r : rem)
                                         if (el.elem == r.elem && el.dir == r.dir &&
                                             el.q1 == r.q1 && el.q2 == r.q2)
                                             return true;
                                     return false;
                                 }),
                  vec.end());
    };

    auto buildActiveModel = [&](Model2D& mc) {
        mc.elements.clear();
        for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
            if (active[ei]) mc.elements.push_back(m.elements[ei]);
        mc.loads = cumLoads;
        mc.elemLoads = cumElemLoads;
        mc.prestresses.clear();
        mc.creepSpecs.clear();
        mc.shrinkageSpecs.clear();
        mc.tempProfiles.clear();
        mc.uniformTemp = 0.0;
        mc.elemTemps.clear();
        mc.selfWeights.clear();

        std::map<int, bool> connected;
        for (const auto& e : mc.elements) {
            connected[e.conn[0]] = true;
            connected[e.conn[1]] = true;
        }
        for (const auto& [id, c] : m.nodes) {
            if (connected[id]) continue;
            for (int dof = 1; dof <= 3; ++dof)
                mc.constraints.push_back({id, dof, 0.0});
        }
    };

    DenseVec prevU;
    bool hasPrevU = false;
    std::vector<char> prevActive;
    bool hasPrevActive = false;
    std::vector<StageResult> results;
    // 同正装：结构不变时跨阶段复用符号分解（P1-7）。
    ReusableSolver reuse(type);

    for (auto it = stages.rbegin(); it != stages.rend(); ++it) {
        const Model2D::Stage& st = *it;
        Model2D mc = m;
        buildActiveModel(mc);

        const DenseVec* U0 = (m.analysisType == "NONLINEAR" && hasPrevU) ? &prevU : nullptr;
        const std::vector<char>* A0 = (m.analysisType == "NONLINEAR" && hasPrevActive) ? &prevActive : nullptr;
        const ModelResult r = solveModel2D(mc, type, U0, A0, nullptr, &reuse);
        prevU = r.U;
        hasPrevU = true;

        StageResult sr;
        sr.stage = st.id;
        sr.U = r.U;
        sr.R = r.R;
        sr.forces = elementForces2D(mc, r.U);
        if (m.analysisType == "NONLINEAR") {
            prevActive.assign(mc.elements.size(), 1);
            for (std::size_t ei = 0; ei < mc.elements.size(); ++ei) {
                const auto& e = mc.elements[ei];
                if (e.type == "CABLE") prevActive[ei] = (sr.forces[ei].N1 > 0.0) ? 1 : 0;
                else if (e.type == "STRUT") prevActive[ei] = (sr.forces[ei].N1 < 0.0) ? 1 : 0;
                else prevActive[ei] = 1;
            }
            hasPrevActive = true;
        }
        try {
            sr.stresses = elementStress2D(mc, sr.forces);
        } catch (const std::exception&) {
            sr.stresses.clear();
        }
        results.push_back(sr);

        removeLoads(cumLoads, st.loads);
        removeElemLoads(cumElemLoads, st.elemLoads);
        for (int id : st.activate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 0;
        }
        for (int id : st.deactivate) {
            for (std::size_t ei = 0; ei < m.elements.size(); ++ei)
                if (m.elements[ei].id == id) active[ei] = 1;
        }
    }
    return results;
}

double forwardBackwardStageResidual(const Model2D& m,
                                           SolverType type) {
    const auto forward = solveModel2DStages(m, type);
    const auto backward = solveModel2DStagesBackward(m, type);
    if (forward.empty() || backward.empty())
        return 0.0;
    return (forward.back().U - backward.front().U).norm();
}

StageIterationResult iterateStagePrestress(const Model2D& m,
                                                  SolverType type,
                                                  double tol,
                                                  int maxIter) {
    if (m.stageControls.empty() || m.stageTargets.empty())
        throw std::runtime_error("多变量正装-倒拆迭代需要 *STAGE_CONTROL 和 *STAGE_TARGET");
    const int nv = static_cast<int>(m.stageControls.size());
    const int nt = static_cast<int>(m.stageTargets.size());
    const int nc = static_cast<int>(m.stageConstraints.size());
    DenseVec x = DenseVec::Ones(nv);

    std::map<int, int> nid;
    int nn = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = nn++;

    auto eval = [&](const DenseVec& s) -> DenseVec {
        Model2D mc = m;
        for (int k = 0; k < nv; ++k) {
            const auto& ctrl = m.stageControls[k];
            bool found = false;
            for (auto& st : mc.stages) {
                if (st.id != ctrl.stage) continue;
                if (ctrl.prestressIndex < 1 ||
                    ctrl.prestressIndex > static_cast<int>(st.prestresses.size()))
                    throw std::runtime_error("*STAGE_CONTROL 的 PRESTRESS 序号不存在");
                st.prestresses[ctrl.prestressIndex - 1].P *= s(k);
                found = true;
                break;
            }
            if (!found) throw std::runtime_error("*STAGE_CONTROL 的阶段号不存在");
        }

        const auto stages = solveModel2DStages(mc, type);
        if (stages.empty()) throw std::runtime_error("施工阶段结果为空");
        const DenseVec U = stages.back().U;
        DenseVec y(nt);
        for (int j = 0; j < nt; ++j) {
            const auto& tg = m.stageTargets[j];
            if (tg.kind == "DISP" || tg.kind == "REACTION") {
                const int d = 3 * nid.at(tg.node) + (tg.dof - 1);
                y(j) = (tg.kind == "REACTION") ? stages.back().R(d) : U(d);
            } else if (tg.kind == "FORCE") {
                double val = 0.0;
                bool found = false;
                for (const auto& f : stages.back().forces) {
                    if (f.id != tg.elem) continue;
                    found = true;
                    if (tg.component == "N") val = (tg.end == 1) ? f.N1 : f.N2;
                    else if (tg.component == "V") val = (tg.end == 1) ? f.V1 : f.V2;
                    else if (tg.component == "M") val = (tg.end == 1) ? f.M1 : f.M2;
                    else throw std::runtime_error("*STAGE_TARGET FORCE 未知分量: " + tg.component);
                    break;
                }
                if (!found) throw std::runtime_error("*STAGE_TARGET FORCE 找不到单元");
                y(j) = val;
            } else if (tg.kind == "STRESS") {
                double val = 0.0;
                bool found = false;
                for (const auto& s : stages.back().stresses) {
                    if (s.id != tg.elem) continue;
                    found = true;
                    if (tg.component == "AXIAL") val = (tg.end == 1) ? s.s1 : s.s2;
                    else if (tg.component == "TOP") val = (tg.end == 1) ? s.stop1 : s.stop2;
                    else if (tg.component == "BOT") val = (tg.end == 1) ? s.sbot1 : s.sbot2;
                    else throw std::runtime_error("*STAGE_TARGET STRESS 未知分量: " + tg.component);
                    break;
                }
                if (!found) throw std::runtime_error("*STAGE_TARGET STRESS 找不到单元");
                y(j) = val;
            } else {
                throw std::runtime_error("*STAGE_TARGET 未知类型: " + tg.kind);
            }
        }
        return y;
    };

    auto evalConstraints = [&](const DenseVec& s) -> DenseVec {
        Model2D mc = m;
        for (int k = 0; k < nv; ++k) {
            const auto& ctrl = m.stageControls[k];
            for (auto& st : mc.stages) {
                if (st.id != ctrl.stage) continue;
                st.prestresses[ctrl.prestressIndex - 1].P *= s(k);
                break;
            }
        }
        const auto stages = solveModel2DStages(mc, type);
        DenseVec c(nc);
        for (int k = 0; k < nc; ++k) {
            const auto& sc = m.stageConstraints[k];
            const auto& tg = sc.target;
            double resp = 0.0;
            if (tg.kind == "DISP" || tg.kind == "REACTION") {
                const int d = 3 * nid.at(tg.node) + (tg.dof - 1);
                resp = (tg.kind == "REACTION") ? stages.back().R(d) : stages.back().U(d);
            } else if (tg.kind == "FORCE") {
                for (const auto& f : stages.back().forces) {
                    if (f.id != tg.elem) continue;
                    if (tg.component == "N") resp = (tg.end == 1) ? f.N1 : f.N2;
                    else if (tg.component == "V") resp = (tg.end == 1) ? f.V1 : f.V2;
                    else if (tg.component == "M") resp = (tg.end == 1) ? f.M1 : f.M2;
                    break;
                }
            } else if (tg.kind == "STRESS") {
                for (const auto& sts : stages.back().stresses) {
                    if (sts.id != tg.elem) continue;
                    if (tg.component == "AXIAL") resp = (tg.end == 1) ? sts.s1 : sts.s2;
                    else if (tg.component == "TOP") resp = (tg.end == 1) ? sts.stop1 : sts.stop2;
                    else if (tg.component == "BOT") resp = (tg.end == 1) ? sts.sbot1 : sts.sbot2;
                    break;
                }
            }
            double violation = 0.0;
            if (resp < sc.lower) violation = resp - sc.lower;
            else if (resp > sc.upper) violation = resp - sc.upper;
            c(k) = sc.target.weight * violation;
        }
        return c;
    };

    for (int iter = 0; iter < maxIter; ++iter) {
        const DenseVec y = eval(x);
        DenseVec r(nt + nc);
        for (int j = 0; j < nt; ++j)
            r(j) = m.stageTargets[j].weight * (y(j) - m.stageTargets[j].value);
        if (nc > 0) r.tail(nc) = evalConstraints(x);
        if (r.norm() < tol) {
            StageIterationResult res;
            res.scales = x;
            res.finalValues = y;
            res.iterations = iter;
            Model2D mc = m;
            for (int k = 0; k < nv; ++k) {
                for (auto& st : mc.stages) {
                    if (st.id == m.stageControls[k].stage) {
                        st.prestresses[m.stageControls[k].prestressIndex - 1].P *= x(k);
                        break;
                    }
                }
            }
            res.stages = solveModel2DStages(mc, type);
            return res;
        }

        Eigen::MatrixXd J(nt + nc, nv);
        for (int j = 0; j < nv; ++j) {
            const double h = 1e-6 * std::max(1.0, std::abs(x(j)));
            DenseVec xp = x, xm = x;
            xp(j) += h;
            xm(j) -= h;
            const DenseVec yp = eval(xp);
            const DenseVec ym = eval(xm);
            DenseVec gp = (yp - ym) / (2.0 * h);
            for (int j = 0; j < nt; ++j) gp(j) *= m.stageTargets[j].weight;
            if (nc > 0) {
                const DenseVec cp = evalConstraints(xp);
                const DenseVec cm = evalConstraints(xm);
                const DenseVec gpc = (cp - cm) / (2.0 * h);
                DenseVec gfull(nt + nc);
                gfull.head(nt) = gp;
                gfull.tail(nc) = gpc;
                J.col(j) = gfull;
            } else {
                J.col(j) = gp;
            }
        }
        const double lambda = 1e-8;
        Eigen::MatrixXd A = J.transpose() * J;
        for (int j = 0; j < nv; ++j) A(j, j) += lambda * std::max(1.0, A(j, j));
        const Eigen::VectorXd rhs = -(J.transpose() * r);
        const Eigen::VectorXd dx = A.fullPivLu().solve(rhs);
        if (!dx.allFinite())
            throw std::runtime_error("多变量正装-倒拆迭代线性求解失败");
        x += dx;
        for (int k = 0; k < nv; ++k) {
            x(k) = std::max(m.stageControls[k].lower, std::min(m.stageControls[k].upper, x(k)));
        }
    }
    throw std::runtime_error("多变量正装-倒拆迭代未收敛");
}

double evaluateInfluenceResponse(const Model2D& m, const ModelResult& r,
                                        const Model2D::InfluenceRequest& spec) {
    std::map<int, int> nid;
    int n = 0;
    for (const auto& [id, c] : m.nodes) nid[id] = n++;

    if (spec.kind == "DISP") {
        const int d = 3 * nid.at(spec.node) + (spec.dof - 1);
        return r.U(d);
    }
    if (spec.kind == "REACTION") {
        const int d = 3 * nid.at(spec.node) + (spec.dof - 1);
        return r.R(d);
    }

    const std::vector<ElementForce> forces = elementForces2D(m, r.U);
    if (spec.kind == "FORCE") {
        for (const auto& f : forces) {
            if (f.id != spec.elem) continue;
            if (spec.component == "N") return spec.end == 1 ? f.N1 : f.N2;
            if (spec.component == "V") return spec.end == 1 ? f.V1 : f.V2;
            if (spec.component == "M") return spec.end == 1 ? f.M1 : f.M2;
            throw std::runtime_error("INFLUENCE FORCE 未知分量: " + spec.component);
        }
        throw std::runtime_error("INFLUENCE FORCE 找不到单元: " + std::to_string(spec.elem));
    }

    if (spec.kind == "STRESS") {
        const std::vector<ElementStress> stresses = elementStress2D(m, forces);
        for (const auto& s : stresses) {
            if (s.id != spec.elem) continue;
            if (spec.component == "AXIAL") return spec.end == 1 ? s.s1 : s.s2;
            if (spec.component == "TOP") return spec.end == 1 ? s.stop1 : s.stop2;
            if (spec.component == "BOT") return spec.end == 1 ? s.sbot1 : s.sbot2;
            throw std::runtime_error("INFLUENCE STRESS 未知分量: " + spec.component);
        }
        throw std::runtime_error("INFLUENCE STRESS 找不到单元: " + std::to_string(spec.elem));
    }

    throw std::runtime_error("INFLUENCE 未知类型: " + spec.kind);
}

InfluenceCurve computeInfluenceCurve(const Model2D& m,
                                             const Model2D::LaneLoad& lane,
                                             const Model2D::InfluenceRequest& spec,
                                             SolverType type) {
    if (lane.nodes.empty())
        throw std::runtime_error("LANE_LOAD 节点列表为空");

    InfluenceCurve out;
    out.label = spec.kind + " id=" + std::to_string(spec.id);
    out.s.reserve(lane.nodes.size());
    out.eta.reserve(lane.nodes.size());

    const int loadDof = (lane.dir == "GX") ? 1 : 2;

    // 车道纵向里程（与求解方式无关，先算出来）。
    double s = 0.0;
    for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
        if (k > 0) {
            const Eigen::Vector2d pa = m.nodes.at(lane.nodes[k - 1]);
            const Eigen::Vector2d pb = m.nodes.at(lane.nodes[k]);
            s += (pb - pa).norm();
        }
        out.s.push_back(s);
    }

    // 索/压杆的激活状态依赖位移，且非线性分析每步都要重算切线，
    // 这两种情况只能逐位置重新求解。
    bool needsReSolve = (m.analysisType == "NONLINEAR");
    for (const auto& e : m.elements)
        if (e.type == "CABLE" || e.type == "STRUT") { needsReSolve = true; break; }

    if (needsReSolve) {
        const ModelResult base = solveModel2D(m, type);
        const double baseValue = evaluateInfluenceResponse(m, base, spec);
        for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
            Model2D mc = m;
            mc.loads.push_back({lane.nodes[k], loadDof, 1.0});
            const ModelResult r = solveModel2D(mc, type);
            out.eta.push_back(evaluateInfluenceResponse(m, r, spec) - baseValue);
        }
    } else {
        // 线性模型：K 与加载位置无关，装配一次、分解一次，
        // 之后每个位置只做一次三角回代（旧实现每个位置都重新装配+重新分解）。
        const Assembled2D a = assembleModel2DLinear(m);

        // 旧实现在这里硬编码 SimplicialLDLT，直接无视调用方传入的 type，
        // 于是 CHOLMOD 等后端在影响线扫描里根本用不上。改用可复用句柄：
        // 既尊重 type，也把 analyze/factorize 分离的能力真正用起来。
        ReusableSolver solver(type);
        solver.factorizeReuse(a.K);

        const DenseVec Ubase = solver.solve(a.F);
        const ModelResult base = makeResult2D(m, a, Ubase, a.F0);
        const double baseValue = evaluateInfluenceResponse(m, base, spec);

        DenseVec F0k(a.nDofs), Fk(a.nDofs);
        for (std::size_t k = 0; k < lane.nodes.size(); ++k) {
            F0k = a.F0;
            Fk  = a.F;
            const int d = 3 * a.index.at(lane.nodes[k]) + (loadDof - 1);
            // 加载点若被约束，单位荷载由支座直接承担，对结构无影响。
            if (!a.isConstrained[d]) { F0k(d) += 1.0; Fk(d) += 1.0; }
            const DenseVec U = solver.solve(Fk);
            const ModelResult r = makeResult2D(m, a, U, F0k);
            out.eta.push_back(evaluateInfluenceResponse(m, r, spec) - baseValue);
        }
    }

    out.vmax = out.eta[0];
    out.vmin = out.eta[0];
    out.smax = out.s[0];
    out.smin = out.s[0];
    for (std::size_t k = 0; k < out.eta.size(); ++k) {
        if (out.eta[k] > out.vmax) { out.vmax = out.eta[k]; out.smax = out.s[k]; }
        if (out.eta[k] < out.vmin) { out.vmin = out.eta[k]; out.smin = out.s[k]; }
    }
    return out;
}

std::vector<double> laneTributaryWeights(const Model2D& m,
                                                const Model2D::LaneLoad& lane) {
    const std::size_t n = lane.nodes.size();
    std::vector<double> w(n, 0.0);
    if (n == 1) {
        w[0] = 1.0;
        return w;
    }

    std::vector<double> s(n, 0.0);
    for (std::size_t k = 1; k < n; ++k) {
        const Eigen::Vector2d a = m.nodes.at(lane.nodes[k - 1]);
        const Eigen::Vector2d b = m.nodes.at(lane.nodes[k]);
        s[k] = s[k - 1] + (b - a).norm();
    }

    w[0] = 0.5 * (s[1] - s[0]);
    for (std::size_t k = 1; k + 1 < n; ++k)
        w[k] = 0.5 * (s[k + 1] - s[k - 1]);
    w[n - 1] = 0.5 * (s[n - 1] - s[n - 2]);
    return w;
}

EnvelopeValue envelopeFromInfluence(const Model2D& m,
                                           const Model2D::LaneLoad& lane,
                                           const InfluenceCurve& curve) {
    const std::vector<double> w = laneTributaryWeights(m, lane);

    double posArea = 0.0, negArea = 0.0;
    for (std::size_t k = 0; k < curve.eta.size(); ++k) {
        if (curve.eta[k] > 0.0) posArea += curve.eta[k] * w[k];
        else negArea += curve.eta[k] * w[k];
    }

    const double qMaxArea = (lane.q >= 0.0) ? posArea : negArea;
    const double qMinArea = (lane.q >= 0.0) ? negArea : posArea;
    const double pMaxOrd = (lane.P >= 0.0) ? curve.vmax : curve.vmin;
    const double pMinOrd = (lane.P >= 0.0) ? curve.vmin : curve.vmax;
    const double pSMax = (lane.P >= 0.0) ? curve.smax : curve.smin;
    const double pSMin = (lane.P >= 0.0) ? curve.smin : curve.smax;

    EnvelopeValue ev;
    ev.label = curve.label;
    ev.maxVal = lane.q * qMaxArea + lane.P * pMaxOrd;
    ev.minVal = lane.q * qMinArea + lane.P * pMinOrd;
    ev.smax = pSMax;
    ev.smin = pSMin;
    return ev;
}

EnvelopeValue nonlinearEnvelopeFromLane(const Model2D& m,
                                               const Model2D::LaneLoad& lane,
                                               const Model2D::InfluenceRequest& spec,
                                               SolverType type) {
    // 先用线性影响线确定 q 的正负加载区，然后用真实 q、P 做非线性工况扫描。
    const InfluenceCurve lin = computeInfluenceCurve(m, lane, spec, type);
    const std::vector<double> w = laneTributaryWeights(m, lane);
    const int loadDof = (lane.dir == "GX") ? 1 : 2;

    EnvelopeValue ev;
    ev.label = "NONLINEAR " + lin.label;
    ev.maxVal = -std::numeric_limits<double>::max();
    ev.minVal = std::numeric_limits<double>::max();
    ev.smax = 0.0;
    ev.smin = 0.0;

    for (std::size_t p = 0; p < lane.nodes.size(); ++p) {
        Model2D mcMax = m;
        Model2D mcMin = m;

        for (std::size_t j = 0; j < lane.nodes.size(); ++j) {
            const double qj = lane.q * w[j];
            const double prod = lane.q * lin.eta[j];
            if (prod > 0.0)
                mcMax.loads.push_back({lane.nodes[j], loadDof, qj});
            else if (prod < 0.0)
                mcMin.loads.push_back({lane.nodes[j], loadDof, qj});
        }

        mcMax.loads.push_back({lane.nodes[p], loadDof, lane.P});
        mcMin.loads.push_back({lane.nodes[p], loadDof, lane.P});

        const ModelResult rMax = solveModel2D(mcMax, type);
        const double vMax = evaluateInfluenceResponse(m, rMax, spec);
        if (vMax > ev.maxVal) { ev.maxVal = vMax; ev.smax = lin.s[p]; }

        const ModelResult rMin = solveModel2D(mcMin, type);
        const double vMin = evaluateInfluenceResponse(m, rMin, spec);
        if (vMin < ev.minVal) { ev.minVal = vMin; ev.smin = lin.s[p]; }
    }
    return ev;
}

std::vector<EnvelopeSpec> defaultEnvelopeSpecs(const Model2D& m) {
    std::vector<EnvelopeSpec> specs;
    int id = 0;
    for (const auto& e : m.elements) {
        if (isRigid(e.type)) continue;

        for (int end = 1; end <= 2; ++end) {
            Model2D::InfluenceRequest req;
            req.id = ++id;
            req.kind = "FORCE";
            req.elem = e.id;
            req.end = end;
            req.component = "N";
            specs.push_back({req, "E" + std::to_string(e.id) + "_N" + std::to_string(end)});

            if (e.type == "BEAM") {
                req.id = ++id;
                req.component = "V";
                specs.push_back({req, "E" + std::to_string(e.id) + "_V" + std::to_string(end)});
                req.id = ++id;
                req.component = "M";
                specs.push_back({req, "E" + std::to_string(e.id) + "_M" + std::to_string(end)});
            }
        }

        if (isAxial(e.type)) {
            for (int end = 1; end <= 2; ++end) {
                Model2D::InfluenceRequest req;
                req.id = ++id;
                req.kind = "STRESS";
                req.elem = e.id;
                req.end = end;
                req.component = "AXIAL";
                specs.push_back({req, "E" + std::to_string(e.id) + "_S" + std::to_string(end)});
            }
        } else if (e.type == "BEAM" && e.Wt > 0.0 && e.Wb > 0.0) {
            for (int end = 1; end <= 2; ++end) {
                const std::vector<std::string> comps = {"TOP", "BOT"};
                for (const std::string& comp : comps) {
                    Model2D::InfluenceRequest req;
                    req.id = ++id;
                    req.kind = "STRESS";
                    req.elem = e.id;
                    req.end = end;
                    req.component = comp;
                    specs.push_back({req, "E" + std::to_string(e.id) + "_" + comp + std::to_string(end)});
                }
            }
        }
    }
    return specs;
}

void writeResult(const std::string& path, const Model2D& model,
                        const ModelResult& result,
                        const std::vector<ElementForce>& forces) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法写入结果文件: " + path);

    out << std::setprecision(8);
    out << "# FEM 计算结果\n";
    out << "# 模型: " << (model.title.empty() ? "未命名" : model.title) << "\n\n";

    out << "节点位移 (m, rad)\n";
    out << "节点\tu(x)\t\tv(y)\t\ttheta\n";
    for (std::size_t i = 0; i < result.nodeIds.size(); ++i) {
        const int d = 3 * static_cast<int>(i);
        out << result.nodeIds[i] << '\t'
            << result.U(d + 0) << '\t' << result.U(d + 1) << '\t'
            << result.U(d + 2) << '\n';
    }

    out << "\n单元内力 (N, N*m)\n";
    out << "单元\t类型\tN1(轴力i)\tN2(轴力j)\tV1(剪力)\tM1(弯矩)\tV2(剪力)\tM2(弯矩)\n";
    for (const auto& f : forces) {
        if (isRigid(f.type)) continue;
        out << f.id << '\t' << f.type << '\t' << f.N1 << '\t' << f.N2;
        if (f.type == "BEAM") {
            out << '\t' << f.V1 << '\t' << f.M1 << '\t' << f.V2 << '\t' << f.M2;
        } else {
            out << "\t0\t0\t0\t0";
        }
        out << '\n';
    }

    try {
        const std::vector<ElementStress> stresses = elementStress2D(model, forces);
        out << "\n单元应力 (Pa)\n";
        out << "单元\t类型\tσ1\t\tσ2\t\tσtop1\t\tσbot1\t\tσtop2\t\tσbot2\n";
        for (const auto& s : stresses) {
            if (isRigid(s.type)) continue;
            if (isAxial(s.type)) {
                out << s.id << "\t" << s.type << "\t" << s.s1 << '\t' << s.s2 << "\t-\t-\t-\t-\n";
            } else {
                out << s.id << "\tBEAM\t-\t-\t" << s.stop1 << '\t' << s.sbot1
                    << '\t' << s.stop2 << '\t' << s.sbot2 << '\n';
            }
        }
    } catch (const std::runtime_error& e) {
        out << "\n单元应力: (跳过，原因: " << e.what() << ")\n";
    }

    out << "\n支座反力 (N, N*m)\n";
    out << "节点\tdof\t反力\n";
    double maxAbsR = 0.0;
    for (int d = 0; d < result.R.size(); ++d) maxAbsR = std::max(maxAbsR, std::abs(result.R(d)));
    const double rThresh = 1e-9 * std::max(1.0, maxAbsR); // 相对阈值，避免大模型漏噪声、小模型滤真值
    for (std::size_t i = 0; i < result.nodeIds.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            const double v = result.R(3 * static_cast<int>(i) + d);
            if (std::abs(v) > rThresh)
                out << result.nodeIds[i] << '\t' << (d + 1) << '\t' << v << '\n';
        }
    }
}

} // namespace fem
