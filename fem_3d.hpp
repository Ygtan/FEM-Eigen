#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fem {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using SparseMat3D = Eigen::SparseMatrix<double>;

struct Model3D;
inline std::vector<double> solveModel3D(const Model3D& m);

// 三维欧拉-伯努利梁局部刚度矩阵，12 自由度：
// [u1,v1,w1,θx1,θy1,θz1,u2,v2,w2,θx2,θy2,θz2]
inline Eigen::Matrix<double, 12, 12> beam3DLocalStiffness(
    double E, double A, double Iy, double Iz, double J, double G, double L) {
    if (L <= 0.0) throw std::runtime_error("三维梁单元长度为零");
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    const double EAL = E * A / L;
    const double kz = E * Iz / (L * L * L);
    const double ky = E * Iy / (L * L * L);
    const double GJL = G * J / L;

    K(0, 0) = EAL;  K(0, 6) = -EAL;
    K(6, 0) = -EAL; K(6, 6) = EAL;

    K(1, 1) = 12 * kz; K(1, 5) = 6 * kz * L; K(1, 7) = -12 * kz; K(1, 11) = 6 * kz * L;
    K(5, 1) = 6 * kz * L; K(5, 5) = 4 * kz * L * L; K(5, 7) = -6 * kz * L; K(5, 11) = 2 * kz * L * L;
    K(7, 1) = -12 * kz; K(7, 5) = -6 * kz * L; K(7, 7) = 12 * kz; K(7, 11) = -6 * kz * L;
    K(11, 1) = 6 * kz * L; K(11, 5) = 2 * kz * L * L; K(11, 7) = -6 * kz * L; K(11, 11) = 4 * kz * L * L;

    K(2, 2) = 12 * ky; K(2, 4) = -6 * ky * L; K(2, 8) = -12 * ky; K(2, 10) = -6 * ky * L;
    K(4, 2) = -6 * ky * L; K(4, 4) = 4 * ky * L * L; K(4, 8) = 6 * ky * L; K(4, 10) = 2 * ky * L * L;
    K(8, 2) = -12 * ky; K(8, 4) = 6 * ky * L; K(8, 8) = 12 * ky; K(8, 10) = 6 * ky * L;
    K(10, 2) = -6 * ky * L; K(10, 4) = 2 * ky * L * L; K(10, 8) = 6 * ky * L; K(10, 10) = 4 * ky * L * L;

    K(3, 3) = GJL; K(3, 9) = -GJL;
    K(9, 3) = -GJL; K(9, 9) = GJL;
    return K;
}

inline Eigen::Matrix<double, 12, 12> beam3DTransformation(const Vec3& p1, const Vec3& p2) {
    Vec3 ex = p2 - p1;
    const double L = std::max(ex.norm(), 1e-12);
    ex /= L;
    Vec3 ez(0.0, 0.0, 1.0);
    if (std::abs(ex.dot(ez)) > 0.999) ez = Vec3(0.0, 1.0, 0.0);
    Vec3 ey = ez.cross(ex); ey.normalize();
    ez = ex.cross(ey); ez.normalize();
    Eigen::Matrix3d R;
    R.col(0)=ex; R.col(1)=ey; R.col(2)=ez;
    Eigen::Matrix<double, 12, 12> T = Eigen::Matrix<double, 12, 12>::Zero();
    for (int b=0; b<4; ++b) T.block<3,3>(3*b,3*b)=R;
    return T;
}

inline Eigen::Matrix<double, 6, 6> truss3DStiffness(double E, double A, const Vec3& p1, const Vec3& p2) {
    const Vec3 d = p2 - p1;
    const double L = d.norm();
    if (L <= 0.0) throw std::runtime_error("三维桁架单元长度为零");
    const Vec3 c = d / L;
    const Mat3 K = (E * A / L) * (c * c.transpose());
    Eigen::Matrix<double, 6, 6> Ke = Eigen::Matrix<double, 6, 6>::Zero();
    Ke.block<3, 3>(0, 0) = K;
    Ke.block<3, 3>(0, 3) = -K;
    Ke.block<3, 3>(3, 0) = -K;
    Ke.block<3, 3>(3, 3) = K;
    return Ke;
}

struct Node3D {
    int id = 0;
    Vec3 x = Vec3::Zero();
};

struct Element3D {
    int id = 0;
    std::string type = "TRUSS";
    double E = 0.0, A = 0.0;
    double Iy = 0.0, Iz = 0.0, J = 0.0, G = 0.0;
    double istrn = 0.0;
    double Wtop = 0.0, Wbot = 0.0;
    int n1 = 0, n2 = 0;
};

struct Model3D {
    std::map<int, Vec3> nodes;
    std::vector<Element3D> elements;
    std::vector<std::pair<int, int>> constraints; // node, dof 1..3
    std::vector<std::tuple<int, int, double>> loads;

    struct Stage3D {
        int id = 0;
        std::vector<int> activate;
        std::vector<int> deactivate;
        std::vector<std::tuple<int, int, double>> loads;
    };
    std::vector<Stage3D> stages;
    double creepPhi = 0.0;
    double creepChi = 0.8;
    double shrinkageEcs = 0.0;

    struct Influence3D {
        int node = 0;
        int dof = 1;
    };
    std::vector<Influence3D> influences;

    struct ElementInfluence3D {
        int elem = 0;
    };
    std::vector<ElementInfluence3D> elementInfluences;

    struct Lane3D {
        std::vector<int> nodes;
        double q = 0.0;
        double p = 0.0;
        int dof = 2;
        double reduction = 1.0;
        double impact = 1.0;
    };
    std::vector<Lane3D> lanes;

    struct LaneSurface3D {
        std::vector<int> xNodes;
        std::vector<int> yNodes;
        double q = 0.0;
        double p = 0.0;
        int dof = 2;
        double reduction = 1.0;
        int laneCount = 1;
        double impact = 1.0;
    };
    LaneSurface3D laneSurface;
};

struct StageResult3D {
    int stage = 0;
    std::vector<double> U;
};

inline std::vector<StageResult3D> solveModel3DStages(const Model3D& m) {
    std::vector<char> active(m.elements.size(), 0);
    std::vector<StageResult3D> out;
    for (const auto& st : m.stages) {
        for (int id : st.activate) {
            bool found=false;
            for (std::size_t i=0;i<m.elements.size();++i) if(m.elements[i].id==id){active[i]=1;found=true;break;}
            if(!found) throw std::runtime_error("三维阶段激活了不存在单元");
        }
        for (int id : st.deactivate) {
            bool found=false;
            for (std::size_t i=0;i<m.elements.size();++i) if(m.elements[i].id==id){active[i]=0;found=true;break;}
            if(!found) throw std::runtime_error("三维阶段钝化了不存在单元");
        }
        Model3D mc;
        mc.nodes = m.nodes;
        mc.constraints = m.constraints;
        mc.loads = st.loads;
        for (std::size_t i=0;i<m.elements.size();++i) if(active[i]) mc.elements.push_back(m.elements[i]);
        std::map<int,bool> connected;
        for (const auto& e : mc.elements) { connected[e.n1]=true; connected[e.n2]=true; }
        bool hasBeamStage = false;
        for (const auto& e : mc.elements) if (e.type=="BEAM") { hasBeamStage=true; break; }
        const int dofPerNodeStage = hasBeamStage ? 6 : 3;
        for (const auto& [id,x] : m.nodes)
            if (!connected[id])
                for (int dof=1; dof<=dofPerNodeStage; ++dof)
                    mc.constraints.push_back({id,dof});
        StageResult3D sr;
        sr.stage = st.id;
        sr.U = solveModel3D(mc);
        out.push_back(sr);
    }
    return out;
}

struct Envelope3D {
    int id = 0;
    double Nmin = 0.0, Nmax = 0.0;
};

inline std::vector<Envelope3D> compute3DStageEnvelope(const Model3D& m,
                                                      const std::vector<StageResult3D>& stages) {
    bool hasBeam = false;
    for (const auto& e : m.elements) if (e.type == "BEAM") { hasBeam = true; break; }
    const int dofPerNode = hasBeam ? 6 : 3;
    std::map<int, Envelope3D> env;
    for (const auto& sr : stages) {
        std::map<int, int> idx;
        int n = 0;
        for (const auto& [id, x] : m.nodes) idx[id] = n++;
        Model3D mc;
        mc.nodes = m.nodes;
        for (const auto& e : m.elements) {
            // 用当前阶段位移计算轴力
            const int i1=idx.at(e.n1), i2=idx.at(e.n2);
            const Vec3 d = m.nodes.at(e.n2)-m.nodes.at(e.n1);
            const double L = std::max(d.norm(),1e-12);
            const Vec3 c=d/L;
            const double dL = (sr.U[dofPerNode*i2]-sr.U[dofPerNode*i1])*c(0)
                            + (sr.U[dofPerNode*i2+1]-sr.U[dofPerNode*i1+1])*c(1)
                            + (sr.U[dofPerNode*i2+2]-sr.U[dofPerNode*i1+2])*c(2);
            const double N = e.E*e.A/L*dL;
            auto& r=env[e.id]; r.id=e.id;
            r.Nmin=std::min(r.Nmin,N); r.Nmax=std::max(r.Nmax,N);
        }
    }
    std::vector<Envelope3D> out;
    for (const auto& [id,r] : env) out.push_back(r);
    return out;
}

struct Camber3D {
    int node = 0, dof = 0;
    double camber = 0.0;
};

struct DisplacementEnvelope3D {
    int node = 0, dof = 0;
    double minVal = 0.0, maxVal = 0.0;
};

struct BeamEnvelope3D {
    int id = 0;
    double Nmin=0,Nmax=0,VyMin=0,VyMax=0,VzMin=0,VzMax=0;
    double MyMin=0,MyMax=0,MzMin=0,MzMax=0;
};

inline std::vector<BeamEnvelope3D> compute3DBeamEnvelope(
    const Model3D& m, const std::vector<StageResult3D>& stages) {
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    std::map<int, BeamEnvelope3D> env;
    for (const auto& sr : stages) {
        for (const auto& e : m.elements) {
            if (e.type != "BEAM") continue;
            const int i1=idx.at(e.n1), i2=idx.at(e.n2);
            Eigen::VectorXd ug(12);
            for (int k=0;k<6;++k){ ug(k)=sr.U[6*i1+k]; ug(6+k)=sr.U[6*i2+k]; }
            const Vec3 p1=m.nodes.at(e.n1), p2=m.nodes.at(e.n2);
            const auto T=beam3DTransformation(p1,p2);
            const Eigen::VectorXd ul=T*ug;
            const double L=(p2-p1).norm();
            const auto Kl=beam3DLocalStiffness(e.E,e.A,e.Iy,e.Iz,e.J,e.G,L);
            const Eigen::VectorXd fl=Kl*ul;
            auto& r=env[e.id]; r.id=e.id;
            double N=-fl(0), Vy=fl(1), Vz=fl(2), My=fl(4), Mz=fl(5);
            r.Nmin=std::min(r.Nmin,N); r.Nmax=std::max(r.Nmax,N);
            r.VyMin=std::min(r.VyMin,Vy); r.VyMax=std::max(r.VyMax,Vy);
            r.VzMin=std::min(r.VzMin,Vz); r.VzMax=std::max(r.VzMax,Vz);
            r.MyMin=std::min(r.MyMin,My); r.MyMax=std::max(r.MyMax,My);
            r.MzMin=std::min(r.MzMin,Mz); r.MzMax=std::max(r.MzMax,Mz);
        }
    }
    std::vector<BeamEnvelope3D> out;
    for (auto& [id,r] : env) out.push_back(r);
    return out;
}

inline std::vector<DisplacementEnvelope3D> compute3DDisplacementEnvelope(
    const Model3D& m, const std::vector<StageResult3D>& stages) {
    bool hasBeam = false;
    for (const auto& e : m.elements) if (e.type == "BEAM") { hasBeam = true; break; }
    const int dofPerNode = hasBeam ? 6 : 3;
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    std::vector<DisplacementEnvelope3D> out;
    for (const auto& [id,x] : m.nodes) {
        for (int dof=1; dof<=dofPerNode; ++dof) {
            double mn=1e300, mx=-1e300;
            for (const auto& s : stages) {
                double v = s.U[dofPerNode*idx[id]+dof-1];
                mn=std::min(mn,v); mx=std::max(mx,v);
            }
            out.push_back({id,dof,mn,mx});
        }
    }
    return out;
}

inline std::vector<double> compute3DInfluence(const Model3D& m,
                                               const Model3D::Influence3D& spec,
                                               const Model3D::Lane3D& lane) {
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    bool hasBeam=false;
    for (const auto& e : m.elements) if (e.type=="BEAM") { hasBeam=true; break; }
    const int dofPerNode = hasBeam ? 6 : 3;
    std::vector<double> eta;
    Model3D base = m;
    base.loads.clear();
    auto baseU = solveModel3D(base);
    const double baseResp = baseU[dofPerNode*idx.at(spec.node)+spec.dof-1];
    for (int node : lane.nodes) {
        Model3D mc = m;
        mc.loads.clear();
        mc.loads.push_back({node,lane.dof,1.0});
        auto U = solveModel3D(mc);
        eta.push_back(U[dofPerNode*idx.at(spec.node)+spec.dof-1] - baseResp);
    }
    return eta;
}

inline std::vector<double> compute3DBeamAxialInfluence(const Model3D& m,
                                                       int elemId,
                                                       const Model3D::Lane3D& lane) {
    const Element3D* target=nullptr;
    for (const auto& e : m.elements) if (e.id==elemId) { target=&e; break; }
    if (!target) throw std::runtime_error("三维梁单元不存在");
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    bool hasBeam=false; for(const auto& e:m.elements) if(e.type=="BEAM"){hasBeam=true;break;}
    const int dofPerNode=hasBeam?6:3;
    Model3D base=m; base.loads.clear();
    auto baseU=solveModel3D(base);
    auto axial=[&](const std::vector<double>& U){
        const int i1=idx.at(target->n1), i2=idx.at(target->n2);
        const Vec3 d=m.nodes.at(target->n2)-m.nodes.at(target->n1);
        const double L=std::max(d.norm(),1e-12); const Vec3 c=d/L;
        double dl=0;
        for(int k=0;k<3;++k) dl += (U[dofPerNode*i2+k]-U[dofPerNode*i1+k])*c(k);
        return target->E*target->A/L*dl;
    };
    const double baseN=axial(baseU);
    std::vector<double> eta;
    for(int node:lane.nodes){
        Model3D mc=m; mc.loads.clear(); mc.loads.push_back({node,lane.dof,1.0});
        auto U=solveModel3D(mc);
        eta.push_back(axial(U)-baseN);
    }
    return eta;
}

inline std::vector<double> compute3DBeamMomentInfluence(const Model3D& m,
                                                        int elemId,
                                                        const std::string& comp,
                                                        int end,
                                                        const Model3D::Lane3D& lane) {
    const Element3D* target=nullptr;
    for (const auto& e : m.elements) if (e.id==elemId) { target=&e; break; }
    if (!target || target->type!="BEAM") throw std::runtime_error("三维梁单元不存在");
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    auto localForce=[&](const std::vector<double>& U)->Eigen::VectorXd{
        const int i1=idx.at(target->n1), i2=idx.at(target->n2);
        Eigen::VectorXd ug(12);
        for(int k=0;k<6;++k){ ug(k)=U[6*i1+k]; ug(6+k)=U[6*i2+k]; }
        const Vec3 p1=m.nodes.at(target->n1), p2=m.nodes.at(target->n2);
        const auto T=beam3DTransformation(p1,p2);
        const double L=(p2-p1).norm();
        const auto Kl=beam3DLocalStiffness(target->E,target->A,target->Iy,target->Iz,target->J,target->G,L);
        return Kl*(T*ug);
    };
    auto baseU=solveModel3D(m);
    const auto baseF=localForce(baseU);
    auto moment=[&](const Eigen::VectorXd& fl){
        int idxLocal = (end==1)?4:10;
        if (comp=="MY") idxLocal += 0; else if (comp=="MZ") idxLocal += 1; else throw std::runtime_error("未知弯矩分量");
        return fl(idxLocal);
    };
    const double baseM=moment(baseF);
    std::vector<double> eta;
    for(int node:lane.nodes){
        Model3D mc=m; mc.loads.clear(); mc.loads.push_back({node,lane.dof,1.0});
        auto U=solveModel3D(mc);
        eta.push_back(moment(localForce(U))-baseM);
    }
    return eta;
}

inline std::vector<double> compute3DBeamShearInfluence(const Model3D& m,
                                                       int elemId,
                                                       const std::string& comp,
                                                       int end,
                                                       const Model3D::Lane3D& lane) {
    const Element3D* target=nullptr;
    for (const auto& e : m.elements) if (e.id==elemId) { target=&e; break; }
    if (!target || target->type!="BEAM") throw std::runtime_error("三维梁单元不存在");
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    auto localForce=[&](const std::vector<double>& U)->Eigen::VectorXd{
        const int i1=idx.at(target->n1), i2=idx.at(target->n2);
        Eigen::VectorXd ug(12);
        for(int k=0;k<6;++k){ ug(k)=U[6*i1+k]; ug(6+k)=U[6*i2+k]; }
        const Vec3 p1=m.nodes.at(target->n1), p2=m.nodes.at(target->n2);
        const auto T=beam3DTransformation(p1,p2);
        const double L=(p2-p1).norm();
        const auto Kl=beam3DLocalStiffness(target->E,target->A,target->Iy,target->Iz,target->J,target->G,L);
        return Kl*(T*ug);
    };
    auto baseU=solveModel3D(m);
    const auto baseF=localForce(baseU);
    auto shear=[&](const Eigen::VectorXd& fl){
        if (comp=="VY") return (end==1)?fl(1):fl(7);
        if (comp=="VZ") return (end==1)?fl(2):fl(8);
        throw std::runtime_error("未知剪力分量");
    };
    const double baseV=shear(baseF);
    std::vector<double> eta;
    for(int node:lane.nodes){
        Model3D mc=m; mc.loads.clear(); mc.loads.push_back({node,lane.dof,1.0});
        auto U=solveModel3D(mc);
        eta.push_back(shear(localForce(U))-baseV);
    }
    return eta;
}

inline std::vector<double> compute3DBeamStressInfluence(const Model3D& m,
                                                        int elemId,
                                                        const std::string& comp,
                                                        int end,
                                                        const Model3D::Lane3D& lane) {
    const Element3D* target=nullptr;
    for (const auto& e : m.elements) if (e.id==elemId) { target=&e; break; }
    if (!target || target->type!="BEAM") throw std::runtime_error("三维梁单元不存在");
    if (target->Wtop<=0 || target->Wbot<=0) throw std::runtime_error("梁应力需要 WT/WB");
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    auto localForce=[&](const std::vector<double>& U)->Eigen::VectorXd{
        const int i1=idx.at(target->n1), i2=idx.at(target->n2);
        Eigen::VectorXd ug(12);
        for(int k=0;k<6;++k){ ug(k)=U[6*i1+k]; ug(6+k)=U[6*i2+k]; }
        const Vec3 p1=m.nodes.at(target->n1), p2=m.nodes.at(target->n2);
        const auto T=beam3DTransformation(p1,p2);
        const double L=(p2-p1).norm();
        const auto Kl=beam3DLocalStiffness(target->E,target->A,target->Iy,target->Iz,target->J,target->G,L);
        return Kl*(T*ug);
    };
    auto baseU=solveModel3D(m);
    const auto baseF=localForce(baseU);
    auto stress=[&](const Eigen::VectorXd& fl){
        double N = (end==1)? -fl(0) : fl(6);
        double Mz = (end==1)? fl(5) : fl(11);
        if (comp=="TOP") return N/target->A - Mz/target->Wtop;
        if (comp=="BOT") return N/target->A + Mz/target->Wbot;
        throw std::runtime_error("未知应力分量");
    };
    const double baseS=stress(baseF);
    std::vector<double> eta;
    for(int node:lane.nodes){
        Model3D mc=m; mc.loads.clear(); mc.loads.push_back({node,lane.dof,1.0});
        auto U=solveModel3D(mc);
        eta.push_back(stress(localForce(U))-baseS);
    }
    return eta;
}

inline std::vector<std::tuple<int,int,double>> compute3DInfluenceSurface(
    const Model3D& m, const Model3D::Influence3D& spec,
    const Model3D::LaneSurface3D& surface) {
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    bool hasBeam=false;
    for (const auto& e : m.elements) if (e.type=="BEAM") { hasBeam=true; break; }
    const int dofPerNode = hasBeam ? 6 : 3;
    Model3D base = m; base.loads.clear();
    auto baseU = solveModel3D(base);
    const double baseResp = baseU[dofPerNode*idx.at(spec.node)+spec.dof-1];
    std::vector<std::tuple<int,int,double>> out;
    for (int yn : surface.yNodes) {
        for (int xn : surface.xNodes) {
            Model3D mc = m; mc.loads.clear(); mc.loads.push_back({xn,surface.dof,1.0});
            auto U = solveModel3D(mc);
            double resp = U[dofPerNode*idx.at(spec.node)+spec.dof-1] - baseResp;
            out.push_back({xn,yn,resp});
        }
    }
    return out;
}

inline std::pair<double,double> envelope3DFromSurface(
    const Model3D::LaneSurface3D& surface,
    const std::vector<std::tuple<int,int,double>>& vals) {
    std::map<int,double> bestMax, bestMin;
    double globalMax=-1e300, globalMin=1e300;
    for (auto [xn,yn,v] : vals) {
        bestMax[xn] = std::max(bestMax[xn], v);
        bestMin[xn] = std::min(bestMin[xn], v);
        globalMax = std::max(globalMax, v);
        globalMin = std::min(globalMin, v);
    }
    double pos=0, neg=0;
    for (auto [xn,v] : bestMax) {
        if (v>0) pos += v; else neg += v;
    }
    const double qMaxArea = surface.q>=0 ? pos : neg;
    const double qMinArea = surface.q>=0 ? neg : pos;
    const double pMaxOrd = surface.p>=0 ? globalMax : globalMin;
    const double pMinOrd = surface.p>=0 ? globalMin : globalMax;
    double factor = surface.reduction;
    if (surface.laneCount == 2) factor = 1.0;
    else if (surface.laneCount == 3) factor = 0.9;
    else if (surface.laneCount >= 4) factor = 0.8;
    return {surface.impact*factor*(surface.q*qMaxArea + surface.p*pMaxOrd),
            surface.impact*factor*(surface.q*qMinArea + surface.p*pMinOrd)};
}

inline std::pair<double,double> envelope3DFromInfluence(const Model3D::Lane3D& lane,
                                                        const std::vector<double>& eta) {
    double pos=0,neg=0,maxv=-1e300,minv=1e300;
    for (std::size_t i=0;i<eta.size();++i) {
        if (eta[i]>0) pos += eta[i]; else neg += eta[i];
        maxv=std::max(maxv,eta[i]); minv=std::min(minv,eta[i]);
    }
    const double qMaxArea = lane.q>=0 ? pos : neg;
    const double qMinArea = lane.q>=0 ? neg : pos;
    const double pMaxOrd = lane.p>=0 ? maxv : minv;
    const double pMinOrd = lane.p>=0 ? minv : maxv;
    return {lane.impact*lane.reduction*(lane.q*qMaxArea + lane.p*pMaxOrd),
            lane.impact*lane.reduction*(lane.q*qMinArea + lane.p*pMinOrd)};
}

struct StageLiveEnvelope3D {
    int stage = 0;
    double maxVal = 0.0, minVal = 0.0;
};

inline std::vector<StageLiveEnvelope3D> compute3DStageLiveEnvelope(
    const Model3D& m, const Model3D::Influence3D& spec) {
    std::vector<StageLiveEnvelope3D> out;
    std::vector<char> active(m.elements.size(),0);
    for (const auto& st : m.stages) {
        for (int id : st.activate)
            for (std::size_t i=0;i<m.elements.size();++i) if(m.elements[i].id==id) active[i]=1;
        for (int id : st.deactivate)
            for (std::size_t i=0;i<m.elements.size();++i) if(m.elements[i].id==id) active[i]=0;
        Model3D mc;
        mc.nodes=m.nodes; mc.constraints=m.constraints;
        mc.creepPhi=m.creepPhi; mc.creepChi=m.creepChi; mc.shrinkageEcs=m.shrinkageEcs;
        for (std::size_t i=0;i<m.elements.size();++i) if(active[i]) mc.elements.push_back(m.elements[i]);
        std::map<int,bool> connected;
        for(const auto& e:mc.elements){connected[e.n1]=true;connected[e.n2]=true;}
        bool hasBeam=false; for(const auto& e:mc.elements) if(e.type=="BEAM"){hasBeam=true;break;}
        const int dofPerNode=hasBeam?6:3;
        for(const auto& [id,x]:m.nodes) if(!connected[id]) for(int dof=1;dof<=dofPerNode;++dof) mc.constraints.push_back({id,dof});
        double maxAll=-1e300,minAll=1e300;
        for(const auto& lane:m.lanes){
            auto eta=compute3DInfluence(mc,spec,lane);
            auto env=envelope3DFromInfluence(lane,eta);
            maxAll=std::max(maxAll,env.first); minAll=std::min(minAll,env.second);
        }
        out.push_back({st.id,maxAll,minAll});
    }
    return out;
}

inline std::vector<Camber3D> compute3DCamber(const Model3D& m,
                                             const std::vector<StageResult3D>& stages) {
    if (stages.empty()) return {};
    bool hasBeam = false;
    for (const auto& e : m.elements) if (e.type == "BEAM") { hasBeam = true; break; }
    const int dofPerNode = hasBeam ? 6 : 3;
    std::vector<Camber3D> out;
    std::map<int,int> idx; int n=0;
    for (const auto& [id,x] : m.nodes) idx[id]=n++;
    for (const auto& [id,x] : m.nodes) {
        for (int dof=1; dof<=dofPerNode; ++dof) {
            out.push_back({id,dof,-stages.back().U[dofPerNode*idx[id]+dof-1]});
        }
    }
    return out;
}

inline Model3D readModel3D(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("无法打开三维输入文件: " + path);
    Model3D m;
    std::map<int, std::string> types;
    std::map<int, double> mats;
    std::map<int, double> areas;
    std::map<int, double> realIy, realIz, realJ, realG;
    std::map<int, double> realIstrn;
    std::map<int,double> realWtop, realWbot;
    int activeType = 0, activeMat = 0, activeReal = 0;
    std::string section;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        if (line[0] == '*') {
            const std::size_t sep = line.find_first_of(", \t");
            const std::string kw = (sep == std::string::npos) ? line : line.substr(0, sep);
            if (kw == "*TYPE" || kw == "*MAT" || kw == "*REAL" || kw == "*EN") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                if (kw == "*TYPE") ss >> activeType;
                else if (kw == "*MAT") ss >> activeMat;
                else if (kw == "*REAL") ss >> activeReal;
                else {
                    int id,n1,n2; if (ss >> id >> n1 >> n2) {
                        Element3D e; e.id=id; e.type=types[activeType]; e.E=mats[activeMat]; e.A=areas[activeReal]; e.n1=n1; e.n2=n2;
                        e.Iy=realIy[activeReal]; e.Iz=realIz[activeReal]; e.J=realJ[activeReal]; e.G=realG[activeReal];
                        e.istrn=realIstrn[activeReal];
                        e.Wtop=realWtop[activeReal]; e.Wbot=realWbot[activeReal];
                        m.elements.push_back(e);
                    }
                }
                continue;
            }
            if (kw == "*CREEP" || kw == "*SHRINKAGE") {
                std::string args = (sep == std::string::npos) ? "" : line.substr(sep + 1);
                std::replace(args.begin(), args.end(), ',', ' ');
                std::istringstream ss(args);
                std::string kv;
                while (ss >> kv) {
                    if (kw == "*CREEP") {
                        if (kv.rfind("PHI=",0)==0) m.creepPhi=std::stod(kv.substr(4));
                        else if (kv.rfind("CHI=",0)==0) m.creepChi=std::stod(kv.substr(4));
                    } else {
                        if (kv.rfind("ECS=",0)==0) m.shrinkageEcs=std::stod(kv.substr(4));
                    }
                }
                continue;
            }
            if (kw == "*INFLUENCE_3D" || kw == "*LANE_3D" || kw == "*LANE_SURFACE_3D" || kw == "*ELEMENT_INFLUENCE_3D") {
                section = line;
                continue;
            }
            section = line;
            continue;
        }
        std::istringstream ss(line);
        if (section == "*NODE") {
            int id; double x,y,z;
            if (ss >> id >> x >> y >> z) m.nodes[id] = Vec3(x,y,z);
        } else if (section == "*ET") {
            int id; std::string t;
            if (ss >> id >> t) types[id] = t;
        } else if (section == "*MP") {
            int id; std::string kv;
            ss >> id;
            double E=0;
            while (ss >> kv) {
                if (kv.rfind("E=",0)==0) E=std::stod(kv.substr(2));
            }
            mats[id]=E;
        } else if (section == "*R") {
            int id; std::string kv;
            ss >> id;
            double A=0, Iy=0, Iz=0, J=0, G=0;
            while (ss >> kv) {
                if (kv.rfind("A=",0)==0) A=std::stod(kv.substr(2));
                else if (kv.rfind("IY=",0)==0) Iy=std::stod(kv.substr(3));
                else if (kv.rfind("IZ=",0)==0) Iz=std::stod(kv.substr(3));
                else if (kv.rfind("J=",0)==0) J=std::stod(kv.substr(2));
                else if (kv.rfind("G=",0)==0) G=std::stod(kv.substr(2));
                else if (kv.rfind("ISTRN=",0)==0) realIstrn[id]=std::stod(kv.substr(6));
                else if (kv.rfind("WT=",0)==0) realWtop[id]=std::stod(kv.substr(3));
                else if (kv.rfind("WB=",0)==0) realWbot[id]=std::stod(kv.substr(3));
            }
            areas[id]=A;
            realIy[id]=Iy; realIz[id]=Iz; realJ[id]=J; realG[id]=G;
        } else if (section == "*BOUNDARY") {
            int node; std::string ds; if (ss >> node >> ds) {
                if (ds.find('-') != std::string::npos) {
                    int a=std::stoi(ds.substr(0,ds.find('-'))); int b=std::stoi(ds.substr(ds.find('-')+1));
                    for (int d=a; d<=b; ++d) m.constraints.push_back({node,d});
                } else m.constraints.push_back({node,std::stoi(ds)});
            }
        } else if (section == "*CLOAD") {
            int node,dof; double val; if (ss >> node >> dof >> val) m.loads.push_back({node,dof,val});
        } else if (section == "*STAGE") {
            int id; if (ss >> id) {
                Model3D::Stage3D st; st.id=id;
                std::string kv;
                while (ss >> kv) {
                    if (kv.rfind("ACTIVE=",0)==0) {
                        std::string list=kv.substr(7);
                        std::stringstream ls(list); std::string tok;
                        while(std::getline(ls,tok,',')) st.activate.push_back(std::stoi(tok));
                    } else if (kv.rfind("DEACTIVE=",0)==0) {
                        std::string list=kv.substr(9);
                        std::stringstream ls(list); std::string tok;
                        while(std::getline(ls,tok,',')) st.deactivate.push_back(std::stoi(tok));
                    }
                }
                m.stages.push_back(st);
            }
        } else if (section == "*STAGE_LOAD") {
            int stageId; std::string kv; ss >> stageId;
            int node=0,dof=0; double val=0;
            while (ss >> kv) {
                if (kv.rfind("NODE=",0)==0) node=std::stoi(kv.substr(5));
                else if (kv.rfind("DOF=",0)==0) dof=std::stoi(kv.substr(4));
                else if (kv.rfind("VALUE=",0)==0) val=std::stod(kv.substr(6));
            }
            if (stageId>=1 && stageId<=static_cast<int>(m.stages.size()))
                m.stages[stageId-1].loads.push_back({node,dof,val});
        } else if (section == "*INFLUENCE_3D") {
            int node,dof; if (ss >> node >> dof) m.influences.push_back({node,dof});
        } else if (section == "*ELEMENT_INFLUENCE_3D") {
            int elem; if (ss >> elem) m.elementInfluences.push_back({elem});
        } else if (section == "*LANE_3D") {
            std::string kv;
            std::string nodes;
            Model3D::Lane3D lane;
            while (ss >> kv) {
                if (kv.rfind("NODES=",0)==0) nodes=kv.substr(6);
                else if (kv.rfind("Q=",0)==0) lane.q=std::stod(kv.substr(2));
                else if (kv.rfind("P=",0)==0) lane.p=std::stod(kv.substr(2));
                else if (kv.rfind("DOF=",0)==0) lane.dof=std::stoi(kv.substr(4));
                else if (kv.rfind("REDUCTION=",0)==0) lane.reduction=std::stod(kv.substr(10));
                else if (kv.rfind("IMPACT=",0)==0) lane.impact=std::stod(kv.substr(7));
            }
            std::stringstream ns(nodes); std::string tok;
            while(std::getline(ns,tok,',')) if(!tok.empty()) lane.nodes.push_back(std::stoi(tok));
            m.lanes.push_back(lane);
        } else if (section == "*LANE_SURFACE_3D") {
            std::string kv, xn, yn;
            while (ss >> kv) {
                if (kv.rfind("X_NODES=",0)==0) xn=kv.substr(8);
                else if (kv.rfind("Y_NODES=",0)==0) yn=kv.substr(8);
                else if (kv.rfind("Q=",0)==0) m.laneSurface.q=std::stod(kv.substr(2));
                else if (kv.rfind("P=",0)==0) m.laneSurface.p=std::stod(kv.substr(2));
                else if (kv.rfind("DOF=",0)==0) m.laneSurface.dof=std::stoi(kv.substr(4));
                else if (kv.rfind("REDUCTION=",0)==0) m.laneSurface.reduction=std::stod(kv.substr(10));
                else if (kv.rfind("LANE_COUNT=",0)==0) m.laneSurface.laneCount=std::stoi(kv.substr(11));
                else if (kv.rfind("IMPACT=",0)==0) m.laneSurface.impact=std::stod(kv.substr(7));
            }
            auto fill=[&](const std::string& s, std::vector<int>& out){
                std::stringstream ls(s); std::string tok;
                while(std::getline(ls,tok,',')) if(!tok.empty()) out.push_back(std::stoi(tok));
            };
            fill(xn,m.laneSurface.xNodes); fill(yn,m.laneSurface.yNodes);
        }
    }
    return m;
}

inline std::vector<double> solveModel3D(const Model3D& m) {
    std::map<int, int> idx;
    int n = 0;
    for (const auto& [id, x] : m.nodes) idx[id] = n++;
    bool hasBeam = false;
    for (const auto& e : m.elements) if (e.type == "BEAM") { hasBeam = true; break; }
    const int dofPerNode = hasBeam ? 6 : 3;
    const int nd = dofPerNode * n;
    std::vector<Eigen::Triplet<double>> trips;
    for (const auto& e : m.elements) {
        const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
        const int i1 = idx.at(e.n1), i2 = idx.at(e.n2);
        if (e.type == "BEAM") {
            const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
            const double L = (p2-p1).norm();
            const auto Kl = beam3DLocalStiffness(Eeff,e.A,e.Iy,e.Iz,e.J,e.G,L);
            const auto T = beam3DTransformation(p1,p2);
            const auto Ke = T.transpose() * Kl * T;
            const int g[12] = {dofPerNode*i1+0,dofPerNode*i1+1,dofPerNode*i1+2,dofPerNode*i1+3,dofPerNode*i1+4,dofPerNode*i1+5,
                               dofPerNode*i2+0,dofPerNode*i2+1,dofPerNode*i2+2,dofPerNode*i2+3,dofPerNode*i2+4,dofPerNode*i2+5};
            for (int a=0;a<12;++a)
                for (int b=0;b<12;++b)
                    if (Ke(a,b)!=0.0) trips.emplace_back(g[a],g[b],Ke(a,b));
        } else {
            const Eigen::Matrix<double, 6, 6> Ke = truss3DStiffness(Eeff, e.A, m.nodes.at(e.n1), m.nodes.at(e.n2));
            const int g[6] = {dofPerNode*i1+0,dofPerNode*i1+1,dofPerNode*i1+2,dofPerNode*i2+0,dofPerNode*i2+1,dofPerNode*i2+2};
            for (int a=0;a<6;++a)
                for (int b=0;b<6;++b)
                if (Ke(a,b) != 0.0) trips.emplace_back(g[a], g[b], Ke(a,b));
        }
    }
    SparseMat3D K(nd, nd);
    K.setFromTriplets(trips.begin(), trips.end());
    Eigen::VectorXd F = Eigen::VectorXd::Zero(nd);
    for (const auto& [node, dof, val] : m.loads) F(dofPerNode*idx.at(node)+dof-1) += val;
    for (const auto& e : m.elements) {
        if (e.istrn == 0.0 && m.shrinkageEcs == 0.0) continue;
        const Vec3 p1 = m.nodes.at(e.n1), p2 = m.nodes.at(e.n2);
        const Vec3 d = p2 - p1;
        const double L = std::max(d.norm(), 1e-12);
        const Vec3 c = d / L;
        const double Eeff = e.E / (1.0 + m.creepChi * m.creepPhi);
        const double P = Eeff * e.A * (e.istrn - m.shrinkageEcs);
        const int i1 = idx.at(e.n1), i2 = idx.at(e.n2);
        for (int k=0; k<3; ++k) {
            F(dofPerNode*i1+k) += P*c(k);
            F(dofPerNode*i2+k) -= P*c(k);
        }
    }
    for (const auto& [node, dof] : m.constraints) {
        const int d = dofPerNode*idx.at(node)+dof-1;
        for (int col = 0; col < K.outerSize(); ++col)
            for (SparseMat3D::InnerIterator it(K,col); it; ++it)
                if (it.row()==d || it.col()==d) it.valueRef() = (it.row()==it.col()) ? 1.0 : 0.0;
        K.coeffRef(d, d) = 1.0;
        F(d) = 0.0;
    }
    Eigen::MatrixXd Kd = K;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(Kd);
    if (!lu.isInvertible()) throw std::runtime_error("三维桁架求解失败");
    Eigen::VectorXd U = lu.solve(F);
    return std::vector<double>(U.data(), U.data()+nd);
}

} // namespace fem
