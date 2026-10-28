#include "fem_2d.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <cmath>

using namespace fem;

struct EnvelopeRow {
    int id = 0;
    std::string type;
    double Nmax = -1e300, Nmin = 1e300;
    double Vmax = -1e300, Vmin = 1e300;
    double Mmax = -1e300, Mmin = 1e300;
    double StopMax = -1e300, StopMin = 1e300;
    double SbotMax = -1e300, SbotMin = 1e300;
};

namespace {

const double kW = 900.0, kH = 650.0, kMargin = 60.0;
struct BBox { double xmin=1e300,xmax=-1e300,ymin=1e300,ymax=-1e300; };
void add(BBox& b, double x, double y){ b.xmin=std::min(b.xmin,x); b.xmax=std::max(b.xmax,x); b.ymin=std::min(b.ymin,y); b.ymax=std::max(b.ymax,y); }
void add(BBox& b, const Eigen::Vector2d& p){ add(b,p(0),p(1)); }
struct Canvas { BBox b; double scale=1,ox=0,oy=0; double X(double x)const{return ox+(x-b.xmin)*scale;} double Y(double y)const{return oy+(b.ymax-y)*scale;} };
Canvas makeCanvas(const BBox& b){ double rx=std::max(b.xmax-b.xmin,1e-9), ry=std::max(b.ymax-b.ymin,1e-9); double s=std::min((kW-2*kMargin)/rx,(kH-2*kMargin)/ry); Canvas c; c.b=b; c.scale=s; c.ox=(kW-rx*s)/2.0; c.oy=(kH-ry*s)/2.0; return c; }
Eigen::Vector2d perp(const Eigen::Vector2d& a,const Eigen::Vector2d& b){ Eigen::Vector2d d=b-a; double L=std::max(d.norm(),1e-12); return Eigen::Vector2d(-d(1)/L,d(0)/L); }
struct Elem{ int id; Eigen::Vector2d p1,p2; };
struct Seg{ int id; Eigen::Vector2d p1,p2; double v1,v2; };

void writeEnvelopeFigure(const std::string& path, const std::string& title,
                         const Model2D& model, const std::vector<Elem>& elems,
                         const std::vector<Seg>& segsMax, const std::vector<Seg>& segsMin,
                         const std::string& colorMax, const std::string& colorMin) {
    BBox box;
    for (const auto& [id,c] : model.nodes) add(box,c);
    double maxAbs=1e-12;
    for(const auto& s:segsMax) maxAbs=std::max(maxAbs,std::max(std::abs(s.v1),std::abs(s.v2)));
    for(const auto& s:segsMin) maxAbs=std::max(maxAbs,std::max(std::abs(s.v1),std::abs(s.v2)));
    double span=std::max(box.xmax-box.xmin, box.ymax-box.ymin);
    double vscale=0.22*span/maxAbs;
    for(const auto& s:segsMax){ auto n=perp(s.p1,s.p2); add(box,s.p1+n*(s.v1*vscale)); add(box,s.p2+n*(s.v2*vscale)); }
    for(const auto& s:segsMin){ auto n=perp(s.p1,s.p2); add(box,s.p1+n*(s.v1*vscale)); add(box,s.p2+n*(s.v2*vscale)); }
    Canvas cv=makeCanvas(box);
    std::ofstream out(path);
    out<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<kW<<"\" height=\""<<kH<<"\" viewBox=\"0 0 "<<kW<<" "<<kH<<"\">\n";
    out<<"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out<<"<text x=\"20\" y=\"26\" font-size=\"16\">"<<title<<"</text>\n";
    for(const auto& e:elems){ out<<"<line x1=\""<<cv.X(e.p1(0))<<"\" y1=\""<<cv.Y(e.p1(1))<<"\" x2=\""<<cv.X(e.p2(0))<<"\" y2=\""<<cv.Y(e.p2(1))<<"\" stroke=\"#999\" stroke-width=\"2\"/>\n"; }
    auto draw=[&](const std::vector<Seg>& segs,const std::string& color){
        for(const auto& s:segs){ auto n=perp(s.p1,s.p2); auto o1=s.p1+n*(s.v1*vscale); auto o2=s.p2+n*(s.v2*vscale); out<<"<polygon points=\""<<cv.X(s.p1(0))<<","<<cv.Y(s.p1(1))<<" "<<cv.X(o1(0))<<","<<cv.Y(o1(1))<<" "<<cv.X(o2(0))<<","<<cv.Y(o2(1))<<" "<<cv.X(s.p2(0))<<","<<cv.Y(s.p2(1))<<"\" fill=\""<<color<<"\" fill-opacity=\"0.3\" stroke=\""<<color<<"\" stroke-width=\"1.5\"/>\n"; }
    };
    draw(segsMax,colorMax); draw(segsMin,colorMin);
    out<<"</svg>\n";
}

} // namespace

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "verify_stage.dat";
        if (argc > 1) path = argv[1];
        const Model2D model = readModel2D(path);
        const auto stages = solveModel2DStages(model);

        std::map<int, EnvelopeRow> rows;
        for (const auto& e : model.elements) {
            if (isRigid(e.type)) continue;
            EnvelopeRow r;
            r.id = e.id;
            r.type = e.type;
            rows[e.id] = r;
        }

        for (const auto& sr : stages) {
            for (const auto& f : sr.forces) {
                if (!rows.count(f.id)) continue;
                auto& r = rows[f.id];
                r.Nmax = std::max(r.Nmax, f.N1);
                r.Nmin = std::min(r.Nmin, f.N1);
                if (f.type == "BEAM") {
                    r.Vmax = std::max(r.Vmax, f.V1);
                    r.Vmin = std::min(r.Vmin, f.V1);
                    r.Mmax = std::max(r.Mmax, f.M1);
                    r.Mmin = std::min(r.Mmin, f.M1);
                }
            }
            for (const auto& s : sr.stresses) {
                if (!rows.count(s.id) || s.type != "BEAM") continue;
                auto& r = rows[s.id];
                r.StopMax = std::max(r.StopMax, s.stop1);
                r.StopMin = std::min(r.StopMin, s.stop1);
                r.SbotMax = std::max(r.SbotMax, s.sbot1);
                r.SbotMin = std::min(r.SbotMin, s.sbot1);
            }
        }

        std::string base = path;
        const auto p = base.rfind(".dat");
        if (p != std::string::npos && p + 4 == base.size()) base = base.substr(0, p);
        const std::string outPath = base + "_stage_envelope.txt";
        std::ofstream out(outPath);
        out << std::setprecision(12);
        out << "elem\ttype\tNmax\tNmin\tVmax\tVmin\tMmax\tMmin\tStopMax\tStopMin\tSbotMax\tSbotMin\n";
        for (const auto& [id, r] : rows) {
            out << r.id << '\t' << r.type << '\t'
                << r.Nmax << '\t' << r.Nmin << '\t'
                << r.Vmax << '\t' << r.Vmin << '\t'
                << r.Mmax << '\t' << r.Mmin << '\t'
                << r.StopMax << '\t' << r.StopMin << '\t'
                << r.SbotMax << '\t' << r.SbotMin << '\n';
        }
        std::cout << "已写入施工阶段包络: " << outPath << "\n";

        std::vector<Elem> elems;
        for (const auto& e : model.elements) {
            Elem el;
            el.id = e.id;
            el.p1 = model.nodes.at(e.conn[0]);
            el.p2 = model.nodes.at(e.conn[1]);
            elems.push_back(el);
        }
        std::vector<Seg> Mmax, Mmin, TopMax, TopMin, BotMax, BotMin;
        std::vector<Seg> Vmax, Vmin, NmaxSeg, NminSeg;
        for (const auto& [id, r] : rows) {
            auto it = std::find_if(elems.begin(), elems.end(), [&](const Elem& e){ return e.id == id; });
            if (it == elems.end()) continue;
            NmaxSeg.push_back({id, it->p1, it->p2, r.Nmax, r.Nmax});
            NminSeg.push_back({id, it->p1, it->p2, r.Nmin, r.Nmin});
            if (r.type != "BEAM") continue;
            Mmax.push_back({id, it->p1, it->p2, r.Mmax, r.Mmax});
            Mmin.push_back({id, it->p1, it->p2, r.Mmin, r.Mmin});
            Vmax.push_back({id, it->p1, it->p2, r.Vmax, r.Vmax});
            Vmin.push_back({id, it->p1, it->p2, r.Vmin, r.Vmin});
            TopMax.push_back({id, it->p1, it->p2, r.StopMax, r.StopMax});
            TopMin.push_back({id, it->p1, it->p2, r.StopMin, r.StopMin});
            BotMax.push_back({id, it->p1, it->p2, r.SbotMax, r.SbotMax});
            BotMin.push_back({id, it->p1, it->p2, r.SbotMin, r.SbotMin});
        }
        writeEnvelopeFigure(base + "_stage_envelope_moment.svg", "Moment envelope", model, elems, Mmax, Mmin, "#ff8c00", "#1e90ff");
        writeEnvelopeFigure(base + "_stage_envelope_shear.svg", "Shear envelope", model, elems, Vmax, Vmin, "#8a2be2", "#1e90ff");
        writeEnvelopeFigure(base + "_stage_envelope_axial.svg", "Axial force envelope", model, elems, NmaxSeg, NminSeg, "#2e8b57", "#1e90ff");
        writeEnvelopeFigure(base + "_stage_envelope_stress_top.svg", "Top stress envelope", model, elems, TopMax, TopMin, "#e60000", "#1e90ff");
        writeEnvelopeFigure(base + "_stage_envelope_stress_bottom.svg", "Bottom stress envelope", model, elems, BotMax, BotMin, "#e60000", "#1e90ff");
        std::cout << "已写入施工阶段包络图\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
