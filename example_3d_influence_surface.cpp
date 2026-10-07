#include "fem_3d.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace fem;

namespace {
struct BBox { double xmin=1e300,xmax=-1e300,ymin=1e300,ymax=-1e300; };
void add(BBox& b,double x,double y){ b.xmin=std::min(b.xmin,x); b.xmax=std::max(b.xmax,x); b.ymin=std::min(b.ymin,y); b.ymax=std::max(b.ymax,y); }
struct Canvas { BBox b; double scale=1,ox=0,oy=0; double X(double x)const{return ox+(x-b.xmin)*scale;} double Y(double y)const{return oy+(b.ymax-y)*scale;} };
Canvas makeCanvas(const BBox& b){ double rx=std::max(b.xmax-b.xmin,1e-9), ry=std::max(b.ymax-b.ymin,1e-9); double s=std::min(780.0/rx,580.0/ry); Canvas c; c.b=b; c.scale=s; c.ox=40+(780-rx*s)/2.0; c.oy=40+(580-ry*s)/2.0; return c; }
void writeSurfaceSVG(const std::string& path,const std::vector<std::tuple<int,int,double>>& vals){
    BBox box; double maxAbs=1e-12;
    for(auto [x,y,v]:vals){ add(box,x,y); maxAbs=std::max(maxAbs,std::abs(v)); }
    Canvas cv=makeCanvas(box);
    std::ofstream out(path);
    out<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"650\">\n";
    out<<"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out<<"<text x=\"20\" y=\"26\" font-size=\"16\">3D influence surface</text>\n";
    for(auto [x,y,v]:vals){
        int r=3+std::min(12,static_cast<int>(std::abs(v)/maxAbs*12));
        std::string color=v>=0?"#ff0000":"#0000ff";
        out<<"<circle cx=\""<<cv.X(x)<<"\" cy=\""<<cv.Y(y)<<"\" r=\""<<r<<"\" fill=\""<<color<<"\" fill-opacity=\"0.55\"/>\n";
        out<<"<text x=\""<<(cv.X(x)+6)<<"\" y=\""<<(cv.Y(y)-6)<<"\" font-size=\"10\">"<<v<<"</text>\n";
    }
    out<<"</svg>\n";
}
} // namespace

int main(int argc, char** argv) {
    std::string path = "verify_3d_influence_surface.dat";
    if (argc > 1) path = argv[1];
    Model3D m = readModel3D(path);
    for (const auto& spec : m.influences) {
        auto vals = compute3DInfluenceSurface(m, spec, m.laneSurface);
        auto env = envelope3DFromSurface(m.laneSurface, vals);
        std::cout << "N" << spec.node << " dof" << spec.dof
                  << " max=" << env.first << " min=" << env.second << "\n";
        std::ofstream csv(path.substr(0,path.rfind(".dat"))+"_surface.csv");
        csv << "x_node,y_node,eta\n";
        for (auto [xn,yn,v] : vals) {
            std::cout << "x=" << xn << " y=" << yn << " eta=" << v << "\n";
            csv << xn << "," << yn << "," << v << "\n";
        }
        writeSurfaceSVG(path.substr(0,path.rfind(".dat"))+"_surface.svg", vals);
    }
    return 0;
}
