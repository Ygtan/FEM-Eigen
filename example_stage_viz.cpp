#include "fem_2d.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace fem;

namespace {

const double kW = 900.0, kH = 650.0, kMargin = 60.0;

struct BBox {
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
};

void add(BBox& b, double x, double y) {
    b.xmin = std::min(b.xmin, x); b.xmax = std::max(b.xmax, x);
    b.ymin = std::min(b.ymin, y); b.ymax = std::max(b.ymax, y);
}
void add(BBox& b, const Eigen::Vector2d& p) { add(b, p(0), p(1)); }

struct Canvas {
    BBox b;
    double scale = 1.0, ox = 0.0, oy = 0.0;
    double X(double x) const { return ox + (x - b.xmin) * scale; }
    double Y(double y) const { return oy + (b.ymax - y) * scale; }
};

Eigen::Vector2d perp(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    const Eigen::Vector2d d = b - a;
    const double L = std::max(d.norm(), 1e-12);
    return Eigen::Vector2d(-d(1) / L, d(0) / L);
}

struct Elem {
    int id = 0;
    std::string type;
    int n1 = 0, n2 = 0;
    Eigen::Vector2d p1, p2;
};

struct Seg {
    int id = 0;
    Eigen::Vector2d p1, p2;
    double v1 = 0.0, v2 = 0.0;
};

Canvas makeCanvas(const BBox& b) {
    const double rx = std::max(b.xmax - b.xmin, 1e-9);
    const double ry = std::max(b.ymax - b.ymin, 1e-9);
    const double s = std::min((kW - 2 * kMargin) / rx, (kH - 2 * kMargin) / ry);
    Canvas c;
    c.b = b;
    c.scale = s;
    c.ox = (kW - rx * s) / 2.0;
    c.oy = (kH - ry * s) / 2.0;
    return c;
}

void writeModel(const std::string& path, const std::string& title,
                const Model2D& model, const std::vector<Eigen::Vector2d>& pts) {
    BBox box;
    for (const auto& [id, c] : model.nodes) add(box, c);
    const Canvas cv = makeCanvas(box);
    std::ofstream out(path);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << kW << "\" height=\"" << kH
        << "\" viewBox=\"0 0 " << kW << " " << kH << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out << "<text x=\"20\" y=\"26\" font-size=\"16\">" << title << "</text>\n";
    for (const auto& e : model.elements) {
        const auto& a = model.nodes.at(e.conn[0]);
        const auto& b = model.nodes.at(e.conn[1]);
        out << "<line x1=\"" << cv.X(a(0)) << "\" y1=\"" << cv.Y(a(1))
            << "\" x2=\"" << cv.X(b(0)) << "\" y2=\"" << cv.Y(b(1))
            << "\" stroke=\"#999\" stroke-width=\"2\"/>\n";
    }
    for (const auto& [id, c] : model.nodes) {
        out << "<circle cx=\"" << cv.X(c(0)) << "\" cy=\"" << cv.Y(c(1)) << "\" r=\"3.5\" fill=\"black\"/>\n";
        out << "<text x=\"" << (cv.X(c(0)) + 5) << "\" y=\"" << (cv.Y(c(1)) - 5) << "\">" << id << "</text>\n";
    }
    out << "</svg>\n";
}

void writeDisplacement(const std::string& path, const std::string& title,
                       const Model2D& model, const DenseVec& U,
                       const std::map<int, int>& index) {
    BBox box;
    for (const auto& [id, c] : model.nodes) add(box, c);
    double maxDisp = 1e-12;
    for (const auto& [id, c] : model.nodes) {
        const int i = index.at(id);
        maxDisp = std::max(maxDisp, std::sqrt(U(3 * i) * U(3 * i) + U(3 * i + 1) * U(3 * i + 1)));
    }
    const double span = std::max(box.xmax - box.xmin, box.ymax - box.ymin);
    const double exag = 0.12 * span / maxDisp;
    for (const auto& [id, c] : model.nodes) {
        const int i = index.at(id);
        add(box, Eigen::Vector2d(c(0) + U(3 * i) * exag, c(1) + U(3 * i + 1) * exag));
    }
    const Canvas cv = makeCanvas(box);
    std::ofstream out(path);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << kW << "\" height=\"" << kH
        << "\" viewBox=\"0 0 " << kW << " " << kH << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out << "<text x=\"20\" y=\"26\" font-size=\"16\">" << title << "</text>\n";
    for (const auto& e : model.elements) {
        const auto& a = model.nodes.at(e.conn[0]);
        const auto& b = model.nodes.at(e.conn[1]);
        out << "<line x1=\"" << cv.X(a(0)) << "\" y1=\"" << cv.Y(a(1))
            << "\" x2=\"" << cv.X(b(0)) << "\" y2=\"" << cv.Y(b(1))
            << "\" stroke=\"#999\" stroke-width=\"2\"/>\n";
    }
    for (const auto& e : model.elements) {
        const int i1 = index.at(e.conn[0]);
        const int i2 = index.at(e.conn[1]);
        const auto a = model.nodes.at(e.conn[0]);
        const auto b = model.nodes.at(e.conn[1]);
        out << "<line x1=\"" << cv.X(a(0) + U(3 * i1) * exag) << "\" y1=\"" << cv.Y(a(1) + U(3 * i1 + 1) * exag)
            << "\" x2=\"" << cv.X(b(0) + U(3 * i2) * exag) << "\" y2=\"" << cv.Y(b(1) + U(3 * i2 + 1) * exag)
            << "\" stroke=\"#e60000\" stroke-width=\"2.5\"/>\n";
    }
    out << "</svg>\n";
}

void writeDiagramFigure(const std::string& path, const std::string& title,
                        const Model2D& model, const std::vector<Elem>& elems,
                        const std::vector<Seg>& segs, const std::string& color) {
    BBox box;
    for (const auto& [id, c] : model.nodes) add(box, c);
    double maxAbs = 1e-12;
    for (const auto& s : segs) maxAbs = std::max(maxAbs, std::max(std::abs(s.v1), std::abs(s.v2)));
    const double span = std::max(box.xmax - box.xmin, box.ymax - box.ymin);
    const double vscale = 0.22 * span / maxAbs;
    for (const auto& s : segs) {
        const Eigen::Vector2d n = perp(s.p1, s.p2);
        add(box, s.p1 + n * (s.v1 * vscale));
        add(box, s.p2 + n * (s.v2 * vscale));
    }
    const Canvas cv = makeCanvas(box);
    std::ofstream out(path);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << kW << "\" height=\"" << kH
        << "\" viewBox=\"0 0 " << kW << " " << kH << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out << "<text x=\"20\" y=\"26\" font-size=\"16\">" << title << "</text>\n";
    for (const auto& e : elems) {
        out << "<line x1=\"" << cv.X(e.p1(0)) << "\" y1=\"" << cv.Y(e.p1(1))
            << "\" x2=\"" << cv.X(e.p2(0)) << "\" y2=\"" << cv.Y(e.p2(1))
            << "\" stroke=\"#999\" stroke-width=\"2\"/>\n";
    }
    for (const auto& s : segs) {
        const Eigen::Vector2d n = perp(s.p1, s.p2);
        const Eigen::Vector2d o1 = s.p1 + n * (s.v1 * vscale);
        const Eigen::Vector2d o2 = s.p2 + n * (s.v2 * vscale);
        out << "<polygon points=\"" << cv.X(s.p1(0)) << "," << cv.Y(s.p1(1)) << " "
            << cv.X(o1(0)) << "," << cv.Y(o1(1)) << " "
            << cv.X(o2(0)) << "," << cv.Y(o2(1)) << " "
            << cv.X(s.p2(0)) << "," << cv.Y(s.p2(1))
            << "\" fill=\"" << color << "\" fill-opacity=\"0.3\" stroke=\"" << color
            << "\" stroke-width=\"1.5\"/>\n";
    }
    out << "</svg>\n";
}

} // namespace

int main(int argc, char** argv) {
    setupConsoleUTF8();
    try {
        std::string path = "verify_stage.dat";
        if (argc > 1) path = argv[1];
        const Model2D model = readModel2D(path);
        const auto stages = solveModel2DStages(model);

        std::map<int, int> index;
        int n = 0;
        for (const auto& [id, c] : model.nodes) index[id] = n++;

        std::string base = path;
        const auto p = base.rfind(".dat");
        if (p != std::string::npos && p + 4 == base.size()) base = base.substr(0, p);

        for (const auto& sr : stages) {
            const std::string suffix = "_stage" + std::to_string(sr.stage);
            writeModel(base + suffix + "_model.svg", "Stage " + std::to_string(sr.stage), model, {});
            writeDisplacement(base + suffix + "_disp.svg", "Stage " + std::to_string(sr.stage) + " displacement", model, sr.U, index);

            std::vector<Elem> elems;
            for (const auto& e : model.elements) {
                Elem el;
                el.id = e.id;
                el.type = e.type;
                el.n1 = e.conn[0];
                el.n2 = e.conn[1];
                el.p1 = model.nodes.at(e.conn[0]);
                el.p2 = model.nodes.at(e.conn[1]);
                elems.push_back(el);
            }
            std::set<int> activeIds;
            for (const auto& f : sr.forces) activeIds.insert(f.id);
            auto findF = [&](int id) -> const ElementForce& {
                for (const auto& f : sr.forces) if (f.id == id) return f;
                throw std::runtime_error("internal: element not found");
            };
            auto findS = [&](int id) -> const ElementStress& {
                for (const auto& s : sr.stresses) if (s.id == id) return s;
                throw std::runtime_error("internal: element not found");
            };
            auto findElem = [&](int id) -> const Model2D::Element& {
                for (const auto& e : model.elements) if (e.id == id) return e;
                throw std::runtime_error("internal: element not found");
            };

            {
                std::vector<Seg> segs;
                const int nSub = 4;
                for (const auto& e : elems) {
                    if (e.type != "BEAM" || !activeIds.count(e.id)) continue;
                    const auto& f = findF(e.id);
                    const auto& me = findElem(e.id);
                    const double qd = localTransverseUDL(model, me, e.p1, e.p2);
                    const Eigen::Vector2d dir = e.p2 - e.p1;
                    const double L = dir.norm();
                    for (int k = 0; k < nSub; ++k) {
                        const double x0 = L * k / nSub, x1 = L * (k + 1) / nSub;
                        const double m0 = f.M1 + (f.M2 - f.M1) * (x0 / L) + qd * x0 * (L - x0) / 2.0;
                        const double m1 = f.M1 + (f.M2 - f.M1) * (x1 / L) + qd * x1 * (L - x1) / 2.0;
                        segs.push_back({e.id, e.p1 + dir * (x0 / L), e.p1 + dir * (x1 / L), m0, m1});
                    }
                }
                writeDiagramFigure(base + suffix + "_moment.svg", "Stage " + std::to_string(sr.stage) + " moment", model, elems, segs, "#ff8c00");
            }
            {
                std::vector<Seg> segs;
                for (const auto& e : elems) {
                    if (e.type != "BEAM" || !activeIds.count(e.id)) continue;
                    const auto& f = findF(e.id);
                    segs.push_back({e.id, e.p1, e.p2, f.V1, f.V2});
                }
                writeDiagramFigure(base + suffix + "_shear.svg", "Stage " + std::to_string(sr.stage) + " shear", model, elems, segs, "#8a2be2");
            }
            {
                std::vector<Seg> segs;
                for (const auto& e : elems) {
                    if (isRigid(e.type) || !activeIds.count(e.id)) continue;
                    const auto& f = findF(e.id);
                    segs.push_back({e.id, e.p1, e.p2, f.N1, f.N2});
                }
                writeDiagramFigure(base + suffix + "_axial.svg", "Stage " + std::to_string(sr.stage) + " axial", model, elems, segs, "#2e8b57");
            }
            if (!sr.stresses.empty()) {
                std::vector<Seg> top, bot;
                for (const auto& e : elems) {
                    if (e.type != "BEAM" || !activeIds.count(e.id)) continue;
                    const auto& s = findS(e.id);
                    top.push_back({e.id, e.p1, e.p2, s.stop1, s.stop2});
                    bot.push_back({e.id, e.p1, e.p2, s.sbot1, s.sbot2});
                }
                writeDiagramFigure(base + suffix + "_stress_top.svg", "Stage " + std::to_string(sr.stage) + " stress top", model, elems, top, "#1e90ff");
                writeDiagramFigure(base + suffix + "_stress_bottom.svg", "Stage " + std::to_string(sr.stage) + " stress bottom", model, elems, bot, "#228b22");
            }
            std::cout << base << suffix << "_model.svg\n";
            std::cout << base << suffix << "_disp.svg\n";
            std::cout << base << suffix << "_moment.svg\n";
            std::cout << base << suffix << "_shear.svg\n";
            std::cout << base << suffix << "_axial.svg\n";
            if (!sr.stresses.empty()) {
                std::cout << base << suffix << "_stress_top.svg\n";
                std::cout << base << suffix << "_stress_bottom.svg\n";
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
