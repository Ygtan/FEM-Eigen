#include "fem_io.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace fem;

namespace {

const double kW = 900.0, kH = 650.0, kMargin = 60.0;

struct BBox {
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
};

inline void add(BBox& b, double x, double y) {
    b.xmin = std::min(b.xmin, x); b.xmax = std::max(b.xmax, x);
    b.ymin = std::min(b.ymin, y); b.ymax = std::max(b.ymax, y);
}
inline void add(BBox& b, const Eigen::Vector2d& p) { add(b, p(0), p(1)); }

struct Canvas {
    BBox b;
    double scale = 1.0, ox = 0.0, oy = 0.0;
    double X(double x) const { return ox + (x - b.xmin) * scale; }
    double Y(double y) const { return oy + (b.ymax - y) * scale; }
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

// 写一张图：模型（灰色）+ 节点编号 + 可选偏移图（彩色）。
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
    out << std::setprecision(6);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << kW << "\" height=\"" << kH
        << "\" viewBox=\"0 0 " << kW << " " << kH << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out << "<text x=\"20\" y=\"26\" font-size=\"16\" font-family=\"sans-serif\">" << title << "</text>\n";

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
    for (const auto& [id, c] : model.nodes) {
        out << "<circle cx=\"" << cv.X(c(0)) << "\" cy=\"" << cv.Y(c(1))
            << "\" r=\"3.5\" fill=\"black\"/>\n";
        out << "<text x=\"" << (cv.X(c(0)) + 5) << "\" y=\"" << (cv.Y(c(1)) - 5)
            << "\" font-size=\"11\" font-family=\"sans-serif\">" << id << "</text>\n";
    }
    out << "</svg>\n";
}

// 位移图：灰色模型 + 红色放大变形。
void writeDisplacementFigure(const std::string& path, const std::string& title,
                             const Model2D& model, const std::vector<Elem>& elems,
                             const DenseVec& U, const std::map<int, int>& index) {
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
    out << std::setprecision(6);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << kW << "\" height=\"" << kH
        << "\" viewBox=\"0 0 " << kW << " " << kH << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out << "<text x=\"20\" y=\"26\" font-size=\"16\" font-family=\"sans-serif\">" << title << "</text>\n";

    for (const auto& e : elems) {
        out << "<line x1=\"" << cv.X(e.p1(0)) << "\" y1=\"" << cv.Y(e.p1(1))
            << "\" x2=\"" << cv.X(e.p2(0)) << "\" y2=\"" << cv.Y(e.p2(1))
            << "\" stroke=\"#999\" stroke-width=\"2\"/>\n";
    }
    for (const auto& e : elems) {
        const int i1 = index.at(e.n1), i2 = index.at(e.n2);
        const Eigen::Vector2d d1(e.p1(0) + U(3 * i1) * exag, e.p1(1) + U(3 * i1 + 1) * exag);
        const Eigen::Vector2d d2(e.p2(0) + U(3 * i2) * exag, e.p2(1) + U(3 * i2 + 1) * exag);
        out << "<line x1=\"" << cv.X(d1(0)) << "\" y1=\"" << cv.Y(d1(1))
            << "\" x2=\"" << cv.X(d2(0)) << "\" y2=\"" << cv.Y(d2(1))
            << "\" stroke=\"#e60000\" stroke-width=\"2.5\"/>\n";
    }
    for (const auto& [id, c] : model.nodes) {
        out << "<circle cx=\"" << cv.X(c(0)) << "\" cy=\"" << cv.Y(c(1))
            << "\" r=\"3.5\" fill=\"black\"/>\n";
        out << "<text x=\"" << (cv.X(c(0)) + 5) << "\" y=\"" << (cv.Y(c(1)) - 5)
            << "\" font-size=\"11\" font-family=\"sans-serif\">" << id << "</text>\n";
    }
    out << "</svg>\n";
}

} // namespace

static int runViz(int argc, char** argv) {
    setupConsoleUTF8();
    std::string path = "mixed_truss_beam.dat";
    if (argc > 1) path = argv[1];

    std::string base = path;
    const std::size_t p = base.rfind(".dat");
    if (p != std::string::npos && p + 4 == base.size()) base = base.substr(0, p);

    const Model2D model = readModel2D(path);
    const ModelResult r = solveModel2D(model);
    const std::vector<ElementForce> forces = elementForces2D(model, r.U);
    std::vector<ElementStress> stresses;
    bool hasStress = true;
    try {
        stresses = elementStress2D(model, forces);
    } catch (const std::runtime_error&) {
        hasStress = false;
    }

    std::map<int, int> index;
    for (std::size_t i = 0; i < r.nodeIds.size(); ++i) index[r.nodeIds[i]] = static_cast<int>(i);

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

    auto findF = [&](int id) -> const ElementForce& {
        for (const auto& f : forces) if (f.id == id) return f;
        throw std::runtime_error("internal: element not found");
    };
    auto findS = [&](int id) -> const ElementStress& {
        for (const auto& s : stresses) if (s.id == id) return s;
        throw std::runtime_error("internal: element not found");
    };
    auto findElem = [&](int id) -> const Model2D::Element& {
        for (const auto& e : model.elements) if (e.id == id) return e;
        throw std::runtime_error("internal: element not found");
    };

    const std::string ttl = model.title.empty() ? path : model.title;

    // 模型图
    writeDiagramFigure(base + "_model.svg", ttl + " - 模型图", model, elems, {}, "#000000");

    // 位移图
    writeDisplacementFigure(base + "_disp.svg", ttl + " - 位移图（放大）", model, elems, r.U, index);

    // 弯矩图
    {
        std::vector<Seg> segs;
        const int nSub = 4;
        for (const auto& e : elems) {
            if (e.type != "BEAM") continue;
            const auto& f = findF(e.id);
            const auto& me = findElem(e.id);
            const double qd = localTransverseUDL(model, me, e.p1, e.p2);  // 向下为正
            const Eigen::Vector2d dir = e.p2 - e.p1;
            const double L = dir.norm();
            // 弯矩 = 端弯矩线性插值 + 分布荷载的“简支”抛物线项 qd·x·(L-x)/2
            for (int k = 0; k < nSub; ++k) {
                const double x0 = L * k / nSub, x1 = L * (k + 1) / nSub;
                const double m0 = f.M1 + (f.M2 - f.M1) * (x0 / L) + qd * x0 * (L - x0) / 2.0;
                const double m1 = f.M1 + (f.M2 - f.M1) * (x1 / L) + qd * x1 * (L - x1) / 2.0;
                segs.push_back({e.id, e.p1 + dir * (x0 / L), e.p1 + dir * (x1 / L), m0, m1});
            }
        }
        writeDiagramFigure(base + "_moment.svg", ttl + " - 弯矩图", model, elems, segs, "#ff8c00");
    }

    // 剪力图
    {
        std::vector<Seg> segs;
        for (const auto& e : elems) {
            if (e.type != "BEAM") continue;
            const auto& f = findF(e.id);
            segs.push_back({e.id, e.p1, e.p2, f.V1, f.V2});
        }
        writeDiagramFigure(base + "_shear.svg", ttl + " - 剪力图", model, elems, segs, "#8a2be2");
    }

    // 轴力图（梁 + 桁架）
    {
        std::vector<Seg> segs;
        for (const auto& e : elems) {
            if (isRigid(e.type)) continue;
            const auto& f = findF(e.id);
            segs.push_back({e.id, e.p1, e.p2, f.N1, f.N2});
        }
        writeDiagramFigure(base + "_axial.svg", ttl + " - 轴力图", model, elems, segs, "#2e8b57");
    }

    if (hasStress) {
        std::vector<Seg> top, bot;
        for (const auto& e : elems) {
            if (e.type != "BEAM") continue;
            const auto& s = findS(e.id);
            top.push_back({e.id, e.p1, e.p2, s.stop1, s.stop2});
            bot.push_back({e.id, e.p1, e.p2, s.sbot1, s.sbot2});
        }
        writeDiagramFigure(base + "_stress_top.svg", ttl + " - 上翼缘应力图", model, elems, top, "#1e90ff");
        writeDiagramFigure(base + "_stress_bottom.svg", ttl + " - 下翼缘应力图", model, elems, bot, "#228b22");
    }

    std::cout << "已生成图形（SVG）:\n";
    std::cout << "  " << base << "_model.svg\n";
    std::cout << "  " << base << "_disp.svg\n";
    std::cout << "  " << base << "_moment.svg\n";
    std::cout << "  " << base << "_shear.svg\n";
    std::cout << "  " << base << "_axial.svg\n";
    if (hasStress) {
        std::cout << "  " << base << "_stress_top.svg\n";
        std::cout << "  " << base << "_stress_bottom.svg\n";
    } else {
        std::cout << "  （梁缺少 WT/WB，跳过上下翼缘应力图）\n";
    }
    return 0;
}

int main(int argc, char** argv) {
    try {
        return runViz(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
