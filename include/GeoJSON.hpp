#ifndef GEOJSON_HPP
#define GEOJSON_HPP
#include "stdafx.h"
#include "json/json.h"
#include <fstream>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Segment_2.h>
namespace GeoJSON
{
    struct Point {
        double x, y;
        Point() : x(0), y(0) {}
        Point(double a, double b) : x(a), y(b) {}
    };

    struct Line {
        Point a, b;
        Line(Point p0, Point p1) : a(p0), b(p1) {}
    };

    struct Block {
        std::vector<std::vector<Point>> coords;
        std::string category;
        std::string sw;
        std::string privacy;
    };

    struct parseout {
        std::vector<Block> data;
        double minx = DBL_MAX, maxx = DBL_MIN, miny = DBL_MAX, maxy = DBL_MIN;
    };

    struct Fe {
        Point p;
        std::string type;
        int rspace;
    };

    template <typename K>
    struct valout {
        std::string type;
        std::vector<std::vector<typename K::Line_2>> oS;
    };

    void parse_geojson(const std::string& datastr, parseout& out);
    Json::Value parse_propsjson(const std::string& datastr);
    std::string out2str(std::vector<std::pair<Point, Point>> out);
    std::string out2str(const std::vector<Block>& ucs_blocks, const std::vector<Point>& corr_ucs);
    std::string out2str(const std::vector<Block>& rect_blocks, const std::vector<Block>& centerline_blocks);
    std::string out2str(const Block& block);

	template <typename K, typename Segment_2 = CGAL::Segment_2<K>, typename Point_2 = CGAL::Point_2<K>>
    std::string segments_to_geojson(const std::vector<Segment_2>& segs) {
        std::vector<std::pair<Point, Point>> res;
        for (auto& seg : segs) {
            Point_2 src = seg.source(), dest = seg.target();
            Point p0(CGAL::to_double(src.x()), CGAL::to_double(src.y()));
            Point p1(CGAL::to_double(dest.x()), CGAL::to_double(dest.y()));
            res.push_back(std::make_pair(p0, p1));
        }
        return out2str(res);
    }

    template <typename K, typename Polygon_2 = CGAL::Polygon_2<K>, typename Point_2 = CGAL::Point_2<K>>
    Polygon_2 convert_poly(const std::vector<GeoJSON::Point>& points) {
        std::vector<Point_2> pts;
        Point_2 tmp;
        for (auto& p : points) {
            Point_2 pt(p.x, p.y);
            if (!pts.empty() && tmp == pt) continue;
            pts.push_back(pt);
            tmp = pt;
        }
        if (pts.front() == pts.back()) pts.pop_back();
        return Polygon_2(pts.begin(), pts.end());
    }
    template <typename K, typename Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>, typename Polygon_2 = CGAL::Polygon_2<K>>
    Polygon_with_holes_2 geojson_to_Pwh(const std::string& geojson) {
        parseout res;
        parse_geojson(geojson, res);
        auto& data = res.data[0];
        Polygon_2 poly = convert_poly<K>(data.coords[0]);  // the first group geojson data, which is the outer boundary
        if (poly.is_clockwise_oriented()) poly.reverse_orientation();
        std::vector<Polygon_2> holes;
        for (int i = 1; i < data.coords.size(); ++i) {
            Polygon_2 hole = convert_poly<K>(data.coords[i]);
            if (hole.is_counterclockwise_oriented()) hole.reverse_orientation();
            holes.push_back(hole);
        }
        if (holes.empty()) return Polygon_with_holes_2(poly);
        else return Polygon_with_holes_2(poly, holes.begin(), holes.end());
    }
    template <typename K, typename Polygon_2 = CGAL::Polygon_2<K>, typename Vector_2 = CGAL::Vector_2<K>>
    std::vector<Point> convert_points(const Polygon_2& poly, const Vector_2& offset) {
        std::vector<Point> ans;
        auto it = poly.vertices_begin();
        if (it == poly.vertices_end()) {
            throw("poly is empty");
        }
        for (; it != poly.vertices_end(); ++it)
            ans.emplace_back(CGAL::to_double(it->x() + offset.x()), CGAL::to_double(it->y() + offset.y()));
        ans.emplace_back(CGAL::to_double(poly.vertices_begin()->x() + offset.x()), CGAL::to_double(poly.vertices_begin()->y() + offset.y()));
        return ans;
    }
    template <typename K, typename Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>, typename Vector_2 = CGAL::Vector_2<K>>
    std::string Pwh_to_geojson(const Polygon_with_holes_2& polygon, Vector_2 offset = Vector_2(0, 0)) {
        Block block;
        block.coords.push_back(convert_points<K>(polygon.outer_boundary(), offset));
        for (auto h_it = polygon.holes_begin(); h_it != polygon.holes_end(); ++h_it) {
            block.coords.push_back(convert_points<K>(*h_it, offset));
        }
        return out2str(block);
    }
    template <typename K, typename Polygon_2 = CGAL::Polygon_2<K>, typename Vector_2 = CGAL::Vector_2<K>>
    std::string Poly_to_geojson(const CGAL::Polygon_2<K>& polygon, CGAL::Vector_2<K> offset = CGAL::Vector_2<K>(0, 0)) {
        Block block;
        block.coords.push_back(convert_points<K>(polygon, offset));
        return out2str(block);
    }
}

namespace GeoJSON
{
    void parse_geojson(const std::string& datastr, parseout& out)
    {
        bool res;
        std::string errs;
        Json::Value root;
        Json::CharReaderBuilder readerBuilder;
        std::unique_ptr<Json::CharReader> const jsonReader(readerBuilder.newCharReader());

        res = jsonReader->parse(datastr.c_str(), datastr.c_str() + datastr.length(), &root, &errs);
        if (!res || !errs.empty()) {
            std::cout << "parseJson err. " << errs << std::endl;
        }
        int size = root["features"].size();
        for (auto i = 0; i < size; ++i)
        {
            auto g = root["features"][i]["geometry"];
            auto gtype = g["type"].asString();
            auto coor = g["coordinates"];
            auto p = root["features"][i]["properties"];
            auto cat = p["Category"].asString();
            auto privacy = p["Privacy"];
            auto sw = p["Switch"];
            auto n = p["Name"];
            bool flag = false;
            std::string name = cat;
            std::string ssw, sprivacy;
            if (!sw.isNull())
            {
                ssw = sw.asString();
            }
            if (!privacy.isNull())
            {
                sprivacy = privacy.asString();
            }
            if (!n.isNull())
            {
                name = n.asString();
            }
            if (gtype == "Polygon")
            {
                Block tb;
                std::vector<std::vector<Point>> tpwh;
                for (auto j = 0; j < coor.size(); ++j)
                {
                    std::vector<Point> tpoly;
                    Point last_p;
                    for (auto k = 0; k < coor[j].size(); ++k)
                    {
                        Point p{ coor[j][k][0].asDouble(),coor[j][k][1].asDouble() };
                        if (k && p.x == last_p.x && p.y == last_p.y) continue;
                        if (p.x < out.minx) out.minx = p.x;
                        if (p.x > out.maxx) out.maxx = p.x;
                        if (p.y < out.miny) out.miny = p.y;
                        if (p.y > out.maxy) out.maxy = p.y;
                        tpoly.push_back(p);
                    }
                    if (tpoly.size() > 2 || j == 0) tpwh.push_back(tpoly);
                }
                tb.coords = tpwh;
                tb.category = cat;
                tb.sw = ssw;
                tb.privacy = sprivacy;
                out.data.push_back(tb);
            }
        }
    }

    Json::Value parse_propsjson(const std::string& datastr) {
        Json::Reader reader;
        Json::Value root;

        std::ifstream in;
        in.open(datastr, std::ios::in || std::ios::binary);
        if (!in.is_open())
        {
            std::cout << "Error: open file is failed! " << std::endl;
        }
        if (!reader.parse(in, root))
        {
            std::cout << "parse string error!" << std::endl;
            return -1;
        }

        return root;
    }

    std::string out2str(std::vector<std::pair<Point, Point>> out)
    {
        //auto pt = out.p;
        //auto type = out.type;
        //auto oS = out.oS;
        //double sample_degree = pi / sample_num;
        Json::Value root, rf, lst;
        Json::Value oc, og, op, coords, point, pc;
        for (auto pr : out) {
            coords.clear();
            point.clear();
            point.append(pr.first.x);
            point.append(pr.first.y);
            coords.append(point);
            point.clear();
            point.append(pr.second.x);
            point.append(pr.second.y);
            coords.append(point);
            og["type"] = "LineString";
            og["coordinates"] = coords;
            oc["type"] = "Feature";
            oc["geometry"] = og;
            rf.append(oc);
        }

        root["features"] = rf;
        root["type"] = "FeatureCollection";
        return root.toStyledString();
    }

    Json::Value dump_block(const Block& block) { // polygon_with_holes: [polygon: [point: [x, y] ] ]
        Json::Value polys, poly, point;
        for (auto lst : block.coords) {
            poly.clear();
            for (auto p : lst) {
                point.clear();
                point.append(p.x); point.append(p.y);
                poly.append(point);
            }
            polys.append(poly);
        }
        return polys;
    }
    std::string out2str(const std::vector<Block>& ucs_blocks, const std::vector<Point>& corr_ucs) {
        Json::Value root, features;
        Json::Value feature, geometry, poly, coords, point, pc, v;
        for (size_t i = 0;i < ucs_blocks.size();++i) {
            auto& block = ucs_blocks[i];
            auto& ucs = corr_ucs[i];
            geometry["type"] = "Polygon";
            //coords = dump_block(block);
            coords.clear();
            for (auto lst : block.coords) {
                poly.clear();
                for (auto p : lst) {
                    point.clear();
                    point.append(p.x); point.append(p.y);
                    poly.append(point);
                }
                coords.append(poly);
            }
            geometry["coordinates"] = coords;
            feature["type"] = "Feature";
            feature["geometry"] = geometry;
            v.clear();
            v.append(ucs.x); v.append(ucs.y);
            feature["properties"]["ucs_direction"] = v;
            features.append(feature);
        }
        root["features"] = features;
        root["type"] = "FeatureCollection";
        return root.toStyledString();
    }
    std::string out2str(const std::vector<Block>& rect_blocks, const std::vector<Block>& centerline_blocks) {
        //auto pt = out.p;
        //auto type = out.type;
        //auto oS = out.oS;
        //double sample_degree = pi / sample_num;
        Json::Value root, features;
        Json::Value feature, geometry, poly, coords, point, pc;
        for (const auto& block : rect_blocks) {
            geometry["type"] = "Polygon";
            //coords = dump_block(block);
            coords.clear();
            for (auto lst : block.coords) {
                poly.clear();
                for (auto p : lst) {
                    point.clear();
                    point.append(p.x); point.append(p.y);
                    poly.append(point);
                }
                coords.append(poly);
            }
            geometry["coordinates"] = coords;
            feature["type"] = "Feature";
            feature["geometry"] = geometry;
            feature["properties"]["is_centerline_covered"] = false;
            features.append(feature);
        }
        for (const auto& block : centerline_blocks) {
            geometry["type"] = "Polygon";
            //coords = dump_block(block);
            coords.clear();
            for (auto lst : block.coords) {
                poly.clear();
                for (auto p : lst) {
                    point.clear();
                    point.append(p.x); point.append(p.y);
                    poly.append(point);
                }
                coords.append(poly);
            }
            geometry["coordinates"] = coords;
            feature["type"] = "Feature";
            feature["geometry"] = geometry;
            feature["properties"]["is_centerline_covered"] = true;
            features.append(feature);
        }
        root["features"] = features;
        root["type"] = "FeatureCollection";
        return root.toStyledString();
    }
    std::string out2str(const Block& block) {
        Json::Value root, features;
        Json::Value feature, geometry, poly, coords, point, pc;

        geometry["type"] = "Polygon";
        //coords = dump_block(block);
        coords.clear();
        for (auto lst : block.coords) {
            poly.clear();
            for (auto p : lst) {
                point.clear();
                point.append(p.x); point.append(p.y);
                poly.append(point);
            }
            coords.append(poly);
        }
        geometry["coordinates"] = coords;
        feature["type"] = "Feature";
        feature["geometry"] = geometry;
        feature["properties"]["is_centerline_covered"] = false;
        features.append(feature);

        root["features"] = features;
        root["type"] = "FeatureCollection";
        return root.toStyledString();
    }
}

#endif // GEOJSON_HPP