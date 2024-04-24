#include "Partition.hpp"
#include "Placement.hpp"
#include "PartitionPlacementContext.hpp"
#include "PPGeoJSON.hpp"
#include "stdafx.h"

#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>

namespace ParitionPlacement
{
    using PartitionProps = Partition::PartitionProps;
    using PlacementProps = Placement::PlacementProps;

    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
    using Polygon_2 = CGAL::Polygon_2<K>;
    using Point_2 = CGAL::Point_2<K>;

    class Solver
    {
    public:
        Polygon_with_holes_2 origin_space;
    private:
        Partition::Solver PartitionSolver;
        Placement::Solver PlacementSolver;
        
    public:
        Solver(const std::string& geojson = "", const Context& context = Context());

        static Polygon_2 convert_poly(std::vector<Point>& points);
        static Polygon_with_holes_2 geojson_to_Pwh(std::string geojson);
    };

    /// <summary>
    /// Constructor
    /// </summary>
    /// <param name="geojson">GeoJSON string</param>
    /// <param name="context">Context</param>
    /// <returns></returns>
    Solver::Solver(const std::string& geojson, const Context& context) {
        std::cout << "Hello PartitionPlacement" << std::endl;
        origin_space = geojson_to_Pwh(geojson);
        PartitionSolver = Partition::Solver(origin_space, context.partProps);
        PlacementSolver = Placement::Solver(context.placeProps);
    }

    Polygon_2 Solver::convert_poly(std::vector<Point>& points) {
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

    Polygon_with_holes_2 Solver::geojson_to_Pwh(std::string geojson) {
        parseout res;
        parse_geojson(geojson, res);
        auto& data = res.data[0];
        Polygon_2 poly = convert_poly(data.coords[0]);  // the first group geojson data, which is the outer boundary
        if (poly.is_clockwise_oriented()) poly.reverse_orientation();
        std::vector<Polygon_2> holes;
        for (int i = 1; i < data.coords.size(); ++i) {
            Polygon_2 hole = convert_poly(data.coords[i]);
            if (hole.is_counterclockwise_oriented()) hole.reverse_orientation();
            holes.push_back(hole);
        }
        if (holes.empty()) return Polygon_with_holes_2(poly);
        else return Polygon_with_holes_2(poly, holes.begin(), holes.end());
    }
} // namespace PartitionPlacement
