#ifndef PARTITIONPLACEMENT_HPP
#define PARTITIONPLACEMENT_HPP
#include "Partition.hpp"
#include "Placement.hpp"
#include "PartitionPlacementContext.hpp"
#include "GeoJSON.hpp"
#include "stdafx.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>

namespace ParitionPlacement
{
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
    };
}

namespace ParitionPlacement
{
    /// <summary>
    /// Constructor
    /// </summary>
    /// <param name="geojson">GeoJSON string</param>
    /// <param name="context">Context</param>
    /// <returns></returns>
    Solver::Solver(const std::string& geojson, const Context& context) {
        std::cout << "Hello PartitionPlacement" << std::endl;
        origin_space = GeoJSON::geojson_to_Pwh<K>(geojson);
        PartitionSolver = Partition::Solver(origin_space, context.ppProps.partProps);
        PlacementSolver = Placement::Solver(context.ppProps.placeProps);
    }
} // namespace PartitionPlacement
#endif // PARTITIONPLACEMENT_HPP
