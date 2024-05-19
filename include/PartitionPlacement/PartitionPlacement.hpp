#ifndef PARTITIONPLACEMENT_HPP
#define PARTITIONPLACEMENT_HPP
#include "Partition.hpp"
#include "Placement.hpp"
#include "GeoJSON.hpp"
#include "stdafx.h"

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>

namespace ParitionPlacement
{
    template <typename K>
    class Solver
    {
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;

        using PartitionProps = Partition::PartitionProps<K>;
        using PlacementProps = Placement::PlacementProps<K>;
    public:
        Polygon_with_holes_2 origin_space;
    private:
        Partition::Solver<K> PartitionSolver;
        Placement::Solver<K> PlacementSolver;

    public:
        Solver(const std::string& geojson = "", const PartitionProps& partitionProps = PartitionProps<K>(), const PlacementProps& placementProps  = PlacementProps());
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
    template <typename K>
    Solver<K>::Solver(const std::string& geojson, const PartitionProps& partitionProps, const PlacementProps& placementProps) {
        std::cout << "Hello PartitionPlacement" << std::endl;
        origin_space = GeoJSON::geojson_to_Pwh<K>(geojson);
        PartitionSolver = Partition::Solver<K>(origin_space, partitionProps);
        PlacementSolver = Placement::Solver<K>(placementProps);
    }
} // namespace PartitionPlacement
#endif // PARTITIONPLACEMENT_HPP
