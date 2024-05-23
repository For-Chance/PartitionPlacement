#ifndef PARTITIONPLACEMENT_HPP
#define PARTITIONPLACEMENT_HPP
#include "Partition.hpp"
#include "Placement.hpp"
#include "GeoJSON.hpp"
#include "stdafx.h"

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>

namespace ParitionPlacement
{
    template <typename K>
    class Solver
    {
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;
		using Segment_2 = CGAL::Segment_2<K>;

        using PartitionProps = Partition::PartitionProps<K>;
        using PlacementProps = Placement::PlacementProps<K>;
    public:
        Polygon_with_holes_2 origin_space;
        Polygon_with_holes_2 polygon;
        std::vector<Segment_2> get_skeleton_segmnts() const { return PartitionSolver.skeleton_segments; }
		std::vector<Segment_2> get_skeleton_centerlines() const { return PartitionSolver.skeleton_centerlines; }
		std::vector<Segment_2> get_skeleton_otherlines() const { return PartitionSolver.skeleton_otherlines; }
		std::vector<Polygon_2> get_skeleton_faces() const { return PartitionSolver.skeleton_faces; }
		std::vector<std::vector<Polygon_2>> get_init_partition() const{ return PartitionSolver.init_partition; }
		std::vector<Polygon_2> get_uncertain_parts() const { return PartitionSolver.uncertain_parts; }
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
        if(geojson == "")
			throw std::invalid_argument("\'geojson\' is empty!");
        origin_space = GeoJSON::geojson_to_Pwh<K>(geojson);
        PartitionSolver = Partition::Solver<K>(origin_space, partitionProps);
		polygon = PartitionSolver.polygon;
        PlacementSolver = Placement::Solver<K>(placementProps);
    }
} // namespace PartitionPlacement
#endif // PARTITIONPLACEMENT_HPP
