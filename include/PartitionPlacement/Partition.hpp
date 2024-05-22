#ifndef PARTITION_HPP
#define PARTITION_HPP
#include "stdafx.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <list>
#include <chrono>

#include "SimplifyBoundary.hpp"

namespace Partition
{
    template <typename K>
    struct PartitionProps {
        bool withSimplifyBoundary;
        PartitionProps() {
            withSimplifyBoundary = true;
        }
    };

    template <typename K>
    class Solver
    {
        using FT = typename K::FT;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;
		using Segment_2 = CGAL::Segment_2<K>;
        using Polygon_list = std::list<Polygon_2>;
        using Ss = CGAL::Straight_skeleton_2<K>;
        using Halfedge_iterator = typename Ss::Halfedge_iterator;
        using Halfedge_handle = typename Ss::Halfedge_handle;
        using Vertex_handle = typename Ss::Vertex_handle;
        using Vertex_iterator = typename Ss::Vertex_iterator;
		using Face_handle = typename Ss::Face_handle;
		using Face_iterator = typename Ss::Face_iterator;
        using SsBuilderTraits = CGAL::Straight_skeleton_builder_traits_2<K>;
        using SsBuilder = CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss>;

        using PartitionProps = PartitionProps<K>;
    private:
        /// <summary>
		/// build skeleton of a Polygon_with_holes_2
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        boost::shared_ptr<Ss> build_skeleton(const Polygon_with_holes_2& polygon) {
            std::cout << "start building skeleton..." << std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            SsBuilder ssb;
            ssb.enter_contour(polygon.outer_boundary().vertices_begin(), polygon.outer_boundary().vertices_end());
            for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
                ssb.enter_contour(hole->vertices_begin(), hole->vertices_end());
			boost::shared_ptr<Ss> skeleton = ssb.construct_skeleton();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
			std::cout << "skeleton built! Time elapsed:" << diff.count() << " s, num of skeleton face:" << skeleton->size_of_faces()<< std::endl;
            return skeleton;
        };
    public:
        Solver() {}
        Solver(const Polygon_with_holes_2& space, const PartitionProps& props = PartitionProps());
        
        Polygon_with_holes_2 origin_space;
        Polygon_with_holes_2 polygon;
        boost::shared_ptr<Ss> skeleton;
		std::vector<Segment_2> skeleton_segments;
        std::vector<Polygon_2> skeleton_faces;
    };
}

namespace Partition
{
    template <typename K>
    Solver<K>::Solver(const Polygon_with_holes_2& space, const PartitionProps& props)
    {
		this->origin_space = space;
        if (props.withSimplifyBoundary) {
            // 1. Simplify Boundary
			SimplifyBoundary::Solver<K> sbSolver(this->origin_space);
            this->polygon = sbSolver.simplify_space;
            
            // 2. skeleton 
			this->skeleton = build_skeleton(this->polygon);
            for (Halfedge_iterator it = this->skeleton->halfedges_begin(); it != this->skeleton->halfedges_end(); ++it) {
                Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
				skeleton_segments.push_back(Segment_2(from->point(), to->point()));
            }
            
            // 3. partition
            // 3.1 save to skeleton_faces
            for (Face_iterator it = this->skeleton->faces_begin(); it != this->skeleton->faces_end(); ++it) {
                Polygon_2 poly;
                Face_handle face = it;
                Halfedge_handle h = face->halfedge();
                Halfedge_handle first = h;
                do {
                    Point_2 point = h->vertex()->point();
                    poly.push_back(point);
                    h = h->next();
                } while (h != first);
                skeleton_faces.push_back(poly);
            }
            // 3.2 find center line
        }
		else
			polygon = origin_space;
    }
} // namespace Partition
#endif // PARTITION_HPP