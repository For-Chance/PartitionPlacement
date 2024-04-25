#ifndef PARTITION_HPP
#define PARTITION_HPP
#include "stdafx.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/partition_2.h>
#include <list>

namespace Partition
{
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
    using Polygon_2 = CGAL::Polygon_2<K>;
    using Point_2 = CGAL::Point_2<K>;
    using Polygon_list = std::list<Polygon_2>;

    struct PartitionProps {

        PartitionProps() {

        }
    };

    class Solver
    {
    private:
        /* data */
    public:
        Solver() {}
        Solver(const Polygon_with_holes_2& space, const PartitionProps& props = PartitionProps());
    };
}

namespace Partition
{
    Solver::Solver(const Polygon_with_holes_2& space, const PartitionProps& props)
    {
        
       
    }
} // namespace Partition
#endif // PARTITION_HPP