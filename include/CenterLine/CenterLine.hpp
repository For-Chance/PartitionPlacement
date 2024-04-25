#ifndef CENTERLINE_HPP
#define CENTERLINE_HPP
#include "CenterLineContext.hpp"
#include "GeoJSON.hpp"
#include "stdafx.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>

namespace CenterLine {
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
    using Polygon_2 = CGAL::Polygon_2<K>;
    using Point_2 = CGAL::Point_2<K>;

    class Solver {
    public:
        Polygon_with_holes_2 origin_space;
    private:

    public:
        Solver() {}
        Solver(const std::string& geojson = "", const Context& context = Context());
    };
}

namespace CenterLine {
    Solver::Solver(const std::string& geojson, const Context& context) {

    }
}
#endif // CENTERLINE_HPP