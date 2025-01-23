#ifndef SIMPLIFYBOUNDARYEX_HPP
#define SIMPLIFYBOUNDARYEX_HPP
#include "SimplifyBoundaryContext.hpp"
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace SimplifyBoundaryEx {
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    inline std::string Simplify(const std::string& geojson, const SimplifyBoundary::Context<K>& context) {
        return "";
    }
}
#endif // SIMPLIFYBOUNDARYEX_HPP