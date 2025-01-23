#ifndef CENTERLINEEX_HPP
#define CENTERLINEEX_HPP
#include "CenterLineContext.hpp"
#include <string>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CenterLineEx {
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    inline std::string Simplify(const std::string& geojson, const CenterLine::Context<K>& context) {
        return "";
    }
}
#endif // CENTERLINEEX_HPP