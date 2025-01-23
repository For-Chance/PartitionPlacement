#ifndef PARTITIONPLACEMENTEX_HPP
#define PARTITIONPLACEMENTEX_HPP
#include "PartitionPlacementContext.hpp"
#include <string>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace PartitionPlacementEx {
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    inline std::string API(const std::string& geojson, const ParitionPlacement::Context<K>& context) {
        return "";
    }
}
#endif // PARTITIONPLACEMENTEX_HPP