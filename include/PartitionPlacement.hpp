#include "Partition.hpp"
#include "Placement.hpp"
#include "PartitionPlacementContext.hpp"
#include "stdafx.h"

#include <iostream>
#include <string>

namespace ParitionPlacement
{
    using PartitionProps = Partition::PartitionProps;
    using PlacementProps = Placement::PlacementProps;
    class Solver
    {
    private:
        Partition::Solver PartitionSolver;
        Placement::Solver PlacementSolver;
    public:
        Solver(const std::string& geojson, const Context& context);
    };
    Solver::Solver(const std::string& geojson, const Context& context) {
        std::cout << "Hello PartitionPlacement" << std::endl;
        PartitionSolver = Partition::Solver(context.partProps);
        PlacementSolver = Placement::Solver(context.placeProps);
    }

} // namespace PartitionPlacement
