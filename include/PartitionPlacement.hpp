#include "Partition.hpp"
#include "Placement.hpp"

#include <iostream>
#include <string>

namespace PP
{
    class Solver
    {
    private:
        Partition::Solver PartitionSolver;
        Placement::Solver PlacementSolver;
    public:
        Solver();
        Solver(std::string config_file);
        ~Solver();
    };

    Solver::Solver()
    {
    }

    Solver::Solver(std::string config_file) {
        std::cout << "Hello PartitionPlacement" << std::endl;
    }

    Solver::~Solver()
    {
    }

} // namespace PartitionPlacement
