#include "PartitionPlacement.hpp"

#include <iostream>
#include <string>

namespace PP {
    struct PartitionPlacementTest
    {
        PartitionPlacementTest(const std::string& config_file);

        void showResult(const PP::Solver& ppSolver);
    };

    PartitionPlacementTest::PartitionPlacementTest(const std::string& config_file) {
        std::cout << "PartitionPlacementTest(const std::string& config_file)" << std::endl;
        PP::Solver ppSolver(config_file);
        showResult(ppSolver);
    }

    void PartitionPlacementTest::showResult(const PP::Solver& ppSolver) {
        std::cout << "Show Result:" << std::endl;
    }
}