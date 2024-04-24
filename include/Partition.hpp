#ifndef PARTITION_HPP
#define PARTITION_HPP
#include "stdafx.h"

namespace Partition
{
    struct PartitionProps {

        PartitionProps() {
               
        }
    };

    class Solver
    {
    private:
        /* data */
    public:
        Solver(const PartitionProps& props);
    };

    Solver::Solver(const PartitionProps& props = PartitionProps())
    {
    }
} // namespace Partition
#endif // PARTITION_HPP