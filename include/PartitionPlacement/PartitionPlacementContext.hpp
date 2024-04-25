#ifndef PARTITIONPLACEMENTCONTEXT_HPP
#define PARTITIONPLACEMENTCONTEXT_HPP
#include "Partition.hpp"
#include "Placement.hpp"
#include <string>
namespace ParitionPlacement {
    struct Context {
        using PartitionProps = Partition::PartitionProps;
        using PlacementProps = Placement::PlacementProps;

        PartitionProps partProps;
        PlacementProps placeProps;
        Context() :partProps(PartitionProps()), placeProps(PlacementProps()) {}
        Context(const PartitionProps& partProps, const PlacementProps& placeProps) :partProps(partProps), placeProps(placeProps) {}
    };
}
#endif // !PARTITIONPLACEMENTCONTEXT_HPP