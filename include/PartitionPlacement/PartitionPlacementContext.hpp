#ifndef PARTITIONPLACEMENTCONTEXT_HPP
#define PARTITIONPLACEMENTCONTEXT_HPP
#include "Partition.hpp"
#include "Placement.hpp"
#include <string>

namespace ParitionPlacement {
    using PartitionProps = Partition::PartitionProps;
    using PlacementProps = Placement::PlacementProps;

    struct PartitionPlacementProps {
        PartitionProps partProps;
        PlacementProps placeProps;

        PartitionPlacementProps() :partProps(PartitionProps()), placeProps(PlacementProps()) {}
        PartitionPlacementProps(const PartitionProps& partProps, const PlacementProps& placeProps) :partProps(partProps), placeProps(placeProps) {}
    };

    struct Context {
        PartitionPlacementProps ppProps;

        Context() :ppProps(PartitionPlacementProps()) {}
        Context(const PartitionPlacementProps& ppProps) :ppProps(ppProps) {}
    };
}
#endif // !PARTITIONPLACEMENTCONTEXT_HPP