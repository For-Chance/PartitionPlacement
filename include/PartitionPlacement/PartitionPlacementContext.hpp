#ifndef PARTITIONPLACEMENTCONTEXT_HPP
#define PARTITIONPLACEMENTCONTEXT_HPP
#include "Partition.hpp"
#include "Placement.hpp"
#include "PartitionPlacement.hpp"
#include <string>

namespace ParitionPlacement {
    template <typename K>
    struct PartitionPlacementProps {
        using PartitionProps = Partition::PartitionProps<K>;
        using PlacementProps = Placement::PlacementProps<K>;

        PartitionProps partProps;
        PlacementProps placeProps;

        PartitionPlacementProps() :partProps(PartitionProps()), placeProps(PlacementProps()) {}
        PartitionPlacementProps(const PartitionProps& partProps, const PlacementProps& placeProps) :partProps(partProps), placeProps(placeProps) {}
    };

    template <typename K>
    struct Context {
        using PartitionPlacementProps = PartitionPlacementProps<K>;

        PartitionPlacementProps ppProps;

        Context() :ppProps(PartitionPlacementProps()) {}
        Context(const PartitionPlacementProps& ppProps) :ppProps(ppProps) {}
    };
}
#endif // !PARTITIONPLACEMENTCONTEXT_HPP