#ifndef PLACEMENT_HPP
#define PLACEMENT_HPP
#include "stdafx.h"

namespace Placement
{
    template <typename K>
    struct PlacementProps
    {
        PlacementProps() {

        }
    };
    
    template <typename K>
    class Solver
    {
        using FT = typename K::FT;

        using PlacementProps = PlacementProps<K>;
    private:
        /* data */
    public:
        Solver(const PlacementProps& props = PlacementProps());
    };
}

namespace Placement 
{
    template <typename K>
    Solver<K>::Solver(const PlacementProps& props)
    {
    }

} // namespace Placement
#endif // PLACEMENT_HPP