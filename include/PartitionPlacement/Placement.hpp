#ifndef PLACEMENT_HPP
#define PLACEMENT_HPP
#include "stdafx.h"

namespace Placement
{
    struct PlacementProps
    {
        PlacementProps() {

        }
    };
    class Solver
    {
    private:
        /* data */
    public:
        Solver(const PlacementProps& props = PlacementProps());
    };
}

namespace Placement 
{
    Solver::Solver(const PlacementProps& props)
    {
    }

} // namespace Placement
#endif // PLACEMENT_HPP