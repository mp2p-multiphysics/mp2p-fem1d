#ifndef BOUNDARY_PHYSICSGROUP
#define BOUNDARY_PHYSICSGROUP
#include <vector>
#include "boundary_line2.hpp"

class BoundaryPhysicsGroup
{

    public:

    // vector with boundaries in group
    std::vector<BoundaryLine2*> boundary_ptr_vec;

    // default constructor
    BoundaryPhysicsGroup()
    {

    }

    // constructor
    BoundaryPhysicsGroup(std::vector<BoundaryLine2*> boundary_ptr_vec_in)
    {
        boundary_ptr_vec = boundary_ptr_vec_in;
    }

};

#endif
