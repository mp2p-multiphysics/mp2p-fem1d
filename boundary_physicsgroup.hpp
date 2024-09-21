#ifndef BOUNDARY_PHYSICSGROUP
#define BOUNDARY_PHYSICSGROUP
#include <vector>
#include "boundary_line2.hpp"

class BoundaryPhysicsGroup
{

    public:

    // variables
    std::vector<BoundaryLine2Struct*> boundary_ptr_vec;

    // functions

    // default constructor
    BoundaryPhysicsGroup()
    {

    }

    // constructor
    BoundaryPhysicsGroup(std::vector<BoundaryLine2Struct*> boundary_ptr_vec_in)
    {
        boundary_ptr_vec = boundary_ptr_vec_in;
    }

};

#endif
