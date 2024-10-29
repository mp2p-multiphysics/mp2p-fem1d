#ifndef BOUNDARY_PHYSICSGROUP
#define BOUNDARY_PHYSICSGROUP
#include <vector>
#include "boundary_line2.hpp"

class BoundaryField
{
    /*

    Groups boundary conditions (BC) that are used in the same physics.

    Variables
    =========
    boundary_l2_ptr_vec_in : vector<BoundaryLine2*>
        vector with pointers to BoundaryLine2 objects.
    
    */

    public:

    // vector with boundaries in group
    std::vector<BoundaryLine2*> boundary_l2_ptr_vec;

    // default constructor
    BoundaryField()
    {

    }

    // constructor
    BoundaryField(std::vector<BoundaryLine2*> boundary_l2_ptr_vec_in)
    {
        boundary_l2_ptr_vec = boundary_l2_ptr_vec_in;
    }

};

#endif
