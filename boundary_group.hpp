#ifndef BOUNDARY_GROUP
#define BOUNDARY_GROUP
#include <vector>
#include "boundary_line2.hpp"

class BoundaryGroup
{
    /*

    Groups boundary conditions (BC) that are used in the same physics.

    Variables
    =========
    boundary_l2_ptr_vec_in : vector<BoundaryLine2*>
        vector with pointers to BoundaryLine2 objects.
    
    Functions
    =========
    update_parameter : void
        Recalculates non-constant boundary condition parameters.

    */

    public:

    // vector with boundaries in group
    std::vector<BoundaryLine2*> boundary_l2_ptr_vec;

    // functions
    void set_boundary_type(VectorInt boundarytype_essential_vec, VectorInt boundarytype_natural_vec);
    void update_parameter();

    // default constructor
    BoundaryGroup() {}

    // constructor
    BoundaryGroup(std::vector<BoundaryLine2*> boundary_l2_ptr_vec_in)
    {
        boundary_l2_ptr_vec = boundary_l2_ptr_vec_in;
    }

};


void BoundaryGroup::set_boundary_type(VectorInt boundarytype_essential_vec, VectorInt boundarytype_natural_vec)
{

    // iterate through each boundary
    for (auto boundary_ptr : boundary_l2_ptr_vec)
    {
        boundary_ptr->set_boundary_type(boundarytype_essential_vec, boundarytype_natural_vec);
    }

}

void BoundaryGroup::update_parameter()
{
    /*

    Recalculates non-constant boundary condition parameters.

    Arguments
    =========
    (none)
    
    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (auto boundary_ptr : boundary_l2_ptr_vec)
    {
        boundary_ptr->update_parameter();
    }

}

#endif
