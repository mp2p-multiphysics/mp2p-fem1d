#ifndef BOUNDARY_PHYSICSGROUP_VECTOR
#define BOUNDARY_PHYSICSGROUP_VECTOR
#include <vector>
#include "boundary_physicsgroup.hpp"

class BoundaryPhysicsGroupVector
{

    public:

    // variables
    int num_entry = 0;
    std::vector<BoundaryPhysicsGroup*> boundary_physicsgroup_ptr_vec;

    // functions
    BoundaryPhysicsGroup* get_entry(int vector_row);

    // default constructor
    BoundaryPhysicsGroupVector()
    {

    }

    // constructor
    BoundaryPhysicsGroupVector(std::vector<BoundaryPhysicsGroup*> boundary_physicsgroup_ptr_vec_in)
    {

        // store vector of variables
        boundary_physicsgroup_ptr_vec = boundary_physicsgroup_ptr_vec_in;

        // count number of variables
        num_entry = boundary_physicsgroup_ptr_vec.size();

    }

};

BoundaryPhysicsGroup* BoundaryPhysicsGroupVector::get_entry(int vector_row)
{
    return boundary_physicsgroup_ptr_vec[vector_row];
}

#endif
