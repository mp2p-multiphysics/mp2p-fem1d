#ifndef BOUNDARY_PHYSICSGROUP_VECTOR
#define BOUNDARY_PHYSICSGROUP_VECTOR
#include <vector>
#include "boundary_physicsgroup.hpp"

class BoundaryPhysicsGroupVector
{
    /*

    Groups BoundaryPhysicsGroup objects into a column vector.
    Used in physics with systems of partial differential equations.

    Variables
    =========
    boundary_physicsgroup_ptr_vec : vector<BoundaryPhysicsGroup*>
        vector with pointers to BoundaryPhysicsGroup objects.
    
    Functions
    =========
    get_entry : BoundaryPhysicsGroup*
        Returns the pointer to the BoundaryPhysicsGroup object at a given vector position.

    */

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
    /*

    Returns the pointer to the BoundaryPhysicsGroup object at a vector position.

    Arguments
    =========
    vector_row : int
        Row in column vector
    
    Returns
    =======
    boundary_physicsgroup_ptr : BoundaryPhysicsGroup*
        Pointer to the BoundaryPhysicsGroup object at the specified row.

    */
    
    return boundary_physicsgroup_ptr_vec[vector_row];
}

#endif
