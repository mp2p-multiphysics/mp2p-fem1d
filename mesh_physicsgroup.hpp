#ifndef MESH_PHYSICSGROUP
#define MESH_PHYSICSGROUP
#include <vector>
#include "mesh_line2.hpp"

class MeshPhysicsGroup
{
    /*

    Groups mesh domains that are used in the same physics.

    Variables
    =========
    mesh_ptr_vec_in : vector<MeshLine2*>
        vector with pointers to MeshLine2 objects.
    
    */

    public:

    // vector with meshes in group
    std::vector<MeshLine2Struct*> mesh_ptr_vec;

    // default constructor
    MeshPhysicsGroup()
    {

    }

    // constructor
    MeshPhysicsGroup(std::vector<MeshLine2Struct*> mesh_ptr_vec_in)
    {
        
        // store variables
        mesh_ptr_vec = mesh_ptr_vec_in;

    }

};

#endif
