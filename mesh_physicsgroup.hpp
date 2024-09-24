#ifndef MESH_PHYSICSGROUP
#define MESH_PHYSICSGROUP
#include <vector>
#include "mesh_line2.hpp"

class MeshPhysicsGroup
{

    public:

    // variables
    int num_domain = 0;
    std::vector<MeshLine2Struct*> mesh_ptr_vec;
    int num_points_physics = 0;
    

    // functions

    // default constructor
    MeshPhysicsGroup()
    {

    }

    // constructor
    MeshPhysicsGroup(std::vector<MeshLine2Struct*> mesh_ptr_vec_in)
    {
        
        // store variables
        mesh_ptr_vec = mesh_ptr_vec_in;

        // get number of domains
        num_domain = mesh_ptr_vec.size();

    }

};

#endif
