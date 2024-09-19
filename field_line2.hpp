#ifndef FIELD_LINE2
#define FIELD_LINE2
#include "container_mesh.hpp"

class FieldLine2
{

    public:

    // variables
    MeshLine2Struct *mesh_l2_ptr;
    VectorDouble point_u_vec;

    // default constructor
    FieldLine2()
    {

    }

    // constructor
    FieldLine2(MeshLine2Struct &mesh_l2_in, double u_init_in)
    {

        // store mesh
        mesh_l2_ptr = &mesh_l2_in;

        // populate initial values
        for (int indx_i = 0; indx_i < mesh_l2_ptr->num_point; indx_i++)
        {
            point_u_vec.push_back(u_init_in);
        }

    }

};

#endif
