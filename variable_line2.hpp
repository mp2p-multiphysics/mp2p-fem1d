#ifndef VARIABLE_LINE2
#define VARIABLE_LINE2
#include "container_mesh.hpp"
#include "container_typedef.hpp"

class VariableLine2
{

    public:

    // mesh and data points
    MeshLine2Struct* mesh_l2_ptr;
    VectorDouble point_u_vec;

    // starting column in matrix equation
    int start_col = -1;

    // default constructor
    VariableLine2()
    {

    }

    // constructor
    VariableLine2(MeshLine2Struct &mesh_l2_in, double u_init_in)
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
