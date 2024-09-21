#ifndef SCALAR_LINE2
#define SCALAR_LINE2
#include "mesh_line2.hpp"

class ScalarLine2
{

    public:

    // variables
    MeshLine2Struct *mesh_l2_ptr;
    VectorDouble point_u_vec;

    // default constructor
    ScalarLine2()
    {

    }

    // constructor
    ScalarLine2(MeshLine2Struct &mesh_l2_in, double u_init_in)
    {

        // store mesh
        mesh_l2_ptr = &mesh_l2_in;

        // populate initial values
        for (int indx_i = 0; indx_i < mesh_l2_ptr->num_domain_point; indx_i++)
        {
            point_u_vec.push_back(u_init_in);
        }

    }

};

#endif
