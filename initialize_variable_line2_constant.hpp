#ifndef INITIALIZE_VARIABLE_LINE2_CSV
#define INITIALIZE_VARIABLE_LINE2_CSV
#include <fstream>
#include <sstream>
#include <vector>
#include "container_typedef.hpp"
#include "mesh_line2.hpp"
#include "variable_line2.hpp"

VariableLine2Struct initialize_variable_line2_constant(MeshLine2Struct &mesh_l2s, double u)
{

    // create new variable
    VariableLine2Struct variable_l2s;

    // set mesh
    variable_l2s.mesh_l2s = &mesh_l2s;

    // set initial values
    for (int indx_i = 0; indx_i < mesh_l2s.num_point; indx_i++)
    {
        variable_l2s.point_u_vec.push_back(u);
    }

    return variable_l2s;

}


#endif
