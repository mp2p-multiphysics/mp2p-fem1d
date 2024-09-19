#ifndef VARIABLE_LINE2
#define VARIABLE_LINE2
#include "container_typedef.hpp"
#include "mesh_line2.hpp"

struct VariableLine2Struct
{

    // mesh where variable is present
    MeshLine2Struct* mesh_l2s;

    // data in mesh points
    VectorDouble point_u_vec;

};

#endif
