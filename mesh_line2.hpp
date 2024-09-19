#ifndef MESH_LINE2
#define MESH_LINE2
#include <vector>
#include "container_typedef.hpp"

struct MeshLine2Struct
{

    // point data
    int num_point = 0;
    VectorInt point_id_vec;
    VectorDouble point_pos_x_vec;

    // element data
    int num_element = 0;
    VectorInt element_id_vec;
    VectorInt element_p0_id_vec;
    VectorInt element_p1_id_vec;

};

#endif
