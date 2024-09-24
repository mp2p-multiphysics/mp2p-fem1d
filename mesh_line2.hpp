#ifndef MESH_LINE2
#define MESH_LINE2
#include <unordered_map>
#include "container_typedef.hpp"

struct MeshLine2Struct
{

    // did - domain ID
    // gid - global ID
    // vectors use did as input

    // point data
    int num_point_domain = 0;
    VectorInt point_gid_vec;
    VectorDouble point_position_x_vec;
    MapIntInt point_gid_to_did_map;

    // element data
    int num_element_domain = 0;
    VectorInt element_gid_vec;
    VectorInt element_p0_gid_vec;
    VectorInt element_p1_gid_vec;
    MapIntInt element_gid_to_did_map;

};

#endif
