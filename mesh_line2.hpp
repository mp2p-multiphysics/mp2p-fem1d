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
    int num_domain_point = 0;
    VectorInt point_gid_vec;
    VectorDouble point_position_x_vec;

    // element data
    int num_domain_element = 0;
    VectorInt element_gid_vec;
    VectorInt element_p0_gid_vec;
    VectorInt element_p1_gid_vec;

    // map of global ID to domain ID
    std::unordered_map<int, int> point_gid_to_did_map;
    std::unordered_map<int, int> element_gid_to_did_map;

};

#endif
