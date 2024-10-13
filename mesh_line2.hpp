#ifndef MESH_LINE2
#define MESH_LINE2
#include <unordered_map>
#include "container_typedef.hpp"

struct MeshLine2Struct
{
    /*

    struct with points and elements of a mesh domain made of line2 elements.

    Variables
    =========
    num_point_domain : int
        Number of points in the domain.
    point_gid_vec : VectorInt
        vector with global IDs of points in the domain.
    point_position_x_vec : VectorDouble
        vector with x-coordinates of points in the domain.
    point_gid_to_did_map : MapIntInt
        map that outputs the global ID of a point given its domain ID.
    num_element_domain : int
        Number of elements in the domain.
    element_gid_vec : VectorInt
        vector with global IDs of elements in the domain.
    element_p0_gid_vec : VectorInt
        vector with global IDs of point 0 in the domain.
    element_p1_gid_vec : VectorInt
        vector with global IDs of point 1 in the domain.
    point_gid_to_did_map : MapIntInt
        map that outputs the global ID of an element given its domain ID.

    Notes
    ====
    The global ID is a unique label for each point.
    The domain ID applies only to a domain and is used to iterate through the vectors in this code.
    Point 0 and 1 refer to the left and right points of the each element.

    */

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
