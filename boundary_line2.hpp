#ifndef BOUNDARY_LINE2
#define BOUNDARY_LINE2
#include <vector>
#include "container_typedef.hpp"

struct BoundaryConfigLine2Struct
{

    std::string boundary_type_str;
    VectorDouble boundary_parameter_vec;

};

struct BoundaryLine2Struct
{

    // gid - global ID
    // lid - local ID

    // flux boundary condition data
    int num_element_flux_domain = 0;
    VectorInt element_flux_gid_vec;
    VectorInt element_flux_pa_lid_vec;
    VectorInt element_flux_config_id_vec;

    // value boundary condition data
    int num_element_value_domain = 0;
    VectorInt element_value_gid_vec;
    VectorInt element_value_pa_lid_vec;
    VectorInt element_value_config_id_vec;

    // boundary condition data
    std::vector<BoundaryConfigLine2Struct> boundary_config_vec;

};

#endif
