#ifndef GRID_LINE2
#define GRID_LINE2
#include <vector>

struct GridLine2Struct
{

    // point data
    int num_point = 0;
    std::vector<int> point_id_vec;
    std::vector<double> point_pos_x_vec;

    // element data
    int num_element = 0;
    std::vector<int> element_id_vec;
    std::vector<int> element_p0_id_vec;
    std::vector<int> element_p1_id_vec;

};

struct BoundaryConfigLine2Struct
{

    std::string boundary_type_str;
    std::vector<double> boundary_parameter_vec;

};

struct BoundaryLine2Struct
{

    // flux boundary condition data
    int num_element_flux = 0;
    std::vector<int> element_flux_id_vec;
    std::vector<int> element_flux_pa_loc_vec;
    std::vector<int> element_flux_config_id_vec;
    
    // value boundary condition data
    int num_element_value = 0;
    std::vector<int> element_value_id_vec;
    std::vector<int> element_value_pa_loc_vec;
    std::vector<int> element_value_config_id_vec;

    // boundary condition data
    std::vector<BoundaryConfigLine2Struct> boundary_config_vec;

};

#endif
