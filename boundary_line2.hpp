#ifndef BOUNDARY_LINE2
#define BOUNDARY_LINE2
#include <fstream>
#include <functional>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "container_typedef.hpp"
#include "variable_line2.hpp"

class BoundaryLine2
{
    /*

    Boundary conditions (BC) for line2 mesh elements.

    Variables
    =========
    file_in_flux_str_in : string
        Path to CSV file with data for flux-type BCs.
    file_in_value_str_in : string
        Path to CSV file with data for value-type BCs.

    Functions
    =========
    set_boundarycondition : void
        Assigns a BC type and parameters to a BC configuration ID.
    set_boundarycondition_parameter : void
        Assigns or modifies the parameters to a BC.

    Notes
    ====
    Both CSV files must have the following columns:
        global element ID where BC is applied
        local point ID where BC is applied (0 or 1)
        BC configuration ID
    Flux-type BCs add additional terms to the linearized equations (e.g., Neumann, Robin)
    Value-type BCs completely replace the linearized equations (e.g., Dirichlet)

    */

    public:

    // lid - local ID
    // bid - boundary ID - used to index BC data
    // bfid - boundary flux ID - used to index flux-type BC data
    // bvid - boundary value ID - used to index value-type BC data
    // bcid - boundary config ID - denotes type of boundary conditions
    // vectors use did as input

    // mesh where variable is applied
    MeshLine2* mesh_ptr;

    // file names
    std::string file_in_str;

    // flux boundary condition data
    int num_boundary_domain = 0;
    VectorInt boundary_element_gid_vec;
    VectorInt boundary_pa_lid_vec;
    VectorInt boundary_pa_bcid_vec;

    // flux boundary condition data
    int num_boundary_flux_domain = 0;
    VectorInt boundary_flux_element_gid_vec;
    VectorInt boundary_flux_pa_lid_vec;
    VectorInt boundary_flux_pa_bcid_vec;
    std::vector<std::string> boundary_flux_type_str_vec;
    Vector2D boundary_flux_parameter_vec;

    // used for non-constant flux boundary condition data
    // use bcid to index
    std::unordered_map<int, VectorInt> bcid_to_bfid_map;
    std::unordered_map<int, bool> boundary_flux_is_parameter_constant_map;
    std::unordered_map<int, std::function<VectorDouble(double, VectorDouble)>> boundary_flux_parameter_function_map;
    std::unordered_map<int, std::vector<VariableLine2*>> boundary_flux_variable_ptr_map;

    // value boundary condition data
    int num_boundary_value_domain = 0;
    VectorInt boundary_value_element_gid_vec;
    VectorInt boundary_value_pa_lid_vec;
    VectorInt boundary_value_pa_bcid_vec;
    std::vector<std::string> boundary_value_type_str_vec;
    Vector2D boundary_value_parameter_vec;

    // used for non-constant value boundary condition data
    // use bcid to index
    std::unordered_map<int, VectorInt> bcid_to_bvid_map;
    std::unordered_map<int, bool> boundary_value_is_parameter_constant_map;
    std::unordered_map<int, std::function<VectorDouble(double, VectorDouble)>> boundary_value_parameter_function_map;
    std::unordered_map<int, std::vector<VariableLine2*>> boundary_value_variable_ptr_map;

    // functions
    void set_boundary_flux(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec);
    void set_boundary_value(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec);
    void set_boundary_flux(int boundaryconfig_id, std::string type_str, std::function<VectorDouble(double, VectorDouble)> parameter_function, std::vector<VariableLine2*> variable_ptr_vec);
    void set_boundary_value(int boundaryconfig_id, std::string type_str, std::function<VectorDouble(double, VectorDouble)> parameter_function, std::vector<VariableLine2*> variable_ptr_vec);
    void update_parameter();

    // default constructor
    BoundaryLine2()
    {

    }

    // constructor
    BoundaryLine2(MeshLine2 &mesh_in, std::string file_in_str_in)
    {

        // store variables
        mesh_ptr = &mesh_in;
        file_in_str = file_in_str_in;

        // read input files and store boundary condition data
        read_boundary(file_in_str);

    }

    private:
    void read_boundary(std::string file_in_str);

};

void BoundaryLine2::set_boundary_flux(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec)
{
    /*

    Assigns a BC type and parameters to a BC configuration ID.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    type_str : string
        Type of boundary condition.
    parameter_vec : VectorDouble
        vector with parameters for the BC.

    Returns
    =======
    (none)

    Notes
    ====
    type_str can be "neumann" or "robin" if boundaryconfig_id refers to flux-type BCs.
    type_str can be "dirichlet" if boundaryconfig_id refers to value-type BCs.

    */

    // iterate through each boundary condition
    for (int bid = 0; bid < num_boundary_domain; bid++)
    {

        // skip if different boundary config ID
        if (boundary_pa_bcid_vec[bid] != boundaryconfig_id)
        {
            continue;
        }

        // map boundaryconfig_id to boundaryflux_id
        bcid_to_bfid_map[boundaryconfig_id].push_back(num_boundary_flux_domain);

        // transfer data to flux-type boundary arrays
        num_boundary_flux_domain++;
        boundary_flux_element_gid_vec.push_back(boundary_element_gid_vec[bid]);
        boundary_flux_pa_lid_vec.push_back(boundary_pa_lid_vec[bid]);
        boundary_flux_pa_bcid_vec.push_back(boundary_pa_bcid_vec[bid]);
        boundary_flux_type_str_vec.push_back(type_str);
        boundary_flux_parameter_vec.push_back(parameter_vec);

        // mark as constant BC parameters
        boundary_flux_is_parameter_constant_map[boundaryconfig_id] = true;

    }

}

void BoundaryLine2::set_boundary_flux(int boundaryconfig_id, std::string type_str, std::function<VectorDouble(double, VectorDouble)> parameter_function, std::vector<VariableLine2*> variable_ptr_vec)
{
    /*

    Assigns a BC type and parameters to a BC configuration ID.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    type_str : string
        Type of boundary condition.
    parameter_vec : VectorDouble
        vector with parameters for the BC.

    Returns
    =======
    (none)

    Notes
    ====
    type_str can be "neumann" or "robin" if boundaryconfig_id refers to flux-type BCs.
    type_str can be "dirichlet" if boundaryconfig_id refers to value-type BCs.

    */

    // iterate through each boundary condition
    for (int bid = 0; bid < num_boundary_domain; bid++)
    {

        // skip if different boundary config ID
        if (boundary_pa_bcid_vec[bid] != boundaryconfig_id)
        {
            continue;
        }

        // map boundaryconfig_id to boundaryflux_id
        bcid_to_bfid_map[boundaryconfig_id].push_back(num_boundary_flux_domain);

        // transfer data to flux-type boundary arrays
        num_boundary_flux_domain++;
        boundary_flux_element_gid_vec.push_back(boundary_element_gid_vec[bid]);
        boundary_flux_pa_lid_vec.push_back(boundary_pa_lid_vec[bid]);
        boundary_flux_pa_bcid_vec.push_back(boundary_pa_bcid_vec[bid]);
        boundary_flux_type_str_vec.push_back(type_str);
        boundary_flux_parameter_vec.push_back({});

        // mark as non-constant BC parameters
        boundary_flux_is_parameter_constant_map[boundaryconfig_id] = false;
        boundary_flux_parameter_function_map[boundaryconfig_id] = parameter_function;
        boundary_flux_variable_ptr_map[boundaryconfig_id] = variable_ptr_vec;

    }

}

void BoundaryLine2::set_boundary_value(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec)
{
    /*

    Assigns a BC type and parameters to a BC configuration ID.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    type_str : string
        Type of boundary condition.
    parameter_vec : VectorDouble
        vector with parameters for the BC.

    Returns
    =======
    (none)

    Notes
    ====
    type_str can be "neumann" or "robin" if boundaryconfig_id refers to flux-type BCs.
    type_str can be "dirichlet" if boundaryconfig_id refers to value-type BCs.

    */

    // iterate through each boundary condition
    for (int bid = 0; bid < num_boundary_domain; bid++)
    {

        // skip if different boundary config ID
        if (boundary_pa_bcid_vec[bid] != boundaryconfig_id)
        {
            continue;
        }

        // map boundaryconfig_id to boundaryvalue_id
        bcid_to_bvid_map[boundaryconfig_id].push_back(num_boundary_value_domain);

        // transfer data to flux-type boundary arrays
        num_boundary_value_domain++;
        boundary_value_element_gid_vec.push_back(boundary_element_gid_vec[bid]);
        boundary_value_pa_lid_vec.push_back(boundary_pa_lid_vec[bid]);
        boundary_value_pa_bcid_vec.push_back(boundary_pa_bcid_vec[bid]);
        boundary_value_type_str_vec.push_back(type_str);
        boundary_value_parameter_vec.push_back(parameter_vec);

        // mark as constant BC parameters
        boundary_value_is_parameter_constant_map[boundaryconfig_id] = true;

    }

}

void BoundaryLine2::set_boundary_value(int boundaryconfig_id, std::string type_str, std::function<VectorDouble(double, VectorDouble)> parameter_function, std::vector<VariableLine2*> variable_ptr_vec)
{
    /*

    Assigns a BC type and parameters to a BC configuration ID.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    type_str : string
        Type of boundary condition.
    parameter_vec : VectorDouble
        vector with parameters for the BC.

    Returns
    =======
    (none)

    Notes
    ====
    type_str can be "neumann" or "robin" if boundaryconfig_id refers to flux-type BCs.
    type_str can be "dirichlet" if boundaryconfig_id refers to value-type BCs.

    */

    // iterate through each boundary condition
    for (int bid = 0; bid < num_boundary_domain; bid++)
    {

        // skip if different boundary config ID
        if (boundary_pa_bcid_vec[bid] != boundaryconfig_id)
        {
            continue;
        }

        // map boundaryconfig_id to boundaryvalue_id
        bcid_to_bvid_map[boundaryconfig_id].push_back(num_boundary_value_domain);

        // transfer data to flux-type boundary arrays
        num_boundary_value_domain++;
        boundary_value_element_gid_vec.push_back(boundary_element_gid_vec[bid]);
        boundary_value_pa_lid_vec.push_back(boundary_pa_lid_vec[bid]);
        boundary_value_pa_bcid_vec.push_back(boundary_pa_bcid_vec[bid]);
        boundary_value_type_str_vec.push_back(type_str);
        boundary_value_parameter_vec.push_back({});

        // mark as non-constant BC parameters
        boundary_value_is_parameter_constant_map[boundaryconfig_id] = false;
        boundary_value_parameter_function_map[boundaryconfig_id] = parameter_function;
        boundary_value_variable_ptr_map[boundaryconfig_id] = variable_ptr_vec;

    }

}

void BoundaryLine2::update_parameter()
{

    // iterate through flux-type boundary config
    for (auto boundary_flux_pair : boundary_flux_is_parameter_constant_map)
    {

        // skip if parameters are constant
        if (boundary_flux_pair.second)
        {
            continue;
        }

        // get data associated with boundary config
        int bcid = boundary_flux_pair.first;
        VectorInt bfid_vec = bcid_to_bfid_map[bcid];
        std::function<VectorDouble(double, VectorDouble)> parameter_function = boundary_flux_parameter_function_map[bcid];
        std::vector<VariableLine2*> variable_ptr_vec = boundary_flux_variable_ptr_map[bcid];

        // iterate through each flux-type BC and recompute parameters
        for (auto bfid : bfid_vec)
        {

            // get point domain ID
            int e_gid = boundary_flux_element_gid_vec[bfid];
            int e_did = mesh_ptr->element_gid_to_did_map[e_gid];
            int p0_gid = mesh_ptr->element_p0_gid_vec[e_did];
            int p1_gid = mesh_ptr->element_p1_gid_vec[e_did];
            int pa_lid = boundary_flux_pa_lid_vec[bfid];
            int p_gid_arr[2] = {p0_gid, p1_gid};
            int pa_gid = p_gid_arr[pa_lid];
            int pa_did = mesh_ptr->point_gid_to_did_map[pa_gid];

            // get mesh coordinate
            double position_x = mesh_ptr->point_position_x_vec[pa_did];

            // iterate through each variable that scalar depends on
            VectorDouble value_vec;
            for (auto variable_ptr : variable_ptr_vec)
            {
                value_vec.push_back(variable_ptr->point_value_vec[pa_did]);
            }

            // calculate parameter value
            boundary_flux_parameter_vec[bfid] = parameter_function(position_x, value_vec);

        }
        
    }

    // iterate through value-type boundary config
    for (auto boundary_value_pair : boundary_value_is_parameter_constant_map)
    {

        // skip if parameters are constant
        if (boundary_value_pair.second)
        {
            continue;
        }

        // get data associated with boundary config
        int bcid = boundary_value_pair.first;
        VectorInt bvid_vec = bcid_to_bvid_map[bcid];
        std::function<VectorDouble(double, VectorDouble)> parameter_function = boundary_value_parameter_function_map[bcid];
        std::vector<VariableLine2*> variable_ptr_vec = boundary_value_variable_ptr_map[bcid];

        // iterate through each value-type BC and recompute parameters
        for (auto bvid : bvid_vec)
        {

            // get point domain ID
            int e_gid = boundary_value_element_gid_vec[bvid];
            int e_did = mesh_ptr->element_gid_to_did_map[e_gid];
            int p0_gid = mesh_ptr->element_p0_gid_vec[e_did];
            int p1_gid = mesh_ptr->element_p1_gid_vec[e_did];
            int pa_lid = boundary_value_pa_lid_vec[bvid];
            int p_gid_arr[2] = {p0_gid, p1_gid};
            int pa_gid = p_gid_arr[pa_lid];
            int pa_did = mesh_ptr->point_gid_to_did_map[pa_gid];

            // get mesh coordinate
            double position_x = mesh_ptr->point_position_x_vec[pa_did];

            // iterate through each variable that scalar depends on
            VectorDouble value_vec;
            for (auto variable_ptr : variable_ptr_vec)
            {
                value_vec.push_back(variable_ptr->point_value_vec[pa_did]);
            }

            // calculate parameter value
            boundary_value_parameter_vec[bvid] = parameter_function(position_x, value_vec);

        }
        
    }

}

void BoundaryLine2::read_boundary(std::string file_in_str)
{

    // read file with flux BC data
    std::ifstream file_in_stream(file_in_str);

    // initialize for iteration
    bool is_header = true;  // true while reading header
    std::string line_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_stream, line_str))
    {

        // skip header
        if (is_header)
        {
            is_header = false; // not reading header
            continue;
        }

        // count number of particles
        num_boundary_domain++;

        // convert line string into stringstream
        std::stringstream line_stream(line_str);

        // initialize for iteration
        int value_num = 0;  // counts position of value
        std::string value_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_stream, value_str, ','))
        {

            // store values in appropriate vector
            switch (value_num)
            {
                case 0: boundary_element_gid_vec.push_back(std::stoi(value_str)); break;
                case 1: boundary_pa_lid_vec.push_back(std::stoi(value_str)); break;
                case 2: boundary_pa_bcid_vec.push_back(std::stoi(value_str)); break;
            }

            // increment value count
            value_num++;

        }

    }

    // close point file
    file_in_stream.close();    

}

#endif
