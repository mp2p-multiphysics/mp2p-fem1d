#ifndef BOUNDARY_LINE2
#define BOUNDARY_LINE2
#include <fstream>
#include <sstream>
#include <vector>
#include "container_typedef.hpp"

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
    // bcid - boundary config ID
    // vectors use did as input

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

    // value boundary condition data
    int num_boundary_value_domain = 0;
    VectorInt boundary_value_element_gid_vec;
    VectorInt boundary_value_pa_lid_vec;
    VectorInt boundary_value_pa_bcid_vec;
    std::vector<std::string> boundary_value_type_str_vec;
    Vector2D boundary_value_parameter_vec; 

    // functions
    void set_boundary_flux(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec);
    void set_boundary_value(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec);

    // default constructor
    BoundaryLine2()
    {

    }

    // constructor
    BoundaryLine2(std::string file_in_str_in)
    {

        // store variables
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

        // transfer data to flux-type boundary arrays
        num_boundary_flux_domain++;
        boundary_flux_element_gid_vec.push_back(boundary_element_gid_vec[bid]);
        boundary_flux_pa_lid_vec.push_back(boundary_pa_lid_vec[bid]);
        boundary_flux_pa_bcid_vec.push_back(boundary_pa_bcid_vec[bid]);
        boundary_flux_type_str_vec.push_back(type_str);
        boundary_flux_parameter_vec.push_back(parameter_vec);

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

        // transfer data to flux-type boundary arrays
        num_boundary_value_domain++;
        boundary_value_element_gid_vec.push_back(boundary_element_gid_vec[bid]);
        boundary_value_pa_lid_vec.push_back(boundary_pa_lid_vec[bid]);
        boundary_value_pa_bcid_vec.push_back(boundary_pa_bcid_vec[bid]);
        boundary_value_type_str_vec.push_back(type_str);
        boundary_value_parameter_vec.push_back(parameter_vec);

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
