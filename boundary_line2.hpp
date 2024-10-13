#ifndef BOUNDARY_LINE2
#define BOUNDARY_LINE2
#include <fstream>
#include <sstream>
#include <vector>
#include "container_typedef.hpp"

struct BoundaryConfigLine2Struct
{

    std::string type_str;
    VectorDouble parameter_vec;

};

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
        local element ID where BC is applied
        BC configuration ID
    Flux-type BCs add additional terms to the linearized equations (e.g., Neumann, Robin)
    Value-type BCs completely replace the linearized equations (e.g., Dirichlet)

    */

    public:

    // file names
    std::string file_in_flux_str;
    std::string file_in_value_str;

    // flux boundary condition data
    int num_element_flux_domain = 0;
    VectorInt element_flux_gid_vec;
    VectorInt element_flux_pa_lid_vec;
    VectorInt element_flux_boundaryconfig_id_vec;

    // value boundary condition data
    int num_element_value_domain = 0;
    VectorInt element_value_gid_vec;
    VectorInt element_value_pa_lid_vec;
    VectorInt element_value_boundaryconfig_id_vec;

    // boundary condition type
    int num_boundaryconfig = 0;
    std::vector<BoundaryConfigLine2Struct> boundaryconfig_vec;

    // functions
    void set_boundarycondition(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec);
    void set_boundarycondition_parameter(int boundaryconfig_id, VectorDouble parameter_vec);

    // default constructor
    BoundaryLine2()
    {

    }

    // constructor
    BoundaryLine2(std::string file_in_flux_str_in, std::string file_in_value_str_in)
    {

        // store variables
        file_in_flux_str = file_in_flux_str_in;
        file_in_value_str = file_in_value_str_in;

        // read input files and store boundary condition data
        read_boundary_flux(file_in_flux_str);
        read_boundary_value(file_in_value_str);

        // get largest boundary config id from flux BC
        for (auto boundaryconfig_id : element_flux_boundaryconfig_id_vec)
        {
            if (boundaryconfig_id > num_boundaryconfig)
            {
                num_boundaryconfig = boundaryconfig_id;
            }
        }

        // get largest boundary config id from value BC
        for (auto boundaryconfig_id : element_value_boundaryconfig_id_vec)
        {
            if (boundaryconfig_id > num_boundaryconfig)
            {
                num_boundaryconfig = boundaryconfig_id;
            }
        }

        // assume that first boundary config id is zero
        // add one to get number of boundary config
        num_boundaryconfig += 1;

        // initialize boundary config id vector
        boundaryconfig_vec = std::vector<BoundaryConfigLine2Struct>(num_boundaryconfig);

        // initialize boundary config id vector with zero flux
        BoundaryConfigLine2Struct boundaryconfig_zeroflux;
        boundaryconfig_zeroflux.type_str = "neumann";
        boundaryconfig_zeroflux.parameter_vec = {0};
        for (auto boundaryconfig_id : element_flux_boundaryconfig_id_vec)
        {
            boundaryconfig_vec[boundaryconfig_id] = boundaryconfig_zeroflux;
        }

        // initialize boundary config id vector with zero values
        BoundaryConfigLine2Struct boundaryconfig_zerovalue;
        boundaryconfig_zerovalue.type_str = "dirichlet";
        boundaryconfig_zerovalue.parameter_vec = {0};
        for (auto boundaryconfig_id : element_value_boundaryconfig_id_vec)
        {
            boundaryconfig_vec[boundaryconfig_id] = boundaryconfig_zerovalue;
        }

    }

    private:
    void read_boundary_flux(std::string file_in_flux_str);
    void read_boundary_value(std::string file_in_value_str);

};

void BoundaryLine2::set_boundarycondition(int boundaryconfig_id, std::string type_str, VectorDouble parameter_vec)
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

    // modify struct properties
    boundaryconfig_vec[boundaryconfig_id].type_str = type_str;
    boundaryconfig_vec[boundaryconfig_id].parameter_vec = parameter_vec;

}

void BoundaryLine2::set_boundarycondition_parameter(int boundaryconfig_id, VectorDouble parameter_vec)
{
    /*

    Assigns or modifies the parameters to a BC.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    parameter_vec : VectorDouble
        vector with parameters for the BC.

    Returns
    =======
    (none)

    */

    // modify struct properties
    boundaryconfig_vec[boundaryconfig_id].parameter_vec = parameter_vec;

}

void BoundaryLine2::read_boundary_flux(std::string file_in_flux_str)
{

    // read file with flux BC data
    std::ifstream file_in_flux_stream(file_in_flux_str);

    // initialize for iteration
    bool is_flux_header = true;  // true while reading header
    std::string line_flux_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_flux_stream, line_flux_str))
    {

        // skip header
        if (is_flux_header)
        {
            is_flux_header = false; // not reading header
            continue;
        }

        // count number of particles
        num_element_flux_domain++;

        // convert line string into stringstream
        std::stringstream line_flux_stream(line_flux_str);

        // initialize for iteration
        int value_flux_num = 0;  // counts position of value
        std::string value_flux_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_flux_stream, value_flux_str, ','))
        {

            // store values in appropriate vector
            switch (value_flux_num)
            {
                case 0: element_flux_gid_vec.push_back(std::stoi(value_flux_str)); break;
                case 1: element_flux_pa_lid_vec.push_back(std::stoi(value_flux_str)); break;
                case 2: element_flux_boundaryconfig_id_vec.push_back(std::stoi(value_flux_str)); break;
            }

            // increment value count
            value_flux_num++;

        }

    }

    // close point file
    file_in_flux_stream.close();    

}

void BoundaryLine2::read_boundary_value(std::string file_in_value_str)
{

    // read file with value BC data
    std::ifstream file_in_value_stream(file_in_value_str);

    // initialize for iteration
    bool is_value_header = true;  // true while reading header
    std::string line_value_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_value_stream, line_value_str))
    {

        // skip header
        if (is_value_header)
        {
            is_value_header = false; // not reading header
            continue;
        }

        // count number of particles
        num_element_value_domain++;

        // convert line string into stringstream
        std::stringstream line_value_stream(line_value_str);

        // initialize for iteration
        int value_value_num = 0;  // counts position of value
        std::string value_value_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_value_stream, value_value_str, ','))
        {

            // store values in appropriate vector
            switch (value_value_num)
            {
                case 0: element_value_gid_vec.push_back(std::stoi(value_value_str)); break;
                case 1: element_value_pa_lid_vec.push_back(std::stoi(value_value_str)); break;
                case 2: element_value_boundaryconfig_id_vec.push_back(std::stoi(value_value_str)); break;
            }

            // increment value count
            value_value_num++;

        }

    }

    // close point file
    file_in_value_stream.close();

}

#endif
