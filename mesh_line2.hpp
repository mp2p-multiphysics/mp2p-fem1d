#ifndef MESH_LINE2
#define MESH_LINE2
#include <unordered_map>
#include "container_typedef.hpp"

class MeshLine2
{
    /*

    Mesh domain made of line2 elements.

    Variables
    =========
    file_in_point_str : string
        Path to CSV file with data for mesh points.
    file_in_element_str : string
        Path to CSV file with data for mesh elements.

    Notes
    ====
    The CSV file with point data must have the following columns:
        global point ID
        x-coordinate of point
    The CSV file with element data must have the following columns:
        global element ID
        global point ID of local point 0
        global point ID of local point 1
    Point 0 and 1 refer to the left and right points of the each element.

    */

    // did - domain ID
    // gid - global ID
    // vectors use did as input

    public:

    // file names
    std::string file_in_point_str;
    std::string file_in_element_str;

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

    // default constructor
    MeshLine2 ()
    {

    }

    // construcotr
    MeshLine2 (std::string file_in_point_str_in, std::string file_in_element_str_in)
    {

        // store variables
        file_in_point_str = file_in_point_str_in;
        file_in_element_str = file_in_element_str_in;

        // read csv files
        read_mesh_point(file_in_point_str);
        read_mesh_element(file_in_element_str);

    }
    
    private:

    // functions
    void read_mesh_point(std::string file_in_point_str);
    void read_mesh_element(std::string file_in_element_str);

};

void MeshLine2::read_mesh_point(std::string file_in_point_str)
{

    // read file with points
    std::ifstream file_in_point_stream(file_in_point_str);

    // initialize for iteration
    bool is_point_header = true;  // true while reading header
    std::string line_point_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_point_stream, line_point_str))
    {

        // skip header
        if (is_point_header)
        {
            is_point_header = false; // not reading header
            continue;
        }

        // count number of particles
        num_point_domain++;

        // convert line string into stringstream
        std::stringstream line_point_stream(line_point_str);

        // initialize for iteration
        int value_point_num = 0;  // counts position of value
        std::string value_point_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_point_stream, value_point_str, ','))
        {

            // store values in appropriate vector
            switch (value_point_num)
            {
                case 0: point_gid_vec.push_back(std::stoi(value_point_str)); break;
                case 1: point_position_x_vec.push_back(std::stod(value_point_str)); break;
            }

            // increment value count
            value_point_num++;

        }

    }

    // close point file
    file_in_point_stream.close();

    // generate map of global to domain ID for points
    for (int point_did = 0; point_did < num_point_domain; point_did++)
    {
        int point_gid = point_gid_vec[point_did];
        point_gid_to_did_map[point_gid] = point_did;
    }

}

void MeshLine2::read_mesh_element(std::string file_in_element_str)
{

    // read file with elements
    std::ifstream file_in_element_stream(file_in_element_str);  

    // initialize for iteration
    bool is_element_header = true;  // true while reading header
    std::string line_element_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_element_stream, line_element_str))
    {

        // skip header
        if (is_element_header)
        {
            is_element_header = false; // not reading header
            continue;
        }

        // count number of particles
        num_element_domain++;

        // convert line string into stringstream
        std::stringstream line_element_stream(line_element_str);

        // initialize for iteration
        int value_element_num = 0;  // counts position of value
        std::string value_element_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_element_stream, value_element_str, ','))
        {

            // store values in appropriate vector
            switch (value_element_num)
            {
                case 0: element_gid_vec.push_back(std::stoi(value_element_str)); break;
                case 1: element_p0_gid_vec.push_back(std::stod(value_element_str)); break;
                case 2: element_p1_gid_vec.push_back(std::stod(value_element_str)); break;
            }

            // increment value count
            value_element_num++;

        }

    }

    // close element file
    file_in_element_stream.close();

    // generate map of global to domain ID for elements
    for (int element_did = 0; element_did < num_element_domain; element_did++)
    {
        int element_gid = element_gid_vec[element_did];
        element_gid_to_did_map[element_gid] = element_did;
    }

}

#endif
