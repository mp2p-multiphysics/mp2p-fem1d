#ifndef INITIALIZE_MESH_LINE2_CSV
#define INITIALIZE_MESH_LINE2_CSV
#include <fstream>
#include <sstream>
#include <vector>
#include "mesh_line2.hpp"

MeshLine2Struct initialize_mesh_line2_csv(std::string file_in_point_str, std::string file_in_element_str)
{

    // read file with points
    std::ifstream file_in_point_stream(file_in_point_str);

    // initialize struct with mesh data
    MeshLine2Struct mesh_l2;

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
        mesh_l2.num_domain_point++;

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
                case 0: mesh_l2.point_gid_vec.push_back(std::stoi(value_point_str)); break;
                case 1: mesh_l2.point_position_x_vec.push_back(std::stod(value_point_str)); break;
            }

            // increment value count
            value_point_num++;

        }

    }

    // close point file
    file_in_point_stream.close();

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
        mesh_l2.num_domain_element++;

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
                case 0: mesh_l2.element_gid_vec.push_back(std::stoi(value_element_str)); break;
                case 1: mesh_l2.element_p0_gid_vec.push_back(std::stod(value_element_str)); break;
                case 2: mesh_l2.element_p1_gid_vec.push_back(std::stod(value_element_str)); break;
            }

            // increment value count
            value_element_num++;

        }

    }

    // close element file
    file_in_element_stream.close();

    // generate map of global to domain ID for points
    for (int p_did = 0; p_did < mesh_l2.num_domain_point; p_did++)
    {
        int p_gid = mesh_l2.point_gid_vec[p_did];
        mesh_l2.point_gid_to_did_map[p_gid] = p_did;
    }

    // generate map of global to domain ID for elements
    for (int e_did = 0; e_did < mesh_l2.num_domain_point; e_did++)
    {
        int e_gid = mesh_l2.element_gid_vec[e_did];
        mesh_l2.element_gid_to_did_map[e_gid] = e_did;
    }

    return mesh_l2;

}

#endif
