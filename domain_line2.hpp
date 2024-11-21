#ifndef DOMAIN_LINE2
#define DOMAIN_LINE2
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "container_typedef.hpp"

class DomainLine2
{
    /*

    Domain made of line2 elements.

    Variables
    =========
    file_in_point_str : string
        Path to CSV file with data for domain points.
    file_in_element_str : string
        Path to CSV file with data for domain elements.

    Notes
    ====
    The CSV file with point data must have the following columns:
        global point ID
        x-coordinate of point
    The CSV file with element data must have the following columns:
        global element ID
        global point ID of local point 0
        global point ID of local point 1
    Point 0 and 1 refer to the left and right points of each element.

    */

    // pdid - point domain ID
    // pgid - point global ID
    // vectors use pdid as input

    public:

    // file names
    std::string file_in_point_str;
    std::string file_in_element_str;

    // point data
    int num_point = 0;
    VectorInt point_pdid_to_pgid_vec;
    MapIntInt point_pgid_to_pdid_map;
    VectorDouble point_position_x_vec;

    // element data
    int num_element = 0;
    VectorInt element_edid_to_egid_vec;
    MapIntInt element_egid_to_edid_map;
    VectorInt element_p0_pgid_vec;
    VectorInt element_p1_pgid_vec;

    // default constructor
    DomainLine2() {}

    // constructor
    DomainLine2(std::string file_in_point_str_in, std::string file_in_element_str_in)
    {

        // store variables
        file_in_point_str = file_in_point_str_in;
        file_in_element_str = file_in_element_str_in;

        // read csv files
        read_domain_point(file_in_point_str);
        read_domain_element(file_in_element_str);

    }
    
    private:

    // functions
    void read_domain_point(std::string file_in_point_str);
    void read_domain_element(std::string file_in_element_str);

};

void DomainLine2::read_domain_point(std::string file_in_point_str)
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

        // count number of points
        num_point++;

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
                case 0: point_pdid_to_pgid_vec.push_back(std::stoi(value_point_str)); break;
                case 1: point_position_x_vec.push_back(std::stod(value_point_str)); break;
            }

            // increment value count
            value_point_num++;

        }

    }

    // close point file
    file_in_point_stream.close();

    // generate map of global to domain ID for points
    for (int pdid = 0; pdid < num_point; pdid++)
    {
        int pgid = point_pdid_to_pgid_vec[pdid];
        point_pgid_to_pdid_map[pgid] = pdid;
    }

}

void DomainLine2::read_domain_element(std::string file_in_element_str)
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

        // count number of elements
        num_element++;

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
                case 0: element_edid_to_egid_vec.push_back(std::stoi(value_element_str)); break;
                case 1: element_p0_pgid_vec.push_back(std::stod(value_element_str)); break;
                case 2: element_p1_pgid_vec.push_back(std::stod(value_element_str)); break;
            }

            // increment value count
            value_element_num++;

        }

    }

    // close element file
    file_in_element_stream.close();

    // generate map of global to domain ID for elements
    for (int edid = 0; edid < num_element; edid++)
    {
        int egid = element_edid_to_egid_vec[edid];
        element_egid_to_edid_map[egid] = edid;
    }

}

#endif
