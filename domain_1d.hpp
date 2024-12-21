#ifndef DOMAIN_1D
#define DOMAIN_1D
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "container_typedef.hpp"

namespace FEM1D
{

class Domain1D
{
    /*

    Representation of a 1D domain.

    Variables
    =========
    file_in_point_str_in : string
        Path to CSV file with data for domain points.
    file_in_element_str_in : string
        Path to CSV file with data for domain elements.
    type_element_in : int
        Denotes type of element.
        Use 0 or 1 for line2 or line3.

    Notes
    ====
    The CSV file with element data must have the following columns:
        global element ID
        global point ID of local point 0
        global point ID of local point 1
        ...

    Point 0 and 1 refer to the left and right points of a line2 element.
    Point 0, 1, and 2 refer to the left, middle, and right points of a line3 element.

    */

    public:

    // type of element
    // 0 - line2; 1 - line3
    int type_element = 0;
    int num_neighbor = 0;

    // point data
    int num_point = 0;
    VectorInt point_pdid_to_pgid_vec;
    MapIntInt point_pgid_to_pdid_map;
    VectorDouble point_position_x_vec;

    // element data
    int num_element = 0;
    VectorInt element_edid_to_egid_vec;
    MapIntInt element_egid_to_edid_map;
    VectorInt2D element_edid_plid_to_pgid_vec;  // [edid][plid] -> pgid

    // file names
    std::string file_in_point_str;
    std::string file_in_element_str;

    // default constructor
    Domain1D () {}

    // construcotr
    Domain1D(std::string file_in_point_str_in, std::string file_in_element_str_in, int type_element_in)
    {

        // store variables
        file_in_point_str = file_in_point_str_in;
        file_in_element_str = file_in_element_str_in;
        type_element = type_element_in;

        // get number of neighboring points per element
        switch (type_element)
        {
            case 0: num_neighbor = 2; break;  // line2
            case 1: num_neighbor = 3; break;  // line3
        }

        // read csv files
        read_domain_point();
        read_domain_element();

    }
    
    private:

    // functions
    void read_domain_point();
    void read_domain_element();

};

void Domain1D::read_domain_point()
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

void Domain1D::read_domain_element()
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
        VectorInt pgid_vec_sub;  // vector of global point IDs

        // iterate through each value
        while (std::getline(line_element_stream, value_element_str, ','))
        {

            // store values in appropriate vector
            if (value_element_num == 0)
            {
                element_edid_to_egid_vec.push_back(std::stoi(value_element_str));
            }
            else if (value_element_num < num_neighbor + 1)
            {
                pgid_vec_sub.push_back(std::stoi(value_element_str));
            }

            // store vector of global point IDs
            if (value_element_num == num_neighbor)
            {
                element_edid_plid_to_pgid_vec.push_back(pgid_vec_sub);
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

}

#endif
