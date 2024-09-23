#ifndef VARIABLE_FIELDGROUP
#define VARIABLE_FIELDGROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "variable_line2.hpp"

class VariableFieldGroup
{

    public:

    // variables
    int num_field_point = 0;
    std::set<int> point_gid_set;  // set of global IDs
    std::vector<VariableLine2*> variable_ptr_vec;  // vector of variables
    std::unordered_map<MeshLine2Struct*, VariableLine2*> variable_ptr_map;
    std::unordered_map<int, int> point_gid_to_fid_map;  // key: global point id; value: field point id
    VectorDouble point_value_vec;

    // starting column of variables in matrix equation
    int start_col = -1;

    // functions
    void output_csv(std::string file_out_str);

    // default constructor
    VariableFieldGroup()
    {

    }

    // constructor
    VariableFieldGroup(std::vector<VariableLine2*> variable_ptr_vec_in)
    {
        
        // store vector of variables
        variable_ptr_vec = variable_ptr_vec_in;

        // iterate through each variable and get set of global IDs
        for (auto variable_ptr : variable_ptr_vec)
        {
            for (auto &point_gid : variable_ptr->mesh_l2_ptr->point_gid_vec)
            {
                point_gid_set.insert(point_gid);
            }
        }

        // initialize field ID
        int point_field_id = 0;

        // iterate through each global ID and assign a field ID
        for (auto point_gid : point_gid_set)
        {

            // skip if global ID is already recorded
            if (point_gid_to_fid_map.count(point_gid))
            {
                continue;
            }

            // assign field ID and increment ID counter
            point_gid_to_fid_map[point_gid] = point_field_id;
            point_field_id++;

        }

        // number of field points = highest field id
        num_field_point = point_field_id;

        // iterate through each vector
        // map mesh to each vector
        for (auto variable_ptr : variable_ptr_vec)
        {
            variable_ptr_map[variable_ptr->mesh_l2_ptr] = variable_ptr;
        }

    }

};

void VariableFieldGroup::output_csv(std::string file_out_str)
{

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "id,value\n";
    for (auto &point_gid : point_gid_set)
    {
        file_out_stream << point_gid << ",";
        file_out_stream << point_value_vec[point_gid_to_fid_map[point_gid]] << "\n";
    }

}

#endif
