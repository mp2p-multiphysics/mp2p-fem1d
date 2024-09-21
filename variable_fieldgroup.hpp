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
    std::vector<VariableLine2*> variable_ptr_vec;  // vector of variables
    std::unordered_map<int, int> point_global_field_id_map;  // key: global point id; value: field point id

    // starting column in matrix equation
    int start_col = -1;

    // functions

    // default constructor
    VariableFieldGroup()
    {

    }

    // constructor
    VariableFieldGroup(std::vector<VariableLine2*> variable_ptr_vec_in)
    {
        
        // store vector of variables
        variable_ptr_vec = variable_ptr_vec_in;

        // initialize set of global IDs
        std::set<int> point_global_id_set;

        // iterate through each variable and get set of global IDs
        for (auto variable_ptr : variable_ptr_vec)
        {
            for (auto &point_global_id : variable_ptr->mesh_l2_ptr->point_global_id_vec)
            {
                point_global_id_set.insert(point_global_id);
            }
        }

        // initialize field ID
        int point_field_id = 0;

        // iterate through each global ID and assign a field ID
        for (auto point_global_id : point_global_id_set)
        {

            // skip if global ID is already recorded
            if (point_global_field_id_map.count(point_global_id))
            {
                continue;
            }

            // assign field ID and increment ID counter
            point_global_field_id_map[point_global_id] = point_field_id;
            point_field_id++;

        }

        // number of field points = highest field id
        num_field_point = point_field_id;

    }

};

#endif
