#ifndef VARIABLE_FIELDGROUP
#define VARIABLE_FIELDGROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "variable_line2.hpp"

class VariableFieldGroup
{
    /*

    Groups variables that are applied to the same field.

    Variables
    =========
    variable_ptr_vec_in : vector<VariableLine2*>
        vector with pointers to VariableLine2 objects.
    
    */

    public:

    // number of unique points in field
    int num_point_field = 0;

    // point IDs
    VectorInt point_gid_vec;  // key: field ID; value: global ID
    MapIntInt point_gid_to_fid_map;  // key: global ID; value: field ID

    // variables and meshes
    std::vector<VariableLine2*> variable_ptr_vec;  // vector of variables
    std::unordered_map<MeshLine2Struct*, VariableLine2*> mesh_to_variable_ptr_map;  // key: mesh; value: variable
   
    // starting column of variables in matrix equation
    int start_col = -1;

    // default constructor
    VariableFieldGroup()
    {

    }

    // constructor
    VariableFieldGroup(std::vector<VariableLine2*> variable_ptr_vec_in)
    {
        
        // store vector of variables
        variable_ptr_vec = variable_ptr_vec_in;

        // map mesh to variables
        for (auto variable_ptr : variable_ptr_vec)
        {
            mesh_to_variable_ptr_map[variable_ptr->mesh_l2_ptr] = variable_ptr;
        }

        // get set of global IDs
        // map global IDs and field IDs

        // initialize set of global IDs
        std::set<int> point_gid_set;  

        // iterate through each variable and get set of global IDs
        for (auto variable_ptr : variable_ptr_vec)
        {
            for (auto &point_gid : variable_ptr->mesh_l2_ptr->point_gid_vec)
            {
                point_gid_set.insert(point_gid);
            }
        }

        // initialize field ID
        int point_fid = 0;

        // iterate through each global ID and assign a field ID
        for (auto point_gid : point_gid_set)
        {

            // skip if global ID is already recorded
            if (point_gid_to_fid_map.count(point_gid))
            {
                continue;
            }

            // map global ID to field ID and vice versa
            point_gid_to_fid_map[point_gid] = point_fid;
            point_gid_vec.push_back(point_gid);
            
            // increment field ID
            point_fid++;

        }

        // total number of field points
        num_point_field = point_fid;

    }

};

#endif
