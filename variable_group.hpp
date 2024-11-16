#ifndef VARIABLE_FIELDGROUP
#define VARIABLE_FIELDGROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "variable_line2.hpp"

class VariableGroup
{
    /*

    Groups variables that are applied to the same group.

    Variables
    =========
    variable_l2_ptr_vec_in : vector<VariableLine2*>
        vector with pointers to VariableLine2 objects.
    
    */

    public:

    // number of unique points in group
    int num_point_group = 0;

    // point IDs
    VectorInt point_gid_vec;  // key: group ID; value: global ID
    MapIntInt point_gid_to_fid_map;  // key: global ID; value: group ID

    // variables and domains
    std::vector<VariableLine2*> variable_l2_ptr_vec;  // vector of variables
    std::unordered_map<DomainLine2*, VariableLine2*> domain_to_variable_ptr_map;  // key: domain; value: variable
   
    // starting column of variables in matrix equation
    int start_col = -1;

    // default constructor
    VariableGroup()
    {

    }

    // constructor
    VariableGroup(std::vector<VariableLine2*> variable_l2_ptr_vec_in)
    {
        
        // store vector of variables
        variable_l2_ptr_vec = variable_l2_ptr_vec_in;

        // map domain to variables
        for (auto variable_ptr : variable_l2_ptr_vec)
        {
            domain_to_variable_ptr_map[variable_ptr->domain_ptr] = variable_ptr;
        }

        // get set of global IDs
        // map global IDs and group IDs

        // initialize set of global IDs
        std::set<int> point_gid_set;  

        // iterate through each variable and get set of global IDs
        for (auto variable_ptr : variable_l2_ptr_vec)
        {
            for (auto &point_gid : variable_ptr->domain_ptr->point_gid_vec)
            {
                point_gid_set.insert(point_gid);
            }
        }

        // initialize group ID
        int point_fid = 0;

        // iterate through each global ID and assign a group ID
        for (auto point_gid : point_gid_set)
        {

            // skip if global ID is already recorded
            if (point_gid_to_fid_map.count(point_gid))
            {
                continue;
            }

            // map global ID to group ID and vice versa
            point_gid_to_fid_map[point_gid] = point_fid;
            point_gid_vec.push_back(point_gid);
            
            // increment group ID
            point_fid++;

        }

        // total number of group points
        num_point_group = point_fid;

    }

};

#endif
