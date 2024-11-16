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
    VectorInt point_pgid_vec;  // key: group ID; value: global ID
    MapIntInt point_pgid_to_pfid_map;  // key: global ID; value: group ID

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
            for (auto &pgid : variable_ptr->domain_ptr->point_pgid_vec)
            {
                point_gid_set.insert(pgid);
            }
        }

        // initialize group ID
        int pfid = 0;

        // iterate through each global ID and assign a group ID
        for (auto pgid : point_gid_set)
        {

            // skip if global ID is already recorded
            if (point_pgid_to_pfid_map.count(pgid))
            {
                continue;
            }

            // map global ID to group ID and vice versa
            point_pgid_to_pfid_map[pgid] = pfid;
            point_pgid_vec.push_back(pgid);
            
            // increment group ID
            pfid++;

        }

        // total number of group points
        num_point_group = pfid;

    }

};

#endif
