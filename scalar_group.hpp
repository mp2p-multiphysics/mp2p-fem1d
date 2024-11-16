#ifndef SCALAR_FIELDGROUP
#define SCALAR_FIELDGROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "domain_line2.hpp"
#include "scalar_line2.hpp"

class ScalarGroup
{
    /*

    Groups scalars that are applied to the same group.

    Variables
    =========
    scalar_l2_ptr_vec_in : vector<ScalarLine2*>
        vector with pointers to ScalarLine2 objects.

    Functions
    =========
    update_value : void
        Recalculates non-constant values.

    */

    public:

    // number of unique points in group
    int num_point_group = 0;

    // point IDs
    VectorInt point_gid_vec;  // key: group ID; value: global ID
    MapIntInt point_gid_to_fid_map;  // key: global ID; value: group ID

    // scalars and domains
    std::vector<ScalarLine2*> scalar_l2_ptr_vec;  // vector of scalars
    std::unordered_map<DomainLine2*, ScalarLine2*> scalar_ptr_map;  // key: domain; value: scalar

    // functions
    void update_value();

    // default constructor
    ScalarGroup()
    {

    }

    // constructor
    ScalarGroup(std::vector<ScalarLine2*> scalar_l2_ptr_vec_in)
    {
        
        // store vector of scalars
        scalar_l2_ptr_vec = scalar_l2_ptr_vec_in;

        // map domain to scalars
        for (auto scalar_ptr : scalar_l2_ptr_vec)
        {
            scalar_ptr_map[scalar_ptr->domain_ptr] = scalar_ptr;
        }

        // get set of global IDs
        // map global IDs and group IDs

        // initialize set of global IDs
        std::set<int> point_gid_set;  

        // iterate through each variable and get set of global IDs
        for (auto scalar_ptr : scalar_l2_ptr_vec)
        {
            for (auto &point_gid : scalar_ptr->domain_ptr->point_gid_vec)
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
            point_fid++;

        }

        // total number of group points
        num_point_group = point_fid;

    }

};

void ScalarGroup::update_value()
{
    /*

    Recalculates non-constant values.

    Arguments
    =========
    (none)
    
    Returns
    =========
    (none)

    */

    // iterate through each scalar
    for (auto scalar_ptr : scalar_l2_ptr_vec)
    {
        scalar_ptr->update_value();
    }

}

#endif
