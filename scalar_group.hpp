#ifndef SCALAR_GROUP
#define SCALAR_GROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "container_typedef.hpp"
#include "domain_line2.hpp"
#include "scalar_line2.hpp"

namespace FEM1D
{

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
    output_csv : void
        Outputs a CSV file with the values of the scalar.
    update_value : void
        Recalculates non-constant values.

    */

    public:

    // number of unique points in group
    int num_point_group = 0;

    // point IDs
    VectorInt point_pfid_to_pgid_vec;  // key: group ID; value: global ID
    MapIntInt point_pgid_to_pfid_map;  // key: global ID; value: group ID

    // scalars and domains
    std::vector<ScalarLine2*> scalar_l2_ptr_vec;  // vector of scalars

    // functions
    void output_csv();
    void output_csv(int ts);
    void update_value();

    // default constructor
    ScalarGroup() {}

    // constructor
    ScalarGroup(std::vector<ScalarLine2*> scalar_l2_ptr_vec_in)
    {
        
        // store vector of scalars
        scalar_l2_ptr_vec = scalar_l2_ptr_vec_in;

        // get set of global IDs
        // map global IDs and group IDs

        // initialize set of global IDs
        std::set<int> point_pgid_set;  

        // iterate through each variable and get set of global IDs
        for (auto scalar_ptr : scalar_l2_ptr_vec)
        {
            for (auto &pgid : scalar_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {
                point_pgid_set.insert(pgid);
            }
        }

        // initialize group ID
        int pfid = 0;

        // iterate through each global ID and assign a group ID
        for (auto pgid : point_pgid_set)
        {

            // skip if global ID is already recorded
            if (point_pgid_to_pfid_map.count(pgid))
            {
                continue;
            }

            // map global ID to group ID and vice versa
            point_pgid_to_pfid_map[pgid] = pfid;
            point_pfid_to_pgid_vec.push_back(pgid);
            pfid++;

        }

        // total number of group points
        num_point_group = pfid;

    }

};

void ScalarGroup::output_csv()
{
    /*

    Outputs a CSV file with the values of the scalar.

    Arguments
    =========
    (none)

    Returns
    =======
    (none)

    */

    // iterate through each scalar
    for (auto scalar_ptr : scalar_l2_ptr_vec)
    {
        scalar_ptr->output_csv();
    }

}

void ScalarGroup::output_csv(int ts)
{
    /*

    Outputs a CSV file with the values of the scalar.

    Arguments
    =========
    ts : int
        Timestep number.

    Returns
    =======
    (none)

    */

    // iterate through each scalar
    for (auto scalar_ptr : scalar_l2_ptr_vec)
    {
        scalar_ptr->output_csv(ts);
    }

}

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

}

#endif
