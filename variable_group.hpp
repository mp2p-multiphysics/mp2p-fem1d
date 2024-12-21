#ifndef VARIABLE_GROUP
#define VARIABLE_GROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "container_typedef.hpp"
#include "domain_0d.hpp"
#include "domain_1d.hpp"
#include "variable_1d.hpp"

namespace FEM1D
{

class VariableGroup
{
    /*

    Groups variables that are applied to the same group.

    Variables
    =========
    variable_ptr_vec_in : vector<Variable1D*>
        vector with pointers to Variable1D objects.

    Functions
    =========
    output_csv : void
        Outputs a CSV file with the values of the variable.
    get_neighbor_pfid : VectorInt
        Outputs a vector with the group IDs of the points surrounding an element.

    */

    public:

    // number of unique points in group
    int num_point = 0;

    // point IDs
    VectorInt point_pfid_to_pgid_vec;  // key: group ID; value: global ID
    MapIntInt point_pgid_to_pfid_map;  // key: global ID; value: group ID

    // variables and domains
    std::vector<Variable1D*> variable_ptr_vec;  // vector of variables
    std::unordered_map<Domain1D*, Variable1D*> domain_to_variable_ptr_map;

    // starting column of variables in matrix equation
    int start_col = -1;

    // functions
    void output_csv();
    void output_csv(int ts);
    VectorInt get_neighbor_pfid(Domain0D* domain_ptr, int edid);
    VectorInt get_neighbor_pfid(Domain1D* domain_ptr, int edid);
    VectorDouble get_neighbor_value(Domain1D* domain_ptr, int edid);

    // default constructor
    VariableGroup() {}

    // constructor
    VariableGroup(std::vector<Variable1D*> variable_ptr_vec_in)
    {
        
        // store vector of variables
        variable_ptr_vec = variable_ptr_vec_in;

        // get set of global IDs
        // map global IDs and group IDs

        // initialize set of global IDs
        std::set<int> point_pgid_set;  

        // iterate through each variable
        for (auto variable_ptr : variable_ptr_vec)
        {
            
            // map domain to variables
            domain_to_variable_ptr_map[variable_ptr->domain_ptr] = variable_ptr;

            // get set of global IDs
            for (auto &pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
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
            
            // increment group ID
            pfid++;

        }

        // total number of group points
        num_point = pfid;

    }

};

void VariableGroup::output_csv()
{
    /*

    Outputs a CSV file with the values of the variable.

    Arguments
    =========
    (none)

    Returns
    =======
    (none)

    */

    // iterate through each variable
    for (auto variable_ptr : variable_ptr_vec)
    {
        variable_ptr->output_csv();
    }

}

void VariableGroup::output_csv(int ts)
{
    /*

    Outputs a CSV file with the values of the variable.

    Arguments
    =========
    ts : int
        Timestep number.

    Returns
    =======
    (none)

    */

    // iterate through each variable
    for (auto variable_ptr : variable_ptr_vec)
    {
        variable_ptr->output_csv(ts);
    }

}

VectorInt VariableGroup::get_neighbor_pfid(Domain0D* domain_ptr, int edid)
{
    /*
    
    Outputs a vector with the group IDs of the points surrounding an element.

    Arguments
    =========
    domain_ptr : Domain0D*
        Pointer to Domain0D object with element.
    edid : int
        Domain ID of the element.

    Returns
    =======
    pfid_vec : VectorInt
        vector with group IDs of the points surrounding the element.

    */

    // get surrounding points
    VectorInt pgid_vec = domain_ptr->element_edid_plid_to_pgid_vec[edid];

    // initialize pfid vector
    VectorInt pfid_vec;

    // iterate for each point and get pfid
    for (int pgid : pgid_vec)
    {
        int pfid = point_pgid_to_pfid_map[pgid];
        pfid_vec.push_back(pfid);
    }
    
    return pfid_vec;

}

VectorInt VariableGroup::get_neighbor_pfid(Domain1D* domain_ptr, int edid)
{
    /*
    
    Outputs a vector with the group IDs of the points surrounding an element.

    Arguments
    =========
    domain_ptr : Domain1D*
        Pointer to Domain1D object with element.
    edid : int
        Domain ID of the element.

    Returns
    =======
    pfid_vec : VectorInt
        vector with group IDs of the points surrounding the element.

    */

    // get surrounding points
    VectorInt pgid_vec = domain_ptr->element_edid_plid_to_pgid_vec[edid];

    // initialize pfid vector
    VectorInt pfid_vec;

    // iterate for each point and get pfid
    for (int pgid : pgid_vec)
    {
        int pfid = point_pgid_to_pfid_map[pgid];
        pfid_vec.push_back(pfid);
    }
    
    return pfid_vec;

}

VectorDouble VariableGroup::get_neighbor_value(Domain1D* domain_ptr, int edid)
{
    /*
    
    Outputs a vector with the variable values at points surrounding an element.

    Arguments
    =========
    domain_ptr : Domain1D*
        Pointer to Domain1D object with element.
    edid : int
        Domain ID of the element.

    Returns
    =======
    value_vec : VectorInt
        vector with variable values at points surrounding an element.

    */

    // get variable
    Variable1D *variable_ptr = domain_to_variable_ptr_map[domain_ptr];

    // get surrounding points
    VectorInt pgid_vec = domain_ptr->element_edid_plid_to_pgid_vec[edid];

    // initialize vector with values
    VectorDouble value_vec;

    // iterate for each point and get pfid
    for (int pgid : pgid_vec)
    {
        int pdid = domain_ptr->point_pgid_to_pdid_map[pgid];
        double value_sub = variable_ptr->point_value_vec[pdid];
        value_vec.push_back(value_sub);
    }
    
    return value_vec;

}

}

#endif
