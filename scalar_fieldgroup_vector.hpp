#ifndef SCALAR_FIELDGROUP_VECTOR
#define SCALAR_FIELDGROUP_VECTOR
#include <vector>
#include "scalar_fieldgroup.hpp"

class ScalarFieldGroupVector
{
    /*

    Groups ScalarFieldGroup objects into a column vector.
    Used in physics with systems of partial differential equations.

    Variables
    =========
    scalar_fieldgroup_ptr_vec_in : vector<ScalarFieldGroup*>
        vector of pointers to ScalarFieldGroup objects.
    
    Functions
    =========
    get_entry : ScalarFieldGroup*
        Returns the pointer to the ScalarFieldGroup object at a given vector position.

    */

    public:

    // variables
    int num_entry = 0;
    std::vector<ScalarFieldGroup*> scalar_fieldgroup_ptr_vec;

    // functions
    ScalarFieldGroup* get_entry(int vector_row);

    // default constructor
    ScalarFieldGroupVector()
    {

    }

    // constructor
    ScalarFieldGroupVector(std::vector<ScalarFieldGroup*> scalar_fieldgroup_ptr_vec_in)
    {

        // store vector of variables
        scalar_fieldgroup_ptr_vec = scalar_fieldgroup_ptr_vec_in;

        // count number of variables
        num_entry = scalar_fieldgroup_ptr_vec.size();

    }

};

ScalarFieldGroup* ScalarFieldGroupVector::get_entry(int vector_row)
{
    /*

    Returns the pointer to the ScalarFieldGroup object at a given vector position.

    Arguments
    =========
    vector_row : int
        Row in matrix.
    
    Returns
    =========
    scalar_fieldgroup_ptr : ScalarFieldGroup*
        Pointer to the ScalarFieldGroup object.

    */
    
    return scalar_fieldgroup_ptr_vec[vector_row];

}

#endif
