#ifndef VARIABLE_FIELDGROUP_VECTOR
#define VARIABLE_FIELDGROUP_VECTOR
#include <vector>
#include "variable_fieldgroup.hpp"

class VariableFieldGroupVector
{
    /*

    Groups VariableFieldGroup objects into a column vector.
    Used in physics with systems of partial differential equations.

    Variables
    =========
    variable_fieldgroup_ptr_vec_in : vector<VariableFieldGroup*>
        vector of pointers to VariableFieldGroup objects.
    
    Functions
    =========
    get_entry : VariableFieldGroup*
        Returns the pointer to the VariableFieldGroup object at a given vector position.
    get_vector : vector<VariableFieldGroup*>
        Returns the vector with the VariableFieldGroup objects.

    */

    public:

    // variables
    int num_entry = 0;
    std::vector<VariableFieldGroup*> variable_fieldgroup_ptr_vec;

    // functions
    VariableFieldGroup* get_entry(int vector_row);
    std::vector<VariableFieldGroup*> get_vector();

    // default constructor
    VariableFieldGroupVector()
    {

    }

    // constructor
    VariableFieldGroupVector(std::vector<VariableFieldGroup*> variable_fieldgroup_ptr_vec_in)
    {

        // store vector of variables
        variable_fieldgroup_ptr_vec = variable_fieldgroup_ptr_vec_in;

        // count number of variables
        num_entry = variable_fieldgroup_ptr_vec.size();

    }

};

VariableFieldGroup* VariableFieldGroupVector::get_entry(int vector_row)
{
    /*

    Returns the pointer to the VariableFieldGroup object at a given vector position.

    Arguments
    =========
    vector_row : int
        Row in matrix.
    
    Returns
    =========
    variable_fieldgroup_ptr : VariableFieldGroup*
        Pointer to the VariableFieldGroup object.

    */

    return variable_fieldgroup_ptr_vec[vector_row];

}

std::vector<VariableFieldGroup*> VariableFieldGroupVector::get_vector()
{
    /*

    Returns the vector with the VariableFieldGroup objects.

    Arguments
    =========
    (none)
    
    Returns
    =========
    variable_fieldgroup_ptr_vec : vector<VariableFieldGroup*>
        Vector with the VariableFieldGroup objects.

    */
    
    return variable_fieldgroup_ptr_vec;

}

#endif
