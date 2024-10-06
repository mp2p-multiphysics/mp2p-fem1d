#ifndef VARIABLE_FIELDGROUP_VECTOR
#define VARIABLE_FIELDGROUP_VECTOR
#include <vector>
#include "variable_fieldgroup.hpp"

class VariableFieldGroupVector
{

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
    return variable_fieldgroup_ptr_vec[vector_row];
}

std::vector<VariableFieldGroup*> VariableFieldGroupVector::get_vector()
{
    return variable_fieldgroup_ptr_vec;
}

#endif
