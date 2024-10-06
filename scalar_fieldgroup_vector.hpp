#ifndef SCALAR_FIELDGROUP_VECTOR
#define SCALAR_FIELDGROUP_VECTOR
#include <vector>
#include "scalar_fieldgroup.hpp"

class ScalarFieldGroupVector
{

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
    return scalar_fieldgroup_ptr_vec[vector_row];
}

#endif
