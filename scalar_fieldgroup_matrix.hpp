#ifndef SCALAR_FIELDGROUP_MATRIX
#define SCALAR_FIELDGROUP_MATRIX
#include <unordered_map>
#include <vector>
#include "scalar_fieldgroup.hpp"

class ScalarFieldGroupMatrix
{

    public:

    // variables
    int num_entry = 0;
    std::vector<std::unordered_map<int, ScalarFieldGroup*>> scalar_fieldgroup_ptr_mat;

    // function
    ScalarFieldGroup* get_entry(int matrix_row, int matrix_column);
    std::unordered_map<int, ScalarFieldGroup*> get_row(int matrix_row);
    void set_entry(int matrix_row, int matrix_column, ScalarFieldGroup *scalar_fieldgroup_ptr);

    // default constructor
    ScalarFieldGroupMatrix()
    {

    }

    // constructor
    ScalarFieldGroupMatrix(std::vector<std::unordered_map<int, ScalarFieldGroup*>> scalar_fieldgroup_ptr_mat_in)
    {

        // store vector of variables
        scalar_fieldgroup_ptr_mat = scalar_fieldgroup_ptr_mat_in;

        // count number of variables
        num_entry = scalar_fieldgroup_ptr_mat.size();

    }

    ScalarFieldGroupMatrix(int num_entry_in)
    {

        // store number of variables
        num_entry = num_entry_in;

        // allocate vector
        scalar_fieldgroup_ptr_mat = std::vector<std::unordered_map<int, ScalarFieldGroup*>> (num_entry);

    }

};

ScalarFieldGroup* ScalarFieldGroupMatrix::get_entry(int matrix_row, int matrix_column)
{
    return scalar_fieldgroup_ptr_mat[matrix_row][matrix_column];
}

std::unordered_map<int, ScalarFieldGroup*> ScalarFieldGroupMatrix::get_row(int matrix_row)
{
    return scalar_fieldgroup_ptr_mat[matrix_row];
}

void ScalarFieldGroupMatrix::set_entry(int matrix_row, int matrix_column, ScalarFieldGroup *scalar_fieldgroup_ptr)
{
    scalar_fieldgroup_ptr_mat[matrix_row][matrix_column] = scalar_fieldgroup_ptr;
}

#endif
