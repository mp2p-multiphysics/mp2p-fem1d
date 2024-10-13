#ifndef SCALAR_FIELDGROUP_MATRIX
#define SCALAR_FIELDGROUP_MATRIX
#include <unordered_map>
#include <vector>
#include "scalar_fieldgroup.hpp"

class ScalarFieldGroupMatrix
{
    /*

    Groups ScalarFieldGroup objects into a square matrix.
    Used in physics with systems of partial differential equations.

    Variables
    =========
    num_entry_in : int
        Number of rows or columns in the matrix.
    
    Functions
    =========
    get_entry : ScalarFieldGroup*
        Returns the pointer to the ScalarFieldGroup object at a given matrix position.
    get_row : std::unordered_map<int, ScalarFieldGroup*>
        Returns the entries in a given row.
    set_entry : void
        Assigns a ScalarFieldGroup object to a given matrix position.

    */

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
    /*

    Returns the pointer to the ScalarFieldGroup object at a given matrix position.

    Arguments
    =========
    matrix_row : int
        Row in matrix.
    matrix_column : int
        Column in matrix
    
    Returns
    =========
    scalar_fieldgroup_ptr : ScalarFieldGroup*
        Pointer to the ScalarFieldGroup object.

    */
    
    return scalar_fieldgroup_ptr_mat[matrix_row][matrix_column];

}

std::unordered_map<int, ScalarFieldGroup*> ScalarFieldGroupMatrix::get_row(int matrix_row)
{
    /*

    Returns the entries in a given row.

    Arguments
    =========
    matrix_row : int
        Row in matrix.
    
    Returns
    =========
    scalar_fieldgroup_ptr_mat : map<int, ScalarFieldGroup*>
        map of pointers to ScalarFieldGroup objects in the specified row.
        The keys in the map are the columns in the matrix.

    */
    
    return scalar_fieldgroup_ptr_mat[matrix_row];

}

void ScalarFieldGroupMatrix::set_entry(int matrix_row, int matrix_column, ScalarFieldGroup *scalar_fieldgroup_ptr)
{
    /*

    Assigns a ScalarFieldGroup object to a given matrix position.

    Arguments
    =========
    matrix_row : int
        Row in matrix.
    matrix_column : int
        Column in matrix
    scalar_fieldgroup_ptr : ScalarFieldGroup*
        Pointer to the ScalarFieldGroup object.

    Returns
    =========
    (none)

    */
    
    scalar_fieldgroup_ptr_mat[matrix_row][matrix_column] = scalar_fieldgroup_ptr;

}

#endif
