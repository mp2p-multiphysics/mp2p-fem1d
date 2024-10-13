#ifndef MATRIXEQUATION_STEADY
#define MATRIXEQUATION_STEADY
#include <set>
#include <vector>
#include "Eigen/Eigen"
#include "physicssteady_base.hpp"
#include "variable_fieldgroup.hpp"

class MatrixEquationSteady
{
    /*

    Represents the matrix equation (system of linear equations) Ax = b for use in steady-state problems.

    Variables
    =========
    physics_ptr_vec_in : vector<PhysicsSteadyBase*>
        vector with steady-state physics classes.

    Functions
    =========
    iterate_solution : void
        Solves for x in Ax = b.
        Uses the currently stored values in A and b.
    store_solution : void
        Transfers the solutions in x into variable objects.

    Notes
    =====
    In the code; A, x, and b are referred to as a_mat, x_vec, and b_vec respectively.

    */

    public:

    // vector of physics
    std::vector<PhysicsSteadyBase*> physics_ptr_vec;
    std::vector<VariableFieldGroup*> variable_field_ptr_vec;

    // matrix equation variables
    Eigen::SparseMatrix<double> a_mat;
    Eigen::VectorXd b_vec;
    Eigen::VectorXd x_vec;
    int num_equation = 0;

    // functions
    void iterate_solution();
    void store_solution();

    // default constructor
    MatrixEquationSteady()
    {

    }

    // constructor
    MatrixEquationSteady(std::vector<PhysicsSteadyBase*> physics_ptr_vec_in)
    {

        // store vector of pointers to physics
        physics_ptr_vec = physics_ptr_vec_in;

        // generate starting rows and columns
        // assign starting row in matrix to test functions (physics)
        // assign starting column in matrix to variables

        // initialize starting rows and columns
        int assign_start_row = 0;
        int assign_start_col = 0;

        // iterate through each physics
        for (auto physics_ptr : physics_ptr_vec)
        {

            // iterate through each variable field
            for (auto variable_field_ptr : physics_ptr->get_variable_field_ptr_vec())
            {
                
                // assign starting column to variable if none yet
                // increment assign_start_col by number of mesh points
                if (variable_field_ptr->start_col == -1)
                {
                    variable_field_ptr->start_col = assign_start_col;
                    assign_start_col += variable_field_ptr->num_point_field;
                }

                // assign starting row to physics if none yet
                // increment assign_start_row by number of new mesh points
                if (physics_ptr->get_start_row() == -1)
                {
                    physics_ptr->set_start_row(assign_start_row);
                    assign_start_row = assign_start_col;
                }

            }

        }

        // get number of linear equations (total number of mesh points)
        num_equation = assign_start_col;

        // get vector of variable fields
        
        // initialize set of variable fields
        std::set<VariableFieldGroup*> variable_field_ptr_set;

        // iterate through each physics
        for (auto physics_ptr : physics_ptr_vec)
        {

            // iterate through each variable field
            for (auto variable_field_ptr : physics_ptr->get_variable_field_ptr_vec())
            {
                
                // store variable field
                variable_field_ptr_set.insert(variable_field_ptr);

            }

        }

        // convert to vector
        variable_field_ptr_vec = std::vector<VariableFieldGroup*>(variable_field_ptr_set.begin(), variable_field_ptr_set.end());

        // initialize matrix equation variables
        a_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
        b_vec = Eigen::VectorXd::Zero(num_equation);
        x_vec = Eigen::VectorXd::Zero(num_equation);
        
        // populate x_vec with initial values

        // iterate through each variable field
        for (auto variable_field_ptr : variable_field_ptr_vec)
        {

            // get starting row
            // note: column in a_mat = row in x_vec
            int start_row = variable_field_ptr->start_col;

            // iterate through each variable
            for (auto variable_ptr : variable_field_ptr->variable_ptr_vec)
            {

                // iterate through each global ID
                for (auto point_gid : variable_ptr->mesh_l2_ptr->point_gid_vec)
                {

                    // get domain and field IDs
                    int point_fid = variable_field_ptr->point_gid_to_fid_map[point_gid];
                    int point_did = variable_ptr->mesh_l2_ptr->point_gid_to_did_map[point_gid];

                    // get value from variable
                    double value = variable_ptr->point_value_vec[point_did];

                    // store value in x_vec
                    int vec_row = start_row + point_fid;
                    x_vec.coeffRef(vec_row) = value;

                }

            }

        }

    }

};

void MatrixEquationSteady::iterate_solution()
{
    /*
    
    Solves for x in Ax = b.
    Uses the currently stored values in A and b.

    Arguments
    =========
    (none)

    Returns
    =======
    (none)

    */

    // fill up a_mat and b_vec with each physics
    for (auto physics_ptr : physics_ptr_vec)
    {
        physics_ptr->matrix_fill(a_mat, b_vec, x_vec);

    }

    // solve the matrix equation
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(a_mat);
    solver.factorize(a_mat);
    x_vec = solver.solve(b_vec);

}

void MatrixEquationSteady::store_solution()
{
    /*
    
    Transfers the solutions in x into variable objects.

    Arguments
    =========
    (none)

    Returns
    =======
    (none)

    */

    // iterate through each variable field
    for (auto variable_field_ptr : variable_field_ptr_vec)
    {

        // get starting row
        // note: column in a_mat = row in x_vec
        int start_row = variable_field_ptr->start_col;

        // iterate through each variable
        for (auto variable_ptr : variable_field_ptr->variable_ptr_vec)
        {

            // iterate through each global ID
            for (auto point_gid : variable_ptr->mesh_l2_ptr->point_gid_vec)
            {

                // get domain and field IDs
                int point_fid = variable_field_ptr->point_gid_to_fid_map[point_gid];
                int point_did = variable_ptr->mesh_l2_ptr->point_gid_to_did_map[point_gid];

                // get value from x_vec
                int vec_row = start_row + point_fid;
                double value = x_vec.coeffRef(vec_row);

                // store value in variable
                variable_ptr->point_value_vec[point_did] = value;

            }

        }

    }

}

#endif
