#ifndef MATRIXEQUATION_TRANSIENT
#define MATRIXEQUATION_TRANSIENT
#include <set>
#include <vector>
#include "Eigen/Eigen"
#include "physicstransient_base.hpp"
#include "variable_group.hpp"

namespace FEM1D
{

class MatrixEquationTransient
{
    /*

    Represents the matrix equation (system of linear equations) Ax = b for use in transient problems.

    Variables
    =========
    physics_ptr_vec_in : vector<PhysicsTransientBase*>
        vector with pointers to PhysicsTransientBase objects.

    Functions
    =========
    set_last_timestep_solution : void
        Sets the current value of x(t+1) to x(t) for use in the next timestep.
    iterate_solution : void
        Solves for x in Ax = b.
        Uses the currently stored values in A and b.
    store_solution : void
        Transfers the solutions in x into variable objects.

    Notes
    =====
    The equation Ax(t+1) = b is expanded into Ax(t+1) = Cx(t) + d for convenience.
    In the code; A, x, C, and d are referred to as a_mat, x_vec, c_mat, and d_vec respectively.

    */

    public:

    // vector of physics
    std::vector<PhysicsTransientBase*> physics_ptr_vec;
    std::vector<BoundaryGroup*> boundary_group_ptr_vec;
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;

    // matrix equation variables
    Eigen::SparseMatrix<double> a_mat;
    Eigen::SparseMatrix<double> c_mat;
    Eigen::VectorXd d_vec;
    Eigen::VectorXd x_vec;
    Eigen::VectorXd x_last_timestep_vec;
    int num_equation = 0;

    // functions
    void set_last_timestep_solution();
    void iterate_solution(double dt);
    void store_solution();

    // default constructor
    MatrixEquationTransient()
    {

    }

    // constructor
    MatrixEquationTransient(std::vector<PhysicsTransientBase*> physics_ptr_vec_in)
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

            // iterate through each variable group
            for (auto variable_group_ptr : physics_ptr->get_variable_group_ptr_vec())
            {
                
                // assign starting column to variable if none yet
                // increment assign_start_col by number of domain points
                if (variable_group_ptr->start_col == -1)
                {
                    variable_group_ptr->start_col = assign_start_col;
                    assign_start_col += variable_group_ptr->num_point;
                }

                // assign starting row to physics
                // increment assign_start_row by number of new domain points
                if (physics_ptr->get_start_row() == -1)
                {
                    physics_ptr->set_start_row(assign_start_row);
                    assign_start_row = assign_start_col;
                }

            }

        }

        // get number of linear equations (total number of domain points)
        num_equation = assign_start_col;

        // get vector of boundary groups
        // iterate through each physics and store boundary groups
        std::set<BoundaryGroup*> boundary_group_ptr_set;
        for (auto physics_ptr : physics_ptr_vec)
        {
            boundary_group_ptr_set.insert(physics_ptr->get_boundary_group_ptr());
        }
        boundary_group_ptr_vec = std::vector<BoundaryGroup*>(boundary_group_ptr_set.begin(), boundary_group_ptr_set.end());

        // get vector of scalar groups
        // iterate through each physics and store scalar groups
        std::set<ScalarGroup*> scalar_group_ptr_set;
        for (auto physics_ptr : physics_ptr_vec) {
        for (auto scalar_group_ptr : physics_ptr->get_scalar_group_ptr_vec()) {
            scalar_group_ptr_set.insert(scalar_group_ptr);
        }}
        scalar_group_ptr_vec = std::vector<ScalarGroup*>(scalar_group_ptr_set.begin(), scalar_group_ptr_set.end());

        // get vector of variable groups
        // iterate through each physics and store variable groups
        std::set<VariableGroup*> variable_group_ptr_set;
        for (auto physics_ptr : physics_ptr_vec) {
        for (auto variable_group_ptr : physics_ptr->get_variable_group_ptr_vec()) {
            variable_group_ptr_set.insert(variable_group_ptr);
        }}
        variable_group_ptr_vec = std::vector<VariableGroup*>(variable_group_ptr_set.begin(), variable_group_ptr_set.end());

        // initialize matrix equation variables
        a_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
        c_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
        d_vec = Eigen::VectorXd::Zero(num_equation);
        x_vec = Eigen::VectorXd::Zero(num_equation);
        
        // populate x_vec with initial values

        // iterate through each variable group
        for (auto variable_group_ptr : variable_group_ptr_vec)
        {

            // get starting row
            // note: column in a_mat = row in x_vec
            int start_row = variable_group_ptr->start_col;

            // iterate through each variable
            for (auto variable_ptr : variable_group_ptr->variable_l2_ptr_vec)
            {

                // iterate through each global ID
                for (auto pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
                {

                    // get domain and group IDs
                    int pfid = variable_group_ptr->point_pgid_to_pfid_map[pgid];
                    int pdid = variable_ptr->domain_ptr->point_pgid_to_pdid_map[pgid];

                    // get value from variable
                    double value = variable_ptr->point_value_vec[pdid];

                    // store value in x_vec
                    int vec_row = start_row + pfid;
                    x_vec.coeffRef(vec_row) = value;

                }

            }

        }

        // use initial values as previous values
        x_last_timestep_vec = x_vec;

    }

};

void MatrixEquationTransient::set_last_timestep_solution()
{
    /*
    
    Sets the current value of x(t+1) to x(t) for use in the next timestep.

    Arguments
    =========
    (none)

    Returns
    =======
    (none)

    */

    x_last_timestep_vec = x_vec;

}

void MatrixEquationTransient::iterate_solution(double dt)
{
    /*
    
    Solves for x in Ax = b.
    Uses the currently stored values in A and b.

    Arguments
    =========
    dt : double
        Length of the timestep.

    Returns
    =======
    (none)

    */

    // update boundary parameters using most recent variable values
    for (auto boundary_group_ptr : boundary_group_ptr_vec)
    {
        boundary_group_ptr->update_parameter();
    }

    // update scalar values using most recent variable values
    for (auto scalar_group_ptr : scalar_group_ptr_vec)
    {
        scalar_group_ptr->update_value();
    }

    // reset matrices
    a_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
    c_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
    d_vec = Eigen::VectorXd::Zero(num_equation);
    x_vec = Eigen::VectorXd::Zero(num_equation);

    // fill up a_mat, c_mat, and d_vec with each physics
    for (auto physics_ptr : physics_ptr_vec)
    {
        physics_ptr->matrix_fill(a_mat, c_mat, d_vec, x_vec, x_last_timestep_vec, dt);
    }

    // solve the matrix equation
    // b_vec = c_mat*x_last_timestep_vec + d_vec
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(a_mat);
    solver.factorize(a_mat);
    x_vec = solver.solve(c_mat*x_last_timestep_vec + d_vec);

}

void MatrixEquationTransient::store_solution()
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

    // iterate through each variable group
    for (auto variable_group_ptr : variable_group_ptr_vec)
    {

        // get starting row
        // note: column in a_mat = row in x_vec
        int start_row = variable_group_ptr->start_col;

        // iterate through each variable
        for (auto variable_ptr : variable_group_ptr->variable_l2_ptr_vec)
        {

            // iterate through each global ID
            for (auto pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {

                // get domain and group IDs
                int pfid = variable_group_ptr->point_pgid_to_pfid_map[pgid];
                int pdid = variable_ptr->domain_ptr->point_pgid_to_pdid_map[pgid];

                // get value from x_vec
                int vec_row = start_row + pfid;
                double value = x_vec.coeffRef(vec_row);

                // store value in variable
                variable_ptr->point_value_vec[pdid] = value;

            }

        }

    }

}

}

#endif
