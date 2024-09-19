#ifndef MATRIXEQUATION_STEADY
#define MATRIXEQUATION_STEADY
#include <set>
#include <vector>
#include "Eigen/Eigen"
#include "physics_solidheattransfer_steady_line2.hpp"

class MatrixEquationSteady
{

    public:

    // physics
    std::vector<PhysicsSolidHeatTransferSteadyLine2*> physics_ptr_vec;

    // matrix equation variables
    Eigen::SparseMatrix<double> a_mat;
    Eigen::VectorXd b_vec;
    Eigen::VectorXd x_last_iteration_vec;
    int num_equation = 0;

    // functions

    // default constructor
    MatrixEquationSteady()
    {

    }

    // constructor
    MatrixEquationSteady(std::vector<PhysicsSolidHeatTransferSteadyLine2*> physics_ptr_vec_in)
    {

        // store vector of pointers to physics
        physics_ptr_vec = physics_ptr_vec_in;

        // initialize starting rows and columns
        int assign_start_row = 0;
        int assign_start_col = 0;

        // iterate through each physics
        // assign starting row in matrix to physics
        // assign starting column in matrix to variables
        for (auto physics_ptr : physics_ptr_vec)
        {

            // iterate through each variable
            for (auto variable_ptr : physics_ptr->variable_ptr_vec)
            {
                
                // assign starting column to variable if none yet
                // increment assign_start_col by number of mesh points
                if (variable_ptr->start_col == -1)
                {
                    variable_ptr->start_col = assign_start_col;
                    assign_start_col += variable_ptr->mesh_l2_ptr->num_point;
                }

                // assign starting row to physics
                // increment assign_start_row by number of new mesh points
                physics_ptr->start_row = assign_start_row;
                assign_start_row = assign_start_col;

            }

        }

        // get number of linear equations (total number of mesh points)
        num_equation = assign_start_row;

        // initialize matrix equation variables
        a_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
        b_vec = Eigen::VectorXd::Zero(num_equation);
        x_last_iteration_vec = Eigen::VectorXd::Zero(num_equation);

    }

};




#endif
