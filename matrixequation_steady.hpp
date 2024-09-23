#ifndef MATRIXEQUATION_STEADY
#define MATRIXEQUATION_STEADY
#include <set>
#include <vector>
#include "Eigen/Eigen"
#include "physics_solidheattransfer_steady.hpp"

class MatrixEquationSteady
{

    public:

    // physics
    std::vector<PhysicsSolidHeatTransferSteady*> physics_ptr_vec;

    // matrix equation variables
    Eigen::SparseMatrix<double> a_mat;
    Eigen::VectorXd b_vec;
    Eigen::VectorXd x_vec;
    int num_equation = 0;

    // functions
    void iterate_x_vec();
    void store_x_vec();

    // default constructor
    MatrixEquationSteady()
    {

    }

    // constructor
    MatrixEquationSteady(std::vector<PhysicsSolidHeatTransferSteady*> physics_ptr_vec_in)
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
            for (auto variable_field_ptr : physics_ptr->variable_field_ptr_vec)
            {
                
                // assign starting column to variable if none yet
                // increment assign_start_col by number of mesh points
                if (variable_field_ptr->start_col == -1)
                {
                    variable_field_ptr->start_col = assign_start_col;
                    assign_start_col += variable_field_ptr->num_field_point;
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
        x_vec = Eigen::VectorXd::Zero(num_equation);  // edit later to use initial values

    }

};

void MatrixEquationSteady::iterate_x_vec()
{

    // fill up a_mat and b_vec
    for (auto physics_ptr : physics_ptr_vec)
    {
        physics_ptr->matrix_fill(a_mat, b_vec, x_vec);

    }

//    // DEBUG - PRINT A
//    std::ofstream file_amat_out("a_mat.csv");
//    for (int i = 0; i < num_equation; i++)
//    {
//        for (int j = 0; j < num_equation; j++)
//        {
//            
//            // last x value for given y
//            if (j == num_equation-1)
//            {
//                file_amat_out << a_mat.coeffRef(i, j) << "\n";
//                continue;
//            }
//
//            // output content of a matrix
//            file_amat_out << a_mat.coeffRef(i, j) << ",";
//
//        }
//    }

    // DEBUG - PRINT B
    std::ofstream file_bvec_out("b_vec.csv");
    for (int i = 0; i < num_equation; i++)
    {
        file_bvec_out << b_vec.coeffRef(i) << "\n";
    }

    // solve the matrix equation
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(a_mat);
    solver.factorize(a_mat);
    x_vec = solver.solve(b_vec);

}


void MatrixEquationSteady::store_x_vec()
{

    // transfer solutions from x_vec to variables
    for (auto physics_ptr : physics_ptr_vec)
    {
        physics_ptr->matrix_store(x_vec);
    }

}

#endif
