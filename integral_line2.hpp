#ifndef INTEGRAL_LINE2
#define INTEGRAL_LINE2
#include <vector>
#include "grid_line2.hpp"
#include "Eigen/Eigen"

class IntegralLine2Class
{

    public:
    
    // variables
    GridLine2Struct gl2s;

    // vectors with test functions and derivatives
    std::vector<std::vector<double>> jacobian_determinant_vec;
    std::vector<std::vector<std::vector<double>>> N_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_x_vec;

    // vectors with integrals
    std::vector<std::vector<double>> integral_Ni_line2_vec;
    std::vector<std::vector<double>> integral_derivative_Ni_line2_x_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_line2_Nj_line2_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_line2_derivative_Nj_line2_x_vec;
    std::vector<std::vector<std::vector<double>>> integral_div_Ni_line2_dot_div_Nj_line2_vec;

    // functions for computing integrals
    void evaluate_test_functions_derivatives();
    void evaluate_integral_Ni_line2();
    void evaluate_integral_derivative_Ni_line2_x();
    void evaluate_integral_Ni_line2_Nj_line2();
    void evaluate_integral_Ni_line2_derivative_Nj_line2_x();
    void evaluate_integral_div_Ni_line2_dot_div_Nj_line2();

    // constructor
    IntegralLine2Class()
    {

    }
    IntegralLine2Class(GridLine2Struct &gl2s_in)
    {
        gl2s = gl2s_in;
    }

};

void IntegralLine2Class::evaluate_test_functions_derivatives()
{

    // integration points
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[2] = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate for each element
    for (int m = 0; m < gl2s.num_element; m++)
    {

        // initialize
        std::vector<double> jacobian_determinant_part_ml_vec;
        std::vector<std::vector<double>> N_part_ml_vec;
        std::vector<std::vector<double>> derivative_N_x_part_ml_vec;

        // get id of points around element
        int n0 = gl2s.element_p0_id_vec[m];
        int n1 = gl2s.element_p1_id_vec[m];

        // unpack x values
        double x0 = gl2s.point_pos_x_vec[n0];
        double x1 = gl2s.point_pos_x_vec[n1];

        // iterate for each integration point
        for (int l = 0; l < 2; l++)
        {

            // initialize
            std::vector<double> N_part_mli_vec;
            std::vector<double> derivative_N_x_part_mli_vec;

            // get a values where function is evaluated
            double a = a_arr[l];

            // get derivatives of x with respect to a
            double derivative_x_a = 0.5*(x1 - x0);

            // get jacobian and its inverse and determinant
            double jacobian_inverse = 1/derivative_x_a;
            double jacobian_determinant = derivative_x_a;

            // iterate for each test function
            for (int i = 0; i < 2; i++)
            {
        
                // get test function N
                double N = 0.;
                switch (i)
                {
                    case 0: N = 0.5*(1 - a); break;
                    case 1: N = 0.5*(1 + a); break;
                }

                // get derivatives of test function N
                double derivative_N_a = 0.;
                switch (i)
                {
                    case 0: derivative_N_a = -0.5; break;
                    case 1: derivative_N_a = +0.5; break;
                }


                // get derivatives of test functions wrt x
                double derivative_N_x = derivative_N_a*jacobian_inverse;

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_N_x_part_mli_vec.push_back(derivative_N_x);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            derivative_N_x_part_ml_vec.push_back(derivative_N_x_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        N_vec.push_back(N_part_ml_vec);
        derivative_N_x_vec.push_back(derivative_N_x_part_ml_vec);
        
    }

}

void IntegralLine2Class::evaluate_integral_Ni_line2()
{
    
    // iterate for each element
    for (int m = 0; m < gl2s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<double> integral_part_i_vec;
    for (int i = 0; i < 2; i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 2; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Ni_line2_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2Class::evaluate_integral_derivative_Ni_line2_x()
{
    
    // iterate for each element
    for (int m = 0; m < gl2s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<double> integral_part_i_vec;
    for (int i = 0; i < 2; i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 2; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * derivative_N_x_vec[m][l][i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_derivative_Ni_line2_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2Class::evaluate_integral_Ni_line2_Nj_line2()
{
    
    // iterate for each element
    for (int m = 0; m < gl2s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 2; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 2; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 2; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i] * N_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_line2_Nj_line2_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2Class::evaluate_integral_Ni_line2_derivative_Nj_line2_x()
{
    
    // iterate for each element
    for (int m = 0; m < gl2s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 2; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 2; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 2; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i] * derivative_N_x_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_line2_derivative_Nj_line2_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2Class::evaluate_integral_div_Ni_line2_dot_div_Nj_line2()
{
    
    // iterate for each element
    for (int m = 0; m < gl2s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 2; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 2; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 2; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * derivative_N_x_vec[m][l][i] * derivative_N_x_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Ni_line2_dot_div_Nj_line2_vec.push_back(integral_part_i_vec);

    }

}

#endif
