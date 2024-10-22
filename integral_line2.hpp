#ifndef INTEGRAL_LINE2
#define INTEGRAL_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "mesh_line2.hpp"
#include "container_typedef.hpp"

class IntegralLine2
{
    /*

    Test function (N) integrals for line2 mesh elements.

    Variables
    =========
    mesh_in : MeshLine2
        struct with mesh data.

    Functions
    =========
    evaluate_Ni_derivative : void
        Calculates test functions and their derivatives.
        Must be called before integrals are evaluated.
    evaluate_integral_Ni : void
        Calculates the integral of Ni.
    evaluate_integral_derivative_Ni_x : void
        Calculates the integral of d(Ni)/dx.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj.
    evaluate_integral_Ni_derivative_Nj_x : void
        Calculates the integral of Ni * d(Nj)/dx.
    evaluate_integral_div_Ni_dot_div_Nj : void
        Calculates the integral of div(Ni) dot div(Nj).
    evaluate_integral_Ni_Nj_derivative_Nk_x : void
        Calculates the integral of Ni * Nj * d(Nk)/dx.

    Notes
    ====
    The calculated integrals are stored in nested vectors.
    Values can be accessed from each vector using the following pattern:
        integral_vec[element_did][i][j]...
    wherein element_gid is the domain element ID and i, j, ... are indices.

    */

    public:
    
    // mesh
    MeshLine2 *mesh_ptr;

    // vectors with test functions and derivatives
    Vector2D jacobian_determinant_vec;
    Vector3D N_vec;
    Vector3D derivative_N_x_vec;

    // vectors with integrals
    Vector2D integral_Ni_vec;
    Vector2D integral_derivative_Ni_x_vec;
    Vector3D integral_Ni_Nj_vec;
    Vector3D integral_Ni_derivative_Nj_x_vec;
    Vector3D integral_div_Ni_dot_div_Nj_vec;
    Vector4D integral_Ni_Nj_derivative_Nk_x_vec;

    // functions for computing integrals
    void evaluate_Ni_derivative();
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();

    // default constructor
    IntegralLine2()
    {

    }

    // constructor
    IntegralLine2(MeshLine2 &mesh_in)
    {
        mesh_ptr = &mesh_in;
    }

};

void IntegralLine2::evaluate_Ni_derivative()
{
    /*

    Calculates test functions (N) and their derivatives.
    Must be called before integrals are evaluated.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1]
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[2] = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++)
    {

        // initialize
        Vector1D jacobian_determinant_part_ml_vec;
        Vector2D N_part_ml_vec;
        Vector2D derivative_N_x_part_ml_vec;

        // get global ID of points around element
        int p0_gid = mesh_ptr->element_p0_gid_vec[element_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[element_did];

        // get domain ID of points
        int p0_did = mesh_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_ptr->point_gid_to_did_map[p1_gid];

        // get x values of points
        double x0 = mesh_ptr->point_position_x_vec[p0_did];
        double x1 = mesh_ptr->point_position_x_vec[p1_did];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 2; indx_l++)
        {

            // initialize
            Vector1D N_part_mli_vec;
            Vector1D derivative_N_x_part_mli_vec;

            // get a values where function is evaluated
            double a = a_arr[indx_l];

            // get derivatives of x with respect to a
            double derivative_x_a = 0.5*(x1 - x0);

            // get jacobian and its inverse and determinant
            double jacobian_inverse = 1./derivative_x_a;
            double jacobian_determinant = derivative_x_a;

            // iterate for each test function
            for (int indx_i = 0; indx_i < 2; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 0.5*(1 - a); break;
                    case 1: N = 0.5*(1 + a); break;
                }

                // get derivatives of test function N
                double derivative_N_a = 0.;
                switch (indx_i)
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

void IntegralLine2::evaluate_integral_Ni()
{
    /*

    Calculates the integral of Ni.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Ni_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2::evaluate_integral_derivative_Ni_x()
{
    /*

    Calculates the integral of d(Ni)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */
    
    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * derivative_N_x_vec[element_did][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_derivative_Ni_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2::evaluate_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * N_vec[element_did][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2::evaluate_integral_Ni_derivative_Nj_x()
{
    /*

    Calculates the integral of Ni * d(Nj)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * derivative_N_x_vec[element_did][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2::evaluate_integral_div_Ni_dot_div_Nj()
{
    /*

    Calculates the integral of div(Ni) dot div(Nj).

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * derivative_N_x_vec[element_did][indx_l][indx_i] * derivative_N_x_vec[element_did][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Ni_dot_div_Nj_vec.push_back(integral_part_i_vec);

    }

}

void IntegralLine2::evaluate_integral_Ni_Nj_derivative_Nk_x()
{
    /*

    Calculates the integral of Ni * Nj * d(Nk)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
    Vector1D integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < 2; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * N_vec[element_did][indx_l][indx_j] * derivative_N_x_vec[element_did][indx_l][indx_k];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_derivative_Nk_x_vec.push_back(integral_part_i_vec);

    }

}

#endif
