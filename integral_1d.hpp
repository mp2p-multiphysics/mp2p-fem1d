#ifndef INTEGRAL_1D
#define INTEGRAL_1D
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_1d.hpp"

namespace FEM1D
{

class Integral1D
{
    /*

    Test function (N) integrals over a 1D domain.

    Variables
    =========
    domain_in : Domain1D
        Domain where element integrals are calculated.

    Functions
    =========
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
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.

    */

    public:
    
    // domain where integral is applied
    Domain1D *domain_ptr;

    // vectors with domain test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_vec;
    VectorDouble3D Ni_vec;
    VectorDouble3D derivative_Ni_x_vec;
    VectorDouble2D jacobian_determinant_vec;

    // vectors with domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble2D integral_Ni_vec;
    VectorDouble2D integral_derivative_Ni_x_vec;
    VectorDouble3D integral_Ni_Nj_vec;
    VectorDouble3D integral_Ni_derivative_Nj_x_vec;
    VectorDouble3D integral_div_Ni_dot_div_Nj_vec;
    VectorDouble4D integral_Ni_Nj_derivative_Nk_x_vec;

    // functions for computing domain integrals
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();

    // default constructor
    Integral1D() {}

    // constructor
    Integral1D(Domain1D &domain_in)
    {
        
        // store domain and boundaries
        domain_ptr = &domain_in;

        // evaluate test functions
        switch (domain_ptr->type_element)
        {
            case 0:
                evaluate_Ni_line2();
            break;
        }

    }

    private:
    void evaluate_Ni_line2();

};

void Integral1D::evaluate_integral_Ni()
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Ni_vec.push_back(integral_part_i_vec);

    }

}

void Integral1D::evaluate_integral_derivative_Ni_x()
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * derivative_Ni_x_vec[edid][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_derivative_Ni_x_vec.push_back(integral_part_i_vec);

    }

}

void Integral1D::evaluate_integral_Ni_Nj()
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * Ni_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_vec.push_back(integral_part_i_vec);

    }

}

void Integral1D::evaluate_integral_Ni_derivative_Nj_x()
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * derivative_Ni_x_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_x_vec.push_back(integral_part_i_vec);

    }

}

void Integral1D::evaluate_integral_div_Ni_dot_div_Nj()
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * derivative_Ni_x_vec[edid][indx_l][indx_i] * derivative_Ni_x_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Ni_dot_div_Nj_vec.push_back(integral_part_i_vec);

    }

}

void Integral1D::evaluate_integral_Ni_Nj_derivative_Nk_x()
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    VectorDouble2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 2; indx_j++){
    VectorDouble integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < 2; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * Ni_vec[edid][indx_l][indx_j] * derivative_Ni_x_vec[edid][indx_l][indx_k];
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

void Integral1D::evaluate_Ni_line2()
{

    // weights for integration
    weight_vec = {1., 1.};

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1]
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[2] = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D derivative_N_x_part_ml_vec;

        // get global ID of points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];

        // get domain ID of points
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];

        // get coordinates of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 2; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble derivative_N_x_part_mli_vec;

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
        Ni_vec.push_back(N_part_ml_vec);
        derivative_Ni_x_vec.push_back(derivative_N_x_part_ml_vec);

    }

}

}

#endif
