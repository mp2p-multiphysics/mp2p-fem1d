#ifndef INTEGRAL_LINE2
#define INTEGRAL_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "boundary_line2.hpp"
#include "container_typedef.hpp"
#include "domain_line2.hpp"

namespace FEM1D
{

class IntegralLine2
{
    /*

    Test function (N) integrals for line2 elements.

    Variables
    =========
    domain_in : DomainLine2
        Domain where element integrals are calculated.
    boundary_in : BoundaryLine2
        Boundaries where element integrals are calculated.

    Functions
    =========
    evaluate_Ni : void
        Calculates test function values and other properties.
        Must be called before domain integrals are evaluated.
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
    evaluate_boundary_Ni : void
        Calculates test functions values and other properties at the boundary.
        Must be called before boundary integrals are evaluated.
    evaluate_integral_boundary_Ni : void
        Calculates the integral of Ni at the boundary.
    evaluate_integral_boundary_Ni_Nj : void
        Calculates the integral of Ni * Nj at the boundary.

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.
    Values of the boundary integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][boundary_key][i][j]...
        wherein boundary_key is an int denoting the location of the boundary.
        For line2 elements, the boundary_key is just the local ID of the point (0 or 1).

    */

    public:
    
    // domain and boundary
    DomainLine2 *domain_ptr;
    BoundaryLine2 *boundary_ptr;

    // vectors with domain test functions
    // index as follows: [edid][integration_point][i]
    Vector3D Ni_vec;
    Vector3D derivative_Ni_x_vec;
    Vector2D jacobian_determinant_vec;

    // vectors with domain integrals
    // index as follows: [edid][i][j][k]
    Vector2D integral_Ni_vec;
    Vector2D integral_derivative_Ni_x_vec;
    Vector3D integral_Ni_Nj_vec;
    Vector3D integral_Ni_derivative_Nj_x_vec;
    Vector3D integral_div_Ni_dot_div_Nj_vec;
    Vector4D integral_Ni_Nj_derivative_Nk_x_vec;

    // vectors with boundary test functions
    // index as follows: [edid][boundary_key][i]
    MapVector3D boundary_Ni_vec;
    MapVector2D boundary_normal_x_vec;

    // vectors with boundary integrals
    // index as follows: [edid][boundary_key][i][j]
    MapVector3D integral_boundary_Ni_vec;
    MapVector4D integral_boundary_Ni_Nj_vec;

    // functions for computing domain integrals
    void evaluate_Ni();
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();

    // functions for computing boundary integrals
    void evaluate_boundary_Ni();
    void evaluate_integral_boundary_Ni();
    void evaluate_integral_boundary_Ni_Nj();

    // default constructor
    IntegralLine2() {}

    // constructor
    IntegralLine2(DomainLine2 &domain_in, BoundaryLine2 &boundary_in)
    {
        
        // store domain and boundaries
        domain_ptr = &domain_in;
        boundary_ptr = &boundary_in;

        // evaluate test functions
        evaluate_Ni();
        evaluate_boundary_Ni();

    }

};

void IntegralLine2::evaluate_Ni()
{
    /*

    Calculates test function values and other properties.
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        Vector1D jacobian_determinant_part_ml_vec;
        Vector2D N_part_ml_vec;
        Vector2D derivative_N_x_part_ml_vec;

        // get global ID of points around element
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];

        // get domain ID of points
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];

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
        Ni_vec.push_back(N_part_ml_vec);
        derivative_Ni_x_vec.push_back(derivative_N_x_part_ml_vec);
        
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 2; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector1D integral_part_ij_vec;
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector1D integral_part_ij_vec;
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 2; indx_i++){  
    Vector1D integral_part_ij_vec;
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
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

void IntegralLine2::evaluate_boundary_Ni()
{
    /*

    Calculates test function values and other properties at the boundary.
    Must be called before bounary integrals are evaluated.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // for 1D elements, just use the local point
        int boundary_key = boundary_ptr->boundary_pa_plid_vec[bid];

        // get dimensionless coordinates
        double a = 0;
        switch (boundary_key)
        {
            case 0: a = -1.; break;
            case 1: a = +1.; break;
        }

        // calculate normal vectors
        double normal_x = 0;
        switch (boundary_key)
        {
            case 0: normal_x = -1.; break;
            case 1: normal_x = +1.; break;
        }
        boundary_normal_x_vec[edid][boundary_key] = normal_x;

        // calculate test function values
        for (int indx_i = 0; indx_i < 2; indx_i++)
        {

            // get test function N
            double N = 0.;
            switch (indx_i)
            {
                case 0: N = 0.5*(1 - a); break;
                case 1: N = 0.5*(1 + a); break;
            }
            boundary_Ni_vec[edid][boundary_key].push_back(N);

        }

    }

}

void IntegralLine2::evaluate_integral_boundary_Ni()
{
    /*

    Calculates the integral of Ni at the boundary.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // for 1D elements, just use the local point
        int boundary_key = boundary_ptr->boundary_pa_plid_vec[bid];

        // iterate for each test function combination
        Vector1D integral_part_i_vec;
        for (int indx_i = 0; indx_i < 2; indx_i++){

            // calculate integral value
            double integral_value = boundary_Ni_vec[edid][boundary_key][indx_i];
            integral_part_i_vec.push_back(integral_value);

        }
        integral_boundary_Ni_vec[edid][boundary_key] = integral_part_i_vec;

    }

}

void IntegralLine2::evaluate_integral_boundary_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj at the boundary.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // for 1D elements, just use the local point
        int boundary_key = boundary_ptr->boundary_pa_plid_vec[bid];

        // iterate for each test function combination
        Vector2D integral_part_i_vec;
        for (int indx_i = 0; indx_i < 2; indx_i++){
        Vector1D integral_part_ij_vec;
        for (int indx_j = 0; indx_j < 2; indx_j++){

            // calculate integral value
            double integral_value = boundary_Ni_vec[edid][boundary_key][indx_i];
            integral_part_ij_vec.push_back(integral_value);

        }
        integral_part_i_vec.push_back(integral_part_ij_vec);
        }
        integral_boundary_Ni_Nj_vec[edid][boundary_key] = integral_part_i_vec;

    }

}

}

#endif
