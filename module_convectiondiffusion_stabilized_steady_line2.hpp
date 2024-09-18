#ifndef MODULE_CONVECTIONDIFFUSION_STABILIZED_STEADY_LINE2
#define MODULE_CONVECTIONDIFFUSION_STABILIZED_STEADY_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "integral_line2.hpp"
#include "grid_line2.hpp"
#include "scalar_line2.hpp"

class ConvectionDiffusionStabilizedSteadyLine2Class
{

    public:

    // variables
    GridLine2Struct gl2s;
    BoundaryLine2Struct bl2s;
    IntegralLine2Class il2c;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
        ScalarLine2Class &diffusion_coefficient_sl2c, ScalarLine2Class &velocity_sl2c, ScalarLine2Class &species_generation_sl2c,
        int start_id
    );

    // constructor
    ConvectionDiffusionStabilizedSteadyLine2Class(GridLine2Struct &gl2s_in, BoundaryLine2Struct &bl2s_in, IntegralLine2Class &il2c_in)
    {
        
        // variables
        gl2s = gl2s_in;
        bl2s = bl2s_in;
        il2c = il2c_in;

        // integrals
        il2c.evaluate_test_functions_derivatives();
        il2c.evaluate_integral_Ni_line2();
        il2c.evaluate_integral_derivative_Ni_line2_x();
        il2c.evaluate_integral_div_Ni_line2_dot_div_Nj_line2();
        il2c.evaluate_integral_Ni_line2_derivative_Nj_line2_x();
        
    }

};

void ConvectionDiffusionStabilizedSteadyLine2Class::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
    ScalarLine2Class &diffusion_coefficient_sl2c, ScalarLine2Class &velocity_sl2c, ScalarLine2Class &species_generation_sl2c,
    int start_id
)
{

    // iterate for each grid element
    for (int m = 0; m < gl2s.num_element; m++)
    {

        // get id of points around element
        int n0 = gl2s.element_p0_id_vec[m];
        int n1 = gl2s.element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get diffusion coefficient of points around element
        double diffcoeff0 = diffusion_coefficient_sl2c.scalar_vec[n0];
        double diffcoeff1 = diffusion_coefficient_sl2c.scalar_vec[n1];
        double diffcoeff_arr[2] = {diffcoeff0, diffcoeff1};

        // get x coordinates of points around element
        double x0 = gl2s.point_pos_x_vec[n0];
        double x1 = gl2s.point_pos_x_vec[n1];

        // get velocity of points around element
        double v0 = velocity_sl2c.scalar_vec[n0];
        double v1 = velocity_sl2c.scalar_vec[n1];
        double v_arr[2] = {v0, v1};

        // get species generation of points around element
        double specgen0 = species_generation_sl2c.scalar_vec[n0];
        double specgen1 = species_generation_sl2c.scalar_vec[n1];
        double specgen_arr[2] = {specgen0, specgen1};

        // get stabilization parameter
        double tau = (x1 - x0)/(v0 + v1);

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable

        // calculate a_mat coefficients
        for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            a_mat.coeffRef(n_arr[i], n_arr[j]) += (
                diffcoeff_arr[i]*il2c.integral_div_Ni_line2_dot_div_Nj_line2_vec[m][i][j]
                + v_arr[i]*il2c.integral_Ni_line2_derivative_Nj_line2_x_vec[m][i][j]
                + tau*v_arr[i]*v_arr[i]*il2c.integral_div_Ni_line2_dot_div_Nj_line2_vec[m][i][j]
            );
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 2; i++)
        {
            b_vec.coeffRef(n_arr[i]) += (
                specgen_arr[i]*il2c.integral_Ni_line2_vec[m][i]
                + tau*v_arr[i]*specgen_arr[i]*il2c.integral_derivative_Ni_line2_x_vec[m][i]
            );
        }

    }

    // iterate for each flux boundary element
    for (int k = 0; k < bl2s.num_element_flux; k++)
    {

        // get id of element
        int m = bl2s.element_flux_id_vec[k];

        // get id of points around element
        int n0 = gl2s.element_p0_id_vec[m];
        int n1 = gl2s.element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get x coordinates of points around element
        double x0 = gl2s.point_pos_x_vec[n0];
        double x1 = gl2s.point_pos_x_vec[n1];

        // get velocity of points around element
        double v0 = velocity_sl2c.scalar_vec[n0];
        double v1 = velocity_sl2c.scalar_vec[n1];
        double v_arr[2] = {v0, v1};

        // get stabilization parameter
        double tau = (x1 - x0)/(v0 + v1);

        // get inverse jacobian
        double inverse_jacobian = 2./(x1 - x0);

        // get points where the boundary is applied
        int a = bl2s.element_flux_pa_loc_vec[k];  // 0 or 1

        // identify boundary type
        int config_id = bl2s.element_flux_config_id_vec[k];
        BoundaryConfigLine2Struct bcl2s = bl2s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2s.boundary_type_str == "neumann")
        {

            // add to b_vec
            b_vec.coeffRef(n_arr[a]) += bcl2s.boundary_parameter_vec[0] * (1. - 0.5*tau*v_arr[a]*inverse_jacobian);

        }

    }

    // clear rows with value boundary elements
    for (int k = 0; k < bl2s.num_element_value; k++)
    {

        // get id of element
        int m = bl2s.element_value_id_vec[k];

        // get id of points around element
        int n0 = gl2s.element_p0_id_vec[m];
        int n1 = gl2s.element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = bl2s.element_value_pa_loc_vec[k];  // 0 or 1

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(n_arr[a]) *= 0.;
            b_vec.coeffRef(n_arr[a]) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int k = 0; k < bl2s.num_element_value; k++)
    {

        // get id of element
        int m = bl2s.element_value_id_vec[k];

        // get id of points around element
        int n0 = gl2s.element_p0_id_vec[m];
        int n1 = gl2s.element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = bl2s.element_value_pa_loc_vec[k];  // 0 or 1
        
        // identify boundary type
        int config_id = bl2s.element_value_config_id_vec[k];
        BoundaryConfigLine2Struct bcl2s = bl2s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(n_arr[a], n_arr[a]) += 1.;
                b_vec.coeffRef(n_arr[a]) += bcl2s.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
