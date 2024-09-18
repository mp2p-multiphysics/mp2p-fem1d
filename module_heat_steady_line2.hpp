#ifndef MODULE_HEAT_STEADY_LINE2
#define MODULE_HEAT_STEADY_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "integral_line2.hpp"
#include "grid_line2.hpp"
#include "scalar_line2.hpp"

class HeatSteadyLine2Class
{

    public:

    // variables
    GridLine2Struct gl2s;
    BoundaryLine2Struct bl2s;
    IntegralLine2Class il2c;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
        ScalarLine2Class &thermal_conductivity_sl2c, ScalarLine2Class &heat_generation_sl2c,
        int start_id
    );

    // constructor
    HeatSteadyLine2Class(GridLine2Struct &gl2s_in, BoundaryLine2Struct &bl2s_in, IntegralLine2Class &il2c_in)
    {
        
        // variables
        gl2s = gl2s_in;
        bl2s = bl2s_in;
        il2c = il2c_in;

        // integrals
        il2c.evaluate_test_functions_derivatives();
        il2c.evaluate_integral_div_Ni_line2_dot_div_Nj_line2();
        il2c.evaluate_integral_Ni_line2();
        
    }

};

void HeatSteadyLine2Class::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
    ScalarLine2Class &thermal_conductivity_sl2c, ScalarLine2Class &heat_generation_sl2c,
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

        // get thermal conductivity of points around element
        double thermcond0 = thermal_conductivity_sl2c.scalar_vec[n0];
        double thermcond1 = thermal_conductivity_sl2c.scalar_vec[n1];
        double thermcond_arr[2] = {thermcond0, thermcond1};

        // get heat generation of points around element
        double heatgen0 = heat_generation_sl2c.scalar_vec[n0];
        double heatgen1 = heat_generation_sl2c.scalar_vec[n1];
        double heatgen_arr[2] = {heatgen0, heatgen1};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable

        // calculate a_mat coefficients
        for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            a_mat.coeffRef(n_arr[i], n_arr[j]) += thermcond_arr[i]*il2c.integral_div_Ni_line2_dot_div_Nj_line2_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 2; i++)
        {
            b_vec.coeffRef(n_arr[i]) += heatgen_arr[i]*il2c.integral_Ni_line2_vec[m][i];
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

        // get points where the boundary is applied
        int a = bl2s.element_flux_pa_loc_vec[k];  // 0 or 1

        // identify boundary type
        int config_id = bl2s.element_flux_config_id_vec[k];
        BoundaryConfigLine2Struct bcl2s = bl2s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2s.boundary_type_str == "neumann")
        {

            // add to b_vec
            b_vec.coeffRef(n_arr[a]) += bcl2s.boundary_parameter_vec[0];

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
