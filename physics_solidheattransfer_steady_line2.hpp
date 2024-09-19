#ifndef PHYSICS_SOLIDHEATTRANSFER_STEADY_LINE2
#define PHYSICS_SOLIDHEATTRANSFER_STEADY_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "container_boundary.hpp"
#include "container_mesh.hpp"
#include "container_typedef.hpp"
#include "field_line2.hpp"
#include "integral_line2.hpp"
#include "variable_line2.hpp"

class PhysicsSolidHeatTransferSteadyLine2
{

    public:

    // mesh and boundary conditions
    MeshLine2Struct *mesh_l2_ptr;
    BoundaryLine2Struct *boundary_l2_ptr;

    // integrals
    IntegralLine2 *integral_l2_ptr;

    // variables
    VariableLine2 *temperature_l2_ptr;
    std::vector<VariableLine2*> variable_ptr_vec;  // list of variables
    
    // fields
    FieldLine2 *thermalconductivity_l2_ptr;
    FieldLine2 *heatgeneration_l2_ptr;

    // starting row in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec);

    // default constructor
    PhysicsSolidHeatTransferSteadyLine2()
    {

    }

    // constructor
    PhysicsSolidHeatTransferSteadyLine2
    (
        MeshLine2Struct &mesh_l2_in, BoundaryLine2Struct &boundary_l2_in, IntegralLine2 &integral_l2_in,
        VariableLine2 &temperature_l2_in, FieldLine2 &thermalconductivity_l2_in, FieldLine2 &heatgeneration_l2_in
    )
    {
        
        // store mesh and boundary conditions
        mesh_l2_ptr = &mesh_l2_in;
        boundary_l2_ptr = &boundary_l2_in;
        
        // store variables
        temperature_l2_ptr = &temperature_l2_in;
        variable_ptr_vec.push_back(&temperature_l2_in);

        // store fields
        thermalconductivity_l2_ptr = &thermalconductivity_l2_in;
        heatgeneration_l2_ptr = &heatgeneration_l2_in;

        // calculate integrals
        integral_l2_ptr = &integral_l2_in;
        integral_l2_ptr->evaluate_Ni_derivative();
        integral_l2_ptr->evaluate_integral_div_Ni_line2_dot_div_Nj_line2();
        integral_l2_ptr->evaluate_integral_Ni_line2();

    }

};

void PhysicsSolidHeatTransferSteadyLine2::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec
)
{

    // iterate for each grid element
    for (int m = 0; m < mesh_l2_ptr->num_element; m++)
    {

        // get id of points around element
        int n0 = mesh_l2_ptr->element_p0_id_vec[m];
        int n1 = mesh_l2_ptr->element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get thermal conductivity of points around element
        double thermcond0 = thermalconductivity_l2_ptr->point_u_vec[n0];
        double thermcond1 = thermalconductivity_l2_ptr->point_u_vec[n1];
        double thermcond_arr[2] = {thermcond0, thermcond1};

        // get heat generation of points around element
        double heatgen0 = heatgeneration_l2_ptr->point_u_vec[n0];
        double heatgen1 = heatgeneration_l2_ptr->point_u_vec[n1];
        double heatgen_arr[2] = {heatgen0, heatgen1};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable

        // calculate a_mat coefficients
        for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            a_mat.coeffRef(start_row + n_arr[i], temperature_l2_ptr->start_col + n_arr[j]) += thermcond_arr[i]*integral_l2_ptr->integral_div_Ni_line2_dot_div_Nj_line2_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 2; i++)
        {
            b_vec.coeffRef(start_row + n_arr[i]) += heatgen_arr[i]*integral_l2_ptr->integral_Ni_line2_vec[m][i];
        }

    }

    // iterate for each flux boundary element
    for (int k = 0; k < boundary_l2_ptr->num_element_flux; k++)
    {

        // get id of element
        int m = boundary_l2_ptr->element_flux_id_vec[k];

        // get id of points around element
        int n0 = mesh_l2_ptr->element_p0_id_vec[m];
        int n1 = mesh_l2_ptr->element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = boundary_l2_ptr->element_flux_pa_loc_vec[k];  // 0 or 1

        // identify boundary type
        int config_id = boundary_l2_ptr->element_flux_config_id_vec[k];
        BoundaryConfigLine2Struct bcl2 = boundary_l2_ptr->boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2.boundary_type_str == "neumann")
        {

            // add to b_vec
            b_vec.coeffRef(start_row + n_arr[a]) += bcl2.boundary_parameter_vec[0];

        }

    }

    // clear rows with value boundary elements
    for (int k = 0; k < boundary_l2_ptr->num_element_value; k++)
    {

        // get id of element
        int m = boundary_l2_ptr->element_value_id_vec[k];

        // get id of points around element
        int n0 = mesh_l2_ptr->element_p0_id_vec[m];
        int n1 = mesh_l2_ptr->element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = boundary_l2_ptr->element_value_pa_loc_vec[k];  // 0 or 1

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_row + n_arr[a]) *= 0.;
            b_vec.coeffRef(start_row + n_arr[a]) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int k = 0; k < boundary_l2_ptr->num_element_value; k++)
    {

        // get id of element
        int m = boundary_l2_ptr->element_value_id_vec[k];

        // get id of points around element
        int n0 = mesh_l2_ptr->element_p0_id_vec[m];
        int n1 = mesh_l2_ptr->element_p1_id_vec[m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = boundary_l2_ptr->element_value_pa_loc_vec[k];  // 0 or 1
        
        // identify boundary type
        int config_id = boundary_l2_ptr->element_value_config_id_vec[k];
        BoundaryConfigLine2Struct bcl2 = boundary_l2_ptr->boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_row + n_arr[a], temperature_l2_ptr->start_col + n_arr[a]) += 1.;
                b_vec.coeffRef(start_row + n_arr[a]) += bcl2.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
