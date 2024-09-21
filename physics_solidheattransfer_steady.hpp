#ifndef PHYSICS_SOLIDHEATTRANSFER_STEADY_LINE2
#define PHYSICS_SOLIDHEATTRANSFER_STEADY_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "boundary_physicsgroup.hpp"
#include "container_typedef.hpp"
#include "mesh_physicsgroup.hpp"
#include "scalar_fieldgroup.hpp"
#include "integral_physicsgroup.hpp"
#include "variable_fieldgroup.hpp"

class PhysicsSolidHeatTransferSteady
{

    public:

    // physics groups
    MeshPhysicsGroup *mesh_physics_ptr;
    BoundaryPhysicsGroup *boundary_physics_ptr;
    IntegralPhysicsGroup *integral_physics_ptr;

    // field groups
    VariableFieldGroup *temperature_field_ptr;
    ScalarFieldGroup *thermalconductivity_field_ptr;
    ScalarFieldGroup *heatgeneration_field_ptr;

    // vector of variable fields
    std::vector<VariableFieldGroup*> variable_field_ptr_vec;

    // starting row in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);

    // default constructor
    PhysicsSolidHeatTransferSteady()
    {

    }

    // constructor
    PhysicsSolidHeatTransferSteady
    (
        MeshPhysicsGroup &mesh_physics_in, BoundaryPhysicsGroup &boundary_physics_in, IntegralPhysicsGroup &integral_physics_in,
        VariableFieldGroup &temperature_field_in, ScalarFieldGroup &thermalconductivity_field_in, ScalarFieldGroup &heatgeneration_field_in
    )
    {
        
        // store physics groups
        mesh_physics_ptr = &mesh_physics_in;
        boundary_physics_ptr = &boundary_physics_in;
        integral_physics_ptr = &integral_physics_in;

        // store field groups
        temperature_field_ptr = &temperature_field_in;
        thermalconductivity_field_ptr = &thermalconductivity_field_in;
        heatgeneration_field_ptr = &heatgeneration_field_in;

        // vector of variable fields 
        variable_field_ptr_vec = {temperature_field_ptr};

        // calculate integrals
        integral_physics_ptr->evaluate_Ni_derivative();
        integral_physics_ptr->evaluate_integral_div_Ni_line2_dot_div_Nj_line2();
        integral_physics_ptr->evaluate_integral_Ni_line2();

    }

    private:
    void matrix_fill_domain
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        MeshLine2Struct *mesh_ptr, BoundaryLine2Struct *boundary_ptr, IntegralLine2 *integral_ptr
    );


};

void PhysicsSolidHeatTransferSteady::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec
)
{

    // iterate through each domain covered by the mesh
    for (int indx_d = 0; indx_d < mesh_physics_ptr->num_domain; indx_d++)
    {

        // subset the mesh, boundary, and intergrals
        MeshLine2Struct *mesh_ptr = mesh_physics_ptr->mesh_ptr_vec[indx_d];
        BoundaryLine2Struct *boundary_ptr = boundary_physics_ptr->boundary_ptr_vec[indx_d];
        IntegralLine2 *integral_ptr = integral_physics_ptr->integral_ptr_vec[indx_d];

        // determine matrix coefficients for the domain
        matrix_fill_domain(a_mat, b_vec, x_vec, mesh_ptr, boundary_ptr, integral_ptr);

    }

}

void PhysicsSolidHeatTransferSteady::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    MeshLine2Struct *mesh_ptr, BoundaryLine2Struct *boundary_ptr, IntegralLine2 *integral_ptr
)
{

    // iterate for each domain element
    for (int indx_m = 0; indx_m < mesh_ptr->num_domain_element; indx_m++)
    {

        // get global id of points around element
        int n0 = mesh_ptr->element_p0_id_vec[indx_m];
        int n1 = mesh_ptr->element_p1_id_vec[indx_m];
        int n_arr[2] = {n0, n1};

        // get thermal conductivity of points around element
        double thermcond0 = thermalconductivity_field_ptr->point_u_map[n0];
        double thermcond1 = thermalconductivity_field_ptr->point_u_map[n1];
        double thermcond_arr[2] = {thermcond0, thermcond1};

        // get heat generation of points around element
        double heatgen0 = heatgeneration_field_ptr->point_u_map[n0];
        double heatgen1 = heatgeneration_field_ptr->point_u_map[n1];
        double heatgen_arr[2] = {heatgen0, heatgen1};

        // calculate a_mat coefficients
        // matrix row = start_row of physics + field ID of variable
        // matrix column = start_column of variable + field ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++){
        for (int indx_j = 0; indx_j < 2; indx_j++){
            int mat_row = start_row + temperature_field_ptr->point_global_field_id_map[n_arr[indx_i]];
            int mat_col = temperature_field_ptr->start_col + temperature_field_ptr->point_global_field_id_map[n_arr[indx_j]];
            a_mat.coeffRef(mat_row, mat_col) += thermcond_arr[indx_i]*integral_ptr->integral_div_Ni_line2_dot_div_Nj_line2_vec[indx_m][indx_i][indx_j];
        }}

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++)
        {
            int mat_row = start_row + temperature_field_ptr->point_global_field_id_map[n_arr[indx_i]];
            b_vec.coeffRef(mat_row) += heatgen_arr[indx_i]*integral_ptr->integral_Ni_line2_vec[indx_m][indx_i];
        }

    }

    // iterate for each flux boundary element
    for (int indx_k = 0; indx_k < boundary_ptr->num_domain_element_flux; indx_k++)
    {

        // get global id of element
        int m = boundary_ptr->element_flux_global_id_vec[indx_k];

        // get index of the element
        auto iter = std::find(mesh_ptr->element_global_id_vec.begin(), mesh_ptr->element_global_id_vec.end(), m); 
        int indx_m = iter - mesh_ptr->element_global_id_vec.begin();

        // get id of points around element
        int n0 = mesh_ptr->element_p0_id_vec[indx_m];
        int n1 = mesh_ptr->element_p1_id_vec[indx_m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = boundary_ptr->element_flux_pa_local_id_vec[indx_k];  // 0 or 1

        // identify boundary type
        int config_id = boundary_ptr->element_flux_config_id_vec[indx_k];
        BoundaryConfigLine2Struct bcl2 = boundary_ptr->boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2.boundary_type_str == "neumann")
        {
            // add to b_vec
            int mat_row = start_row + temperature_field_ptr->point_global_field_id_map[n_arr[a]];
            b_vec.coeffRef(mat_row) += bcl2.boundary_parameter_vec[0];
        }

    }

    // clear rows with value boundary elements
    for (int indx_k = 0; indx_k < boundary_ptr->num_domain_element_value; indx_k++)
    {

        // get global id of element
        int m = boundary_ptr->element_value_global_id_vec[indx_k];

        // get index of the element
        auto iter = std::find(mesh_ptr->element_global_id_vec.begin(), mesh_ptr->element_global_id_vec.end(), m); 
        int indx_m = iter - mesh_ptr->element_global_id_vec.begin();

        // get id of points around element
        int n0 = mesh_ptr->element_p0_id_vec[indx_m];
        int n1 = mesh_ptr->element_p1_id_vec[indx_m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = boundary_ptr->element_value_pa_local_id_vec[indx_k];  // 0 or 1

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            int mat_row = start_row + temperature_field_ptr->point_global_field_id_map[n_arr[a]];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int indx_k = 0; indx_k < boundary_ptr->num_domain_element_value; indx_k++)
    {

        // get global id of element
        int m = boundary_ptr->element_value_global_id_vec[indx_k];

        // get index of the element
        auto iter = std::find(mesh_ptr->element_global_id_vec.begin(), mesh_ptr->element_global_id_vec.end(), m); 
        int indx_m = iter - mesh_ptr->element_global_id_vec.begin();

        // get id of points around element
        int n0 = mesh_ptr->element_p0_id_vec[indx_m];
        int n1 = mesh_ptr->element_p1_id_vec[indx_m];
        int n_arr[2] = {n0, n1};

        // get points where the boundary is applied
        int a = boundary_ptr->element_value_pa_local_id_vec[indx_k];  // 0 or 1
        
        // identify boundary type
        int config_id = boundary_ptr->element_value_config_id_vec[indx_k];
        BoundaryConfigLine2Struct bcl2 = boundary_ptr->boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            int mat_row = start_row + temperature_field_ptr->point_global_field_id_map[n_arr[a]];
            int mat_col = temperature_field_ptr->start_col + temperature_field_ptr->point_global_field_id_map[n_arr[a]];
            if (a != -1)
            {
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                b_vec.coeffRef(mat_row) += bcl2.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
