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

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);
    void matrix_store(Eigen::VectorXd &x_vec);

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
        MeshLine2Struct *mesh_ptr, BoundaryLine2Struct *boundary_ptr, IntegralLine2 *integral_ptr,
        ScalarLine2 *thermalconductivity_ptr, ScalarLine2 *heatgeneration_ptr
    );
    void matrix_store_domain
    (
        Eigen::VectorXd &x_vec, MeshLine2Struct *mesh_ptr, VariableLine2 *temperature_ptr
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

        // get scalar fields
        ScalarLine2 *thermalconductivity_ptr = thermalconductivity_field_ptr->scalar_ptr_map[mesh_ptr];
        ScalarLine2 *heatgeneration_ptr = heatgeneration_field_ptr->scalar_ptr_map[mesh_ptr];

        // determine matrix coefficients for the domain
        matrix_fill_domain(a_mat, b_vec, x_vec, mesh_ptr, boundary_ptr, integral_ptr, thermalconductivity_ptr, heatgeneration_ptr);

    }

}

void PhysicsSolidHeatTransferSteady::matrix_store
(
    Eigen::VectorXd &x_vec
)
{

    // iterate through each domain covered by the mesh
    for (int indx_d = 0; indx_d < mesh_physics_ptr->num_domain; indx_d++)
    {

        // subset the mesh
        MeshLine2Struct *mesh_ptr = mesh_physics_ptr->mesh_ptr_vec[indx_d];

        // get variable fields
        VariableLine2 *temperature_ptr = temperature_field_ptr->variable_ptr_map[mesh_ptr];

        // determine matrix coefficients for the domain
        matrix_store_domain(x_vec, mesh_ptr, temperature_ptr);

    }

}

void PhysicsSolidHeatTransferSteady::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    MeshLine2Struct *mesh_ptr, BoundaryLine2Struct *boundary_ptr, IntegralLine2 *integral_ptr,
    ScalarLine2 *thermalconductivity_ptr, ScalarLine2 *heatgeneration_ptr
)
{

    // iterate for each domain element
    for (int gid_e = 0; gid_e < mesh_ptr->num_domain_element; gid_e++)
    {

        // get global ID of points around element
        int gid_p0 = mesh_ptr->element_p0_gid_vec[gid_e];
        int gid_p1 = mesh_ptr->element_p1_gid_vec[gid_e];

        // get domain ID of points
        // used for getting properties and integrals
        int did_p0 = mesh_ptr->point_gid_to_did_map[gid_p0];
        int did_p1 = mesh_ptr->point_gid_to_did_map[gid_p1];
        int did_arr[2] = {did_p0, did_p1};

        // get field ID of points
        // used for getting matrix rows and columns
        int fid_p0 = temperature_field_ptr->point_gid_to_fid_map[gid_p0];
        int fid_p1 = temperature_field_ptr->point_gid_to_fid_map[gid_p1];
        int fid_arr[2] = {fid_p0, fid_p1};

        // get thermal conductivity of points around element
        double thermcond_p0 = thermalconductivity_ptr->point_value_vec[did_p0];
        double thermcond_p1 = thermalconductivity_ptr->point_value_vec[did_p1];
        double thermcond_arr[2] = {thermcond_p0, thermcond_p1};

        // get heat generation of points around element
        double heatgen_p0 = heatgeneration_ptr->point_value_vec[did_p0];
        double heatgen_p1 = heatgeneration_ptr->point_value_vec[did_p1];
        double heatgen_arr[2] = {heatgen_p0, heatgen_p1};

        // calculate a_mat coefficients
        // matrix row = start_row of test function (physics) + field ID of variable
        // matrix column = start_column of variable + field ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++){
        for (int indx_j = 0; indx_j < 2; indx_j++){
            int mat_row = start_row + fid_arr[indx_i];
            int mat_col = temperature_field_ptr->start_col + fid_arr[indx_j];
            a_mat.coeffRef(mat_row, mat_col) += thermcond_arr[indx_i]*integral_ptr->integral_div_Ni_line2_dot_div_Nj_line2_vec[gid_e][indx_i][indx_j];
        }}

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++)
        {
            int mat_row = start_row + fid_arr[indx_i];
            b_vec.coeffRef(mat_row) += heatgen_arr[indx_i]*integral_ptr->integral_Ni_line2_vec[gid_e][indx_i];
        }

    }

    // iterate for each flux boundary element
    for (int indx_k = 0; indx_k < boundary_ptr->num_domain_element_flux; indx_k++)
    {

        // get global ID of element
        int gid_ea = boundary_ptr->element_flux_gid_vec[indx_k];

        // get domain ID of element
        // used for getting global ID of points
        int did_ea = mesh_ptr->element_gid_to_did_map[gid_ea];

        // get global ID of points
        int gid_p0 = mesh_ptr->element_p0_gid_vec[did_ea];
        int gid_p1 = mesh_ptr->element_p1_gid_vec[did_ea];

        // get field ID of points
        // used for getting matrix rows and columns
        int fid_p0 = temperature_field_ptr->point_gid_to_fid_map[gid_p0];
        int fid_p1 = temperature_field_ptr->point_gid_to_fid_map[gid_p1];
        int fid_arr[2] = {fid_p0, fid_p1};

        // get local ID of point where boundary is applied
        int lid_pa = boundary_ptr->element_flux_pa_lid_vec[indx_k];  // 0 or 1

        // identify boundary type
        int config_id = boundary_ptr->element_flux_config_id_vec[indx_k];
        BoundaryConfigLine2Struct bcl2 = boundary_ptr->boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2.boundary_type_str == "neumann")
        {
            // add to b_vec
            int mat_row = start_row + fid_arr[lid_pa];
            b_vec.coeffRef(mat_row) += bcl2.boundary_parameter_vec[0];
        }

    }

    // clear rows with value boundary elements
    for (int indx_k = 0; indx_k < boundary_ptr->num_domain_element_value; indx_k++)
    {

        // get global ID of element
        int gid_ea = boundary_ptr->element_value_gid_vec[indx_k];

        // get domain ID of element
        // used for getting global ID of points
        int did_ea = mesh_ptr->element_gid_to_did_map[gid_ea];

        // get global ID of points
        int gid_p0 = mesh_ptr->element_p0_gid_vec[did_ea];
        int gid_p1 = mesh_ptr->element_p1_gid_vec[did_ea];

        // get field ID of points
        // used for getting matrix rows and columns
        int fid_p0 = temperature_field_ptr->point_gid_to_fid_map[gid_p0];
        int fid_p1 = temperature_field_ptr->point_gid_to_fid_map[gid_p1];
        int fid_arr[2] = {fid_p0, fid_p1};

        // get local ID of point where boundary is applied
        int lid_pa = boundary_ptr->element_value_pa_lid_vec[indx_k];  // 0 or 1

        // erase entire row
        // -1 values indicate invalid points
        if (lid_pa != -1)
        {
            int mat_row = start_row + fid_arr[lid_pa];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int indx_k = 0; indx_k < boundary_ptr->num_domain_element_value; indx_k++)
    {

        // get global ID of element
        int gid_ea = boundary_ptr->element_value_gid_vec[indx_k];

        // get domain ID of element
        // used for getting global ID of points
        int did_ea = mesh_ptr->element_gid_to_did_map[gid_ea];

        // get global ID of points
        int gid_p0 = mesh_ptr->element_p0_gid_vec[did_ea];
        int gid_p1 = mesh_ptr->element_p1_gid_vec[did_ea];

        // get field ID of points
        // used for getting matrix rows and columns
        int fid_p0 = temperature_field_ptr->point_gid_to_fid_map[gid_p0];
        int fid_p1 = temperature_field_ptr->point_gid_to_fid_map[gid_p1];
        int fid_arr[2] = {fid_p0, fid_p1};

        // get local ID of point where boundary is applied
        int lid_pa = boundary_ptr->element_value_pa_lid_vec[indx_k];  // 0 or 1
        
        // identify boundary type
        int config_id = boundary_ptr->element_value_config_id_vec[indx_k];
        BoundaryConfigLine2Struct bcl2 = boundary_ptr->boundary_config_vec[config_id];

        // apply boundary condition
        if (bcl2.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            int mat_row = start_row + fid_arr[lid_pa];
            int mat_col = temperature_field_ptr->start_col + fid_arr[lid_pa];
            if (lid_pa != -1)
            {
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                b_vec.coeffRef(mat_row) += bcl2.boundary_parameter_vec[0];
            }

        }

    }

}

void PhysicsSolidHeatTransferSteady::matrix_store_domain
(
    Eigen::VectorXd &x_vec, MeshLine2Struct *mesh_ptr, VariableLine2 *temperature_ptr
)
{

    // iterate for each temperature value
    for (int fid_p = 0; fid_p < temperature_field_ptr->num_field_point; fid_p++)
    {
        
        // get temperature from x vector
        int mat_row = temperature_field_ptr->start_col + fid_p;
        double temperature_value = x_vec.coeffRef(mat_row);

        // store in variable
        temperature_field_ptr->point_value_vec.push_back(temperature_value);

    }

}

#endif
