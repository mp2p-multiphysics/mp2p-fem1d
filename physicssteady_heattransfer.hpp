#ifndef PHYSICSSTEADY_HEATTRANSFER
#define PHYSICSSTEADY_HEATTRANSFER
#include <vector>
#include "Eigen/Eigen"
#include "boundary_physicsgroup.hpp"
#include "container_typedef.hpp"
#include "integral_physicsgroup.hpp"
#include "mesh_physicsgroup.hpp"
#include "physicssteady_base.hpp"
#include "scalar_fieldgroup.hpp"
#include "variable_fieldgroup.hpp"

class PhysicsSteadyHeatTransfer : public PhysicsSteadyBase
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
    void set_start_row(int start_row_in);
    std::vector<VariableFieldGroup*> get_variable_field_ptr_vec();

    // default constructor
    PhysicsSteadyHeatTransfer()
    {

    }

    // constructor
    PhysicsSteadyHeatTransfer
    (
        MeshPhysicsGroup &mesh_physics_in, BoundaryPhysicsGroup &boundary_physics_in, IntegralPhysicsGroup &integral_physics_in,
        VariableFieldGroup &temperature_field_in, ScalarFieldGroup &thermalconductivity_field_in, ScalarFieldGroup &heatgeneration_field_in
    )
    {
        
        // store physics groups
        mesh_physics_ptr = &mesh_physics_in;
        boundary_physics_ptr = &boundary_physics_in;
        integral_physics_ptr = &integral_physics_in;

        // store field 
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

};

void PhysicsSteadyHeatTransfer::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec
)
{

    // iterate through each domain covered by the mesh
    for (int indx_d = 0; indx_d < mesh_physics_ptr->mesh_ptr_vec.size(); indx_d++)
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

void PhysicsSteadyHeatTransfer::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    MeshLine2Struct *mesh_ptr, BoundaryLine2Struct *boundary_ptr, IntegralLine2 *integral_ptr,
    ScalarLine2 *thermalconductivity_ptr, ScalarLine2 *heatgeneration_ptr
)
{

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++)
    {

        // get global ID of points around element
        int p0_gid = mesh_ptr->element_p0_gid_vec[element_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[element_did];

        // get domain ID of points
        // used for getting properties and integrals
        int p0_did = mesh_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_ptr->point_gid_to_did_map[p1_gid];
        int did_arr[2] = {p0_did, p1_did};

        // get thermal conductivity of points around element
        double thermcond_p0 = thermalconductivity_ptr->point_value_vec[p0_did];
        double thermcond_p1 = thermalconductivity_ptr->point_value_vec[p1_did];
        double thermcond_arr[2] = {thermcond_p0, thermcond_p1};

        // get heat generation of points around element
        double heatgen_p0 = heatgeneration_ptr->point_value_vec[p0_did];
        double heatgen_p1 = heatgeneration_ptr->point_value_vec[p1_did];
        double heatgen_arr[2] = {heatgen_p0, heatgen_p1};

        // calculate a_mat coefficients
        // matrix row = start_row of test function (physics) + field ID of variable
        // matrix column = start_column of variable + field ID of variable

        // get field ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = temperature_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = temperature_field_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++){
        for (int indx_j = 0; indx_j < 2; indx_j++){
            int mat_row = start_row + fid_arr[indx_i];
            int mat_col = temperature_field_ptr->start_col + fid_arr[indx_j];
            a_mat.coeffRef(mat_row, mat_col) += thermcond_arr[indx_i]*integral_ptr->integral_div_Ni_line2_dot_div_Nj_line2_vec[element_did][indx_i][indx_j];
        }}

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++)
        {
            int mat_row = start_row + fid_arr[indx_i];
            b_vec.coeffRef(mat_row) += heatgen_arr[indx_i]*integral_ptr->integral_Ni_line2_vec[element_did][indx_i];
        }

    }

    // iterate for each flux boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_flux_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_flux_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int ea_lid = boundary_ptr->element_flux_pa_lid_vec[boundary_id];  // 0 or 1

        // identify boundary type
        int config_id = boundary_ptr->element_flux_config_id_vec[boundary_id];
        BoundaryConfigLine2Struct bcl2 = boundary_ptr->boundary_config_vec[config_id];

        // get field ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = temperature_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = temperature_field_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // apply boundary condition
        if (bcl2.boundary_type_str == "neumann")
        {
            // add to b_vec
            int mat_row = start_row + fid_arr[ea_lid];
            b_vec.coeffRef(mat_row) += bcl2.boundary_parameter_vec[0];
        }

    }

    // clear rows with value boundary elements
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_value_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_value_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int ea_lid = boundary_ptr->element_value_pa_lid_vec[boundary_id];  // 0 or 1

        // get field ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = temperature_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = temperature_field_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // erase entire row
        // -1 values indicate invalid points
        if (ea_lid != -1)
        {
            int mat_row = start_row + fid_arr[ea_lid];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_value_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_value_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int ea_lid = boundary_ptr->element_value_pa_lid_vec[boundary_id];  // 0 or 1
        
        // identify boundary type
        int config_id = boundary_ptr->element_value_config_id_vec[boundary_id];
        BoundaryConfigLine2Struct bcl2 = boundary_ptr->boundary_config_vec[config_id];

        // get field ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = temperature_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = temperature_field_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // apply boundary condition
        if (bcl2.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            int mat_row = start_row + fid_arr[ea_lid];
            int mat_col = temperature_field_ptr->start_col + fid_arr[ea_lid];
            if (ea_lid != -1)
            {
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                b_vec.coeffRef(mat_row) += bcl2.boundary_parameter_vec[0];
            }

        }

    }

}

void PhysicsSteadyHeatTransfer::set_start_row(int start_row_in)
{
    start_row = start_row_in;
}

std::vector<VariableFieldGroup*> PhysicsSteadyHeatTransfer::get_variable_field_ptr_vec()
{
    return variable_field_ptr_vec;
}

#endif
