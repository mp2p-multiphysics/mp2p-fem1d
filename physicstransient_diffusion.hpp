#ifndef PHYSICSTRANSIENT_DIFFUSION
#define PHYSICSTRANSIENT_DIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "boundary_group.hpp"
#include "container_typedef.hpp"
#include "integral_group.hpp"
#include "domain_group.hpp"
#include "physicstransient_base.hpp"
#include "scalar_group.hpp"
#include "variable_group.hpp"

class PhysicsTransientDiffusion : public PhysicsTransientBase
{
    /*

    Single-component transient diffusion equation.    
    
    a * du/dt = -div(-b * grad(u)) + c

    Variables
    =========
    domain_group_in : DomainGroup
        Domaines where this physics is applied to.
    boundary_group_in : BoundaryGroup
        Boundary conditions pertinent to this physics.
    integral_group_in : IntegralGroup
        Test function integrals of the domains.
    value_group_in : VariableGroup
        u in a * du/dt = -div(-b * grad(u)) + c.
        This will be solved for by the matrix equation.
    derivativecoefficient_group_ptr : ScalarGroup
        a in a * du/dt = -div(-b * grad(u)) + c.
    diffusioncoefficient_group_in : ScalarGroup
        b in a * du/dt = -div(-b * grad(u)) + c.
    generationcoefficient_group_in : ScalarGroup
        c in a * du/dt = -div(-b * grad(u)) + c.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_start_row : void
        Sets the starting row in A and b where entries are filled up.
    get_start_row : int
        Returns the starting row.
    get_boundary_group_ptr_vec() : vector<BoundaryGroup*>
        Returns the vector containing pointers to BoundaryGroup objects tied to this physics.
    get_scalar_group_ptr_vec() : vector<ScalarGroup*>
        Returns the vector containing pointers to ScalarGroup objects tied to this physics.
    get_variable_group_ptr_vec() : vector<VariableGroup*>
        Returns the vector containing pointers to VariableGroup objects tied to this physics.

    */

    public:

    // variables
    DomainGroup *domain_group_ptr;
    BoundaryGroup *boundary_group_ptr;
    IntegralGroup *integral_group_ptr;
    VariableGroup *value_group_ptr;
    ScalarGroup *derivativecoefficient_group_ptr;
    ScalarGroup *diffusioncoefficient_group_ptr;
    ScalarGroup *generationcoefficient_group_ptr;

    // vector of scalar and variable groups
    BoundaryGroup* get_boundary_group_ptr();
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt
    );
    void set_start_row(int start_row_in);
    int get_start_row();
    std::vector<ScalarGroup*> get_scalar_group_ptr_vec();
    std::vector<VariableGroup*> get_variable_group_ptr_vec();

    // default constructor
    PhysicsTransientDiffusion()
    {

    }

    // constructor
    PhysicsTransientDiffusion
    (
        DomainGroup &domain_group_in, BoundaryGroup &boundary_group_in, IntegralGroup &integral_group_in,
        VariableGroup &value_group_in,
        ScalarGroup &derivativecoefficient_group_in, ScalarGroup &diffusioncoefficient_group_in, ScalarGroup &generationcoefficient_group_in
    )
    {
        
        // store variables
        domain_group_ptr = &domain_group_in;
        boundary_group_ptr = &boundary_group_in;
        integral_group_ptr = &integral_group_in;
        value_group_ptr = &value_group_in;
        derivativecoefficient_group_ptr = &derivativecoefficient_group_in;
        diffusioncoefficient_group_ptr = &diffusioncoefficient_group_in;
        generationcoefficient_group_ptr = &generationcoefficient_group_in;

        // vector of variable groups
        scalar_group_ptr_vec = {derivativecoefficient_group_ptr, diffusioncoefficient_group_ptr, generationcoefficient_group_ptr};
        variable_group_ptr_vec = {value_group_ptr};

        // calculate integrals
        integral_group_ptr->evaluate_Ni_derivative();
        integral_group_ptr->evaluate_integral_div_Ni_dot_div_Nj();
        integral_group_ptr->evaluate_integral_Ni_Nj();
        integral_group_ptr->evaluate_integral_Ni();

    }

    private:
    void matrix_fill_domain
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr, IntegralLine2 *integral_ptr,
        ScalarLine2 *derivativecoefficient_ptr, ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *generationcoefficient_ptr
    );

};

void PhysicsTransientDiffusion::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt
)
{
    /*

    Fill up the matrix equation Ax(t+1) = Cx(t) + d with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : Eigen::SparseMatrix<double>
        A in Ax(t+1) = Cx(t) + d.
    c_mat : Eigen::SparseMatrix<double>
        C in Ax(t+1) = Cx(t) + d.
    d_vec : Eigen::VectorXd
        d in Ax(t+1) = Cx(t) + d.
    x_vec : Eigen::VectorXd
        x(t+1) in Ax(t+1) = Cx(t) + d.
    x_last_timestep_vec : Eigen::VectorXd
        x(t) in Ax(t+1) = Cx(t) + d.
    dt : double
        Length of the timestep.

    Returns
    =======
    (none)

    */

    // iterate through each domain covered by the domain
    for (int indx_d = 0; indx_d < domain_group_ptr->domain_l2_ptr_vec.size(); indx_d++)
    {

        // subset the domain, boundary, and intergrals
        DomainLine2 *domain_ptr = domain_group_ptr->domain_l2_ptr_vec[indx_d];
        BoundaryLine2 *boundary_ptr = boundary_group_ptr->boundary_l2_ptr_vec[indx_d];
        IntegralLine2 *integral_ptr = integral_group_ptr->integral_l2_ptr_vec[indx_d];

        // get scalar groups
        ScalarLine2 *derivativecoefficient_ptr = derivativecoefficient_group_ptr->scalar_ptr_map[domain_ptr];
        ScalarLine2 *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_map[domain_ptr];
        ScalarLine2 *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_map[domain_ptr];

        // determine matrix coefficients for the domain
        matrix_fill_domain(
            a_mat, c_mat, d_vec,
            x_vec, x_last_timestep_vec, dt,
            domain_ptr, boundary_ptr, integral_ptr,
            derivativecoefficient_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr
        );

    }

}

void PhysicsTransientDiffusion::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr, IntegralLine2 *integral_ptr,
    ScalarLine2 *derivativecoefficient_ptr, ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int element_did = 0; element_did < domain_ptr->num_element_domain; element_did++)
    {

        // get global ID of points around element
        int p0_gid = domain_ptr->element_p0_gid_vec[element_did];
        int p1_gid = domain_ptr->element_p1_gid_vec[element_did];

        // get domain ID of points
        // used for getting properties and integrals
        int p0_did = domain_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = domain_ptr->point_gid_to_did_map[p1_gid];
        int did_arr[2] = {p0_did, p1_did};

        // get derivative coefficient of points around element
        double dervcoeff_p0 = derivativecoefficient_ptr->point_value_vec[p0_did];
        double dervcoeff_p1 = derivativecoefficient_ptr->point_value_vec[p1_did];
        double dervcoeff_arr[2] = {dervcoeff_p0, dervcoeff_p1};

        // get diffusion coefficient of points around element
        double diffcoeff_p0 = diffusioncoefficient_ptr->point_value_vec[p0_did];
        double diffcoeff_p1 = diffusioncoefficient_ptr->point_value_vec[p1_did];
        double diffcoeff_arr[2] = {diffcoeff_p0, diffcoeff_p1};

        // get generation coefficient of points around element
        double gencoeff_p0 = generationcoefficient_ptr->point_value_vec[p0_did];
        double gencoeff_p1 = generationcoefficient_ptr->point_value_vec[p1_did];
        double gencoeff_arr[2] = {gencoeff_p0, gencoeff_p1};

        // calculate a_mat coefficients
        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // get group ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = value_group_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_group_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // calculate a_mat and c_mat coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++){
        for (int indx_j = 0; indx_j < 2; indx_j++){
            
            // calculate matrix indices
            int mat_row = start_row + fid_arr[indx_i];
            int mat_col = value_group_ptr->start_col + fid_arr[indx_j];

            // calculate a_mat coefficients
            a_mat.coeffRef(mat_row, mat_col) += (
                (dervcoeff_arr[indx_i]/dt)*integral_ptr->integral_Ni_Nj_vec[element_did][indx_i][indx_j] +
                diffcoeff_arr[indx_i]*integral_ptr->integral_div_Ni_dot_div_Nj_vec[element_did][indx_i][indx_j]
            );

            // calculate c_mat coefficients
            c_mat.coeffRef(mat_row, mat_col) += (dervcoeff_arr[indx_i]/dt)*integral_ptr->integral_Ni_Nj_vec[element_did][indx_i][indx_j];

        }}

        // calculate d_vec coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++)
        {
            int mat_row = start_row + fid_arr[indx_i];
            d_vec.coeffRef(mat_row) += gencoeff_arr[indx_i]*integral_ptr->integral_Ni_vec[element_did][indx_i];
        }

    }

    // iterate for each flux boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_boundary_flux_domain; boundary_id++)
    {

        // get global ID of element
        int e_gid = boundary_ptr->boundary_flux_element_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int e_did = domain_ptr->element_gid_to_did_map[e_gid];

        // get global ID of points
        int p0_gid = domain_ptr->element_p0_gid_vec[e_did];
        int p1_gid = domain_ptr->element_p1_gid_vec[e_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->boundary_flux_pa_lid_vec[boundary_id];  // 0 or 1

        // get group ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = value_group_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_group_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // identify boundary type and parameters
        std::string type_str = boundary_ptr->boundary_flux_type_str_vec[boundary_id];
        VectorDouble parameter_vec = boundary_ptr->boundary_flux_parameter_vec[boundary_id];

        // apply boundary condition
        if (type_str == "neumann")
        {
            int mat_row = start_row + fid_arr[pa_lid];
            d_vec.coeffRef(mat_row) += parameter_vec[0];
        }
        else if (type_str == "robin")
        {
            int mat_row = start_row + fid_arr[pa_lid];
            int mat_col = value_group_ptr->start_col + fid_arr[pa_lid];
            d_vec.coeffRef(mat_row) += parameter_vec[0];
            a_mat.coeffRef(mat_row, mat_col) += -parameter_vec[1];
        }

    }

    // clear rows with value boundary elements
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_boundary_value_domain; boundary_id++)
    {

        // get global ID of element
        int e_gid = boundary_ptr->boundary_value_element_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int e_did = domain_ptr->element_gid_to_did_map[e_gid];

        // get global ID of points
        int p0_gid = domain_ptr->element_p0_gid_vec[e_did];
        int p1_gid = domain_ptr->element_p1_gid_vec[e_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->boundary_value_pa_lid_vec[boundary_id];  // 0 or 1

        // get group ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = value_group_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_group_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // erase entire row
        // -1 values indicate invalid points
        if (pa_lid != -1)
        {
            int mat_row = start_row + fid_arr[pa_lid];
            a_mat.row(mat_row) *= 0.;
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_boundary_value_domain; boundary_id++)
    {

        // get global ID of element
        int e_gid = boundary_ptr->boundary_value_element_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int e_did = domain_ptr->element_gid_to_did_map[e_gid];

        // get global ID of points
        int p0_gid = domain_ptr->element_p0_gid_vec[e_did];
        int p1_gid = domain_ptr->element_p1_gid_vec[e_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->boundary_value_pa_lid_vec[boundary_id];  // 0 or 1
        
        // identify boundary type and parameters
        std::string type_str = boundary_ptr->boundary_value_type_str_vec[boundary_id];
        VectorDouble parameter_vec = boundary_ptr->boundary_value_parameter_vec[boundary_id];

        // get group ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = value_group_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_group_ptr->point_gid_to_fid_map[p1_gid];
        int fid_arr[2] = {p0_fid, p1_fid};

        // apply boundary condition
        if (type_str == "dirichlet")
        {
            int mat_row = start_row + fid_arr[pa_lid];
            int mat_col = value_group_ptr->start_col + fid_arr[pa_lid];
            if (pa_lid != -1)  // -1 values indicate invalid points
            {
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                d_vec.coeffRef(mat_row) += parameter_vec[0];
            }
        }

    }

}

void PhysicsTransientDiffusion::set_start_row(int start_row_in)
{
    /*

    Sets the starting row in A and b where entries are filled up.

    Arguments
    =========
    start_row_in : int
        Starting row in A and b.

    Returns
    =======
    (none)

    */
    
    start_row = start_row_in;

}

int PhysicsTransientDiffusion::get_start_row()
{
    /*

    Returns the starting row.

    Arguments
    =========
    (none)

    Returns
    =======
    start_row : int
        Starting row in A and b.

    */

    return start_row;

}

BoundaryGroup* PhysicsTransientDiffusion::get_boundary_group_ptr()
{
    /*

    Returns the pointer to the BoundaryGroup object tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    boundary_group_ptr : BoundaryGroup*
        Pointer to BoundaryGroup object.

    */
    
    return boundary_group_ptr;

}

std::vector<ScalarGroup*> PhysicsTransientDiffusion::get_scalar_group_ptr_vec()
{
    /*

    Returns the vector containing pointers to ScalarGroup objects tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    scalar_group_ptr : vector<ScalarGroup*>
        Vector containing pointers to ScalarGroup objects.

    */
    
    return scalar_group_ptr_vec;

}

std::vector<VariableGroup*> PhysicsTransientDiffusion::get_variable_group_ptr_vec()
{
    /*

    Returns the vector containing pointers to VariableGroup objects tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    variable_group_ptr : vector<VariableGroup*>
        Vector containing pointers to VariableGroup objects.

    */
    
    return variable_group_ptr_vec;

}

#endif
