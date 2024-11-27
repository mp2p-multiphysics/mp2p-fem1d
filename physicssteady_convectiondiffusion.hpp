#ifndef PHYSICSSTEADY_CONVECTIONDIFFUSION
#define PHYSICSSTEADY_CONVECTIONDIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "boundary_group.hpp"
#include "container_typedef.hpp"
#include "integral_group.hpp"
#include "domain_group.hpp"
#include "physicssteady_base.hpp"
#include "scalar_group.hpp"
#include "variable_group.hpp"

class PhysicsSteadyConvectionDiffusion : public PhysicsSteadyBase
{
    /*

    Single-component steady-state convection-diffusion equation.    
    
    0 = -div(-b * grad(u)) - v * grad(u) + c

    Variables
    =========
    domain_group_in : DomainGroup
        Domains where this physics is applied to.
    boundary_group_in : BoundaryGroup
        Boundaries where this physics is applied to.
    integral_group_in : IntegralGroup
        Test function integrals that this physics uses.
    value_group_in : VariableGroup
        u in 0 = -div(-b * grad(u)) - v * grad(u) + c.
        This will be solved for by the matrix equation.
    diffusioncoefficient_group_in : ScalarGroup
        b in 0 = -div(-b * grad(u)) - v * grad(u) + c.
    velocity_x_group_in : ScalarGroup
        v in 0 = -div(-b * grad(u)) - v * grad(u) + c.
    generationcoefficient_group_in : ScalarGroup
        c in 0 = -div(-b * grad(u)) - v * grad(u) + c.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_start_row : void
        Sets the starting row in A and b where entries are filled up.
    get_start_row : int
        Returns the starting row.
    get_boundary_group_ptr_vec : vector<BoundaryGroup*>
        Returns the vector containing pointers to BoundaryGroup objects tied to this physics.
    get_scalar_group_ptr_vec : vector<ScalarGroup*>
        Returns the vector containing pointers to ScalarGroup objects tied to this physics.
    get_variable_group_ptr_vec : vector<VariableGroup*>
        Returns the vector containing pointers to VariableGroup objects tied to this physics.

    */

    public:

    // variables
    DomainGroup *domain_group_ptr;
    BoundaryGroup *boundary_group_ptr;
    IntegralGroup *integral_group_ptr;
    VariableGroup *value_group_ptr;
    ScalarGroup *diffusioncoefficient_group_ptr;
    ScalarGroup *velocity_x_group_ptr;
    ScalarGroup *generationcoefficient_group_ptr;

    // vector of scalar and variable groups
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);
    void set_start_row(int start_row_in);
    virtual int get_start_row();
    BoundaryGroup* get_boundary_group_ptr();
    std::vector<ScalarGroup*> get_scalar_group_ptr_vec();
    std::vector<VariableGroup*> get_variable_group_ptr_vec();

    // default constructor
    PhysicsSteadyConvectionDiffusion() {}

    // constructor
    PhysicsSteadyConvectionDiffusion
    (
        DomainGroup &domain_group_in, BoundaryGroup &boundary_group_in, IntegralGroup &integral_group_in,
        VariableGroup &value_group_in,
        ScalarGroup &diffusioncoefficient_group_in, ScalarGroup &velocity_x_group_in, ScalarGroup &generationcoefficient_group_in
    )
    {
        
        // store variables
        domain_group_ptr = &domain_group_in;
        boundary_group_ptr = &boundary_group_in;
        integral_group_ptr = &integral_group_in;
        value_group_ptr = &value_group_in;
        diffusioncoefficient_group_ptr = &diffusioncoefficient_group_in;
        velocity_x_group_ptr = &velocity_x_group_in;
        generationcoefficient_group_ptr = &generationcoefficient_group_in;

        // set boundary conditions
        boundary_group_ptr->set_boundary_type({0}, {1, 2});

        // vector of scalar and variable groups
        scalar_group_ptr_vec = {diffusioncoefficient_group_ptr, velocity_x_group_ptr, generationcoefficient_group_ptr};
        variable_group_ptr_vec = {value_group_ptr};

        // calculate integrals
        integral_group_ptr->evaluate_integral_div_Ni_dot_div_Nj();
        integral_group_ptr->evaluate_integral_Ni_derivative_Nj_x();
        integral_group_ptr->evaluate_integral_Ni();
        integral_group_ptr->evaluate_integral_boundary_Ni();
        integral_group_ptr->evaluate_integral_boundary_Ni_Nj();

    }

    private:
    void matrix_fill_domain
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainLine2 *domain_ptr, IntegralLine2 *integral_ptr,
        ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *velocity_x_ptr, ScalarLine2 *generationcoefficient_ptr
    );
    void matrix_fill_natural
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr,  IntegralLine2 *integral_ptr,
        ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *velocity_x_ptr, ScalarLine2 *generationcoefficient_ptr
    );
    void matrix_fill_essential_clear
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr
    );
    void matrix_fill_essential
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr,  IntegralLine2 *integral_ptr,
        ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *velocity_x_ptr, ScalarLine2 *generationcoefficient_ptr
    );

};

void PhysicsSteadyConvectionDiffusion::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec
)
{
    /*

    Fill up the matrix equation Ax = b with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : Eigen::SparseMatrix<double>
        A in Ax = b.
    b_vec : Eigen::VectorXd
        b in Ax = b.
    x_vec : Eigen::VectorXd
        x in Ax = b.

    Returns
    =======
    (none)

    */

   // iterate through each domain
   for (int indx_d = 0; indx_d < domain_group_ptr->domain_l2_ptr_vec.size(); indx_d++)
   {

        // subset the domain and integrals
        DomainLine2 *domain_ptr = domain_group_ptr->domain_l2_ptr_vec[indx_d];
        IntegralLine2 *integral_ptr = integral_group_ptr->integral_l2_ptr_vec[indx_d];

        // subset the scalars
        ScalarLine2 *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_l2_ptr_vec[indx_d];
        ScalarLine2 *velocity_x_ptr = velocity_x_group_ptr->scalar_l2_ptr_vec[indx_d];
        ScalarLine2 *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_l2_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(a_mat, b_vec, x_vec, domain_ptr, integral_ptr, diffusioncoefficient_ptr, velocity_x_ptr, generationcoefficient_ptr);

   }

   // iterate through each natural boundary
   for (int indx_b = 0; indx_b < boundary_group_ptr->boundary_l2_ptr_vec.size(); indx_b++)
   {

        // subset the domain and integrals
        DomainLine2 *domain_ptr = domain_group_ptr->domain_l2_ptr_vec[indx_b];
        BoundaryLine2 *boundary_ptr = boundary_group_ptr->boundary_l2_ptr_vec[indx_b];
        IntegralLine2 *integral_ptr = integral_group_ptr->integral_l2_ptr_vec[indx_b];

        // subset the scalars
        ScalarLine2 *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_l2_ptr_vec[indx_b];
        ScalarLine2 *velocity_x_ptr = velocity_x_group_ptr->scalar_l2_ptr_vec[indx_b];
        ScalarLine2 *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_l2_ptr_vec[indx_b];

        // fill up matrix with boundary conditions
        matrix_fill_natural(a_mat, b_vec, x_vec, domain_ptr, boundary_ptr, integral_ptr, diffusioncoefficient_ptr, velocity_x_ptr, generationcoefficient_ptr);

   }

   // clear equations with essential boundary conditions
   for (int indx_b = 0; indx_b < boundary_group_ptr->boundary_l2_ptr_vec.size(); indx_b++)
   {

        // subset the domain and integrals
        DomainLine2 *domain_ptr = domain_group_ptr->domain_l2_ptr_vec[indx_b];
        BoundaryLine2 *boundary_ptr = boundary_group_ptr->boundary_l2_ptr_vec[indx_b];

        // fill up matrix with boundary conditions
        matrix_fill_essential_clear(a_mat, b_vec, x_vec, domain_ptr, boundary_ptr);

   }

   // iterate through each essential boundary
   for (int indx_b = 0; indx_b < boundary_group_ptr->boundary_l2_ptr_vec.size(); indx_b++)
   {

        // subset the domain and integrals
        DomainLine2 *domain_ptr = domain_group_ptr->domain_l2_ptr_vec[indx_b];
        BoundaryLine2 *boundary_ptr = boundary_group_ptr->boundary_l2_ptr_vec[indx_b];
        IntegralLine2 *integral_ptr = integral_group_ptr->integral_l2_ptr_vec[indx_b];

        // subset the scalars
        ScalarLine2 *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_l2_ptr_vec[indx_b];
        ScalarLine2 *velocity_x_ptr = velocity_x_group_ptr->scalar_l2_ptr_vec[indx_b];
        ScalarLine2 *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_l2_ptr_vec[indx_b];

        // fill up matrix with boundary conditions
        matrix_fill_essential(a_mat, b_vec, x_vec, domain_ptr, boundary_ptr, integral_ptr, diffusioncoefficient_ptr, velocity_x_ptr, generationcoefficient_ptr);

   }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainLine2 *domain_ptr, IntegralLine2 *integral_ptr,
    ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *velocity_x_ptr, ScalarLine2 *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get points around element
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];

        // get diffusion coefficient of points around element
        double diffcoeff_p0 = diffusioncoefficient_ptr->point_value_vec[p0_pdid];
        double diffcoeff_p1 = diffusioncoefficient_ptr->point_value_vec[p1_pdid];
        double diffcoeff_arr[2] = {diffcoeff_p0, diffcoeff_p1};

        // get velocity of points around element
        double velx_p0 = velocity_x_ptr->point_value_vec[p0_pdid];
        double velx_p1 = velocity_x_ptr->point_value_vec[p1_pdid];
        double velx_arr[2] = {velx_p0, velx_p1};

        // get generation coefficient of points around element
        double gencoeff_p0 = generationcoefficient_ptr->point_value_vec[p0_pdid];
        double gencoeff_p1 = generationcoefficient_ptr->point_value_vec[p1_pdid];
        double gencoeff_arr[2] = {gencoeff_p0, gencoeff_p1};

        // calculate a_mat coefficients
        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // get group ID of value points
        // used for getting matrix rows and columns
        int p0_pfid = value_group_ptr->point_pgid_to_pfid_map[p0_pgid];
        int p1_pfid = value_group_ptr->point_pgid_to_pfid_map[p1_pgid];
        int pfid_arr[2] = {p0_pfid, p1_pfid};

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++){
        for (int indx_j = 0; indx_j < 2; indx_j++){
            int mat_row = start_row + pfid_arr[indx_i];
            int mat_col = value_group_ptr->start_col + pfid_arr[indx_j];
            a_mat.coeffRef(mat_row, mat_col) += (
                diffcoeff_arr[indx_i]*integral_ptr->integral_div_Ni_dot_div_Nj_vec[edid][indx_i][indx_j] +
                velx_arr[indx_i]*integral_ptr->integral_Ni_derivative_Nj_x_vec[edid][indx_i][indx_j]
            );
        }}

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < 2; indx_i++)
        {
            int mat_row = start_row + pfid_arr[indx_i];
            b_vec.coeffRef(mat_row) += gencoeff_arr[indx_i]*integral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_natural
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr, IntegralLine2 *integral_ptr,
    ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *velocity_x_ptr, ScalarLine2 *generationcoefficient_ptr
)
{

    // iterate for each natural boundary element
    for (int bnid = 0; bnid < boundary_ptr->num_natural; bnid++)
    {

        // get ID of element
        int egid = boundary_ptr->natural_egid_vec[bnid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get global ID of points
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];
        int pgid_arr[2] = {p0_pgid, p1_pgid};

        // get point where boundary is applied
        int pa_plid = boundary_ptr->natural_pa_plid_vec[bnid];
        int pa_pgid = pgid_arr[pa_plid];
        int pa_pfid = value_group_ptr->point_pgid_to_pfid_map[pa_pgid];

        // identify boundary location, type, and parameters
        int boundary_key = pa_plid;
        int btid = boundary_ptr->natural_btid_vec[bnid];
        VectorDouble parameter_vec = boundary_ptr->natural_parameter_vec[bnid];

        // apply boundary condition
        int mat_row_pa = start_row + pa_pfid;
        int mat_col_pa = value_group_ptr->start_col + pa_pfid;
        switch (btid)
        {
            case 1:  // neumann
                b_vec.coeffRef(mat_row_pa) += parameter_vec[0];  // * integral_ptr->integral_boundary_Ni_vec[edid][boundary_key][pa_plid];
            break;
            case 2:  // robin
                a_mat.coeffRef(mat_row_pa, mat_col_pa) += -parameter_vec[1];  // * integral_ptr->integral_boundary_Ni_Nj_vec[edid][boundary_key][pa_plid][...];
                b_vec.coeffRef(mat_row_pa) += parameter_vec[0];  // * integral_ptr->integral_boundary_Ni_vec[edid][boundary_key][pa_plid];
            break;
        }
        
    }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_essential_clear
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr
)
{

    // iterate for each essential boundary element
    for (int beid = 0; beid < boundary_ptr->num_essential; beid++)
    {

        // get ID of element
        int egid = boundary_ptr->essential_egid_vec[beid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get global ID of points
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];
        int pgid_arr[2] = {p0_pgid, p1_pgid};

        // get point where boundary is applied
        int pa_plid = boundary_ptr->essential_pa_plid_vec[beid];
        int pa_pgid = pgid_arr[pa_plid];
        int pa_pfid = value_group_ptr->point_pgid_to_pfid_map[pa_pgid];

        // erase row
        int mat_row_pa = start_row + pa_pfid;
        a_mat.row(mat_row_pa) *= 0.;
        b_vec.coeffRef(mat_row_pa) = 0.;

    }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_essential
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainLine2 *domain_ptr, BoundaryLine2 *boundary_ptr, IntegralLine2 *integral_ptr,
    ScalarLine2 *diffusioncoefficient_ptr, ScalarLine2 *velocity_x_ptr, ScalarLine2 *generationcoefficient_ptr
)
{

    // iterate for each essential boundary element
    for (int beid = 0; beid < boundary_ptr->num_essential; beid++)
    {

        // get ID of element
        int egid = boundary_ptr->essential_egid_vec[beid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get global ID of points
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];
        int pgid_arr[2] = {p0_pgid, p1_pgid};

        // get point where boundary is applied
        int pa_plid = boundary_ptr->essential_pa_plid_vec[beid];
        int pa_pgid = pgid_arr[pa_plid];
        int pa_pfid = value_group_ptr->point_pgid_to_pfid_map[pa_pgid];

        // identify boundary location, type, and parameters
        int boundary_key = pa_plid;
        int btid = boundary_ptr->essential_btid_vec[beid];
        VectorDouble parameter_vec = boundary_ptr->essential_parameter_vec[beid];

        // apply boundary condition
        int mat_row_pa = start_row + pa_pfid;
        int mat_col_pa = value_group_ptr->start_col + pa_pfid;
        switch (btid)
        {
            case 0:  // dirichlet
                a_mat.coeffRef(mat_row_pa, mat_col_pa) += 1.;
                b_vec.coeffRef(mat_row_pa) += parameter_vec[0];
            break;
        }

    }

}

void PhysicsSteadyConvectionDiffusion::set_start_row(int start_row_in)
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

int PhysicsSteadyConvectionDiffusion::get_start_row()
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

BoundaryGroup* PhysicsSteadyConvectionDiffusion::get_boundary_group_ptr()
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

std::vector<ScalarGroup*> PhysicsSteadyConvectionDiffusion::get_scalar_group_ptr_vec()
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

std::vector<VariableGroup*> PhysicsSteadyConvectionDiffusion::get_variable_group_ptr_vec()
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
