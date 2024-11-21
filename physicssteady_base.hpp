#ifndef PHYSICSSTEADY_BASE
#define PHYSICSSTEADY_BASE
#include <vector>
#include "Eigen/Eigen"
#include "boundary_group.hpp"
#include "container_typedef.hpp"
#include "domain_group.hpp"
#include "scalar_group.hpp"
#include "integral_group.hpp"
#include "variable_group.hpp"

class PhysicsSteadyBase
{
    /*

    Base class for steady-state physics.

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

    // vector of scalar and variable groups
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    virtual void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);
    virtual void set_start_row(int start_row_in);
    virtual int get_start_row();
    virtual BoundaryGroup* get_boundary_group_ptr();
    virtual std::vector<ScalarGroup*> get_scalar_group_ptr_vec();
    virtual std::vector<VariableGroup*> get_variable_group_ptr_vec();

    // default constructor
    PhysicsSteadyBase() {}

};

void PhysicsSteadyBase::matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec)
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

}

void PhysicsSteadyBase::set_start_row(int start_row_in)
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

int PhysicsSteadyBase::get_start_row()
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

BoundaryGroup* PhysicsSteadyBase::get_boundary_group_ptr()
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

std::vector<ScalarGroup*> PhysicsSteadyBase::get_scalar_group_ptr_vec()
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

std::vector<VariableGroup*> PhysicsSteadyBase::get_variable_group_ptr_vec()
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
