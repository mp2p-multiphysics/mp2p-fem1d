#ifndef PHYSICSSTEADY_BASE
#define PHYSICSSTEADY_BASE
#include <vector>
#include "Eigen/Eigen"
#include "boundary_physicsgroup.hpp"
#include "container_typedef.hpp"
#include "mesh_physicsgroup.hpp"
#include "scalar_fieldgroup.hpp"
#include "integral_physicsgroup.hpp"
#include "variable_fieldgroup.hpp"

class PhysicsSteadyBase
{

    public:

    // physics groups
    MeshPhysicsGroup *mesh_physics_ptr;
    BoundaryPhysicsGroup *boundary_physics_ptr;
    IntegralPhysicsGroup *integral_physics_ptr;

    // vector of variable fields
    std::vector<VariableFieldGroup*> variable_field_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    virtual void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);
    virtual void set_start_row(int start_row_in);
    virtual int get_start_row();
    virtual std::vector<VariableFieldGroup*> get_variable_field_ptr_vec();

    // default constructor
    PhysicsSteadyBase()
    {

    }

};

void PhysicsSteadyBase::matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec)
{

}

void PhysicsSteadyBase::set_start_row(int start_row_in)
{
    start_row = start_row_in;
}

int PhysicsSteadyBase::get_start_row()
{
    return start_row;
}

std::vector<VariableFieldGroup*> PhysicsSteadyBase::get_variable_field_ptr_vec()
{
    return variable_field_ptr_vec;
}

#endif
