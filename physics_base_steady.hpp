#ifndef PHYSICS_BASE_STEADY
#define PHYSICS_BASE_STEADY
#include <vector>
#include "Eigen/Eigen"
#include "boundary_physicsgroup.hpp"
#include "container_typedef.hpp"
#include "mesh_physicsgroup.hpp"
#include "scalar_fieldgroup.hpp"
#include "integral_physicsgroup.hpp"
#include "variable_fieldgroup.hpp"

class PhysicsBaseSteady
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
    virtual std::vector<VariableFieldGroup*> get_variable_field_ptr_vec();

    // default constructor
    PhysicsBaseSteady()
    {

    }

};

void PhysicsBaseSteady::matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec)
{

}

void PhysicsBaseSteady::set_start_row(int start_row_in)
{
    start_row = start_row_in;
}

std::vector<VariableFieldGroup*> PhysicsBaseSteady::get_variable_field_ptr_vec()
{
    return variable_field_ptr_vec;
}

#endif
