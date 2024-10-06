#ifndef INTEGRAL_PHYSICSGROUP
#define INTEGRAL_PHYSICSGROUP
#include <vector>
#include "integral_line2.hpp"

class IntegralPhysicsGroup
{

    public:

    // variables
    std::vector<IntegralLine2*> integral_ptr_vec;

    // functions
    void evaluate_Ni_derivative();
    void evaluate_integral_Ni_line2();
    void evaluate_integral_derivative_Ni_line2_x();
    void evaluate_integral_Ni_line2_Nj_line2();
    void evaluate_integral_Ni_line2_derivative_Nj_line2_x();
    void evaluate_integral_div_Ni_line2_dot_div_Nj_line2();
    void evaluate_integral_Ni_line2_Nj_line2_derivative_Nk_line2_x();

    // default constructor
    IntegralPhysicsGroup()
    {

    }

    // constructor
    IntegralPhysicsGroup(std::vector<IntegralLine2*> integral_ptr_vec_in)
    {
        integral_ptr_vec = integral_ptr_vec_in;
    }

};

void IntegralPhysicsGroup::evaluate_Ni_derivative()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_Ni_derivative();
    }

}

void IntegralPhysicsGroup::evaluate_integral_Ni_line2()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_line2();
    }

}

void IntegralPhysicsGroup::evaluate_integral_derivative_Ni_line2_x()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_derivative_Ni_line2_x();
    }

}

void IntegralPhysicsGroup::evaluate_integral_Ni_line2_Nj_line2()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_line2_Nj_line2();
    }

}

void IntegralPhysicsGroup::evaluate_integral_Ni_line2_derivative_Nj_line2_x()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_line2_derivative_Nj_line2_x();
    }

}

void IntegralPhysicsGroup::evaluate_integral_div_Ni_line2_dot_div_Nj_line2()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_div_Ni_line2_dot_div_Nj_line2();
    }

}

void IntegralPhysicsGroup::evaluate_integral_Ni_line2_Nj_line2_derivative_Nk_line2_x()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_line2_Nj_line2_derivative_Nk_line2_x();
    }

}

#endif
