#ifndef INTEGRAL_GROUP
#define INTEGRAL_GROUP
#include <vector>
#include "integral_line2.hpp"

class IntegralGroup
{
    /*

    Groups test function integrals (N) that are used in the same physics.

    Variables
    =========
    integral_l2_ptr_vec_in : vector<IntegralLine2*>
        vector with pointers to BoundaryLine2 objects.

    Functions
    =========
    evaluate_Ni_derivative : void
        Calculates test functions (N) and their derivatives.
        Must be called before integrals are evaluated.
    evaluate_integral_Ni : void
        Calculates the integral of Ni.
    evaluate_integral_derivative_Ni_x : void
        Calculates the integral of d(Ni)/dx.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj.
    evaluate_integral_Ni_derivative_Nj_x : void
        Calculates the integral of Ni * d(Nj)/dx.
    evaluate_integral_div_Ni_dot_div_Nj : void
        Calculates the integral of div(Ni) dot div(Nj).
    evaluate_integral_Ni_Nj_derivative_Nk_x : void
        Calculates the integral of Ni * Nj * d(Nk)/dx.

    */


    public:

    // variables
    std::vector<IntegralLine2*> integral_l2_ptr_vec;

    // functions for computing domain integrals
    void evaluate_Ni();
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();

    // functions for computing boundary integrals
    void evaluate_boundary_Ni();
    void evaluate_integral_boundary_Ni();
    void evaluate_integral_boundary_Ni_Nj();

    // default constructor
    IntegralGroup() {}

    // constructor
    IntegralGroup(std::vector<IntegralLine2*> integral_l2_ptr_vec_in)
    {
        integral_l2_ptr_vec = integral_l2_ptr_vec_in;
    }

};

void IntegralGroup::evaluate_Ni()
{
    /*

    Calculates test functions (N) and their derivatives.
    Must be called before integrals are evaluated.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_Ni();
    }

}

void IntegralGroup::evaluate_integral_Ni()
{
    /*

    Calculates the integral of Ni.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni();
    }

}

void IntegralGroup::evaluate_integral_derivative_Ni_x()
{
    /*

    Calculates the integral of d(Ni)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_derivative_Ni_x();
    }

}

void IntegralGroup::evaluate_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj();
    }

}

void IntegralGroup::evaluate_integral_Ni_derivative_Nj_x()
{
    /*

    Calculates the integral of Ni * d(Nj)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_derivative_Nj_x();
    }

}

void IntegralGroup::evaluate_integral_div_Ni_dot_div_Nj()
{
    /*

    Calculates the integral of div(Ni) dot div(Nj).

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_div_Ni_dot_div_Nj();
    }

}

void IntegralGroup::evaluate_integral_Ni_Nj_derivative_Nk_x()
{
    /*

    Calculates the integral of Ni * Nj * d(Nk)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_x();
    }

}

void IntegralGroup::evaluate_boundary_Ni()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_boundary_Ni();
    }

}

void IntegralGroup::evaluate_integral_boundary_Ni()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_boundary_Ni();
    }

}

void IntegralGroup::evaluate_integral_boundary_Ni_Nj()
{

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_l2_ptr_vec)
    {
        integral_ptr->evaluate_integral_boundary_Ni_Nj();
    }

}

#endif
