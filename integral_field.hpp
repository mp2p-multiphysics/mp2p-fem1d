#ifndef INTEGRAL_PHYSICSGROUP
#define INTEGRAL_PHYSICSGROUP
#include <vector>
#include "integral_line2.hpp"

class IntegralField
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

    // functions
    void evaluate_Ni_derivative();
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();

    // default constructor
    IntegralField()
    {

    }

    // constructor
    IntegralField(std::vector<IntegralLine2*> integral_l2_ptr_vec_in)
    {
        integral_l2_ptr_vec = integral_l2_ptr_vec_in;
    }

};

void IntegralField::evaluate_Ni_derivative()
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
        integral_ptr->evaluate_Ni_derivative();
    }

}

void IntegralField::evaluate_integral_Ni()
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

void IntegralField::evaluate_integral_derivative_Ni_x()
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

void IntegralField::evaluate_integral_Ni_Nj()
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

void IntegralField::evaluate_integral_Ni_derivative_Nj_x()
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

void IntegralField::evaluate_integral_div_Ni_dot_div_Nj()
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

void IntegralField::evaluate_integral_Ni_Nj_derivative_Nk_x()
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

#endif
