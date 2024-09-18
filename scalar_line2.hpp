#ifndef SCALAR_LINE2
#define SCALAR_LINE2
#include <vector>
#include "Eigen/Eigen"
#include "grid_line2.hpp"

class ScalarLine2Class
{

    public:
    
    // variables
    GridLine2Struct gl2s;
    std::vector<double> scalar_vec;

    // functions
    void set_from_constant(double value);
    void set_from_vector(std::vector<double> value_vec, int start_id);
    void set_from_evector(Eigen::VectorXd value_evec, int start_id);

    // constructor
    ScalarLine2Class(GridLine2Struct &gl2s_in)
    {
        gl2s = gl2s_in;
    }

    ScalarLine2Class(GridLine2Struct &gl2s_in, double value)
    {
        gl2s = gl2s_in;
        set_from_constant(value);
    }

    ScalarLine2Class(GridLine2Struct &gl2s_in, std::vector<double> value_vec, int start_id)
    {
        gl2s = gl2s_in;
        set_from_vector(value_vec, start_id);
    }

    ScalarLine2Class(GridLine2Struct &gl2s_in, Eigen::VectorXd value_evec, int start_id)
    {
        gl2s = gl2s_in;
        set_from_evector(value_evec, start_id);
    }

};

void ScalarLine2Class::set_from_constant(double value)
{

    // set each scalar_vec value to the constant
    for (int n = 0; n < gl2s.num_point; n++)
    {
        scalar_vec.push_back(value);
    }

}

void ScalarLine2Class::set_from_vector(std::vector<double> value_vec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gl2s.num_point; n++)
    {
        scalar_vec.push_back(value_vec[start_id + n]);
    }

}

void ScalarLine2Class::set_from_evector(Eigen::VectorXd value_evec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gl2s.num_point; n++)
    {
        scalar_vec.push_back(value_evec.coeffRef(start_id + n));
    }

}

#endif
