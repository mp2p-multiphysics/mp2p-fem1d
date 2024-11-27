#ifndef SCALAR_LINE2
#define SCALAR_LINE2
#include <fstream>
#include <sstream>
#include <functional>
#include "container_typedef.hpp"
#include "domain_line2.hpp"
#include "variable_line2.hpp"

class ScalarLine2
{
    /*

    Scalar applied over line2 domain elements.

    Variables (for constant values)
    =========
    domain_in : DomainLine2
        Domain where scalar value is applied.
    value_constant_in : double
        Value of the scalar.

    Variables (for non-constant values)
    =========
    domain_in : DomainLine2
        Domain where scalar value is applied.
    value_function_in : function(double, VectorDouble) -> double
        Function used to compute scalar values based on variable values.
    variable_ptr_vec_in : vector<VariableLine2*>
        vector of pointers to variable objects needed to compute scalar values.

    Functions
    =========
    output_csv : void
        Outputs a CSV file with the values of the scalar.
    update_value : void
        Recalculates non-constant values.

    Notes
    ====
    The inputs to the value function are the x-coordinate (double) and vector of variable values (VectorDouble) at a specified point.
        The variable values are in the same order as the variables in variable_ptr_vec.
    The output of the value function is the value of the scalar.

    */

    public:

    // domain where variable is applied
    DomainLine2* domain_ptr;

    // values in scalar
    VectorDouble point_value_vec;  // key: point domain ID; value: value
    Vector2D element_value_vec;  // key: element domain ID, point local ID; value: value

    // use for non-constant scalars
    bool is_value_constant = true;
    double value_constant = 0;  // used if value is constant
    std::function<double(double, VectorDouble)> value_function;  // used if value is non-constant
    std::vector<VariableLine2*> variable_ptr_vec;  // variables that values depend on

    // functions
    void output_csv(std::string file_out_str);
    void output_csv(std::string file_out_base_str, int ts);
    void update_value();

    // default constructor
    ScalarLine2() {}
    
    // constructor for constant values
    ScalarLine2(DomainLine2 &domain_in, double value_constant_in)
    {

        // store domain
        domain_ptr = &domain_in;

        // store values
        is_value_constant = true;
        value_constant = value_constant_in;

        // populate initial values
        point_value_vec = VectorDouble(domain_ptr->num_point, value_constant);

    }

    // constructor for non-constant values
    ScalarLine2(DomainLine2 &domain_in, std::function<double(double, VectorDouble)> value_function_in, std::vector<VariableLine2*> variable_ptr_vec_in)
    {

        // store domain
        domain_ptr = &domain_in;

        // store values
        is_value_constant = false;
        value_function = value_function_in;
        variable_ptr_vec = variable_ptr_vec_in;

        // populate initial values
        point_value_vec = VectorDouble(domain_ptr->num_point, 0.);  // unused placeholder

    }

};

void ScalarLine2::output_csv(std::string file_out_str)
{
    /*

    Outputs a CSV file with the values of the scalar.

    Arguments
    =========
    file_out_str : string
        Path to CSV file.

    Returns
    =======
    (none)

    Notes
    =====
    This function is intended to be used with steady-state simulations.

    */

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "gid,position_x,value\n";
    for (int pdid = 0; pdid < domain_ptr->num_point; pdid++)
    {
        file_out_stream << domain_ptr->point_pdid_to_pgid_vec[pdid] << ",";
        file_out_stream << domain_ptr->point_position_x_vec[pdid] << ",";
        file_out_stream << point_value_vec[pdid] << "\n";
    }

}

void ScalarLine2::output_csv(std::string file_out_base_str, int ts)
{
    /*

    Outputs a CSV file with the values of the scalar.

    Arguments
    =========
    file_out_base_str : string
        Path to CSV file with base file name.
    ts : int
        Timestep number.

    Returns
    =======
    (none)

    Notes
    =====
    file_out_base_str must have an asterisk '*', which will be replaced with ts.
    This function is intended to be used with transient simulations.

    */

    // split filename at '*'
    // will be replaced with timestep later
    std::vector<std::string> file_out_base_vec;
    std::stringstream file_out_base_stream(file_out_base_str);
    std::string string_sub;
    while(std::getline(file_out_base_stream, string_sub, '*'))
    {
        file_out_base_vec.push_back(string_sub);
    }

    // create output filename
    // replace '*' with timestep
    std::string file_out_str = file_out_base_vec[0];
    for (int i = 1; i < file_out_base_vec.size(); i++)
    {
        file_out_str += std::to_string(ts) + file_out_base_vec[i];
    }

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "gid,position_x,value\n";
    for (int pdid = 0; pdid < domain_ptr->num_point; pdid++)
    {
        file_out_stream << domain_ptr->point_pdid_to_pgid_vec[pdid] << ",";
        file_out_stream << domain_ptr->point_position_x_vec[pdid] << ",";
        file_out_stream << point_value_vec[pdid] << "\n";
    }

}

void ScalarLine2::update_value()
{
    /*

    Recalculates non-constant values.

    Arguments
    =========
    (none)
    
    Returns
    =========
    (none)

    */

    // skip if constant value
    if (is_value_constant)
    {
        return;
    }

    // iterate through each point in scalar
    for (int pdid = 0; pdid < domain_ptr->num_point; pdid++)
    {

        // get domain coordinate
        double position_x = domain_ptr->point_position_x_vec[pdid];

        // iterate through each variable that scalar depends on
        VectorDouble value_vec;
        for (auto variable_ptr : variable_ptr_vec)
        {
            value_vec.push_back(variable_ptr->point_value_vec[pdid]);
        }

        // calculate scalar value
        point_value_vec[pdid] = value_function(position_x, value_vec);

    }

}

#endif
