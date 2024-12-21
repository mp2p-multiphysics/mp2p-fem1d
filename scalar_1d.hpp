#ifndef SCALAR_1D
#define SCALAR_1D
#include <fstream>
#include <sstream>
#include <functional>
#include "container_typedef.hpp"
#include "domain_1d.hpp"
#include "variable_1d.hpp"

namespace FEM1D
{

class Scalar1D
{
    /*

    Scalar applied over a domain

    Variables (for constant values)
    =========
    domain_in : Domain1D
        Domain where scalar value is applied.
    value_constant_in : double
        Value of the scalar.

    Variables (for non-constant values)
    =========
    domain_in : Domain1D
        Domain where scalar value is applied.
    value_function_in : function(double, double, VectorDouble) -> double
        Function used to compute scalar values based on variable values.
    variable_ptr_vec_in : vector<Variable1D*>
        vector of pointers to variable objects needed to compute scalar values.

    Functions
    =========
    set_output : void
        Set the output CSV file name with the values of the scalar.
    output_csv : void
        Outputs a CSV file with the values of the scalar.
    update_value : void
        Recalculates non-constant values.
    get_neighbor_value : VectorDouble
        Outputs a vector with the scalar values at points surrounding an element.

    Notes
    ====
    The inputs to the value function are the x-coordinate (double) and vector of variable values (VectorDouble) at a specified point.
        The variable values are in the same order as the variables in variable_ptr_vec.
    The output of the value function is the value of the scalar.

    */

    public:

    // domain where variable is applied
    Domain1D* domain_ptr;

    // values in scalar
    VectorDouble point_value_vec;  // key: point domain ID; value: value

    // use for non-constant scalars
    bool is_value_constant = true;
    double value_constant = 0;  // used if value is constant
    std::function<double(double, VectorDouble)> value_function;  // used if value is non-constant
    std::vector<Variable1D*> variable_ptr_vec;  // variables that values depend on

    // use for generating csv file
    std::string file_out_base_str;

    // functions
    void set_output(std::string file_out_str);
    void output_csv();
    void output_csv(int ts);
    void update_value();
    VectorDouble get_neighbor_value(int edid);

    // default constructor
    Scalar1D() {}
    
    // constructor for constant values
    Scalar1D(Domain1D &domain_in, double value_constant_in)
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
    Scalar1D(Domain1D &domain_in, std::function<double(double, VectorDouble)> value_function_in, std::vector<Variable1D*> variable_ptr_vec_in)
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

void Scalar1D::set_output(std::string file_out_str)
{
    /*

    Set the output CSV file name with the values of the scalar.

    Arguments
    =========
    file_out_str : string
        Path to CSV file.

    Returns
    =======
    (none)

    Notes
    =====
    file_out_str must have an asterisk '*' for transient simulations.
    This will be replaced with the timestep number.

    */

    // set file name
    file_out_base_str = file_out_str;

}

void Scalar1D::output_csv()
{
    /*

    Outputs a CSV file with the values of the scalar.

    Arguments
    =========
    (none)

    Returns
    =======
    (none)

    Notes
    =====
    This function is intended to be used with steady-state simulations.

    */

    // do not make file if filename not set
    if (file_out_base_str.empty())
    {
        return;
    }

    // initialize file stream
    std::ofstream file_out_stream(file_out_base_str);

    // write to file
    file_out_stream << "point_id,position_x,value\n";
    for (int pdid = 0; pdid < domain_ptr->num_point; pdid++)
    {
        file_out_stream << domain_ptr->point_pdid_to_pgid_vec[pdid] << ",";
        file_out_stream << domain_ptr->point_position_x_vec[pdid] << ",";
        file_out_stream << point_value_vec[pdid] << "\n";
    }

}

void Scalar1D::output_csv(int ts)
{
    /*

    Outputs a CSV file with the values of the scalar.

    Arguments
    =========
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

    // do not make file if filename not set
    if (file_out_base_str.empty())
    {
        return;
    }

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
    file_out_stream << "point_id,position_x,value\n";
    for (int pdid = 0; pdid < domain_ptr->num_point; pdid++)
    {
        file_out_stream << domain_ptr->point_pdid_to_pgid_vec[pdid] << ",";
        file_out_stream << domain_ptr->point_position_x_vec[pdid] << ",";
        file_out_stream << point_value_vec[pdid] << "\n";
    }

}

void Scalar1D::update_value()
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

VectorDouble Scalar1D::get_neighbor_value(int edid)
{
    /*
    
    Outputs a vector with the scalar values at points surrounding an element.

    Arguments
    =========
    edid : int
        Domain ID of the element.

    Returns
    =======
    value_vec : VectorInt
        vector with scalar values at points surrounding an element.

    */

    // get point surrounding element
    VectorInt pgid_vec = domain_ptr->element_edid_plid_to_pgid_vec[edid];

    // initialize vector with values
    VectorDouble value_vec;

    // iterate through each point and get value
    for (int pgid : pgid_vec)
    {
        int pdid = domain_ptr->point_pgid_to_pdid_map[pgid];
        double value_sub = point_value_vec[pdid];
        value_vec.push_back(value_sub);
    }

    return value_vec;

}

}

#endif
