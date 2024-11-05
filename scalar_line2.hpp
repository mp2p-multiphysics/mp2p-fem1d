#ifndef SCALAR_LINE2
#define SCALAR_LINE2
#include <functional>
#include "mesh_line2.hpp"
#include "variable_line2.hpp"

class ScalarLine2
{
    /*

    Scalar applied over line2 mesh elements.

    Variables
    =========
    mesh_in : MeshLine2
        Mesh where scalar value is applied.
    u_init_in : double
        Initial value of the scalar.

    Functions
    =========
    output_csv : void
        Outputs a CSV file with the values of the scalar.

    */

    public:

    // mesh where variable is applied
    MeshLine2* mesh_ptr;
    int num_point_domain = 0;  // number of points in domain

    // values in variable
    VectorDouble point_value_vec;  // key: domain ID; value: value
    bool is_value_constant = true;
    double value_constant = 0;  // used if value is constant
    std::function<double(double, VectorDouble)> value_function;  // used if value is non-constant
    std::vector<VariableLine2*> variable_ptr_vec;  // variables that values depend on

    // functions
    void output_csv(std::string file_out_str);
    void output_csv(std::string file_out_base_str, int ts);
    void update_value();

    // default constructor
    ScalarLine2()
    {

    }

    // constructor for constant values
    ScalarLine2(MeshLine2 &mesh_in, double value_constant_in)
    {

        // store mesh
        mesh_ptr = &mesh_in;
        num_point_domain = mesh_ptr->num_point_domain;

        // store values
        is_value_constant = true;
        value_constant = value_constant_in;

        // populate initial values
        point_value_vec = VectorDouble(num_point_domain, value_constant);

    }

    // constructor for non-constant values
    ScalarLine2(MeshLine2 &mesh_in, std::function<double(double, VectorDouble)> value_function_in, std::vector<VariableLine2*> variable_ptr_vec_in)
    {

        // store mesh
        mesh_ptr = &mesh_in;
        num_point_domain = mesh_ptr->num_point_domain;

        // store values
        is_value_constant = false;
        value_function = value_function_in;
        variable_ptr_vec = variable_ptr_vec_in;

        // populate initial values
        point_value_vec = VectorDouble(num_point_domain, 0.);  // unused placeholder

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
    for (int point_did = 0; point_did < num_point_domain; point_did++)
    {
        file_out_stream << mesh_ptr->point_gid_vec[point_did] << ",";
        file_out_stream << mesh_ptr->point_position_x_vec[point_did] << ",";
        file_out_stream << point_value_vec[point_did] << "\n";
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
    for (int point_did = 0; point_did < num_point_domain; point_did++)
    {
        file_out_stream << mesh_ptr->point_gid_vec[point_did] << ",";
        file_out_stream << mesh_ptr->point_position_x_vec[point_did] << ",";
        file_out_stream << point_value_vec[point_did] << "\n";
    }

}

void ScalarLine2::update_value()
{

    // skip if constant value
    if (is_value_constant)
    {
        return;
    }

    // iterate through each point in scalar
    for (int point_did = 0; point_did < mesh_ptr->num_point_domain; point_did++)
    {

        // get mesh coordinate
        double position_x = mesh_ptr->point_position_x_vec[point_did];

        // iterate through each variable that scalar depends on
        VectorDouble value_vec;
        for (auto variable_ptr : variable_ptr_vec)
        {
            value_vec.push_back(variable_ptr->point_value_vec[point_did]);
        }

        // calculate scalar value
        point_value_vec[point_did] = value_function(position_x, value_vec);

    }

}

#endif
