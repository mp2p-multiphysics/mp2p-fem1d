#ifndef VARIABLE_LINE2
#define VARIABLE_LINE2
#include <fstream>
#include <sstream>
#include "domain_line2.hpp"
#include "container_typedef.hpp"

namespace FEM1D
{

class VariableLine2
{
    /*

    Variable applied over line2 domain elements.

    Variables
    =========
    domain_in : DomainLine2
        Domain where variable is applied.
    value_initial_in : double
        Initial value of the variable.

    Functions
    =========
    set_output : void
        Set the output CSV file name with the values of the variable.
    output_csv : void
        Outputs a CSV file with the values of the variable.

    */

    public:

    // domain where variable is applied
    DomainLine2* domain_ptr;  

    // values in variable
    int num_point = 0;  // number of points in domain
    VectorDouble point_value_vec;  // key: domain ID; value: value
    
    // use for generating csv file
    bool is_file_out = false;
    std::string file_out_base_str;
    std::vector<std::string> file_out_base_vec;

    // functions
    void set_output(std::string file_out_str);
    void output_csv();
    void output_csv(int ts);

    // default constructor
    VariableLine2() {}

    // constructor
    VariableLine2(DomainLine2 &domain_in, double value_initial_in)
    {

        // store domain
        domain_ptr = &domain_in;

        // get number of domain points
        num_point = domain_ptr->num_point;

        // populate value vector with initial values
        for (int pdid = 0; pdid < num_point; pdid++)
        {
            point_value_vec.push_back(value_initial_in);
        }

    }

    private:
    
};

void VariableLine2::set_output(std::string file_out_str)
{
    /*

    Set the output CSV file name with the values of the variable.

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

    // generate CSV file when output_csv is called
    is_file_out = true;

    // split filename at '*'
    // will be replaced with timestep later
    std::stringstream file_out_base_stream(file_out_base_str);
    std::string string_sub;
    while(std::getline(file_out_base_stream, string_sub, '*'))
    {
        file_out_base_vec.push_back(string_sub);
    }

}

void VariableLine2::output_csv()
{
    /*

    Outputs a CSV file with the values of the variable.

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
    if (!is_file_out)
    {
        return;
    }

    // initialize file stream
    std::ofstream file_out_stream(file_out_base_str);

    // write to file
    file_out_stream << "gid,position_x,value\n";
    for (int pdid = 0; pdid < num_point; pdid++)
    {
        file_out_stream << domain_ptr->point_pdid_to_pgid_vec[pdid] << ",";
        file_out_stream << domain_ptr->point_position_x_vec[pdid] << ",";
        file_out_stream << point_value_vec[pdid] << "\n";
    }

}

void VariableLine2::output_csv(int ts)
{
    /*

    Outputs a CSV file with the values of the variable.

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
    if (!is_file_out)
    {
        return;
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
    for (int pdid = 0; pdid < num_point; pdid++)
    {
        file_out_stream << domain_ptr->point_pdid_to_pgid_vec[pdid] << ",";
        file_out_stream << domain_ptr->point_position_x_vec[pdid] << ",";
        file_out_stream << point_value_vec[pdid] << "\n";
    }

}

}

#endif
