#ifndef VARIABLE_LINE2
#define VARIABLE_LINE2
#include <fstream>
#include "container_mesh.hpp"
#include "container_typedef.hpp"

class VariableLine2
{

    public:

    // mesh and data points
    MeshLine2Struct* mesh_l2_ptr;
    VectorDouble point_u_vec;

    // starting column in matrix equation
    int start_col = -1;

    // function
    void output_csv(std::string file_out_str);

    // default constructor
    VariableLine2()
    {

    }

    // constructor
    VariableLine2(MeshLine2Struct &mesh_l2_in, double u_init_in)
    {
        
        // store mesh
        mesh_l2_ptr = &mesh_l2_in;

        // populate initial values
        for (int indx_i = 0; indx_i < mesh_l2_ptr->num_point; indx_i++)
        {
            point_u_vec.push_back(u_init_in);
        }

    }

};

void VariableLine2::output_csv(std::string file_out_str)
{

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "id,pos_x,value\n";
    for (int n = 0; n < mesh_l2_ptr->num_point; n++)
    {
        file_out_stream << mesh_l2_ptr->point_id_vec[n] << ",";
        file_out_stream << mesh_l2_ptr->point_pos_x_vec[n] << ",";
        file_out_stream << point_u_vec[n] << "\n";
    }

}

#endif
