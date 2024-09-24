#ifndef SCALAR_LINE2
#define SCALAR_LINE2
#include "mesh_line2.hpp"

class ScalarLine2
{

    public:

    // values in variable
    int num_point_domain = 0;  // number of points in domain
    VectorDouble point_value_vec;  // key: domain ID; value: value
    
    // mesh where variable is applied
    MeshLine2Struct* mesh_l2_ptr;  

    // functions
    void output_csv(std::string file_out_str);

    // default constructor
    ScalarLine2()
    {

    }

    // constructor
    ScalarLine2(MeshLine2Struct &mesh_l2_in, double u_init_in)
    {

        // store mesh
        mesh_l2_ptr = &mesh_l2_in;

        // get number of domain points
        num_point_domain = mesh_l2_ptr->num_point_domain;

        // populate initial values
        for (int point_did = 0; point_did < num_point_domain; point_did++)
        {
            point_value_vec.push_back(u_init_in);
        }

    }

};

void ScalarLine2::output_csv(std::string file_out_str)
{

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "id,pos_x,value\n";
    for (int point_did = 0; point_did < num_point_domain; point_did++)
    {
        file_out_stream << mesh_l2_ptr->point_gid_vec[point_did] << ",";
        file_out_stream << mesh_l2_ptr->point_position_x_vec[point_did] << ",";
        file_out_stream << point_value_vec[point_did] << "\n";
    }

}

#endif
