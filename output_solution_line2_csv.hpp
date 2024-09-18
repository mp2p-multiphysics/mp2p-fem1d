#ifndef OUTPUT_SOLUTION_LINE2_CSV
#define OUTPUT_SOLUTION_LINE2_CSV
#include <vector>
#include <fstream>
#include "Eigen/Eigen"
#include "grid_line2.hpp"
#include "scalar_line2.hpp"

void output_solution_line2_csv(ScalarLine2Class &sl2c, std::string file_out_base_str)
{

        // output file name
        std::string file_out_str = file_out_base_str + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridLine2Struct gl2s = sl2c.gl2s;

        // write to file
        file_out_stream << "id,pos_x,value\n";
        for (int n = 0; n < gl2s.num_point; n++)
        {
            file_out_stream << gl2s.point_id_vec[n] << "," << gl2s.point_pos_x_vec[n] << "," << sl2c.scalar_vec[n] << "\n";
        }

}

void output_solution_line2_csv(ScalarLine2Class &sl2c, std::string file_out_base_str, int ts)
{

        // output file name
        std::string file_out_str = file_out_base_str + std::to_string(ts) + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridLine2Struct gl2s = sl2c.gl2s;

        // write to file
        file_out_stream << "id,pos_x,value\n";
        for (int n = 0; n < gl2s.num_point; n++)
        {
            file_out_stream << gl2s.point_id_vec[n] << "," << gl2s.point_pos_x_vec[n] << "," << sl2c.scalar_vec[n] << "\n";
        }

}

#endif
