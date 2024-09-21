#ifndef SCALAR_FIELDGROUP
#define SCALAR_FIELDGROUP
#include <unordered_map>
#include <vector>
#include "scalar_line2.hpp"

class ScalarFieldGroup
{

    public:

    // variables
    std::vector<ScalarLine2*> scalar_ptr_vec;
    std::unordered_map<int, double> point_u_map;  // key: global point id; value: value of scalar

    // functions

    // default constructor
    ScalarFieldGroup()
    {

    }

    // constructor
    ScalarFieldGroup(std::vector<ScalarLine2*> scalar_ptr_vec_in)
    {
        
        // store variables
        scalar_ptr_vec = scalar_ptr_vec_in;

        // iterate through each scalar
        for (auto scalar_ptr : scalar_ptr_vec)
        {
            for (int indx_i = 0; indx_i < scalar_ptr->mesh_l2_ptr->num_domain_point; indx_i++)
            {

                // get global point ID and value of scalar
                int n = scalar_ptr->mesh_l2_ptr->point_global_id_vec[indx_i];
                int u = scalar_ptr->point_u_vec[indx_i];

                // assign to map
                point_u_map[n] = u;

            }
        }

    }

};

#endif
