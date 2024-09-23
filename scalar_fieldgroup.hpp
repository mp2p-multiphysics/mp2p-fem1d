#ifndef SCALAR_FIELDGROUP
#define SCALAR_FIELDGROUP
#include <unordered_map>
#include <vector>
#include "mesh_line2.hpp"
#include "scalar_line2.hpp"

class ScalarFieldGroup
{

    public:

    // variables
    std::vector<ScalarLine2*> scalar_ptr_vec;
    std::unordered_map<MeshLine2Struct*, ScalarLine2*> scalar_ptr_map;

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
        // map mesh to each scalar
        for (auto scalar_ptr : scalar_ptr_vec)
        {
            scalar_ptr_map[scalar_ptr->mesh_l2_ptr] = scalar_ptr;
        }

    }

};

#endif
