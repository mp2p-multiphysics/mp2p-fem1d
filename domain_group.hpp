#ifndef DOMAIN_GROUP
#define DOMAIN_GROUP
#include <vector>
#include "domain_line2.hpp"

namespace FEM1D
{

class DomainGroup
{
    /*

    Groups domains that are used in the same physics.

    Variables
    =========
    domain_l2_ptr_vec_in : vector<DomainLine2*>
        vector with pointers to DomainLine2 objects.
    
    */

    public:

    // vector with domains in group
    std::vector<DomainLine2*> domain_l2_ptr_vec;

    // default constructor
    DomainGroup() {}

    // constructor
    DomainGroup(std::vector<DomainLine2*> domain_l2_ptr_vec_in)
    {
        
        // store variables
        domain_l2_ptr_vec = domain_l2_ptr_vec_in;

    }

};

}

#endif
