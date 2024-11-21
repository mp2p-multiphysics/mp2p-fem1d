/*
###################################################
###################################################
###   __  __ ____ ____  ____                    ###
###  |  \/  |  _ \___ \|  _ \     Multipurpose  ###
###  | |\/| | |_) |__) | |_) |    Multiphysics  ###
###  | |  | |  __// __/|  __/     Package       ###
###  |_|  |_|_|  |_____|_|        (FEM 1D)      ###
###                                             ###
###################################################
###################################################
*/

#include "boundary_group.hpp"
#include "boundary_line2.hpp"
#include "container_typedef.hpp"
#include "domain_group.hpp"
#include "domain_line2.hpp"
#include "integral_group.hpp"
#include "integral_line2.hpp"
#include "matrixequation_steady.hpp"
#include "matrixequation_transient.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_diffusion.hpp"
#include "physicstransient_base.hpp"
#include "physicstransient_diffusion.hpp"
#include "scalar_group.hpp"
#include "scalar_line2.hpp"
#include "variable_group.hpp"
#include "variable_line2.hpp"

// #include "boundary_group.hpp"
// #include "boundary_line2.hpp"
// #include "domain_group.hpp"
// #include "domain_line2.hpp"
// #include "container_typedef.hpp"
// #include "integral_group.hpp"
// #include "integral_line2.hpp"
// #include "matrixequation_steady.hpp"
// #include "matrixequation_transient.hpp"
// #include "physicssteady_base.hpp"
// #include "physicssteady_convectiondiffusion.hpp"
// #include "physicssteady_diffusion.hpp"
// #include "physicstransient_base.hpp"
// #include "physicstransient_convectiondiffusion.hpp"
// #include "physicstransient_diffusion.hpp"
// #include "scalar_group.hpp"
// #include "scalar_line2.hpp"
// #include "variable_group.hpp"
// #include "variable_line2.hpp"
