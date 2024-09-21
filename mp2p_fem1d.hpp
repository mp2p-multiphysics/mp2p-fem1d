/*
####################################################
####################################################
###   __  __ ____ ____  ____                     ###
###  |  \/  |  _ \___ \|  _ \     Multi-purpose  ###
###  | |\/| | |_) |__) | |_) |    Multiphysics   ###
###  | |  | |  __// __/|  __/     Program        ###
###  |_|  |_|_|  |_____|_|        (FEM 1D)       ###
###                                              ###
####################################################
####################################################
*/

#include "boundary_line2.hpp"
#include "boundary_physicsgroup.hpp"
#include "initialize_boundary_line2_csv.hpp"
#include "initialize_mesh_line2_csv.hpp"
#include "integral_line2.hpp"
#include "integral_physicsgroup.hpp"
#include "matrixequation_steady.hpp"
#include "mesh_line2.hpp"
#include "mesh_physicsgroup.hpp"
#include "physics_solidheattransfer_steady.hpp"
#include "scalar_fieldgroup.hpp"
#include "scalar_line2.hpp"
#include "variable_fieldgroup.hpp"
#include "variable_line2.hpp"

// #include "boundary_group.hpp"
// #include "boundary_line2.hpp"
// #include "container_typedef.hpp"
// #include "initialize_mesh_line2_csv.hpp"
// #include "initialize_boundary_line2_csv.hpp"
// #include "integral_group.hpp"
// #include "integral_line2.hpp"
// #include "matrixequation_steady.hpp"
// #include "mesh_group.hpp"
// #include "mesh_line2.hpp"
// #include "physics_solidheattransfer_steady_line2.hpp"
// #include "scalar_group.hpp"
// #include "scalar_line2.hpp"
// #include "variable_group.hpp"
// #include "variable_line2.hpp"
