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

#include "boundary_line2.hpp"
#include "boundary_field.hpp"
#include "container_boundaryconfig.hpp"
#include "container_typedef.hpp"
#include "integral_line2.hpp"
#include "integral_field.hpp"
#include "matrixequation_steady.hpp"
#include "matrixequation_transient.hpp"
#include "mesh_line2.hpp"
#include "mesh_field.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_convectiondiffusion.hpp"
#include "physicssteady_diffusion.hpp"
#include "physicstransient_base.hpp"
#include "physicstransient_convectiondiffusion.hpp"
#include "physicstransient_diffusion.hpp"
#include "scalar_field.hpp"
#include "scalar_line2.hpp"
#include "variable_field.hpp"
#include "variable_line2.hpp"
