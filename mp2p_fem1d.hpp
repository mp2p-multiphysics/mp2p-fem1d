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

#include "domain_0d.hpp"
#include "domain_1d.hpp"
#include "integral_1d.hpp"
#include "matrixequation_steady.hpp"
#include "matrixequation_transient.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_convectiondiffusion.hpp"
#include "physicssteady_diffusion.hpp"
#include "physicstransient_base.hpp"
#include "physicstransient_diffusion.hpp"
#include "physicstransient_convectiondiffusion.hpp"
#include "scalar_0d.hpp"
#include "scalar_1d.hpp"
#include "variable_1d.hpp"
#include "variable_group.hpp"
