#pragma once

#include <vector>

#include "collision/BGK.hpp"
#include "collision/MRT.hpp"
#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"
/**
 * @brief Initializes distribution function based on an initial velocity field and density by 
 * computing the equilibrium distribution
 * 
 * @param f Vector storing the particle distribution function 
 * @param lattice Constant reference to the lattice object 
 * @param geometry Constant reference to the geometry object
 * @param colParams Constant reference to the collision parameters map specific to the Collision Model
 * @param u Vector storing the domain initial velocity
 * @param rho Vector storing the domain initial density
 */
void initializeFields(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry,
                            const ColParamMap& colParams,
                            const std::vector<double>& u,
                            const std::vector<double>& rho);
