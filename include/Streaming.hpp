#pragma once

#include <vector>
#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"

/**
 * @brief Performs the streaming step of the LBM.
 *
 * This function updates the distribution function by propagating the particle
 * populations along the lattice directions according to the lattice model.
 * The updated values are stored in the output vector `fn`.
 *
 * @param f Constant reference to the vector storing the current distribution function
 * @param fn Reference to the vector where the streamed distribution function will be stored
 * @param lattice Constant reference to the lattice model
 * @param geometry Constant reference to the geometry
 */
void performStreaming(std::vector<double>& f, std::vector<double>& fn, const LatticeModel& lattice, const Geometry& geometry);

