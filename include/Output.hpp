#pragma once

#include <fstream>
#include <iomanip>
#include <filesystem>
#include <vector>

#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"
#include "vtkio/xmlvtk.h"

/**
 * @brief Writes the distribution function to a TSV (tab-separated values) file.
 *
 * @param f Constant reference to the vector storing the particle distribution function
 * @param lattice Constant reference to the lattice model
 * @param geo Constant reference to the geometry of the domain
 * @param path Output folder or file path
 * @param t Current timestep (used to name the file)
 */
void writeTSV(const std::vector<double>& f, const LatticeModel& lattice, const Geometry& geo, const std::string& path, int t);

/**
 * @brief Writes the distribution function to a VTI (VTK ImageData) file. 
 * 
 * This format can be read by Paraview for 3D visualization
 *
 * @param f Constant reference to the vector storing the particle distribution function
 * @param lattice Constant reference to the lattice model
 * @param geo Constant reference to the geometry of the domain
 * @param path Output folder or file path
 * @param t Current timestep (used to name the file)
 */

void writeVTI(const std::vector<double>& f, const LatticeModel& lattice, const Geometry& geo, const std::string& path, int t);
