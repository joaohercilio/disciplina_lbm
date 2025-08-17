#pragma once

#include <fstream>
#include <iomanip>
#include <filesystem>
#include <vector>
#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"

void writeTSV(const std::vector<double>& f, const LatticeModel& lattice, const Geometry& geo, const std::string& path, int t);
