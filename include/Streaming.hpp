#pragma once

#include <vector>
#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"

class Streaming {

public:
    ~Streaming() = default;

    void performStreaming(std::vector<double>& f, std::vector<double>& fn, const LatticeModel& lattice, const Geometry& geometry) const;

};