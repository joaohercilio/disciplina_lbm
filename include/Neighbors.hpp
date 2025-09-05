#pragma once

#include <vector>

/**
 * @file Neighbors.hpp
 * @brief Precomputed neighbor table for LBM streaming.
 *
 * For each lattice site (id) and each velocity direction (k),
 * stores the linear index of the neighbor reached by moving one step.
 */
class Geometry;
class LatticeModel;

class Neighbors {
public:

    Neighbors(const Geometry& geo, const LatticeModel& lat);

    // Get neighbor index for site id in direction k
    int operator()(int id, int k, int numOfVel) const {
        return neighbors_[id * numOfVel + k];
    }

    const std::vector<int>& raw() const { return neighbors_; }

private:
    std::vector<int> neighbors_; // Table [id*Q + k] -> neighbor id
};
