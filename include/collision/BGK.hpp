#pragma once

#include "CollisionModel.hpp"

class BGK : public CollisionModel {
public:
    void computeCollision(std::vector<double>& f, const CollisionParameters& colParams, const LatticeModel& lattice, const Geometry& geometry) const override;
};