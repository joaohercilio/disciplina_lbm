#pragma once

#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"

struct CollisionParameters {

    double alphaEq;
    double alphaNonEq;

};

class CollisionModel {
    
public:
    virtual ~CollisionModel() = default;

    virtual void computeCollision(std::vector<double>& f, const CollisionParameters& colParams, const LatticeModel& lattice, const Geometry& geometry) const = 0;

};