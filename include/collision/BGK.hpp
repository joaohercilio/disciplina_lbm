#pragma once

#include "CollisionModel.hpp"

class BGK : public CollisionModel {

public:
    
    ColParamMap prepareColParams(const ColParamMap& raw) const override;

    void computeCollision(std::vector<double>& f,
                          const LatticeModel& lattice,
                          const Geometry& geometry,
                          const ColParamMap& colParams,
                          const std::vector<double>& force) override;
};