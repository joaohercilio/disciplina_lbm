#pragma once

#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"

#include <map>
#include <string>
#include <stdexcept>

using ColParamMap = std::map<std::string,double>;

class CollisionModel {
    
public:
    virtual ~CollisionModel() = default;

    virtual void computeCollision(std::vector<double>& f, 
                                  const LatticeModel& lattice, 
                                  const Geometry& geometry,
                                  const ColParamMap& colParams,
                                  const std::vector<double>& force) = 0;
                                  
    virtual ColParamMap prepareColParams(const ColParamMap& raw) const { return raw; }
};