#pragma once

#include <memory>
#include <cmath>
#include <iostream>
#include <chrono>

#include "lattice/D2Q5.hpp"
#include "lattice/D2Q9.hpp"
#include "lattice/D3Q19.hpp"

#include "collision/BGK.hpp"

#include "boundary/HalfwayBounceAndBack.hpp"
#include "boundary/Periodic.hpp"

#include "Geometry.hpp"
#include "Streaming.hpp"
#include "Output.hpp"

class Simulation {

protected:
    std::unique_ptr<LatticeModel> lattice_;
    std::unique_ptr<CollisionModel> collision_;
    Geometry geometry_;
    std::unique_ptr<BoundaryModel> boundary_;
    std::vector<double> f1_;
    std::vector<double> f2_;

public:
    Simulation(std::unique_ptr<LatticeModel> lattice, std::unique_ptr<CollisionModel> collision, const Geometry& geometry, std::unique_ptr<BoundaryModel> boundary)
        : lattice_(std::move(lattice)),
          collision_(std::move(collision)),
          geometry_(geometry),
          boundary_(std::move(boundary)),
          f1_(geometry_.getNumOfPoints() * lattice_->getNumOfVel(), 0.0),
          f2_(geometry_.getNumOfPoints() * lattice_->getNumOfVel(), 0.0) {}

    virtual ~Simulation() = default;

    virtual void initialize() = 0;
    virtual void run() = 0;

    const LatticeModel& getLattice() const;
    const Geometry&     getGeometry() const;

    static std::unique_ptr<Simulation> create();
};