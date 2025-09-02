#pragma once

#include <memory>
#include <cmath>

#include "lattice/D2Q5.hpp"
#include "lattice/D2Q9.hpp"
#include "lattice/D3Q19.hpp"
#include "collision/BGK.hpp"
#include "collision/MRT.hpp"
#include "Geometry.hpp"
#include "Initialization.hpp"
#include "Streaming.hpp"
#include "Output.hpp"
#include "Timer.hpp"
#include "Logger.hpp"

/**
 * @brief Abstract base class representing a general LBM simulation.
 *
 * This class manages the lattice model, collision model, geometry, boundary model,
 * and distribution function arrays. It defines the interface for initializing
 * the simulation and running the time loop.
 */
class Simulation {

protected:

    std::unique_ptr<LatticeModel> lattice_;     ///< Pointer to the lattice model
    std::unique_ptr<CollisionModel> collision_; ///< Pointer to the collision model
    Geometry geometry_;                         ///< Geometry of the computational domain
    std::vector<double> f1_;                    ///< Distribution function at current step
    std::vector<double> f2_;                    ///< Distribution function at next step

public:

    /**
     * @brief Constructs a Simulation object with the given models and geometry.
     *
     * Initializes the distribution function vectors to zero.
     *
     * @param lattice Unique pointer to a lattice model
     * @param collision Unique pointer to a collision model
     * @param geometry Geometry object describing the computational domain
     * @param boundary Unique pointer to a boundary model
     */
    Simulation(std::unique_ptr<LatticeModel> lattice,
               std::unique_ptr<CollisionModel> collision,
               Geometry geometry)

        : lattice_(std::move(lattice)),
          collision_(std::move(collision)),
          geometry_(std::move(geometry))
    {
    const int totalPoints = geometry_.getNumOfPoints();
    const int numOfVel    = lattice_->getNumOfVel();
    f1_.resize(totalPoints * numOfVel, 0.0);
    f2_.resize(totalPoints * numOfVel, 0.0);
    }

    /// Virtual destructor
    virtual ~Simulation() = default;

    /**
     * @brief Runs the simulation time loop.
     *
     * Pure virtual method that must be implemented in derived classes.
     */
    virtual void run() = 0;

    /**
     * @brief Returns a constant reference to the lattice model.
     * @return Reference to the lattice model
     */
    const LatticeModel& getLattice() const;
    
    /**
     * @brief Returns a constant reference to the geometry.
     * @return Reference to the geometry object
     */
    const Geometry& getGeometry() const;

    /**
     * @brief Factory method to create a concrete simulation instance.
     *
     * @return Unique pointer to a Simulation object
     */
    static std::unique_ptr<Simulation> create();
};