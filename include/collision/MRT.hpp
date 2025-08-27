#pragma once

#include "CollisionModel.hpp"

/**
 * @brief Class that defines the Multi Relatxation Time (MRT) collision operator methods.
 */
class MRT : public CollisionModel {

public:
    /**
     * @brief Prepares the collision parameters for the MRT operator
     * 
     * The MRT operator requires only the "tau" parameter, wich must be bigger than 0.5 
     * 
     * @param raw Constant reference to a map of raw collision parameters (ColParamMap).
     */
    ColParamMap prepareColParams(const ColParamMap& raw) const override;

    /**
     * @brief Computes the collision all over the domain by calculating the moments of the distribution function 
     * and colliding at the moment space with individual relaxation rates.
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param colParams Constant reference to the collision parameters map specific to the Collision Model
     * @param force Constant vector storing the (x, y, z) components of the external force
     */
    void computeCollision(std::vector<double>& f,
                          const LatticeModel& lattice,
                          const Geometry& geometry,
                          const ColParamMap& colParams,
                          const std::vector<double>& force) override;

    /**
     * @brief Initializes the density field compatible with the imposed velocity
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param colParams Constant reference to the collision parameters map specific to the Collision Model
     * 
     * @note This method is pure virtual and defines the interface for all collision models.
     */
    virtual void initializeDensityField(std::vector<double>& f, 
                                        std::vector<double>& fn,
                                        const LatticeModel& lattice, 
                                        const Geometry& geometry,
                                        const ColParamMap& colParams,
                                        const std::vector<double>& u,
                                        const std::vector<double>& force,
                                        const int numberOfIterations) override;
};