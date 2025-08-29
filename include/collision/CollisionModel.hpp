#pragma once

#include "lattice/LatticeModel.hpp"
#include "Geometry.hpp"

#include <map>
#include <string>
#include <stdexcept>
#include <omp.h>

/**
 * @brief Map that stores collision parameters.
 *
 * The key is the parameter name as a string, and the value is a double representing
 * the parameter value. Used to pass collision parameters to the computeCollision method.
 */
using ColParamMap = std::map<std::string,double>;

/**
 * @brief Abstract base class to Collision models definitions
 * 
 * This class offers basic infrastructure to define different collision models
 * 
 */
class CollisionModel {
    
public:

    // Default destructor
    virtual ~CollisionModel() = default;

    /**
     * @brief Prepares collision parameters for the collision model.
     *
     * This method can be overridden in derived classes to modify or augment
     * the raw collision parameters before using them in the simulation.
     *
     * @param raw Constant reference to a map of raw collision parameters (ColParamMap)
     * @return A ColParamMap containing the prepared collision parameters
     *
     * @note The default implementation simply returns the input parameters unmodified.
     *       Derived classes can override this to apply model-specific adjustments.
     */
    virtual ColParamMap prepareColParams(const ColParamMap& raw) const { return raw; }

    /**
     * @brief Computes the collision of the distribution function all over the Fluid nodes
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param colParams Constant reference to the collision parameters map specific to the Collision Model
     * @param force Constant vector storing the (x, y, z) components of the external force
     * 
     * @note This method is pure virtual and defines the interface for all collision models.
     */
    virtual void computeCollision(std::vector<double>& f, 
                                  const LatticeModel& lattice, 
                                  const Geometry& geometry,
                                  const ColParamMap& colParams,
                                  const std::vector<double>& force) = 0;

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
                                        const int numberOfIterations) = 0;
};