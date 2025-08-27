#pragma once

#include "lattice/LatticeModel.hpp"

/**
 *  @brief Forward declaration that allows methods in this class to use Geometry objects
 */ 
class Geometry;

/**
 * @brief Enumerator that lists in what direction the Fluid node is facing an ajdacent Solid node.
 */
enum class BoundaryType
{ 
    None,
    HalfwayBounceAndBackSouth,
    HalfwayBounceAndBackNorth 
};

/**
 *  @brief Abstract base class to Boundary Models definitions
 * 
 * This class offers basic infrastructure to define diferent boundary conditions
 *  
 */
class BoundaryModel {

public:

    // Default destructor
    virtual ~BoundaryModel() = default;

    /**
     * @brief Applies boundary conditions at a Solid node
     * 
     * This method modifies the distribution function according to the 
     * boundary model at the specific Solid node
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param x X-coordinate of the Solid node
     * @param y Y-coordinate of the Solid node
     * @param z Z-coordinate of the Solid node
     * 
     * @note This method is pure virtual and defines the interface for all boundary models.
     */
    virtual void applyBoundary(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry) = 0;

private:

    /**
     * @brief Applies boundary conditions at a Solid node located at the south boundary of the domain
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * 
     * @note This method is pure virtual and defines the interface for all boundary models.
     */
    virtual void applySouthBoundary(std::vector<double>& f, 
                                    const LatticeModel& lattice, 
                                    const Geometry& geometry,
                                    int id) = 0;

    /**
     * @brief Applies boundary conditions at a Solid node located at the north boundary of the domain
     *      
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param id Index of the node
     * 
     * @note This method is pure virtual and defines the interface for all boundary models.
     */
    virtual void applyNorthBoundary(std::vector<double>& f, 
                                    const LatticeModel& lattice, 
                                    const Geometry& geometry,
                                    int id) = 0;

};