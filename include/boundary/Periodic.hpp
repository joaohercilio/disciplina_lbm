#include "BoundaryModel.hpp"

/**
 * @brief Class that declares methods for Periodic Boundary conditions. 
 * 
 * Since the streaming process already handles the distribution function periodically 
 * on the domain boundaries, the methods declared here have no functionality.
 */
class PeriodicBoundary : public BoundaryModel {

public:

    /**
     * @brief Applies periodic boundary conditions
     * 
     * This method does nothing because streaming already handles periodicity.
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param x X-coordinate of the Solid node
     * @param y Y-coordinate of the Solid node
     * @param z Z-coordinate of the Solid node
     */
    void applyBoundary(std::vector<double>& f, 
                        const LatticeModel& lattice, 
                        const Geometry& geometry,
                        int x, int y, int z) override;

private:

    /**
     * @brief Applies periodic boundary conditions 
     * 
     * This method does nothing
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param x X-coordinate of the Solid node
     * @param y Y-coordinate of the Solid node
     * @param z Z-coordinate of the Solid node
     */
    void applySouthBoundary(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry,
                            int x, int y, int z) override;

    /**
     * @brief Applies periodic boundary conditions 
     *      
     * This method does nothing
     * 
     * @param f Vector storing the particle distribution function 
     * @param lattice Constant reference to the lattice object 
     * @param geometry Constant reference to the geometry object
     * @param x X-coordinate of the Solid node
     * @param y Y-coordinate of the Solid node
     * @param z Z-coordinate of the Solid node
     */
    void applyNorthBoundary(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry,
                            int x, int y, int z) override;

};
