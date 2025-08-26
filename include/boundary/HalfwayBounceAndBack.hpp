#include "BoundaryModel.hpp"
#include "Geometry.hpp"

/**
 * @brief Class that define the Halfway Bounce And Back boundary conditions method
 */
class HalfwayBounceAndBack : public BoundaryModel {

public:

    /**
     * @brief Applies Halfway Bounce And Back boundary conditions at a Solid node
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
     */
    void applyBoundary(std::vector<double>& f, 
                        const LatticeModel& lattice, 
                        const Geometry& geometry,
                        int x, int y, int z) override;

private:

    /**
     * @brief Applies boundary conditions at a Solid node located at the south boundary of the domain
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
     * @brief Applies boundary conditions at a Solid node located at the north boundary of the domain
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