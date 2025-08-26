#pragma once

#include <vector>
#include <cmath>

/**
 * @brief Abstract base class to Lattice Models definitions
 * 
 * This class offers basic infrastructure to define different lattice models
 */
class LatticeModel {

protected:
    int numOfDim_;              ///< Number of spatial dimensions (2 or 3)
    int numOfVel_;              ///< Number of discrete velocities
    double cs_;                 ///< Sound speed in lattice units
    std::vector<int> cx_;       ///< Discrete velocities in X-direction
    std::vector<int> cy_;       ///< Discrete velocities in Y-direction
    std::vector<int> cz_;       ///< Discrete velocities in Z-direction (vector of zeros if numOfDim = 2)
    std::vector<double> w_;     ///< Lattice weights

public:

    /**
     * @brief Constructor of a lattice model with given parameters
     * 
     * @param numOfDim Number of spatial dimensions (2 or 3)
     * @param numOfVel Number of discrete velocities
     * @param cs Sound speed in lattice units
     * @param cx Vector of discrete velocities in X-direction
     * @param cy Vector of discrete velocities in Y-direction
     * @param cz Vector of discrete velocities in Z-direction
     * @param w Vector of lattice weights
     */
    LatticeModel(
        int numOfDim, int numOfVel, double cs,
        std::vector<int> cx, std::vector<int> cy, std::vector<int> cz,
        std::vector<double> w
    ):  numOfDim_(numOfDim),
        numOfVel_(numOfVel),
        cs_(cs),
        cx_(cx),
        cy_(cy),
        cz_(cz),
        w_(w)
    {}
 
    /// Virtual default destructor
    virtual ~LatticeModel() = default;

    // ====== Getters =====

    /// Returns the number of discrete velocities
    int getNumOfVel() const { return numOfVel_; }

    /// Returns the number of spatial dimensions
    int getNumOfDim() const { return numOfDim_; }

    /// Returns a constant reference to the x-component velocities
    const std::vector<int>& getCx() const { return cx_; }

    /// Returns a constant reference to the y-component velocities
    const std::vector<int>& getCy() const { return cy_; }

    /// Returns a constant reference to the z-component velocities
    const std::vector<int>& getCz() const { return cz_; }    

    /// Returns a constant reference to the lattice weights
    const std::vector<double>& getWeights() const { return w_; }

    /// Returns the lattice speed of sound
    double getCs() const { return cs_; }

    // ===== Pure virtual methods =====

    /**
     * @brief Computes macroscopic density and velocity from the distribution function.
     *
     * @param f Pointer to the vector of distribution functions
     * @param rho Reference to store computed density
     * @param vx Reference to store X-component of velocity
     * @param vy Reference to store Y-component of velocity
     * @param vz Reference to store Z-component of velocity
     *
     * @note This is a pure virtual method that must be implemented in derived classes.
     */
    virtual void computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const = 0;

    /**
     * @brief Computes the equilibrium distribution for given macroscopic fields.
     *
     * @param feq Pointer to the array to store the equilibrium distribution
     * @param rho Density at the node
     * @param vx X-component of velocity
     * @param vy Y-component of velocity
     * @param vz Z-component of velocity
     *
     * @note Pure virtual; implementation depends on the lattice type.
     */
    virtual void computeEquilibrium(double* feq, const double rho, const double vx, const double vy, const double vz) const = 0;

    /**
     * @brief Computes moments from the distribution function.
     *
     * @param f Pointer to the array of distribution functions
     * @param m Pointer to the array where moments will be stored
     */
    virtual void computeMoments(const double* f, double* m) const = 0;

    /**
     * @brief Computes equilibrium moments from the given moments.
     *
     * @param meq Pointer to the array where equilibrium moments will be stored
     * @param m Pointer to the array of current moments
     */
    virtual void computeEquilibriumMoments(double* meq, const double* m) const = 0;

    /**
     * @brief Reconstructs the distribution function from the moments.
     *
     * @param f Pointer to the array where the distribution function will be stored
     * @param m Pointer to the array of moments
     */
    virtual void reconstructDistribution(double* f, const double* m) const = 0;

    /**
     * @brief Returns the relaxation matrix for the MRT collision model.
     *
     * @param tau Relaxation time
     * @param numOfVel Number of discrete velocities
     * @return Vector containing the relaxation matrix
     */
    virtual std::vector<double> relaxationMatrix(const double tau, const int numOfVel) const = 0;
};
