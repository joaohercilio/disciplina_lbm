#include "Simulation.hpp"

#include <iostream>

// ----- PROBLEM SETUP (the user should edit this namespace) ----- //

namespace user {
    
    Streaming streaming;

    // Choose Lattice Model
    using LatticeModel = D2Q9; 

    // Choose Collision Model and set collision parameters
    using CollisionModel = BGK;
    double tau = 0.8;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }

    const int N = 8;
    const double nu_latt = 1.0/3.0 * (tau - 0.5);
    const double nu_phy = 0.001;
    const double L = 1.0;
    const double gx = 1.0;
    const double Umax = L*L*gx / (8.0 * nu_phy);
    const double ReMax = L*Umax / nu_phy; 

    inline void print() {
        std::cout << "Re_max = " << ReMax << std::endl;
    }

    // Initial conditions
    std::vector<double> initialVelocity(int x, int y, int z)
    {
        return {0.0, 0.0, 0.0}; 
    }

    double initialDensity (int x, int y, int z) {
        return 1.0; // Uniform density
    }

    // External forces
    std::vector<double> externalForce()
    {
        return {gx, 0.0, 0.0};
    }

    // Geometry definition and boundary conditions
    using BoundaryModel = BounceAndBack;
    Geometry problemGeometry() {
        Geometry geo(1, 8, 1);
        //geo.setSolid(0, 0, 0);
        //geo.setSolid(0, 7, 0);

        // Boundary conditions
        geo.setBoundaryType(0, 0, 0, BoundaryType::BounceAndBackSouth);
        geo.setBoundaryType(0, 7, 0, BoundaryType::BounceAndBackNorth);
        
        return geo;
    }

    // Total number of time steps
    int totalSteps () {
        return 100000;
    }

    // Output options
    int writeInterval() {
        return 100000; 
    }
}