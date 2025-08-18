#include "Simulation.hpp"

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
        return {0.1, 0.0, 0.0};
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
        return 10;
    }

    // Output options
    int writeInterval() {
        return 1; // Write output every 2 steps
    }
}