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
    std::vector<double> initialVelocity(int x, int y, int z);

    double initialDensity (int x, int y, int z) {
        return 1.0; // Uniform density
    }

    // Geometry definition
    Geometry problemGeometry() {
        Geometry geo(10, 1, 1);
        // geo.setSolid(0, 0, 0);
        return geo;
    }

    // Total number of time steps
    int totalSteps () {
        return 10;
    }

    // Output options
    int writeInterval() {
        return 2; // Write output every 2 steps
    }
}