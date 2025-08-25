#include "Simulation.hpp"

// ----- PROBLEM SETUP (the user should edit this namespace) ----- //

namespace user {
    
    const int N = 16;

    using LatticeModel = D3Q19; 

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
        return {0.0, 0.0, 0.0};
    }

    // Geometry definition and boundary conditions
    using BoundaryModel = PeriodicBoundary;
    Geometry problemGeometry() {
        Geometry geo(N, N, N);
        return geo;
    }

    // Total number of time steps
    int totalSteps () {
        return 1e3;
    }

    // Output options
    int writeInterval() {
        return 100; // Write output every 2 steps
    }

    void print() {}
}