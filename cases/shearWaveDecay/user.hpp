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
    
    // Problem variables and functions
    const int N = 64;
    const double rho0 = 1.0;
    const double nu = 1.0/3.0 * (tau - 0.5);
    const double nu_phy = 0.001;
    const double Re = 10.0;
    const double L = 1.0;
    const double U_phy = Re*nu_phy/L;
    const double finalPhysicalTime = L*L/(4*M_PI*M_PI*nu_phy)*std::log(2);
    const double h = L/N;
    const double delta = nu/nu_phy * h*h;
    const double U_latt = U_phy*delta/h;
    const double finalSymTime = finalPhysicalTime / delta; 
    const int maxSteps = ceil(finalSymTime);
    const double Ma = U_latt / (1.0 / std::sqrt(3));

    double vy (int x) {
        return U_latt * std::sin(2.0*M_PI/N * (x+0.5));
    }

    // Initial conditions
    std::vector<double> initialVelocity(int x, int y, int z) {
        return {0, vy(x), 0};
    }

    double initialDensity (int x, int y, int z) {
        return rho0;
    }

    // Geometry definition
    using BoundaryModel = BounceAndBack;
    Geometry problemGeometry() {
        Geometry geo(N, 1, 1);
        // geo.setSolid(0, 0, 0);
        
        return geo;
    }

    // External forces
    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }
    
    // Total number of time steps
    int totalSteps () {
        return maxSteps;
    }

    // Output options
    int writeInterval() {
        return maxSteps / 10;
    }
}