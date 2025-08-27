#include "Simulation.hpp"


// ----- PROBLEM SETUP (the user should edit this namespace) ----- //

namespace user {
   
    using LatticeModel = D2Q9; 

    using CollisionModel = BGK;
    double tau = 0.8;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }

    // Problem parameters
    const int N             = 32 + 2; // 2 solid nodes at the boundaries
    const double L_phy      = 1.0;
    const double Umax_phy   = 0.01;
    const double nu_phy     = 1e-3;
    
    const double nu_latt    = (1.0/3.0)*(tau - 0.5);
    const double L_latt     = N - 2;

    const double h          = L_phy / (N-2);
    const double delta      = nu_latt / nu_phy * h*h;
    
    const double Umax_latt  = Umax_phy * delta / h;
    const double gx_latt    = 8.0 * nu_latt * Umax_latt / (L_latt * L_latt);

    inline void print() {
    }

    // Initial conditions
    std::vector<double> initialVelocity(const Geometry& geo)
    {
        std::vector<double> u (geo.getNumOfPoints(), 0.0);
        return u;
    }

    std::vector<double> initialDensity (const Geometry& geo) {
        std::vector<double> rho (geo.getNumOfPoints(), 1.0);
        return rho;
    }

    // External forces
    std::vector<double> externalForce()
    {
        return {gx_latt, 0.0, 0.0};
    }

    // Geometry definition and boundary conditions
    using BoundaryModel = HalfwayBounceAndBack;
    Geometry problemGeometry() {
        Geometry geo(1, N, 1);

        geo.setSolid(0, 0, 0);
        geo.setSolid(0, N-1, 0);

        // Boundary conditions
        geo.setBoundaryType(0,  0 , 0, BoundaryType::HalfwayBounceAndBackSouth);
        geo.setBoundaryType(0, N-1, 0, BoundaryType::HalfwayBounceAndBackNorth);
        
        return geo;
    }

    // Total number of time steps
    int totalSteps () {
        return 1e3;
    }

    // Output options
    int writeInterval() {
        return 100; 
    }
}