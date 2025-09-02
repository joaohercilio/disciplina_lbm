#include "Simulation.hpp"

// ----- PROBLEM SETUP ----- //

namespace user {
   
    using LatticeModel = D2Q9; 

    using CollisionModel = BGK;
    double tau = 0.8;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }

    const int N             = 64 + 2; // 2 solid nodes at the boundaries
    const int TPHY          = 10000;
    const double LPHY       = 1.0;
    const double UPHY       = 0.001;
    const double NUPHY      = 1e-3;
    const double REPHY      = (UPHY * LPHY)/NUPHY;
    const double FPHY       = (8 * NUPHY * UPHY) / (LPHY * LPHY);
    
    const double NULATT     = (1.0/3.0) * (tau - 0.5);
    const double h          = LPHY  / (N-2);
    const double delta      = NULATT / NUPHY * h*h;
    const double FLATT      = FPHY * (delta * delta) / h;
    const double ULATT      = UPHY * delta / h;
    const double TLATT      = int(TPHY / delta);
    const double RELATT     = N * ULATT / NULATT;

    // User configuration summary
    inline std::vector<std::pair<std::string, double>> userLogParams() {
        return {
            {"N", N},
            {"Viscosity", NULATT},
            {"h", h},
            {"delta", delta},
            {"Umax", ULATT},
            {"F", FLATT},
            {"Re", RELATT},
            {"Ma", ULATT/(1/sqrt(3))}
        };
    }    
    
    std::vector<double> initialVelocity(const Geometry& geo)
    {
        std::vector<double> u (3*geo.getNumOfPoints(), 0.0);
        return u;
    }

    std::vector<double> initialDensity (const Geometry& geo) {
        std::vector<double> rho (geo.getNumOfPoints(), 1.0);
        return rho;
    }

    std::vector<double> externalForce()
    {
        return {FLATT, 0.0, 0.0};
    }

    Geometry problemGeometry() {
        Geometry geo(1, N, 1);
        geo.setSolid(0, 0, 0);
        geo.setSolid(0, N-1, 0);
        return geo;
    }

    OutputType outputType() {
        return OutputType::TSV;  
    }

    int totalSteps () {
        return 1000000;
    }

    int writeInterval() {
        return TLATT; 
    }

    int initializePressureIterations()
    {
        return 0;
    }
}