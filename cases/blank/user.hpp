#include "Simulation.hpp"

// ----- PROBLEM SETUP  ----- //

namespace user {
    
    using LatticeModel = D2Q9; 

    using CollisionModel = BGK;
    double tau = 0.8;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }
    
    int N = 10;

    std::vector<double> initialVelocity(const Geometry& geo) {
        std::vector <double> u (3 * geo.getNumOfPoints(), 0.0);
        return u;
    }

    std::vector<double> initialDensity (const Geometry& geo) {
        std::vector <double> rho (geo.getNumOfPoints(), 1.0);
        return rho;
    }
    
    int initializePressureIterations() {
        return 0;
    }

    Geometry problemGeometry() {
        Geometry geo(N, N, N);
        return geo;
    }

    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }
    
    int totalSteps () {
        return 0;
    }

    OutputType outputType() {
        return OutputType::BOTH;  
    }

    int writeInterval() {
        return 0;
    }

    inline std::vector<std::pair<std::string, double>> userLogParams() {
        return {
            {"N", N},
            {"Viscosity", (1.0/3.0) * (tau - 0.5)}
        };
    }    

}