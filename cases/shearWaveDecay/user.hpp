#include "Simulation.hpp"

// ----- PROBLEM SETUP  ----- //

namespace user {
    
    using LatticeModel = D2Q9; 

    using CollisionModel = BGK;
    double tau = 0.8;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }
    
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

    inline std::vector<std::pair<std::string, double>> userLogParams()
    {
        return {
            {"N", N},
            {"Viscosity", nu}
        };
    }


    double vy (int x) {
        return U_latt * std::sin(2.0*M_PI/N * (x+0.5));
    }

    std::vector<double> initialVelocity(const Geometry& geo) {
        std::vector <double> u (3 * geo.getNumOfPoints(), 0.0);
        for (int id = 0; id < geo.getNumOfPoints(); id++)
        {
            int x,y;
            geo.getCoords(id, x, y);
            u[geo.getVelocityIndex(id, 0)] = 0.0;
            u[geo.getVelocityIndex(id, 1)] = vy(x);
            u[geo.getVelocityIndex(id, 2)] = 0.0;
        }
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
        Geometry geo(N, 1, 1);
        return geo;
    }

    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }
    
    int totalSteps () {
        return maxSteps;
    }

    OutputType outputType() {
        return OutputType::BOTH;  
    }
    
    int writeInterval() {
        return maxSteps / 10;
    }
}