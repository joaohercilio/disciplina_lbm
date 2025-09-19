#include "Simulation.hpp"

using namespace std;

// ----- PROBLEM SETUP ----- //

/*
    Nathen, Patrick & Gaudlitz, Daniel & Krause, Mathias & Adams, Nikolaus. (2018).
    On the Stability and Accuracy of the BGK, MRT and RLB Boltzmann Schemes for the Simulation of 
    Turbulent Flows. Communications in Computational Physics. 23. 846-876. 10.4208/cicp.OA-2016-0229. 
*/

namespace user {
    
    // Lattice units
    int N = 64;
    double uL = 0.1;
    double Re = 800;
    double nu = N*uL / Re;
    double tau = 3*nu + 0.5;

    double Ma = uL / (1.0/sqrt(3.0));
    
    // Physical units
    double u_phy = 1.0;
    double L = 2.0*M_PI;
    double nu_phy = u_phy*L/u_phy;

    double dx = L / N; 
    double dt = nu/nu_phy * dx*dx;

    using LatticeModel = D3Q19; 

    using CollisionModel = BGK;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }
    
    inline std::vector<std::pair<std::string, double>> userLogParams() {
        return {
            {"N", N},
            {"Reynolds", Re},
            {"nu", (1.0/3.0) * (tau - 0.5)},
            {"Ma", Ma},
            {"deltaX", dx},
            {"deltaT", dt}
        };
    }    

    std::vector<double> initialVelocity(const Geometry& geo)
    {
        double L = 2.0*M_PI;

        int nx = geo.nx();
        int ny = geo.ny();
        int nz = geo.nz();

        std::vector <double> u (3 * geo.getNumOfPoints(), 0.0);

        for ( int y = 0; y < geo.ny(); y++ )
        {
            for ( int x = 0; x < geo.nx(); x++ )
            {		
                for ( int z = 0; z < geo.nz(); z ++) 
                {
                    int id = geo.getIndex(x, y, z);

                    double xi = L / geo.nx() * (x + 0.5);
                    double yi = L / geo.ny() * (y + 0.5);
                    double zi = L / geo.nz() * (z + 0.5);

                    double vx  = uL * (- 2.0 / sqrt(3.0) * sin( 2.0/3.0 * M_PI) * sin(xi) * cos(yi) * cos(zi));
                    double vy  = uL * (  2.0 / sqrt(3.0) * sin(-2.0/3.0 * M_PI) * cos(xi) * sin(yi) * cos(zi));
                    
                    u[geo.getVelocityIndex(id, 0)] = vx;
                    u[geo.getVelocityIndex(id, 1)] = vy;
                    u[geo.getVelocityIndex(id, 2)] = 0.0;
                }
            }
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

    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }

    Geometry problemGeometry() {
        Geometry geo(N, N, N);
        return geo;
    }

    int totalSteps () {
        return 0;
    }

    OutputType outputType() {
        return OutputType::VTI;  
    }
    
    int writeInterval() {
        return 1; 
    }
}