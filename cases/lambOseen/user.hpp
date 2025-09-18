#include "Simulation.hpp"

using namespace std;

// ----- PROBLEM SETUP ----- //

namespace user {
    
    const int N = 32;

    using LatticeModel = D2Q9; 

    using CollisionModel = MRT;
    double tau = 0.506;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }
    
    inline std::vector<std::pair<std::string, double>> userLogParams() {
        return {
            {"N", N},
            {"Viscosity", (1.0/3.0) * (tau - 0.5)}
        };
    }    

    std::vector<double> initialVelocity(const Geometry& geo)
    {
    double Lx = 1.0;
    double Ly = 1.0;
    double Gamma = 0.1;     // circulação
    double rc = 0.05;       // raio de núcleo
    int nx = geo.nx();
    int ny = geo.ny();

    double hx = Lx / nx;
    double hy = Ly / ny;

    // centro no meio do domínio
    double x0 = 0.5 * Lx;
    double y0 = 0.5 * Ly;

    std::vector<double> u(3 * geo.getNumOfPoints(), 0.0);

    for (int j = 0; j < ny; ++j)
    {
        double y = (j + 0.5) * hy;
            for (int i = 0; i < nx; ++i)
            {
                double x = (i + 0.5) * hx;

                double dx = x - x0;
                double dy = y - y0;
                double r2 = dx*dx + dy*dy;

                double ux = 0.0;
                double uy = 0.0;

                if (r2 > 1e-12) // evita singularidade no centro
                {
                    double factor = (Gamma / (2.0 * M_PI * r2)) * (1.0 - exp(-r2 / (rc*rc)));
                    ux = -factor * dy;
                    uy =  factor * dx;
                }

                u[geo.getVelocityIndex(geo.getIndex(i, j, 0), 0)] = ux;
                u[geo.getVelocityIndex(geo.getIndex(i, j, 0), 1)] = uy;
                u[geo.getVelocityIndex(geo.getIndex(i, j, 0), 2)] = 0.0;
            }
        }
        return u; 
    }

    std::vector<double> initialDensity (const Geometry& geo) {

        std::vector <double> rho (geo.getNumOfPoints(), 0.0);
        return rho;
    }

    int initializePressureIterations() {
        return 1000;
    }

    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }

    Geometry problemGeometry() {
        Geometry geo(N, N, 1);
        return geo;
    }

    int totalSteps () {
        return 200;
    }

    OutputType outputType() {
        return OutputType::BOTH;  
    }
    
    int writeInterval() {
        return 1; 
    }
}