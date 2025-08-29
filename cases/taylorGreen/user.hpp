#include "Simulation.hpp"

using namespace std;

// ----- PROBLEM SETUP ----- //

namespace user {
    
    const int N = 32;

    using LatticeModel = D2Q9; 

    // Choose Collision Model and set collision parameters
    using CollisionModel = MRT;
    double tau = 0.506;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }
    
    // Initial conditions
    std::vector<double> initialVelocity(const Geometry& geo)
    {
        double Lx = 1.0;
        double Ly = 1.0;
        double U0 = 0.05;

        int nx = geo.nx();
        int ny = geo.ny();
        int nz = geo.nz();

        std::vector <double> u (3 * geo.getNumOfPoints(), 0.0);

        double kx = 2*M_PI/Lx;
        double ky = 2*M_PI/Ly;

        for ( int y = 0; y < geo.ny(); y++ )
        {
            for ( int x = 0; x < geo.nx(); x++ )
            {		
                int id = geo.getIndex(x, y, 0);

                double xi = Lx / geo.nx() * (x + 0.5);
                double yi = Ly / geo.ny() * (y + 0.5);

                double vx  = -sqrt(kx/ky) * U0 * cos(kx * xi) * sin(ky * yi);
                double vy  =  sqrt(kx/ky) * U0 * sin(kx * xi) * cos(ky * yi);
                double p = -0.25 * U0*U0 * ( (ky/kx)*cos(2*kx*xi) + (kx/ky)*cos(2*ky*yi) );
                
                u[geo.getVelocityIndex(id, 0)] = vx;
                u[geo.getVelocityIndex(id, 1)] = vy;
                u[geo.getVelocityIndex(id, 2)] = 0.0;
            }
        }

        return u; 
    }

    std::vector<double> initialDensity (const Geometry& geo) {

        std::vector <double> rho (geo.getNumOfPoints(), 1.0);

        double Lx = 1.0;
        double Ly = 1.0;
        double U0 = 0.05;

        int nx = geo.nx();
        int ny = geo.ny();
        int nz = geo.nz();

        double kx = 2*M_PI/Lx;
        double ky = 2*M_PI/Ly;

        for ( int y = 0; y < geo.ny(); y++ )
        {
            for ( int x = 0; x < geo.nx(); x++ )
            {		
                int id = geo.getIndex(x, y, 0);

                double xi = Lx / geo.nx() * (x + 0.5);
                double yi = Ly / geo.ny() * (y + 0.5);

                double p = -0.5 * U0*U0 * cos( kx*(xi + yi) ) * cos( kx*(xi - yi) );

                //rho[id] = 1.0 + 3.0 * p;
            }
        }
        return rho;
    }

    // Number of iterations on the Density Field Initialization procedure
    int initializePressureIterations() {
        return 10000;
    }

    // External forces
    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }

    // Geometry definition and boundary conditions
    using BoundaryModel = PeriodicBoundary;
    Geometry problemGeometry() {
        Geometry geo(N, N, 1);
        return geo;
    }

    // Total number of time steps
    int totalSteps () {
        return 200;
    }

    // Output options
    enum class OutputType { TSV, VTI, BOTH, NONE };

    OutputType outputType() {
        return OutputType::BOTH;  
    }
    
    int writeInterval() {
        return 1; 
    }

    void print() {}
}