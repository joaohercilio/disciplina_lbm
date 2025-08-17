#include "Simulation.hpp"

int main(int argc, char* argv[]) 
{    
    auto simulation = Simulation::create();

    simulation->run();  
}

/*
void stepLBM(const D2Q9& lattice, double* f, double* fn, double* feq, double* rho, double* v, int NX, int NY, double omega, double FX, double FY, bool print, int printInterval, int t)
{
    // Compute macroscopic quantities
    computeFields(f, rho, v, NX, NY);

    // Compute the equilibrium distribution function
    lattice.computeEquilibrium(feq, rho, v, NX, NY, FX, FY);

    // Perform collision (BGK operator)
    computeCollision(f, feq, rho, v, omega, NX, NY);

    // Perform propagation (with periodic boundary conditions)
    streaming(f, fn, NX, NY);

    // Apply specific boundary conditions
    // boundaryCondition_2Dchannel(fn, NY);

    if (t%printInterval == 0 && print) {
        writeFields(fn, rho, v, NX, NY, t);
        //writeDistributionFunction(fn, NX, NY, t);
    }

    if (t%50 == 0) {
        std::cout << "Step: " << t << std::endl;
        checkMassConservation(rho, NX, NY);
        checkMomentumConservation(rho, v, NX, NY);
    }
}

void checkMassConservation(double* rho, int NX, int NY)
{
    double totalMass = 0.0;

    for (int i = 0; i < NX*NY; i++) {
        totalMass += rho[i];
    }
    std::cout << "Mass: " << totalMass << std::endl; 
}

void checkMomentumConservation(double* rho, double* v, int NX, int NY)
{
    double momentumX = 0.0;
    double momentumY = 0.0;

    for (int i = 0; i < NX*NY; i++) {
        momentumX += rho[i] * v[i*NUM_OF_DIR + 0];
        momentumY += rho[i] * v[i*NUM_OF_DIR + 1];
    }

    std::cout << "X-momentum: " << momentumX << std::endl;
    std::cout << "Y-momentum: " << momentumY << std::endl << std::endl; 
}

void initialCondition_pulse(double* f, int k)
{
    f[k] = 1.0;
}

void initialCondition_2Dchannel(double* f, double* rho, double* v, int NY)
{
    for (int y = 0; y < NY; y++) {
        int idx = y;
        rho[idx] = 1.0;

        double vx = 0;
        double vy = 0;

        // Compute equilibrium distribution functions from initial velocity field
        for (int k = 0; k < NUM_OF_VEL; k++) {
            double c_dot_v = (cx[k] * vx + cy[k] * vy);
            f[idx*NUM_OF_VEL + k] = rho[idx]*w[k] * (1 + 3*c_dot_v + 4.5*c_dot_v*c_dot_v - 1.5*(vx*vx + vy*vy));
        }
    }
}

void boundaryCondition_2Dchannel(double* fn, int NY)
{
    int idxN = NY - 1;
    int idxS = 0;

    //FULL WAY BOUNCE AND BACK
    // North boundary condition
    fn[idxN * NUM_OF_VEL + 4] = fn[idxN * NUM_OF_VEL + 3] ;
    fn[idxN * NUM_OF_VEL + 6] = fn[idxN * NUM_OF_VEL + 5] ;
    fn[idxN * NUM_OF_VEL + 8] = fn[idxN * NUM_OF_VEL + 7] ;

    // South boundary condition
    fn[idxS * NUM_OF_VEL + 3] = fn[idxS * NUM_OF_VEL + 4] ;
    fn[idxS * NUM_OF_VEL + 5] = fn[idxS * NUM_OF_VEL + 6] ;
    fn[idxS * NUM_OF_VEL + 7] = fn[idxS * NUM_OF_VEL + 8] ;
}
*/