#include "lattice/D2Q9.hpp"

void D2Q9::computeEquilibrium(double* feq, const double rho, const double vx, const double vy, const double vz) const {
    
    double vx2 = vx*vx;
    double vy2 = vy*vy;
    double rho0 = w_[0]*rho;
    double rho1 = w_[1]*rho;
    double rho2 = w_[5]*rho;
    double ax = 3.0 * vx;
    double ay = 3.0 * vy;
    double b1 =  3.0 * vx + 3.0 * vy;
    double b2 = -3.0 * vx + 3.0 * vy;

    feq[0] = 1.0 - 1.5 * ( vx2 + vy2 );

    feq[1] = feq[0] + 4.5 * vx2;
    feq[2] = rho1 * ( feq[1] - ax );
    feq[1] = rho1 * ( feq[1] + ax );

    feq[3] = feq[0] + 4.5 * vy2;
    feq[4] = rho1 * ( feq[3] - ay );
    feq[3] = rho1 * ( feq[3] + ay );
    
    feq[5] = feq[0] + 0.5 * b1 * b1;
    feq[6] = rho2 * ( feq[5] - b1 );
    feq[5] = rho2 * ( feq[5] + b1 );

    feq[7] = feq[0] + 0.5 * b2 * b2;
    feq[8] = rho2 * ( feq[7] - b2 );
    feq[7] = rho2 * ( feq[7] + b2 );

    feq[0] = rho0 * feq[0];
}

void D2Q9::computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const {
    
    rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];

    double mx =  f[1] - f[2] + f[5] - f[6] - f[7] + f[8];
    double my =  f[3] - f[4] + f[5] - f[6] + f[7] - f[8];

    double oneOverRho = 1.0 / rho;

    vx = mx * oneOverRho;
    vy = my * oneOverRho;
    vz = 0.0;
}
