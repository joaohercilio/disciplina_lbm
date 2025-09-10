#include "lattice/D2Q9.hpp"

void D2Q9::computeEquilibrium(double* feq, const double drho, const double vx, const double vy, const double vz) const {
    
    double vx2 = vx*vx;
    double vy2 = vy*vy;
    
    double rho0 = w_[0]*(1.0 + drho);
    double rho1 = w_[1]*(1.0 + drho);
    double rho2 = w_[5]*(1.0 + drho);

    double drho0 = w_[0]*drho;
    double drho1 = w_[1]*drho;
    double drho2 = w_[5]*drho;

    double ax =  3.0 * vx;
    double ay =  3.0 * vy;
    double b1 =  3.0 * vx + 3.0 * vy;
    double b2 = -3.0 * vx + 3.0 * vy;

    feq[0] = - 1.5 * ( vx2 + vy2 );

    feq[1] = feq[0] + 4.5 * vx2;
    feq[2] = rho1 * ( feq[1] - ax ) + drho1;
    feq[1] = rho1 * ( feq[1] + ax ) + drho1;

    feq[3] = feq[0] + 4.5 * vy2;
    feq[4] = rho1 * ( feq[3] - ay ) + drho1;
    feq[3] = rho1 * ( feq[3] + ay ) + drho1;
    
    feq[5] = feq[0] + 0.5 * b1 * b1;
    feq[6] = rho2 * ( feq[5] - b1 ) + drho2;
    feq[5] = rho2 * ( feq[5] + b1 ) + drho2;

    feq[7] = feq[0] + 0.5 * b2 * b2;
    feq[8] = rho2 * ( feq[7] - b2 ) + drho2;
    feq[7] = rho2 * ( feq[7] + b2 ) + drho2;

    feq[0] = rho0 * feq[0] + drho0;
}

void D2Q9::computeFields(const double* f, double& drho, double& vx, double& vy, double& vz) const {
    
    drho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
    
    double mx =  f[1] - f[2] + f[5] - f[6] - f[7] + f[8];
    double my =  f[3] - f[4] + f[5] - f[6] + f[7] - f[8];

    double oneOverRho = 1.0 / (1.0 + drho);

    vx = mx * oneOverRho;
    vy = my * oneOverRho;
    vz = 0.0;
}

void D2Q9::computeMoments(const double* f, double* m) const {

    m[0] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
    m[1] = f[1] - f[2] + f[5] - f[6] - f[7] + f[8];
    m[2] = f[3] - f[4] + f[5] - f[6] + f[7] - f[8];
    m[3] = -4*f[0] - f[1] - f[2] - f[3] - f[4] + 2*f[5] + 2*f[6] + 2*f[7] + 2*f[8];
    m[4] = 4*f[0] - 2*f[1] - 2*f[2] - 2*f[3] - 2*f[4] + f[5] + f[6] + f[7] + f[8];
    m[5] = -2*f[1] + 2*f[2] + f[5] - f[6] - f[7] + f[8];
    m[6] = -2*f[3] + 2*f[4] + f[5] - f[6] + f[7] - f[8];
    m[7] = f[1] + f[2] - f[3] - f[4];
    m[8] = f[5] + f[6] - f[7] - f[8];
}

void D2Q9::computeEquilibriumMoments(double* meq, const double* m) const {

    double drho = m[0];

    double jx = m[1];	
    double jy = m[2];
    
    double jx2 = jx*jx;
    double jy2 = jy*jy;

    meq[0] = drho;  // Mass density

    meq[1] = jx;	// Jx Momentum
    meq[2] = jy;    // Jy Momentum
    
    meq[3] = -2*drho + 3*(jx2 + jy2); // Kinetic energy
    meq[4] =    drho - 3*(jx2 + jy2); // Squared Kinetic energy

    meq[5] = -jx;	// qx Energy flux
    meq[6] = -jy;	// qy Energy flux

    meq[7] = jx2 - jy2; // pxx Viscous stress tensor
    meq[8] = jx*jy;		// pxy Viscous stress tensor
}

void D2Q9::reconstructDistribution(double* f, const double* m) const {
	
    f[0] = (1.0/9.0)*m[0] - 1.0/9.0*m[3] + (1.0/9.0)*m[4];
    f[1] = (1.0/9.0)*m[0] + (1.0/6.0)*m[1] - 1.0/36.0*m[3] - 1.0/18.0*m[4] - 1.0/6.0*m[5] + (1.0/4.0)*m[7];
    f[2] = (1.0/9.0)*m[0] - 1.0/6.0*m[1] - 1.0/36.0*m[3] - 1.0/18.0*m[4] + (1.0/6.0)*m[5] + (1.0/4.0)*m[7];
    f[3] = (1.0/9.0)*m[0] + (1.0/6.0)*m[2] - 1.0/36.0*m[3] - 1.0/18.0*m[4] - 1.0/6.0*m[6] - 1.0/4.0*m[7];
    f[4] = (1.0/9.0)*m[0] - 1.0/6.0*m[2] - 1.0/36.0*m[3] - 1.0/18.0*m[4] + (1.0/6.0)*m[6] - 1.0/4.0*m[7];
    f[5] = (1.0/9.0)*m[0] + (1.0/6.0)*m[1] + (1.0/6.0)*m[2] + (1.0/18.0)*m[3] + (1.0/36.0)*m[4] + (1.0/12.0)*m[5] + (1.0/12.0)*m[6] + (1.0/4.0)*m[8];
    f[6] = (1.0/9.0)*m[0] - 1.0/6.0*m[1] - 1.0/6.0*m[2] + (1.0/18.0)*m[3] + (1.0/36.0)*m[4] - 1.0/12.0*m[5] - 1.0/12.0*m[6] + (1.0/4.0)*m[8];
    f[7] = (1.0/9.0)*m[0] - 1.0/6.0*m[1] + (1.0/6.0)*m[2] + (1.0/18.0)*m[3] + (1.0/36.0)*m[4] - 1.0/12.0*m[5] + (1.0/12.0)*m[6] - 1.0/4.0*m[8];
    f[8] = (1.0/9.0)*m[0] + (1.0/6.0)*m[1] - 1.0/6.0*m[2] + (1.0/18.0)*m[3] + (1.0/36.0)*m[4] + (1.0/12.0)*m[5] - 1.0/12.0*m[6] - 1.0/4.0*m[8];
}

std::vector<double> D2Q9::relaxationMatrix(const double tau) const {
	
	std::vector<double> s(9, 0.0);
    s[1] =  1.0;
	s[2] =  1.0;
	s[3] =  1.0; 	 //s_e
	s[4] =  1.4;	 //s_epsilon
	s[5] =  1.7;	 //s_q	
	s[6] =  1.7;	 //s_q
	s[7] =  1./tau;	 //s_nu
	s[8] =  1./tau;	 //s_nu
	return s;
}