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

    double rho = m[0];

    double ux = m[1] / rho;
    double uy = m[2] / rho;
    
    double ux2 = ux * ux;
    double uy2 = uy * uy;

    meq[0] = rho;

    meq[1] = rho * ux;	
    meq[2] = rho * uy;
    
    meq[3] = rho * ( 3*ux2 + 3*uy2 - 2); // e
    meq[4] = rho * (-3*ux2 - 3*uy2 + 1); // epsilon

    meq[5] = -rho * ux;	//qx
    meq[6] = -rho * uy;	//qy

    meq[7] = rho * (ux2 - uy2); //pxx
    meq[8] = rho *   ux*uy;		//pxy
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

std::vector<double> D2Q9::relaxationMatrix(const double tau, const int numOfVel) const {
	
	std::vector<double> s(numOfVel, 0.0);
	s[3] =  1.0; 	 //s_e
	s[4] =  1.4;	 //s_epsilon
	s[5] =  1.7;	 //s_q	
	s[6] =  1.7;	 //s_q
	s[7] =  1./tau;	 //s_nu
	s[8] =  1./tau;	 //s_nu
	return s;
}