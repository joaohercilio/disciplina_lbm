#include "lattice/D3Q19.hpp"

void D3Q19::computeEquilibrium(double* feq, const double rho, const double vx, const double vy, const double vz) const {
    
    double vx2 = vx*vx;
    double vy2 = vy*vy;
    double vz2 = vz*vz;

    double rho0 = w_[0]*rho;
    double rho1 = w_[1]*rho;
    double rho2 = w_[7]*rho;

    double ax = 3.0 * vx;
    double ay = 3.0 * vy;
    double az = 3.0 * vz;

	feq[1] = feq[0] + 4.5 * vx2;
	feq[2] = rho1 * (feq[1] - ax);
	feq[1] = rho1 * (feq[1] + ax);

	feq[3] = feq[0] + 4.5 * vy2;
	feq[4] = rho1 * (feq[3] - ay);
	feq[3] = rho1 * (feq[3] + ay);

	feq[5] = feq[0] + 4.5 * vz2;
	feq[6] = rho1 * (feq[5] - az);
	feq[5] = rho1 * (feq[5] + az);

	double cv = (ax + ay);
	feq[ 7] = feq[ 0] + 0.5*cv*cv;
	feq[ 8] = rho2 * (feq[ 7] - cv);
	feq[ 7] = rho2 * (feq[ 7] + cv);

	cv = (ax - ay);
	feq[ 9]  = feq[ 0] + 0.5*cv*cv;
	feq[10] = rho2 * (feq[ 9] - cv);
	feq[ 9] = rho2 * (feq[ 9] + cv);

	cv = (ax + az);
	feq[11] = feq[ 0] + 0.5*cv*cv;
	feq[12] = rho2 * (feq[11] - cv);
	feq[11] = rho2 * (feq[11] + cv);

	cv = (ax - az);
	feq[13] = (feq[0] + 0.5*cv*cv);
	feq[14] = rho2* (feq[13] - cv);
	feq[13] = rho2* (feq[13] + cv);

	cv = (-ay - az);
	feq[15] = feq[ 0] + 0.5*cv*cv;
	feq[16] = rho2* (feq[15] - cv);
	feq[15] = rho2* (feq[15] + cv);

	cv = (-ay + az);
	feq[17] = feq[ 0] + 0.5*cv*cv;
	feq[18] = rho2 * (feq[17] - cv);
	feq[17] = rho2 * (feq[17] + cv);

    feq[0] = rho0 * feq[0];
}

void D3Q19::computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const {
    
    rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];

    double mx =  f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14];
	double my =  f[3] - f[4] + f[7] - f[8] - f[9] + f[10] - f[15] + f[16] - f[17] + f[18];
	double mz =  f[5] - f[6] +f[11] -f[12]- f[13] + f[14] - f[15] + f[16] + f[17] - f[18];
    
    double oneOverRho = 1.0 / rho;

    vx = mx * oneOverRho;
    vy = my * oneOverRho;
    vz = mz * oneOverRho;
}
