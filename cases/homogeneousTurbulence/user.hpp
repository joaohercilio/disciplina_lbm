#include "Simulation.hpp"
#include <fftw3.h>
#include <complex>
#include <random>

using namespace std;

// ----- PROBLEM SETUP ----- //

void initVelocity(std::vector<double>& u, const Geometry& geo, double kmin, double kmax, double A, double m);
double energySpectrum_(double k, double kmin, double kmax, double A, int m);
void uHat(double kx, double ky, double kz, complex<double> vHat[3], double kmin, double kmax, double A, double m);

namespace user {
    
    const int N = 32;

    using LatticeModel = D3Q19; 

    // Choose Collision Model and set collision parameters
    using CollisionModel = MRT;
    double tau = 0.8;
    ColParamMap colParams() {
        return {{"tau", tau}};
    }
    
    // Initial conditions
    std::vector<double> initialVelocity(const Geometry& geo)
    {
        double kmin = 3;
        double kmax = 8;
        double A = 1.4293e-4;
        double m = 4;

        std::vector <double> u (3 * geo.getNumOfPoints());
        
        initVelocity(u, geo, kmin, kmax, A, m);
        
        double c = 1.0/5.0;
        for (int i = 0; i < geo.getNumOfPoints()*3; i++) { u[i] *= c; }

        return u; 
    }

    std::vector<double> initialDensity (const Geometry& geo) {
        std::vector <double> rho (geo.getNumOfPoints(), 1.0);
        return rho;
    }

    // External forces
    std::vector<double> externalForce()
    {
        return {0.0, 0.0, 0.0};
    }

    // Geometry definition and boundary conditions
    using BoundaryModel = PeriodicBoundary;
    Geometry problemGeometry() {
        Geometry geo(N, N, N);
        return geo;
    }

    // Total number of time steps
    int totalSteps () {
        return 100;
    }

    // Output options
    int writeInterval() {
        return 1; 
    }

    void print() {}
}

double energySpectrum_(double k, double kmin, double kmax, double A, int m)
{
    if (k <= kmin || k >= kmax) return 0.0;
    return A * pow(k,m-2) * exp(-0.14 * k*k) / (2 * M_PI );
}

void uHat(double kx, double ky, double kz, complex<double> vHat[3], double kmin, double kmax, double A, double m)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis_phi(0, 2 * M_PI);
    uniform_real_distribution<> dis_theta(-M_PI, M_PI);

    double kMag = sqrt(kx * kx + ky * ky + kz * kz);
    double M = sqrt(kx * kx + ky * ky);
    double phi = dis_phi(gen);
    double theta1 = dis_theta(gen);
    double theta2 = dis_theta(gen);

    complex<double> alpha = sqrt(energySpectrum_(kMag, kmin, kmax, A, m));
    complex<double> beta  =  alpha * exp(complex<double>(0, theta2)) * sin(phi);
                    alpha =  alpha * exp(complex<double>(0, theta1)) * cos(phi);

    if (M == 0)
    {
        vHat[0] = alpha;
        vHat[1] = beta;
        vHat[2] = 0.0;
    }

    else
    {
        vHat[0] = (alpha * kMag * ky + beta * kx * kz) / (kMag * M);
        vHat[1] = (beta * ky * kz - alpha * kMag * kx) / (kMag * M);
        vHat[2] = -(beta * M) / kMag;
    }
}

void initVelocity(std::vector<double>& u, const Geometry& geo, double kmin, double kmax, double A, double m)
{
    int Nx = geo.nx();
    int Ny = geo.ny();
    int Nz = geo.nz();

    int Nz_ = (Nz/2 + 1);

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 3 * Nx * Ny * Nz_ );
    double *out = u.data();

    fftw_complex *vHatx = in;
    fftw_complex *vHaty = vHatx + Nx * Ny * Nz_;
    fftw_complex *vHatz = vHaty + Nx * Ny * Nz_;

    double* vx = out;
    double* vy = vx + Nx * Ny * Nz;
    double* vz = vy + Nx * Ny * Nz;

    fftw_plan planX = fftw_plan_dft_c2r_3d(Nx,Ny,Nz, vHatx , vx, FFTW_ESTIMATE);
    fftw_plan planY = fftw_plan_dft_c2r_3d(Nx,Ny,Nz, vHaty , vy, FFTW_ESTIMATE);
    fftw_plan planZ = fftw_plan_dft_c2r_3d(Nx,Ny,Nz, vHatz , vz, FFTW_ESTIMATE);

    double K = 0;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            for (int z = 0; z < Nz_; z++)
            {
                double kx =  (x < Nx/2) ? x : x - Nx;
                double ky =  (y < Ny/2) ? y : y - Ny;
                double kz =  z;

                complex<double> vHat[3];
                uHat(kx, ky, kz, vHat, kmin, kmax, A, m);
                vHatx[x * Ny * Nz_ + y * Nz_ + z][0] = vHat[0].real();
                vHatx[x * Ny * Nz_ + y * Nz_ + z][1] = vHat[0].imag();
                vHaty[x * Ny * Nz_ + y * Nz_ + z][0] = vHat[1].real();
                vHaty[x * Ny * Nz_ + y * Nz_ + z][1] = vHat[1].imag();
                vHatz[x * Ny * Nz_ + y * Nz_ + z][0] = vHat[2].real();
                vHatz[x * Ny * Nz_ + y * Nz_ + z][1] = vHat[2].imag();
		        K += (z == 0 ? 0.5 : 1 ) * ( norm(vHat[0]) +  norm(vHat[1]) + norm(vHat[2]) );
            }
        }
    }

    cout << "Mean kineric energy (Fourier Space): " << K << endl;

    fftw_execute(planX);
    fftw_execute(planY);
    fftw_execute(planZ);

    fftw_destroy_plan(planX);
    fftw_destroy_plan(planY);
    fftw_destroy_plan(planZ);
    fftw_free(in);
}
