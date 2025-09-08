#include "Simulation.hpp"
#include <fftw3.h>
#include <complex>
#include <random>
#include <stdio.h>
using namespace std;

// ----- PROBLEM SETUP ----- //

void initVelocity(std::vector<double>& u, const Geometry& geo, double kmin, double kmax, double A, double m);
double energySpectrum_(double k, double kmin, double kmax, double A, int m);
void uHat(double kx, double ky, double kz, complex<double> vHat[3], double kmin, double kmax, double A, double m);
void velocityComponent( const int N,
                        const int dir,
                        const double kp,
                        const double u0,
                        double* const uvel );
void setVelocityField(const int N,
                      std::vector<double>& uvel,
                      std::vector<double>& vvel,
                      std::vector<double>& wvel );
namespace user {
    
    const int N = 32;

    using LatticeModel = D3Q19; 

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
        double kmin = 3;
        double kmax = 8;
        double A = 1.4293e-4;
        double m = 4;

        long N3 = N*N*N;
        int i,j,k;
        double dx = 2*M_PI / ((double)N);
        std::vector<double> uvel(geo.getNumOfPoints(), 0.0);
        std::vector<double> vvel(geo.getNumOfPoints(), 0.0);
        std::vector<double> wvel(geo.getNumOfPoints(), 0.0);        
        setVelocityField(N, uvel, vvel, wvel);

        std::vector<double> u;
        u.insert(u.end(), uvel.begin(), uvel.end());
        u.insert(u.end(), vvel.begin(), vvel.end());
        u.insert(u.end(), wvel.begin(), wvel.end());

        /*
        initVelocity(u, geo, kmin, kmax, A, m);
        
        double c = 1.0/5.0;
        for (int i = 0; i < 3 * geo.getNumOfPoints(); i++) { u[i] *= c; }
        */
        return u; 
    }

    std::vector<double> initialDensity (const Geometry& geo) {
        std::vector <double> rho (geo.getNumOfPoints(), 1.0);
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
        Geometry geo(N, N, N);
        return geo;
    }

    int totalSteps () {
        return 10;
    }

    OutputType outputType() {
        return OutputType::VTI;  
    }
    int writeInterval() {
        return 1; 
    }
}

double energySpectrum_(double k, double kmin, double kmax, double A, int m)
{
    if (k <= kmin || k >= kmax) return 0.0;
    return A * pow(k,m-2) * exp(-0.14 * k*k) / (2 * M_PI );
}

void velocityComponent( const int N,
                        const int dir,
                        const double kp,
                        const double u0,
                        double* const uvel )
{
    long N3 = N*N*N;
    long i,j,k;
    double kk = sqrt(3 * (N/2)*(N/2));
    int kkmax = (int) kk;
    fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));    
    fftw_complex *u = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));   
    fftw_plan inv_trans_u;
    inv_trans_u = fftw_plan_dft_3d(N, N, N, uhat, u, 1, FFTW_MEASURE);
  /* Specifying the velocities in Fourier space */
  for (i=0; i<N3; i++) uhat[i][0] = uhat[i][1] = 0.0;
    for (i = 1; i < N/2; i++){
        for (j = 0; j < N/2; j++){
            for (k = 0; k < N/2; k++){
                double kk   = sqrt(i*i + j*j + k*k);
                double th1  = 2*M_PI * ((double)rand())/((double)RAND_MAX);
                double th2  = 2*M_PI * ((double)rand())/((double)RAND_MAX);
                double phi1 = 2*M_PI * ((double)rand())/((double)RAND_MAX);
                double E = 16.0 * sqrt(2.0/M_PI) * (u0*u0/kp) * pow(kk/kp, 4.0) 
                        * exp(-2.0*(kk/kp)*(kk/kp));
                double alfa_real = sqrt(E/(4*M_PI*kk*kk))*cos(th1)*cos(phi1);
                double alfa_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th1)*cos(phi1);
                double beta_real = sqrt(E/(4*M_PI*kk*kk))*cos(th2)*sin(phi1);
                double beta_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th2)*sin(phi1);
        if (dir == 0) {
                  uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
                  uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));
        } else if (dir == 1) {
                  uhat[k+N*(j+N*i)][0] = (beta_real*j*k-alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
                  uhat[k+N*(j+N*i)][1] = (beta_imag*j*k-alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));
        } else {
                  uhat[k+N*(j+N*i)][0] = -(beta_real*sqrt(i*i+j*j))/kk;
                  uhat[k+N*(j+N*i)][1] = -(beta_imag*sqrt(i*i+j*j))/kk;
        }
            }
        }
    }
    for (i = 0; i < 1; i++){
        for (k = 0; k < N/2; k++){
            for (j = 1; j < N/2; j++){
                double kk   = sqrt(i*i + j*j + k*k);
                double th1  = 2*M_PI * ((double)rand())/((double)RAND_MAX);
                double th2  = 2*M_PI * ((double)rand())/((double)RAND_MAX);
                double phi1 = 2*M_PI * ((double)rand())/((double)RAND_MAX);
                double E = 16.0 * sqrt(2.0/M_PI) * (u0*u0/kp) * pow(kk/kp, 4.0) 
                        * exp(-2.0*(kk/kp)*(kk/kp));
                double alfa_real = sqrt(E/(4*M_PI*kk*kk))*cos(th1)*cos(phi1);
                double alfa_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th1)*cos(phi1);
                double beta_real = sqrt(E/(4*M_PI*kk*kk))*cos(th2)*sin(phi1);
                double beta_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th2)*sin(phi1);
        if (dir == 0) {
                  uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
                  uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));
        } else if (dir == 1) {
                  uhat[k+N*(j+N*i)][0] = (beta_real*j*k-alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
                  uhat[k+N*(j+N*i)][1] = (beta_imag*j*k-alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));
        } else {
                  uhat[k+N*(j+N*i)][0] = -(beta_real*sqrt(i*i+j*j))/kk;
                  uhat[k+N*(j+N*i)][1] = -(beta_imag*sqrt(i*i+j*j))/kk;
        }
            }
        }
    }
    for (i = 0; i < 1; i++){
        for (j = 0; j < 1; j++){
            for (k = 0; k < N/2; k++){
                uhat[k+N*(j+N*i)][0] = 0;
                uhat[k+N*(j+N*i)][1] = 0;
            }
        }
    }
  /* The following is necessary to ensure that the inverse Fourier
     transform yields a real velocity field with zero imaginary
     components */
  for (i=N/2+1; i<N; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=N/2+1; k<N; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*((N-j)+N*(N-i))][0];
                uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*((N-j)+N*(N-i))][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=0; k<1; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*(N-i))][0];
                uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*(N-i))][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=0; j<1; j++) {
      for (k=N/2+1; k<N; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*(N-i))][0];
                uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*(N-i))][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=N/2+1; k<N; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*((N-j)+N*i)][0];
                uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*((N-j)+N*i)][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=0; j<1; j++) {
      for (k=N/2+1; k<N; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*i)][0];
                uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*i)][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=0; k<1; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*i)][0];
                uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*i)][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=0; j<1; j++) {
      for (k=0; k<1; k++) {
                uhat[k+N*(j+N*i)][0] =  uhat[k+N*(j+N*(N-i))][0];
                uhat[k+N*(j+N*i)][1] = -uhat[k+N*(j+N*(N-i))][1];
      }
    }
  }
  /* Inverse Fourier transform */
    fftw_execute(inv_trans_u);
    fftw_free(uhat);
    fftw_destroy_plan(inv_trans_u);
  double imag_velocity = 0;
    for (i = 0; i < N3; i++){
        double uu = u[i][1];
        imag_velocity += (uu*uu);
    }
    imag_velocity = sqrt(imag_velocity / ((double)N3));
    printf("RMS of imaginary component of computed velocity: %1.6e\n",imag_velocity);
    for (i = 0; i < N3; i++){
    uvel[i] = u[i][0];
    }
  fftw_free(u);
  return;
}

void setVelocityField(const int N,
                      std::vector<double>& uvel,
                      std::vector<double>& vvel,
                      std::vector<double>& wvel )
{
    double kp = 4.0;
    double u0 = 0.001;
    long N3 = N*N*N;
    long i,j,k;
    double dx = 2*M_PI / ((double)N);
  velocityComponent( N, 0, kp, u0, uvel.data() ); 
  velocityComponent( N, 1, kp, u0, vvel.data() ); 
  velocityComponent( N, 2, kp, u0, wvel.data() ); 
    double rms_velocity = 0;
    for (i = 0; i < N3; i++){
    double uu, vv, ww;
        uu = uvel[i];
        vv = vvel[i];
        ww = wvel[i];
        rms_velocity += (uu*uu + vv*vv + ww*ww);
    }
    rms_velocity = sqrt(rms_velocity / (3*((double)N3)));
  
  /* scale the velocity components so that rms velocity matches u0 */
  double factor = u0 / rms_velocity;
  printf("Scaling factor = %1.16E\n",factor);
    for (i = 0; i < N3; i++){
        uvel[i] *= factor;
        vvel[i] *= factor;
        wvel[i] *= factor;
    }
    rms_velocity = 0;
    for (i = 0; i < N3; i++){
    double uu, vv, ww;
        uu = uvel[i];
        vv = vvel[i];
        ww = wvel[i];
        rms_velocity += (uu*uu + vv*vv + ww*ww);
    }
    rms_velocity = sqrt(rms_velocity / (3*((double)N3)));
    printf("RMS velocity (component-wise): %1.16E\n",rms_velocity);
  /* calculate the divergence of velocity */
  double DivergenceNorm = 0;
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
        double u1, u2, v1, v2, w1, w2;
        u1 = (i==0   ? uvel[k+N*(j+N*(N-1))] : uvel[k+N*(j+N*(i-1))] );
        u2 = (i==N-1 ? uvel[k+N*(j+N*(0  ))] : uvel[k+N*(j+N*(i+1))] );
        v1 = (j==0   ? vvel[k+N*((N-1)+N*i)] : vvel[k+N*((j-1)+N*i)] );
        v2 = (j==N-1 ? vvel[k+N*((0  )+N*i)] : vvel[k+N*((j+1)+N*i)] );
        w1 = (k==0   ? wvel[(N-1)+N*(j+N*i)] : wvel[(k-1)+N*(j+N*i)] );
        w2 = (k==N-1 ? wvel[(0  )+N*(j+N*i)] : wvel[(k+1)+N*(j+N*i)] );
        double Divergence = ( (u2-u1) + (v2-v1) + (w2-w1) ) / (2.0*dx);
        DivergenceNorm += (Divergence*Divergence);
      }
    }
  }
  DivergenceNorm = sqrt(DivergenceNorm / ((double)N3));
  printf("Velocity divergence: %1.16E\n",DivergenceNorm);
  /* calculate the Taylor microscales */
  double TaylorMicroscale[3];
  double Numerator[3] = {0,0,0};
  double Denominator[3] = {0,0,0};
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
        double u1, u2, uu, v1, v2, vv, w1, w2, ww;
        u1 = (i==0   ? uvel[k+N*(j+N*(N-1))] : uvel[k+N*(j+N*(i-1))] );
        u2 = (i==N-1 ? uvel[k+N*(j+N*(0  ))] : uvel[k+N*(j+N*(i+1))] );
        v1 = (j==0   ? vvel[k+N*((N-1)+N*i)] : vvel[k+N*((j-1)+N*i)] );
        v2 = (j==N-1 ? vvel[k+N*((0  )+N*i)] : vvel[k+N*((j+1)+N*i)] );
        w1 = (k==0   ? wvel[(N-1)+N*(j+N*i)] : wvel[(k-1)+N*(j+N*i)] );
        w2 = (k==N-1 ? wvel[(0  )+N*(j+N*i)] : wvel[(k+1)+N*(j+N*i)] );
        uu  = uvel[k+N*(j+N*i)];
        vv  = vvel[k+N*(j+N*i)];
        ww  = wvel[k+N*(j+N*i)];
        double du, dv, dw;
        du = (u2 - u1) / (2.0*dx);
        dv = (v2 - v1) / (2.0*dx);
        dw = (w2 - w1) / (2.0*dx);
        Numerator[0] += (uu*uu);
        Numerator[1] += (vv*vv);
        Numerator[2] += (ww*ww);
        Denominator[0] += (du*du);
        Denominator[1] += (dv*dv);
        Denominator[2] += (dw*dw);
      }
    }
  }
  Numerator[0] /= (N*N*N); Denominator[0] /= (N*N*N);
  Numerator[1] /= (N*N*N); Denominator[1] /= (N*N*N);
  Numerator[2] /= (N*N*N); Denominator[2] /= (N*N*N);
  TaylorMicroscale[0] = sqrt(Numerator[0]/Denominator[0]);
  TaylorMicroscale[1] = sqrt(Numerator[1]/Denominator[1]);
  TaylorMicroscale[2] = sqrt(Numerator[2]/Denominator[2]);
  printf("Taylor microscales: %1.16E, %1.16E, %1.16E\n",
         TaylorMicroscale[0],TaylorMicroscale[1],TaylorMicroscale[2]);
  return;
}