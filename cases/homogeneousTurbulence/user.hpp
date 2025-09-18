#include "Simulation.hpp"
#include <fftw3.h>
#include <complex>
#include <random>

using std::vector;

double energySpectrum(double k, double kmin = 3.0, double kmax = 8.0, double A = 1.4293E-4, int m = 4);

void initVelocity(vector<double>& u, const Geometry& geo, double kmin, double kmax, double A, double m);
void velocityComponent( const int N, double* const uvel, double* const vvel, double* const wvel );
void setVelocityField(const int N, double u0, vector<double>& uvel, vector<double>& vvel, vector<double>& wvel );


// ----- PROBLEM SETUP ----- //

                      
namespace user {
  
// Domain size
  const int N = 32;

// Geometry creation (default periodic boundary conditions)
  Geometry problemGeometry() {
    Geometry geo(N, N, N);
    return geo;
  }

// Lattice model
  using LatticeModel = D3Q19; 

// Collision Operator and parameters
  using CollisionModel = MRT;
  double tau = 0.506;

// Initial Conditions
  vector<double> initialDensity (const Geometry& geo)
  {
    vector <double> rho (geo.getNumOfPoints(), 0.0);

    return rho;
  }

  vector<double> initialVelocity(const Geometry& geo)
  {
    vector<double> uvel(geo.getNumOfPoints(), 0.0);
    vector<double> vvel(geo.getNumOfPoints(), 0.0);
    vector<double> wvel(geo.getNumOfPoints(), 0.0);

    double u0 = 0.001;
    setVelocityField(N, u0, uvel, vvel, wvel);

    vector<double> u;
    u.insert(u.end(), uvel.begin(), uvel.end());
    u.insert(u.end(), vvel.begin(), vvel.end());
    u.insert(u.end(), wvel.begin(), wvel.end());

    return u; 
  }

// External forces
  vector<double> externalForce()
  {
    return {0.0, 0.0, 0.0};
  }

// Number of iterations on the Pressure Field Initialization Routine
  int initializePressureIterations() 
  {
    return 0;
  }

// Total time steps
  int totalSteps () {
    return 100;
  }

// Time interval to write output files
  int writeInterval() {
    return 1; 
  }

// Output format (VTI, TSV, BOTH)
  OutputType outputType() {
    return OutputType::VTI;  
  }

  ColParamMap colParams() {
    return {{"tau", tau}};
  }
  
  inline vector<std::pair<std::string, double>> userLogParams() {
    return {
        {"N", N},
        {"Viscosity", (1.0/3.0) * (tau - 0.5)}
    };
  }  
}

// Energy spectrum formula
double energySpectrum(double k, double kmin, double kmax, double A, int m)
{
  if (k <= kmin || k >= kmax) return 0.0;
  return A * pow(k,m-2) * exp(-0.14 * k*k) / (2 * M_PI );
}

void velocityComponent( const int N,
                        double* const uvel,
                        double* const vvel,
                        double* const wvel )
{
  long N3 = N*N*N;
  long i,j,k;
  double kk = sqrt(3 * (N/2)*(N/2));
  int kkmax = (int) kk;
  fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));    
  fftw_complex *vhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));    
  fftw_complex *what = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));    

  fftw_complex *u = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  fftw_complex *v = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  fftw_complex *w = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));

  fftw_plan inv_trans_u = fftw_plan_dft_3d(N, N, N, uhat, u, 1, FFTW_MEASURE);
  fftw_plan inv_trans_v = fftw_plan_dft_3d(N, N, N, vhat, v, 1, FFTW_MEASURE);
  fftw_plan inv_trans_w = fftw_plan_dft_3d(N, N, N, what, w, 1, FFTW_MEASURE);

  /* Specifying the velocities in Fourier space */
  for (i=0; i<N3; i++){
    uhat[i][0] = uhat[i][1] = 0.0;
    vhat[i][0] = vhat[i][1] = 0.0;
    what[i][0] = what[i][1] = 0.0;
  }
   
  for (i = 1; i < N/2; i++){
    for (j = 0; j < N/2; j++){
      for (k = 0; k < N/2; k++){
        double kk   = sqrt(i*i + j*j + k*k);
        double th1  = 2*M_PI * ((double)rand())/((double)RAND_MAX);
        double th2  = 2*M_PI * ((double)rand())/((double)RAND_MAX);
        double phi1 = 2*M_PI * ((double)rand())/((double)RAND_MAX);
        
        double E = energySpectrum(kk);

        double alfa_real = sqrt(E/(4*M_PI*kk*kk))*cos(th1)*cos(phi1);
        double alfa_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th1)*cos(phi1);
        double beta_real = sqrt(E/(4*M_PI*kk*kk))*cos(th2)*sin(phi1);
        double beta_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th2)*sin(phi1);

        uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
        uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));

        vhat[k+N*(j+N*i)][0] = (beta_real*j*k-alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
        vhat[k+N*(j+N*i)][1] = (beta_imag*j*k-alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));

        what[k+N*(j+N*i)][0] = -(beta_real*sqrt(i*i+j*j))/kk;
        what[k+N*(j+N*i)][1] = -(beta_imag*sqrt(i*i+j*j))/kk;

        // Project to zero divergence
        /*
        std::complex<double> u(uhat[k+N*(j+N*i)][0], uhat[k+N*(j+N*i)][1]);
        std::complex<double> v(vhat[k+N*(j+N*i)][0], vhat[k+N*(j+N*i)][1]);
        std::complex<double> w(what[k+N*(j+N*i)][0], what[k+N*(j+N*i)][1]);

        double kx = (double) i;
        double ky = (double) j;
        double kz = (double) k;

        std::complex<double> dot = kx * u + ky * v + kz * w;        
        printf("Divergence: %1.6e\n", dot);

        uhat[k+N*(j+N*i)][0] -= dot.real() * kx / kk;
        uhat[k+N*(j+N*i)][1] -= dot.imag() * kx / kk;

        vhat[k+N*(j+N*i)][0] -= dot.real() * ky / kk;
        vhat[k+N*(j+N*i)][1] -= dot.imag() * ky / kk;

        what[k+N*(j+N*i)][0] -= dot.real() * kz / kk;
        what[k+N*(j+N*i)][1] -= dot.imag() * kz / kk;
        */
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

        double E = energySpectrum(kk);

        double alfa_real = sqrt(E/(4*M_PI*kk*kk))*cos(th1)*cos(phi1);
        double alfa_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th1)*cos(phi1);
        double beta_real = sqrt(E/(4*M_PI*kk*kk))*cos(th2)*sin(phi1);
        double beta_imag = sqrt(E/(4*M_PI*kk*kk))*sin(th2)*sin(phi1);

        uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
        uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));

        vhat[k+N*(j+N*i)][0] = (beta_real*j*k-alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
        vhat[k+N*(j+N*i)][1] = (beta_imag*j*k-alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));

        what[k+N*(j+N*i)][0] = -(beta_real*sqrt(i*i+j*j))/kk;
        what[k+N*(j+N*i)][1] = -(beta_imag*sqrt(i*i+j*j))/kk;
  
        // Project to zero divergence
        /*
        std::complex<double> u(uhat[k+N*(j+N*i)][0], uhat[k+N*(j+N*i)][1]);
        std::complex<double> v(vhat[k+N*(j+N*i)][0], vhat[k+N*(j+N*i)][1]);
        std::complex<double> w(what[k+N*(j+N*i)][0], what[k+N*(j+N*i)][1]);

        double kx = (double) i;
        double ky = (double) j;
        double kz = (double) k;

        std::complex<double> dot = kx * u + ky * v + kz * w;        

        uhat[k+N*(j+N*i)][0] -= dot.real() * kx / kk;
        uhat[k+N*(j+N*i)][1] -= dot.imag() * kx / kk;

        vhat[k+N*(j+N*i)][0] -= dot.real() * ky / kk;
        vhat[k+N*(j+N*i)][1] -= dot.imag() * ky / kk;

        what[k+N*(j+N*i)][0] -= dot.real() * kz / kk;
        what[k+N*(j+N*i)][1] -= dot.imag() * kz / kk;
        */
      }
    }
  }

  for (i = 0; i < 1; i++) {
    for (j = 0; j < 1; j++) {
      for (k = 0; k < N/2; k++) {
        uhat[k+N*(j+N*i)][0] = 0;
        uhat[k+N*(j+N*i)][1] = 0;
        vhat[k+N*(j+N*i)][0] = 0;
        vhat[k+N*(j+N*i)][1] = 0;
        what[k+N*(j+N*i)][0] = 0;
        what[k+N*(j+N*i)][1] = 0;
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

        vhat[k+N*(j+N*i)][0] =  vhat[(N-k)+N*((N-j)+N*(N-i))][0];
        vhat[k+N*(j+N*i)][1] = -vhat[(N-k)+N*((N-j)+N*(N-i))][1];

        what[k+N*(j+N*i)][0] =  what[(N-k)+N*((N-j)+N*(N-i))][0];
        what[k+N*(j+N*i)][1] = -what[(N-k)+N*((N-j)+N*(N-i))][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=0; k<1; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*(N-i))][1];

        vhat[k+N*(j+N*i)][0] =  vhat[k+N*((N-j)+N*(N-i))][0];
        vhat[k+N*(j+N*i)][1] = -vhat[k+N*((N-j)+N*(N-i))][1];

        what[k+N*(j+N*i)][0] =  what[k+N*((N-j)+N*(N-i))][0];
        what[k+N*(j+N*i)][1] = -what[k+N*((N-j)+N*(N-i))][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=0; j<1; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*(N-i))][1];

        vhat[k+N*(j+N*i)][0] =  vhat[(N-k)+N*(j+N*(N-i))][0];
        vhat[k+N*(j+N*i)][1] = -vhat[(N-k)+N*(j+N*(N-i))][1];

        what[k+N*(j+N*i)][0] =  what[(N-k)+N*(j+N*(N-i))][0];
        what[k+N*(j+N*i)][1] = -what[(N-k)+N*(j+N*(N-i))][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*((N-j)+N*i)][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*((N-j)+N*i)][1];

        vhat[k+N*(j+N*i)][0] =  vhat[(N-k)+N*((N-j)+N*i)][0];
        vhat[k+N*(j+N*i)][1] = -vhat[(N-k)+N*((N-j)+N*i)][1];

        what[k+N*(j+N*i)][0] =  what[(N-k)+N*((N-j)+N*i)][0];
        what[k+N*(j+N*i)][1] = -what[(N-k)+N*((N-j)+N*i)][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=0; j<1; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*i)][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*i)][1];
  
        vhat[k+N*(j+N*i)][0] =  vhat[(N-k)+N*(j+N*i)][0];
        vhat[k+N*(j+N*i)][1] = -vhat[(N-k)+N*(j+N*i)][1];

        what[k+N*(j+N*i)][0] =  what[(N-k)+N*(j+N*i)][0];
        what[k+N*(j+N*i)][1] = -what[(N-k)+N*(j+N*i)][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=0; k<1; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*i)][0];
        uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*i)][1];

        vhat[k+N*(j+N*i)][0] =  vhat[k+N*((N-j)+N*i)][0];
        vhat[k+N*(j+N*i)][1] = -vhat[k+N*((N-j)+N*i)][1];

        what[k+N*(j+N*i)][0] =  what[k+N*((N-j)+N*i)][0];
        what[k+N*(j+N*i)][1] = -what[k+N*((N-j)+N*i)][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=0; j<1; j++) {
      for (k=0; k<1; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[k+N*(j+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[k+N*(j+N*(N-i))][1];

        vhat[k+N*(j+N*i)][0] =  vhat[k+N*(j+N*(N-i))][0];
        vhat[k+N*(j+N*i)][1] = -vhat[k+N*(j+N*(N-i))][1];

        what[k+N*(j+N*i)][0] =  what[k+N*(j+N*(N-i))][0];
        what[k+N*(j+N*i)][1] = -what[k+N*(j+N*(N-i))][1];
      }
    }
  }
  /* Inverse Fourier transform */
  fftw_execute(inv_trans_u);
  fftw_execute(inv_trans_v);
  fftw_execute(inv_trans_w);

  fftw_free(uhat);
  fftw_free(vhat);
  fftw_free(what);

  fftw_destroy_plan(inv_trans_u);
  fftw_destroy_plan(inv_trans_v);
  fftw_destroy_plan(inv_trans_w);

  double u_imag_velocity = 0;
  double v_imag_velocity = 0; 
  double w_imag_velocity = 0;

  for (i = 0; i < N3; i++){
    double uu = u[i][1];
    u_imag_velocity += (uu*uu);

    double vv = v[i][1];
    v_imag_velocity += (vv*vv);

    double ww = w[i][1];
    w_imag_velocity += (ww*ww);
  }

  u_imag_velocity = sqrt(u_imag_velocity / ((double)N3));
  v_imag_velocity = sqrt(v_imag_velocity / ((double)N3));
  w_imag_velocity = sqrt(w_imag_velocity / ((double)N3));

  Logger logger;
  logger.logMessage("RMS of imaginary computed velocity U: " + logger.to_string(u_imag_velocity));
  logger.logMessage("RMS of imaginary computed velocity V: " + logger.to_string(v_imag_velocity));
  logger.logMessage("RMS of imaginary computed velocity W: " + logger.to_string(w_imag_velocity));

  for (i = 0; i < N3; i++){
    uvel[i] = u[i][0];
    vvel[i] = v[i][0];
    wvel[i] = w[i][0];
  }

  fftw_free(u);
  fftw_free(v);
  fftw_free(w);

  return;
}

void setVelocityField(const int N,
                      double u0,
                      vector<double>& uvel,
                      vector<double>& vvel,
                      vector<double>& wvel )
{
  Logger logger;
  long N3 = N*N*N;
  long i,j,k;
  double dx = 2*M_PI / ((double)N);

  velocityComponent( N, uvel.data(), vvel.data(), wvel.data() ); 
    
  double rms_velocity = 0;

  for (i = 0; i < N3; i++) {

    double uu, vv, ww;
    
    uu = uvel[i];
    vv = vvel[i];
    ww = wvel[i];
    
    rms_velocity += (uu*uu + vv*vv + ww*ww);
  }

  rms_velocity = sqrt(rms_velocity / (3*((double)N3)));
  
  /* scale the velocity components so that rms velocity matches u0 */
  double factor = u0 / rms_velocity;
  logger.logMessage("Scaling factor = " + logger.to_string(factor));

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
  logger.logMessage("RMS velocity (component-wise): " + logger.to_string(rms_velocity));

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