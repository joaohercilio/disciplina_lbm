#include "Initialization.hpp"
#include "Logger.hpp"

void initializeFields(std::vector<double>& f, 
                        const LatticeModel& lattice, 
                        const Geometry& geometry,
                        const ColParamMap& colParams,
                        const std::vector<double>& u0,
                        const std::vector<double>& drho0)
{
    int numOfVel = lattice.getNumOfVel();
    int numOfPoints = geometry.getNumOfPoints();

    #pragma omp parallel for
    for (int id = 0; id < numOfPoints; ++id) {
        if (geometry.getNode(id) == NodeType::Fluid) {
            
            double* mapF = f.data() + id*numOfVel;

            std::vector<double> feq(numOfVel, 0.0);

            lattice.computeEquilibrium(feq.data(), drho0[id],u0[geometry.getVelocityIndex(id, 0)], 
                                                             u0[geometry.getVelocityIndex(id, 1)], 
                                                             u0[geometry.getVelocityIndex(id, 2)]); 
            
            for (int k = 0; k < numOfVel; ++k) {
                mapF[k] = feq[k];
            }

        } else {
            // Initialize Solid Nodes with zero velocity and density
            double* mapF = f.data() + id*numOfVel;
            std::vector<double> feq(numOfVel, 0.0);
            lattice.computeEquilibrium(feq.data(), 0.0, 0.0, 0.0, 0.0);
            for (int k = 0; k < numOfVel; ++k) {
                mapF[k] = feq[k];
            }
        }
    }
}

void computeNonEquilibriumMoments(std::vector<double>& f, 
                const LatticeModel& lattice, 
                const Geometry& geometry,
                const ColParamMap& colParams,
                const std::vector<double>& u0) 
{
    double tau = colParams.find("tau")->second;
    double Snu = 1.0 / tau;
    int numOfPoints = geometry.getNumOfPoints();
    int numOfVel = lattice.getNumOfVel();
    int nx = geometry.nx();
    int ny = geometry.ny();
    
    std::vector<double> jx(numOfPoints, 0.0);
    std::vector<double> jy(numOfPoints, 0.0);
    std::vector<double> dxjx(numOfPoints, 0.0);
    std::vector<double> dxjy(numOfPoints, 0.0);
    std::vector<double> pxx(numOfPoints, 0.0);
    std::vector<double> pxy(numOfPoints, 0.0);
    std::vector<double> e(numOfPoints, 0.0);

    for(int id = 0; id < numOfPoints; id++) {
        jx[id] = u0[geometry.getVelocityIndex(id, 0)];
        jy[id] = u0[geometry.getVelocityIndex(id, 1)];
    }

    int complex_size = ny * (nx/2 + 1);
    fftw_complex *JX = (fftw_complex*) fftw_malloc(complex_size * sizeof(fftw_complex));   
    fftw_complex *JY = (fftw_complex*) fftw_malloc(complex_size * sizeof(fftw_complex));   

    fftw_plan planX_forward = fftw_plan_dft_r2c_2d(nx, ny, jx.data(), JX, FFTW_ESTIMATE);
    fftw_plan planY_forward = fftw_plan_dft_r2c_2d(nx, ny, jy.data(), JY, FFTW_ESTIMATE);

    fftw_execute(planX_forward);
    fftw_execute(planY_forward);

    // ik (Re + iIm) = -kIm + ikRe 
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx/2 + 1; x++) {
            int idx = y * (nx/2 + 1) + x;
            
            double kx = (x < nx/2) ? x : x - nx;
            double ky = (y < ny/2) ? y : y - ny;
            
            double re_x = JX[idx][0];
            double im_x = JX[idx][1];

            JX[idx][0] = -kx * im_x;  
            JX[idx][1] =  kx * re_x;  

            double re_y = JY[idx][0];
            double im_y = JY[idx][1];

            JY[idx][0] = -ky * im_y;  
            JY[idx][1] =  ky * re_y; 
        }
    }

    fftw_plan planX_backward = fftw_plan_dft_c2r_2d(nx, ny, JX, dxjx.data(), FFTW_ESTIMATE);
    fftw_plan planY_backward = fftw_plan_dft_c2r_2d(nx, ny, JY, dxjy.data(), FFTW_ESTIMATE);

    fftw_execute(planX_backward);
    fftw_execute(planY_backward);

    double norm = 1.0 / (numOfPoints);
    for(int id = 0; id < numOfPoints; id++) {
        dxjx[id] *= norm;
        dxjy[id] *= norm;
    }

    for(int id = 0; id < numOfPoints; id++) {
        pxx[id] = -2.0 / (3.0 * Snu) * (dxjx[id] - dxjy[id]);
        pxy[id] = -1.0 / (3.0 * Snu) * (dxjx[id] + dxjy[id]);
        e[id] = -1.0 * (dxjx[id] + dxjy[id]); // Se = 1.0
    }

    fftw_destroy_plan(planX_forward);
    fftw_destroy_plan(planY_forward);
    fftw_destroy_plan(planX_backward);
    fftw_destroy_plan(planY_backward);
    fftw_free(JX);
    fftw_free(JY);

    for (int id = 0; id < numOfPoints; id++) {
        double* mapF = f.data() + id * numOfVel;

        std::vector<double> m(numOfVel);

        lattice.computeMoments(mapF, m.data());

        m[3] += e[id];
        m[7] += pxx[id];  
        m[8] += pxy[id];

        lattice.reconstructDistribution(mapF, m.data());
    }
}