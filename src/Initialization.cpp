#include "Initialization.hpp"

void initializeFields(std::vector<double>& f, 
                        const LatticeModel& lattice, 
                        const Geometry& geometry,
                        const ColParamMap& colParams,
                        const std::vector<double>& u,
                        const std::vector<double>& rho)
{
    int numOfVel = lattice.getNumOfVel();
    int numOfPoints = geometry.getNumOfPoints();

    #pragma omp parallel for
    for (int id = 0; id < numOfPoints; ++id) {
        if (geometry.getNode(id) == NodeType::Fluid) {

            std::vector<double> feq(numOfVel, 0.0);

            lattice.computeEquilibrium(feq.data(), rho[id], u[geometry.getVelocityIndex(id, 0)], 
                                                            u[geometry.getVelocityIndex(id, 1)], 
                                                            u[geometry.getVelocityIndex(id, 2)]); 

            for (int k = 0; k < numOfVel; ++k) {
                f[id*numOfVel + k] = feq[k];
            }

        } else {
            // Initialize Solid Nodes with zero velocity and density

            std::vector<double> feq(numOfVel, 0.0);
            lattice.computeEquilibrium(feq.data(), 1.0, 0.0, 0.0, 0.0);
            for (int k = 0; k < numOfVel; ++k) {
                f[id*numOfVel + k] = feq[k];
            }
        }
    }
}
