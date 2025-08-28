#include "Initialization.hpp"

void initializeFields(std::vector<double>& f, 
                        const LatticeModel& lattice, 
                        const Geometry& geometry,
                        const ColParamMap& colParams,
                        const std::vector<double>& u,
                        const std::vector<double>& rho)
{
    int numOfVel = lattice.getNumOfVel();
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (geometry.getNode(x, y, z) == NodeType::Fluid) {

                    int id = geometry.getIndex(x,y,z);

                    std::vector<double> feq(numOfVel, 0.0);

                    lattice.computeEquilibrium(feq.data(), rho[id], u[3*id], u[3*id+1], u[3*id+2]); 

                    for (int k = 0; k < numOfVel; ++k) {
                        f[id*numOfVel + k] = feq[k];
                    }

                } else {
                    
                    // Initialize Solid Nodes with zero velocity and density

                    std::vector<double> feq(numOfVel, 0.0);
                    lattice.computeEquilibrium(feq.data(), 1.0, 0.0, 0.0, 0.0);
                    int id = geometry.getIndex(x, y, z);
                    for (int k = 0; k < numOfVel; ++k) {
                        f[id*numOfVel + k] = feq[k];
                    }
                }
            }
        }
    }
}
