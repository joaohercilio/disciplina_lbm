#include "Streaming.hpp"
#include "Neighbors.hpp"

void performStreaming(std::vector<double>& f, std::vector<double>& fn, const LatticeModel& lattice, const Geometry& geometry, const Neighbors& neighbors)
{
    int numOfVel = lattice.getNumOfVel();
    
    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();
    
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();

    #pragma omp parallel for
    for (int id = 0; id < geometry.getNumOfPoints(); id++) {
        if (geometry.getNode(id) == NodeType::Fluid) {
            
            int x, y, z;
            
            geometry.getCoords(id, x, y, z);

            for (int k = 0; k < numOfVel; k++) {

                int idn = neighbors(id, k, numOfVel);
                
                // Halfway bounce and back
                if (geometry.getNode(idn) == NodeType::Solid) {
                    int k_opposite = (k == 0) ? 0 : (k % 2 == 1) ? k + 1 : k - 1;
                    fn[id * numOfVel + k_opposite] = f[id * numOfVel + k];
                } else {
                    fn[idn * numOfVel + k] = f[id * numOfVel + k];
                }
            }        
        }
    }
}