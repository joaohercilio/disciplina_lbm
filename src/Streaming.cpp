#include "Streaming.hpp"

void Streaming::performStreaming(std::vector<double>& f, std::vector<double>& fn, const LatticeModel& lattice, const Geometry& geometry) const
{
    int numOfVel = lattice.getNumOfVel();
    
    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();
    
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();

    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                //if (geometry.getNode(x,y,z) == NodeType::Fluid) {
                    
                    for (int k = 0; k < numOfVel; k++) {
                        int xn = (x + cx[k] + nx) % nx;
                        int yn = (y + cy[k] + ny) % ny;
                        int zn = (z + cz[k] + nz) % nz;

                        int idx  = x  + y  * nx + z  * nx * ny;
                        int idxn = xn + yn * nx + zn * nx * ny;

                        fn[idxn * numOfVel + k] = f[idx * numOfVel + k];
                    }
              
                //}
            }
        }
    }
}