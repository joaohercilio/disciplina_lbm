#include "Streaming.hpp"

void performStreaming(std::vector<double>& f, std::vector<double>& fn, const LatticeModel& lattice, const Geometry& geometry)
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
                
                int xn = (x + cx[k] + nx) % nx;
                int yn = (y + cy[k] + ny) % ny;
                int zn = (z + cz[k] + nz) % nz;

                int idn = geometry.getIndex(xn, yn, zn);

                fn[idn * numOfVel + k] = f[id * numOfVel + k];
            }        
        }
    }
}