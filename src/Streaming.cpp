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

    for (int i = 0; i < geometry.getNumOfPoints(); i++) {
        if (geometry.getNode(i) == NodeType::Fluid) {
            
            int x, y, z;
            
            geometry.getCoords(i, x, y, z);

            for (int k = 0; k < numOfVel; k++) {
                
                int xn = (x + cx[k] + nx) % nx;
                int yn = (y + cy[k] + ny) % ny;
                int zn = (z + cz[k] + nz) % nz;

                int idxn = geometry.getIndex(xn, yn, zn);
                int idx  = geometry.getIndex(x, y, z);

                fn[idxn * numOfVel + k] = f[idx * numOfVel + k];
            }        
        }
    }
}