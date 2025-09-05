#include "Neighbors.hpp"
#include "Geometry.hpp"
#include "lattice/LatticeModel.hpp"

Neighbors::Neighbors(const Geometry& geometry, const LatticeModel& lattice) : neighbors_(geometry.getNumOfPoints() * lattice.getNumOfVel())
{
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();

    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();

    for (int id = 0; id < geometry.getNumOfPoints(); id++){

        int x, y, z;
        geometry.getCoords(id, x, y, z);

        for (int k = 0; k < lattice.getNumOfVel(); ++k) {
            int xn = (x + cx[k] + nx) % nx;
            int yn = (y + cy[k] + ny) % ny;
            int zn = (z + cz[k] + nz) % nz;
            neighbors_[id * lattice.getNumOfVel() + k] = geometry.getIndex(xn, yn, zn);
        }
    }
}
