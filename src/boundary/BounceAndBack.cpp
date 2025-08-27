#include "boundary/HalfwayBounceAndBack.hpp"

void HalfwayBounceAndBack::applyBoundary(std::vector<double>& f, 
                                    const LatticeModel& lattice, 
                                    const Geometry& geometry)
{
    for (int id = 0; id < geometry.getNumOfPoints(); ++id) {

        auto type = geometry.getBoundaryType(id);

        if (type == BoundaryType::HalfwayBounceAndBackSouth) 
        {
            applySouthBoundary(f, lattice, geometry, id);
        }

        else if (type == BoundaryType::HalfwayBounceAndBackNorth) 
        {
            applyNorthBoundary(f, lattice, geometry, id);
        }
    }
}                                       

void HalfwayBounceAndBack::applySouthBoundary(std::vector<double>& f, 
                                       const LatticeModel& lattice, 
                                       const Geometry& geometry,
                                       int id)
{    
    int x,y,z;
    geometry.getCoords(id, x, y, z);
    
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();

    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();

    int numOfVel = lattice.getNumOfVel();

    double* mapF = f.data() + id*numOfVel;

    for (int k = 0; k < numOfVel; k++) {
        if (cy[k] == -1) {

            int xn = (x + cx[k-1] + nx) % nx;
            int yn = (y + cy[k-1] + ny) % ny;
            int zn = (z + cz[k-1] + nz) % nz;

            int idn = geometry.getIndex(xn, yn, zn);
            double* mapFn = f.data() + idn*numOfVel;
            mapFn[k-1] = mapF[k];
            mapF[k] = 0.0;
        }
    }
}

void HalfwayBounceAndBack::applyNorthBoundary(std::vector<double>& f, 
                                       const LatticeModel& lattice, 
                                       const Geometry& geometry,
                                       int id)
{
    int x,y,z;
    geometry.getCoords(id, x, y, z);
    
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();

    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();

    int numOfVel = lattice.getNumOfVel();

    double* mapF = f.data() + id*numOfVel;

    for (int k = 0; k < numOfVel; k++) {
        if (cy[k] == 1) {

            int xn = (x + cx[k+1] + nx) % nx;
            int yn = (y + cy[k+1] + ny) % ny;
            int zn = (z + cz[k+1] + nz) % nz;
            int idn = geometry.getIndex(xn, yn, zn);

            double* mapFn = f.data() + idn*numOfVel;

            mapFn[k+1] = mapF[k];
            mapF[k] = 0.0;
        }
    }
}