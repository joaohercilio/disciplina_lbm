#include "boundary/BounceAndBack.hpp"
#include <iostream>

void BounceAndBack::applyBoundary(std::vector<double>& f, 
                                    const LatticeModel& lattice, 
                                    const Geometry& geometry,
                                    int x, int y, int z)
{
    BoundaryType type = geometry.getBoundaryType(x, y, z);

    if (type == BoundaryType::BounceAndBackSouth) 
    {
        applySouthBoundary(f, lattice, geometry, x, y, z);
    }
    else if (type == BoundaryType::BounceAndBackNorth) 
    {
        applyNorthBoundary(f, lattice, geometry, x, y, z);
    }
}                                       

void BounceAndBack::applySouthBoundary(std::vector<double>& f, 
                                       const LatticeModel& lattice, 
                                       const Geometry& geometry,
                                       int x, int y, int z)
{
    int idx = geometry.getIndex(x, y, z);

    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();

    int numOfVel = lattice.getNumOfVel();

    double* mapF = f.data() + idx*numOfVel;

    for (int k = 0; k < numOfVel; k++) {
        if (cy[k] == 1) {
            mapF[k] = mapF[k+1];
        }
    }
}

void BounceAndBack::applyNorthBoundary(std::vector<double>& f, 
                                       const LatticeModel& lattice, 
                                       const Geometry& geometry,
                                       int x, int y, int z)
{
    int idx = geometry.getIndex(x, y, z);

    std::vector<int> cx = lattice.getCx();
    std::vector<int> cy = lattice.getCy();
    std::vector<int> cz = lattice.getCz();

    int numOfVel = lattice.getNumOfVel();

    double* mapF = f.data() + idx*numOfVel;

    for (int k = 0; k < numOfVel; k++) {
        if (cy[k] == -1) {
            mapF[k] = mapF[k-1];
        }
    }
}