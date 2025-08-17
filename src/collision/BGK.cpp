#include "collision/BGK.hpp"

void BGK::computeCollision(std::vector<double>& f, const CollisionParameters& colParams, const LatticeModel& lattice, const Geometry& geometry) const {

    int numOfVel = lattice.getNumOfVel();
    int numPoints = geometry.getNumOfPoints();

    for (int i = 0; i < numPoints; ++i) {
        if (geometry.getNode(i) == NodeType::Fluid) {

            double* mapF = f.data() + i*numOfVel;

            double rho, vx, vy, vz;
            
            lattice.computeFields( mapF, rho, vx, vy, vz );

            double feq[numOfVel];

            lattice.computeEquilibrium(feq, rho, vx, vy, vz);

            for (int k = 0; k < numOfVel; k++) 
            {
                mapF[k] = colParams.alphaEq * mapF[k] + colParams.alphaNonEq * feq[k];
            }
        }
    }
}