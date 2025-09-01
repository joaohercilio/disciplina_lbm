#include "collision/BGK.hpp"

ColParamMap BGK::prepareColParams(const ColParamMap& raw) const {
    auto it = raw.find("tau");
    if (it == raw.end()) throw std::invalid_argument("BGK requires 'tau'");
    double tau = it->second;
    if (tau <= 0.0) throw std::invalid_argument("tau must be > 0");

    const double alphaEq = 1.0 / tau;
    const double alphaNonEq = 1.0 - alphaEq;
    
    ColParamMap out = raw;
    out["alphaEq"] = alphaEq;
    out["alphaNonEq"] = alphaNonEq;
    return out;
}

void BGK::computeCollision(std::vector<double>& f,
                           const LatticeModel& lattice,
                           const Geometry& geometry,
                           const ColParamMap& colParams,
                           const std::vector<double>& force)
{
    const double alphaEq    = colParams.find("alphaEq") -> second;
    const double alphaNonEq = colParams.find("alphaNonEq") -> second;
    const double tau        = colParams.find("tau") -> second;

    int numOfVel = lattice.getNumOfVel();
    int numPoints = geometry.getNumOfPoints();

    #pragma omp parallel for
    for (int i = 0; i < numPoints; ++i) {
        if (geometry.getNode(i) == NodeType::Fluid) {

            double* mapF = f.data() + i*numOfVel;

            double rho, vx, vy, vz;
            
            lattice.computeFields( mapF, rho, vx, vy, vz );

            double feq[numOfVel];

            vx += tau*force[0] ;
            vy += tau*force[1] ;
            vz += tau*force[2] ;

            lattice.computeEquilibrium( feq, rho, vx, vy, vz );

            for (int k = 0; k < numOfVel; k++) 
            {
                mapF[k] = alphaNonEq * mapF[k] + alphaEq * feq[k];
            }
        }
    }
}

void BGK::initializeDensityField(std::vector<double>& f,
                                std::vector<double>& fn,
                                const LatticeModel& lattice,
                                const Geometry& geometry,
                                const ColParamMap& colParams,
                                const std::vector<double>& u,
                                const std::vector<double>& force,
                                const int numberOfIterations)
{
}