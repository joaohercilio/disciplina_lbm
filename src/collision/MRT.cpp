#include "collision/MRT.hpp"

ColParamMap MRT::prepareColParams(const ColParamMap& raw) const {
    auto it = raw.find("tau");
    if (it == raw.end()) throw std::invalid_argument("MRT requires 'tau'");
    double tau = it->second;
    if (tau <= 0.0) throw std::invalid_argument("tau must be > 0");

    ColParamMap out = raw;
    out["tau"] = tau;
    return out;
}

void MRT::computeCollision(std::vector<double>& f,
                           const LatticeModel& lattice,
                           const Geometry& geometry,
                           const ColParamMap& colParams,
                           const std::vector<double>& force)
{
    const double tau = colParams.find("tau") -> second;
   
    int numOfVel = lattice.getNumOfVel();
    int numOfDim = lattice.getNumOfDim();
    int numPoints = geometry.getNumOfPoints();

    std::vector<double> s = lattice.relaxationMatrix(tau, numOfVel);

    for (int i = 0; i < numPoints; ++i) {
        if (geometry.getNode(i) == NodeType::Fluid) {

            double* mapF = f.data() + i*numOfVel;

            double m[numOfVel];
            double meq[numOfVel];

            lattice.computeMoments(mapF, m);
            lattice.computeEquilibriumMoments(meq, m);

            for (int i = 1 + numOfDim; i < numOfVel; i++) {
                m[i] = m[i] - s[i] * (m[i] - meq[i]);
            }

            lattice.reconstructDistribution(mapF, m);

        }
    }
}