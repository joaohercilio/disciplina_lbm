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
    int numOfVel = lattice.getNumOfVel();
    int numOfDim = lattice.getNumOfDim();
    int numPoints = geometry.getNumOfPoints();

    const double tau = colParams.find("tau") -> second;
    std::vector<double> s = lattice.relaxationMatrix(tau, numOfVel);

    #pragma omp parallel for
    for (int id = 0; id < numPoints; ++id) {
        if (geometry.getNode(id) == NodeType::Fluid) {

            double* mapF = f.data() + id*numOfVel;

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

void MRT::initializeDensityField(std::vector<double>& f,
                                std::vector<double>& fn,
                                const LatticeModel& lattice,
                                const Geometry& geometry,
                                const ColParamMap& colParams,
                                const std::vector<double>& u,
                                const std::vector<double>& force,
                                const int numberOfIterations,
                                Logger& logger)
{
    logger.logMessage("\nPressure field initialization steps\n");

    int numOfPoints = geometry.getNumOfPoints();
    int nx = geometry.nx();
    int ny = geometry.ny();
    int nz = geometry.nz();
    int numOfVel = lattice.getNumOfVel();
    int numOfDim = lattice.getNumOfDim();

    const double tau = colParams.find("tau") -> second;
    std::vector<double> s = lattice.relaxationMatrix(tau, numOfVel);

    std::vector<double> drhoOld(numOfPoints, 1.0);

    int Istep = 0;

    double varrho;
    while (Istep < numberOfIterations) {
        varrho = 0.0;

        #pragma omp parallel for reduction(+:varrho)
        for (int id = 0; id < numOfPoints; id++)
        {
            double* mapF = f.data() + id*numOfVel;

            double m[numOfVel];
            double meq[numOfVel];

            lattice.computeMoments(mapF, m);
            lattice.computeEquilibriumMoments(meq, m);

            for (int i = 1 + numOfDim; i < numOfVel; i++) 
            {
                m[i] = m[i] - s[i] * (m[i] - meq[i]);
            }

            lattice.reconstructDistribution(mapF, m);

            varrho += fabs(drhoOld[id] - m[0]);
            drhoOld[id] = m[0];
        }

        performStreaming(f, fn, lattice, geometry);
        std::swap(f, fn);
        
        #pragma omp parallel for
        for (int id = 0; id < numOfPoints; id++)
        {
            double* mapF = f.data() + id*numOfVel;
            double m[numOfVel];
            double meq[numOfVel];

            lattice.computeMoments(mapF, m);
            
            double rho = m[0];
            double ux = u[geometry.getVelocityIndex(id, 0)];
            double uy = u[geometry.getVelocityIndex(id, 1)];

            meq[0] = rho;
            meq[1] = rho * ux;	
            meq[2] = rho * uy;
            
            lattice.reconstructDistribution(mapF, m);
        }

        Istep++;

        logger.logStep(Istep,numberOfIterations);
    } 

    logger.endLine();
    logger.logMessage("Final difference: " + std::to_string(varrho) + "\n\n");

}