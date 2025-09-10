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
    for (int id = 0; id < numPoints; ++id) {
        if (geometry.getNode(id) == NodeType::Fluid) {

            double* mapF = f.data() + id*numOfVel;

            double drho, vx, vy, vz;
            
            lattice.computeFields( mapF, drho, vx, vy, vz );

            std::vector<double> feq(numOfVel);

            vx += tau*force[0] ;
            vy += tau*force[1] ;
            vz += tau*force[2] ;

            lattice.computeEquilibrium( feq.data(), drho, vx, vy, vz );

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
                                const Neighbors& neighbors,
                                const int numberOfIterations,
                                Logger& logger)
{
    logger.logMessage("\nPressure field initialization steps\n");

    int numOfPoints = geometry.getNumOfPoints();
    int numOfVel = lattice.getNumOfVel();
    int numOfDim = lattice.getNumOfDim();

    const double alphaEq    = colParams.find("alphaEq") -> second;
    const double alphaNonEq = colParams.find("alphaNonEq") -> second;

    std::vector<double> drhoOld(numOfPoints, 0.0);

    int Istep = 0;

    double varrho;
    while (Istep < numberOfIterations) {
        varrho = 0.0;

        #pragma omp parallel for reduction(+:varrho)
        for (int id = 0; id < numOfPoints; id++)
        {
            double* mapF = f.data() + id*numOfVel;

            double rho, vx, vy, vz;
            
            lattice.computeFields( mapF, rho, vx, vy, vz );

            std::vector<double> feq(numOfVel);

            double ux = u[geometry.getVelocityIndex(id, 0)];
            double uy = u[geometry.getVelocityIndex(id, 1)];
            double uz = u[geometry.getVelocityIndex(id, 2)];

            lattice.computeEquilibrium( feq.data(), rho, ux, uy, uz );

            for (int k = 0; k < numOfVel; k++) 
            {
                mapF[k] = alphaNonEq * mapF[k] + alphaEq * feq[k];
            }

            varrho += fabs(drhoOld[id] - rho);
            drhoOld[id] = rho;
        }

        performStreaming(f, fn, lattice, geometry, neighbors);
        std::swap(f, fn);
        
        Istep++;

        logger.logStep(Istep,numberOfIterations);
    } 

    #pragma omp parallel for
    for (int id = 0; id < numOfPoints; id++)
    {
        double* mapF = f.data() + id*numOfVel;
        std::vector<double> m(numOfVel);
        double ux = u[geometry.getVelocityIndex(id, 0)];
        double uy = u[geometry.getVelocityIndex(id, 1)];
        double uz = u[geometry.getVelocityIndex(id, 2)];
        lattice.computeMoments(mapF, m.data());
        m[1] = (1.0 + m[0]) * ux;
        m[2] = (1.0 + m[0]) * uy;
        if (numOfDim == 3) m[3] = (1.0 + m[0]) * uz;
        lattice.reconstructDistribution(mapF, m.data());
    }

    logger.endLine();
    logger.logMessage("Final difference: " + logger.to_string(varrho) + "\n\n");
}