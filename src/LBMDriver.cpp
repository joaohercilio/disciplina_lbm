#include "user.hpp"

class LBMSimulation : public Simulation {
public:
    LBMSimulation() : Simulation (
            std::make_unique<user::LatticeModel>(), 
            std::make_unique<user::CollisionModel>(),
            user::problemGeometry()
        ) {}

    void run() override 
    {        
        const auto logParams = user::userLogParams();
        const auto outputType = user::outputType();
        const auto colParams = collision_ -> prepareColParams(user::colParams());
        const auto externalForce = user::externalForce();
        const auto initializePressureIterations = user::initializePressureIterations();
        const auto totalSteps = user::totalSteps();
        const auto writeInterval = user::writeInterval();
        const auto u0 = user::initialVelocity(geometry_);
        const auto rho0 = user::initialDensity(geometry_);

        Logger logger;
        Timer timer(logger);

        logger.logUserConfig("USER CONFIGURATION SUMMARY", logParams);

        auto& f_in  = f1_;
        auto& f_out = f2_;

        timer.start("Initialization");
        initializeFields(f_in, *lattice_, geometry_, colParams, u0, rho0);
        Neighbors neighbors(geometry_, *lattice_);
        //computeNonEquilibriumMoments(f_in, *lattice_, geometry_, colParams, u0);
        if(initializePressureIterations > 0)
            collision_->initializeDensityField(f_in, f_out, *lattice_, geometry_, colParams, u0, neighbors, initializePressureIterations, logger);
        timer.stop("Initialization");

        callOutput(f_in, *lattice_, geometry_, 0 ,outputType);

        logger.logMessage("Time steps\n");

        for (int t = 1; t <= totalSteps; ++t) {

            timer.start("Collision");
            collision_->computeCollision(f_in, *lattice_, geometry_, colParams, externalForce);
            timer.stop("Collision");

            timer.start("Streaming");
            performStreaming(f_in, f_out, *lattice_, geometry_, neighbors);
            timer.stop("Streaming");

            if (t % writeInterval == 0) {
                callOutput(f_out, *lattice_, geometry_, t ,outputType);
            }

            logger.logStep(t,totalSteps);

            std::swap(f_in, f_out);
        }

        logger.endLine();
        timer.report();

    }
};

std::unique_ptr<Simulation> Simulation::create() {
    return std::make_unique<LBMSimulation>();
}