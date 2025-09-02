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
        auto& f_in  = f1_;
        auto& f_out = f2_;

        Logger logger;
        Timer timer(logger);

        auto logParams = user::userLogParams();
        logger.logUserConfig("USER CONFIGURATION SUMMARY", logParams);

        auto outputType = user::outputType();

        const auto colParams = collision_ -> prepareColParams(user::colParams());

        timer.start("Initialization");
        std::vector<double> u0 = user::initialVelocity(geometry_);
        std::vector<double> rho0 = user::initialDensity(geometry_);
        initializeFields(f_in, *lattice_, geometry_, colParams, u0, rho0);
        collision_->initializeDensityField(f_in, f_out, *lattice_, geometry_, colParams, u0, user::externalForce(), user::initializePressureIterations(), logger);
        timer.stop("Initialization");

        logger.logMessage("Time steps\n");

        if (outputType == OutputType::TSV || outputType == OutputType::BOTH)
            writeTSV(f_in, *lattice_, geometry_, "../outputTSV", 0);

        if (outputType == OutputType::VTI || outputType == OutputType::BOTH)
            writeVTI(f_in, *lattice_, geometry_, "../outputVTI", 0);

        for (int t = 1; t <= user::totalSteps(); ++t) {

            timer.start("Collision");
            collision_->computeCollision(f_in, *lattice_, geometry_, colParams, user::externalForce());
            timer.stop("Collision");

            timer.start("Streaming");
            performStreaming(f_in, f_out, *lattice_, geometry_);
            timer.stop("Streaming");

            if (t % user::writeInterval() == 0) {
                if (outputType == OutputType::TSV || outputType == OutputType::BOTH)
                    writeTSV(f_out, *lattice_, geometry_, "../outputTSV", t);

                if (outputType == OutputType::VTI || outputType == OutputType::BOTH)
                    writeVTI(f_out, *lattice_, geometry_, "../outputVTI", t);
            }

            logger.logStep(t,user::totalSteps());

            std::swap(f_in, f_out);
        }

        logger.endLine();
        timer.report();

    }
};

std::unique_ptr<Simulation> Simulation::create() {
    return std::make_unique<LBMSimulation>();
}

const LatticeModel& Simulation::getLattice() const {
    return *lattice_;
}

const Geometry& Simulation::getGeometry() const {
    return geometry_;
}