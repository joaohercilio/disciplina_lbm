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
        user::print();
        
        auto& f_in  = f1_;
        auto& f_out = f2_;

        const auto colParams = collision_ -> prepareColParams(user::colParams());

        Timer timer;

        auto outputType = user::outputType();

        timer.start("Initialization");
        std::vector<double> u0 = user::initialVelocity(geometry_);
        std::vector<double> rho0 = user::initialDensity(geometry_);
        initializeFields(f_in, *lattice_, geometry_, colParams, u0, rho0);
        //writeTSV(f_in, *lattice_, geometry_, "../outputTSV", 66666666);
        //writeVTI(f_in, *lattice_, geometry_, "../outputVTI", 66666666);
        //collision_->initializeDensityField(f_in, f_out, *lattice_, geometry_, colParams, u0, user::externalForce(), user::initializePressureIterations());
        timer.stop("Initialization");

        if (outputType == user::OutputType::TSV || outputType == user::OutputType::BOTH)
            writeTSV(f_in, *lattice_, geometry_, "../outputTSV", 0);

        if (outputType == user::OutputType::VTI || outputType == user::OutputType::BOTH)
            writeVTI(f_in, *lattice_, geometry_, "../outputVTI", 0);

        for (int t = 1; t <= user::totalSteps(); ++t) {

            timer.start("Collision");
            collision_->computeCollision(f_in, *lattice_, geometry_, colParams, user::externalForce());
            timer.stop("Collision");

            timer.start("Streaming");
            performStreaming(f_in, f_out, *lattice_, geometry_);
            timer.stop("Streaming");

            if (t % user::writeInterval() == 0) {
                if (outputType == user::OutputType::TSV || outputType == user::OutputType::BOTH)
                    writeTSV(f_out, *lattice_, geometry_, "../outputTSV", t);

                if (outputType == user::OutputType::VTI || outputType == user::OutputType::BOTH)
                    writeVTI(f_out, *lattice_, geometry_, "../outputVTI", t);
            }
            if ( t % 50 == 0 ) {std::cout << t << std::endl;}

            std::swap(f_in, f_out);
        }

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
