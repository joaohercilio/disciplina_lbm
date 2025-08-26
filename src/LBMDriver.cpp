#include "user.hpp"

class LBMSimulation : public Simulation {
public:
    LBMSimulation() : Simulation (
            std::make_unique<user::LatticeModel>(), 
            std::make_unique<user::CollisionModel>(),
            user::problemGeometry(),
            std::make_unique<user::BoundaryModel>()
        ) {}

    void run() override 
    {
        auto& f_in  = f1_;
        auto& f_out = f2_;

        const auto colParams = collision_ -> prepareColParams(user::colParams());

        Timer timer;

        user::print();

        timer.start("Initialization");
        initializeFields(f_in, *lattice_, geometry_, colParams, user::initialVelocity(geometry_), user::initialDensity(geometry_));
        timer.stop("Initialization");

        writeTSV(f_in, *lattice_, geometry_, "../outputTSV", 0);
        writeVTI(f_in, *lattice_, geometry_, "../outputVTI", 0);

        for (int t = 1; t <= user::totalSteps(); ++t) {

            timer.start("Collision");
            collision_->computeCollision(f_in, *lattice_, geometry_, colParams, user::externalForce());
            timer.stop("Collision");

            timer.start("Streaming");
            performStreaming(f_in, f_out, *lattice_, geometry_);
            timer.stop("Streaming");

            timer.start("Boundary");
            for (int z = 0; z < geometry_.nz(); ++z) {
                for (int y = 0; y < geometry_.ny(); ++y) {
                    for (int x = 0; x < geometry_.nx(); ++x) {
                        if (geometry_.getBoundaryType(x, y, z) != BoundaryType::None) {
                            boundary_->applyBoundary(f_out, *lattice_, geometry_, x, y, z);
                        }
                    }
                }
            }
            timer.stop("Boundary");

            if ( t % user::writeInterval() == 0) writeTSV(f_out, *lattice_, geometry_, "../outputTSV", t);
            if ( t % user::writeInterval() == 0) writeVTI(f_out, *lattice_, geometry_, "../outputVTI", t);

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
