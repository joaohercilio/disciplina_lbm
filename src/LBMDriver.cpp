#include "user.hpp"

class LBMSimulation : public Simulation {
public:
    LBMSimulation() : Simulation (
            std::make_unique<user::LatticeModel>(), 
            std::make_unique<user::CollisionModel>(),
            user::problemGeometry(),
            std::make_unique<user::BoundaryModel>()
        ) {}

    void initialize() override 
    {
        int numOfVel = lattice_->getNumOfVel();
        int nx = geometry_.nx();
        int ny = geometry_.ny();
        int nz = geometry_.nz();

        std::vector<double> u = user::initialVelocity(geometry_);
        std::vector<double> rho = user::initialDensity(geometry_);

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    if (geometry_.getNode(x, y, z) == NodeType::Fluid) {

                        int id = geometry_.getIndex(x,y,z);

                        std::vector<double> feq(numOfVel, 0.0);

                        lattice_->computeEquilibrium(feq.data(), rho[id], u[id], u[id+1], u[id+2]); 

                        for (int k = 0; k < numOfVel; ++k) {
                            f1_[id*numOfVel + k] = feq[k];
                        }

                    } else {
                        // Initialize Solid Nodes with zero velocity and density

                        std::vector<double> feq(numOfVel, 0.0);
                        lattice_->computeEquilibrium(feq.data(), 1.0, 0.0, 0.0, 0.0);
                        int id = geometry_.getIndex(x, y, z);
                        for (int k = 0; k < numOfVel; ++k) {
                            f1_[id*numOfVel + k] = feq[k];
                        }
                    }
                }
            }
        }
    }

    void run() override 
    {
        user::print();

        Timer timer;

        timer.start("Initialization");
        initialize();
        timer.stop("Initialization");

        const auto colParams = collision_ -> prepareColParams(user::colParams());

        auto& f_in  = f1_;
        auto& f_out = f2_;

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
