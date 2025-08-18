#include "user.hpp"

class LBMSimulation : public Simulation {
public:
    LBMSimulation() : Simulation (
            std::make_unique<user::LatticeModel>(), 
            std::make_unique<user::CollisionModel>(),
            user::problemGeometry(),
            std::make_unique<user::BoundaryModel>(),
            user::streaming
        ) {}

    void initialize() override 
    {
        int numOfVel = lattice_->getNumOfVel();
        int nx = geometry_.nx();
        int ny = geometry_.ny();
        int nz = geometry_.nz();

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    //if (geometry_.getNode(x, y, z) == NodeType::Fluid) {

                        std::vector<double> u = user::initialVelocity(x, y, z);

                        double rho = user::initialDensity(x, y, z);

                        std::vector<double> feq(numOfVel, 0.0);
                        lattice_->computeEquilibrium(feq.data(), rho, u[0], u[1], u[2]); 

                        int idx = (z * ny * nx + y * nx + x) * numOfVel;
                        for (int k = 0; k < numOfVel; ++k) {
                            f1_[idx + k] = feq[k];
                        }
                    //}
                }
            }
        }
    }

    void run() override 
    {
        initialize();

        const auto colParams = collision_ -> prepareColParams(user::colParams());

        std::vector<double>* f_in  = &f1_;
        std::vector<double>* f_out = &f2_;

        writeTSV(*f_in, *lattice_, geometry_, "../output", 0);

        for (int t = 0; t < user::totalSteps(); ++t) {

            collision_->computeCollision(*f_in, *lattice_, geometry_, colParams, user::externalForce());

            streaming_.performStreaming(*f_in, *f_out, *lattice_, geometry_);

            for (int z = 0; z < geometry_.nz(); ++z) {
                for (int y = 0; y < geometry_.ny(); ++y) {
                    for (int x = 0; x < geometry_.nx(); ++x) {
                        if (geometry_.getBoundaryType(x, y, z) != BoundaryType::None) {
                            boundary_->applyBoundary(*f_out, *lattice_, geometry_, x, y, z);
                        }
                    }
                }
            }
            
            if (t % user::writeInterval() == 0) writeTSV(*f_out, *lattice_, geometry_, "../output", t+1);

            std::swap(f_in, f_out);
        }

        // Guarantee that the final distribution is in f1_
        if (user::totalSteps() % 2 != 0) {
            f1_ = f2_;
        }

        writeTSV(*f_out, *lattice_, geometry_, "../output", user::totalSteps());

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
