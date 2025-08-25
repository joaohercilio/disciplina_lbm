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

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    if (geometry_.getNode(x, y, z) == NodeType::Fluid) {

                        std::vector<double> u = user::initialVelocity(x, y, z);

                        double rho = user::initialDensity(x, y, z);

                        std::vector<double> feq(numOfVel, 0.0);
                        lattice_->computeEquilibrium(feq.data(), rho, u[0], u[1], u[2]); 

                        int id = geometry_.getIndex(x, y, z);
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
        using std::chrono::high_resolution_clock;
        using std::chrono::duration_cast;
        using std::chrono::microseconds;

        microseconds collisionDuration = microseconds::zero();
        microseconds streamingDuration = microseconds::zero();
        microseconds boundaryDuration  = microseconds::zero();

        user::print();

        initialize();

        const auto colParams = collision_ -> prepareColParams(user::colParams());

        std::vector<double>* f_in  = &f1_;
        std::vector<double>* f_out = &f2_;

        writeTSV(*f_in, *lattice_, geometry_, "../output", 0);

        for (int t = 1; t <= user::totalSteps(); ++t) {

            auto startCol = high_resolution_clock::now();
            collision_->computeCollision(*f_in, *lattice_, geometry_, colParams, user::externalForce());
            auto stopCol = high_resolution_clock::now();
            collisionDuration += duration_cast<microseconds>(stopCol - startCol);

            auto startStr = high_resolution_clock::now();
            performStreaming(*f_in, *f_out, *lattice_, geometry_);
            auto stopStr = high_resolution_clock::now();
            streamingDuration += duration_cast<microseconds>(stopStr - startStr);

            auto startBnd = high_resolution_clock::now();
            for (int z = 0; z < geometry_.nz(); ++z) {
                for (int y = 0; y < geometry_.ny(); ++y) {
                    for (int x = 0; x < geometry_.nx(); ++x) {
                        if (geometry_.getBoundaryType(x, y, z) != BoundaryType::None) {
                            boundary_->applyBoundary(*f_out, *lattice_, geometry_, x, y, z);
                        }
                    }
                }
            }
            auto stopBnd = high_resolution_clock::now();
            boundaryDuration += duration_cast<microseconds>(stopBnd - startBnd);

            if ( t % user::writeInterval() == 0) writeTSV(*f_out, *lattice_, geometry_, "../output", t);

            if ( t % 50 == 0 ) {std::cout << t << std::endl;}

            std::swap(f_in, f_out);
        }

        std::cout << "Collision duration: " << collisionDuration.count() / 1000 << " ms\n"
                  << "Streaming duration: " << streamingDuration.count() / 1000 << " ms\n"
                  << "Boundary duration:  "  << boundaryDuration.count() / 1000 << " ms\n";
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
