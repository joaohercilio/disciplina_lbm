#include "Problem.hpp"

// ----- PROBLEM SETUP (the user should edit this namespace) ----- //

namespace user {

    using Lattice = D2Q9; 
    using Collision = BGK;
    Streaming streaming;

    // Collision operator parameters
    const double tau = 0.8;
    const double omega = 1.0/tau;
    CollisionParameters makeParams() {
        CollisionParameters params;
        params.alphaEq = omega;
        params.alphaNonEq = 1.0 - omega;
        return params;
    }

    // Problem setup
    const int N = 8;
    const double rho0 = 1.0;
    const double nu = 1.0/3.0 * (tau - 0.5);
    const double nu_phy = 0.001;
    const double Re = 10.0;
    const double L = 1.0;
    const double U_phy = Re*nu_phy/L;
    const double finalPhysicalTime = L*L/(4*M_PI*M_PI*nu_phy)*std::log(2);
    const double h = L/N;
    const double delta = nu/nu_phy * h*h;
    const double U_latt = U_phy*delta/h;
    const double finalSymTime = finalPhysicalTime / delta; 
    const int maxSteps = ceil(finalSymTime);
    const double Ma = U_latt / (1.0 / std::sqrt(3));

    std::vector<double> initialVelocity(int x, int y, int z) {
        return {U_latt * std::sin(2.0*M_PI/N * (x+0.5)), 0, 0};
    }

    double initialDensity (int x, int y, int z) {
        return rho0;
    }

    int totalSteps () {
        return maxSteps;
    }

    int writeInterval() {
        return 1;
    }

    // Geometry definition
    Geometry problemGeometry() {
        Geometry geo(N, 1, 1);
        // geo.setSolid(0, 0, 0);
        return geo;
    }
}

class Case : public Problem {

public:
    Case() : Problem(
            std::make_unique<user::Lattice>(), 
            std::make_unique<user::Collision>(),
            user::problemGeometry(),
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
                    if (geometry_.getNode(x, y, z) == NodeType::Fluid) {

                        std::vector<double> u = user::initialVelocity(x, y, z);

                        double rho = user::initialDensity(x, y, z);

                        std::vector<double> feq(numOfVel, 0.0);
                        lattice_->computeEquilibrium(feq.data(), rho, u[0], u[1], u[2]); 

                        int idx = (z * ny * nx + y * nx + x) * numOfVel;
                        for (int k = 0; k < numOfVel; ++k) {
                            f1_[idx + k] = feq[k];
                        }
                    }
                }
            }
        }
    }

    void run() override 
    {
        initialize();

        CollisionParameters params = user::makeParams();

        std::vector<double>* f_in  = &f1_;
        std::vector<double>* f_out = &f2_;

        writeTSV(*f_in, *lattice_, geometry_, "../output", 0);

        for (int t = 0; t < user::totalSteps(); ++t) {

            collision_->computeCollision(*f_in, params, *lattice_, geometry_);

            streaming_.performStreaming(*f_in, *f_out, *lattice_, geometry_);

            if (t % user::writeInterval() == 0) writeTSV(*f_out, *lattice_, geometry_, "../output", t+1);

            std::swap(f_in, f_out);
        }

        // Guarantee that the final distribution is in f1_
        if (user::maxSteps % 2 != 0) {
            f1_ = f2_;
        }
    }
};

const LatticeModel& Problem::getLattice() const {
    return *lattice_;
}

const Geometry& Problem::getGeometry() const {
    return geometry_;
}

std::unique_ptr<Problem> Problem::create() {
    return std::make_unique<Case>();
}
