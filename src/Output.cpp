#include "Output.hpp"

void writeTSV(const std::vector<double>& f, const LatticeModel& lattice, const Geometry& geo, const std::string& path, int t) {
    std::filesystem::create_directories(path);
    std::ofstream out(path + "/out_" + std::to_string(t) + ".tsv");
    
    out << std::setprecision(10);
    out << std::scientific;
    out << "x\ty\tz\tvx\tvy\tvz\trho\n";
    
    int numOfVel = lattice.getNumOfVel();
    int n = geo.getNumOfPoints();
    
    for (int i = 0; i < n; ++i) {
        int x, y, z;
        geo.getCoords(i, x, y, z);
        const double* p = f.data() + i * numOfVel;
        double rho, vx, vy, vz;
        lattice.computeFields(p, rho, vx, vy, vz);
        out << x << "\t" << y << "\t" << z << "\t" << vx << "\t" << vy << "\t" << vz << "\t" << rho << "\n";
    }

    out.close();
}
