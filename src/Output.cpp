#include "Output.hpp"

void writeTSV(const std::vector<double>& f, const LatticeModel& lattice, const Geometry& geo, const std::string& path, int t) {
    std::filesystem::create_directories(path);
    std::ofstream out(path + "/out_" + std::to_string(t) + ".tsv");
    
    out << std::setprecision(10);
    out << std::scientific;
    out << "x\ty\tz\tvx\tvy\tvz\trho\n";
    
    int numOfVel = lattice.getNumOfVel();
    int n = geo.getNumOfPoints();
    
    for (int id = 0; id < n; ++id) {
        int x, y, z;
        geo.getCoords(id, x, y, z);
        const double* mapF = f.data() + id * numOfVel;
        double rho, vx, vy, vz;
        lattice.computeFields(mapF, rho, vx, vy, vz);
        out << x << "\t" << y << "\t" << z << "\t" << vx << "\t" << vy << "\t" << vz << "\t" << rho << "\n";
    }

    out.close();
}

void writeVTI(const std::vector<double>& f, const LatticeModel& lattice, const Geometry& geo, const std::string& path, int t) {
    std::filesystem::create_directories(path);

    char filename[256];
    std::sprintf(filename, "%s/output_%07d.vti", path.c_str(), t);

    int nx = geo.nx();
    int ny = geo.ny();
    int nz = geo.nz();
    int numOfVel = lattice.getNumOfVel();
    int numOfPoints = geo.getNumOfPoints();

    std::vector<double> rho(numOfPoints, 0.0);
    std::vector<double> vel(3 * numOfPoints, 0.0);

    for (int id = 0; id < numOfPoints; ++id) {
        int x, y, z;
        geo.getCoords(id, x, y, z);

        const double* mapF = f.data() + id * numOfVel;
        double vx, vy, vz, density;
        lattice.computeFields(mapF, density, vx, vy, vz);

        int pos = x + y * nx + z * nx * ny;
        rho[pos] = density;
        vel[3*pos] = vx;
        vel[3*pos+1] = vy;
        vel[3*pos+2] = vz;
    }

    VTIWriter vti(filename);
    vti.setCompress();
    vti.setWholeExtent(0, 0, 0, nx-1, ny-1, nz-1);
    vti.setSpacing(1.0, 1.0, 1.0);
    vti.setOrigin(0.0, 0.0, 0.0);

    vti.addPointData("Density", "Float64", "binary", 1, (unsigned char*) rho.data());
    vti.addPointData("Velocity", "Float64", "binary", 3, (unsigned char*) vel.data());

    vti.write();
}