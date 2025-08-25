#include "BoundaryModel.hpp"

class PeriodicBoundary : public BoundaryModel {

public:

    void applyBoundary(std::vector<double>& f, 
                        const LatticeModel& lattice, 
                        const Geometry& geometry,
                        int x, int y, int z) override;

private:

    void applySouthBoundary(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry,
                            int x, int y, int z) override;

    void applyNorthBoundary(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry,
                            int x, int y, int z) override;

};
