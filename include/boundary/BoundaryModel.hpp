#pragma once

#include "lattice/LatticeModel.hpp"

//forward declaration
class Geometry;

enum class BoundaryType
{ 
    None,
    BounceAndBackSouth,
    BounceAndBackNorth 
};

class BoundaryModel {

public:

    virtual ~BoundaryModel() = default;

    virtual void applyBoundary(std::vector<double>& f, 
                            const LatticeModel& lattice, 
                            const Geometry& geometry,
                            int x, int y, int z) = 0;

private:

    virtual void applySouthBoundary(std::vector<double>& f, 
                                    const LatticeModel& lattice, 
                                    const Geometry& geometry,
                                    int x, int y, int z) = 0;
                                    
    virtual void applyNorthBoundary(std::vector<double>& f, 
                                    const LatticeModel& lattice, 
                                    const Geometry& geometry,
                                    int x, int y, int z) = 0;

};