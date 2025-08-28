#include "Geometry.hpp"

int Geometry::getNumOfPoints() const 
{ 
    return nx_ * ny_ * nz_; 
}

int Geometry::getIndex(int x, int y, int z) const 
{
    return x + y * nx_ + z * nx_ * ny_;
}

void Geometry::getCoords(int index, int& x, int& y, int& z) const 
{
    z = index / (nx_ * ny_);
    int rem = index % (nx_ * ny_);
    y = rem / nx_;
    x = rem % nx_;
}

void Geometry::getCoords(int index, int& x, int& y) const 
{
    y = index / nx_;
    x = index % nx_;
}

int Geometry::getVelocityIndex(int id, int dir) const
{
    return id + dir * (nx_ * ny_ * nz_);
}
