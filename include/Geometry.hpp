#pragma once

#include "boundary/BoundaryModel.hpp"

#include <vector>
#include <cstdint>

enum class NodeType : uint8_t { Solid, Fluid };

class Geometry {
    
private:
    int nx_, ny_, nz_;
    std::vector<NodeType> nodes_;
    std::vector<BoundaryType> boundaryTypes_;

public:

    Geometry(int nx, int ny, int nz) : 
        nx_(nx), ny_(ny), nz_(nz),
        nodes_(nx * ny * nz, NodeType::Fluid),
        boundaryTypes_(nx * ny * nz, BoundaryType::None) {}

    void setSolid(int x, int y, int z)   { nodes_[getIndex(x, y, z)] = NodeType::Solid; }
    void setFluid(int x, int y, int z)   { nodes_[getIndex(x, y, z)] = NodeType::Fluid; }

    NodeType getNode(int x, int y, int z) const { return nodes_[getIndex(x, y, z)]; }
    NodeType getNode(int i) const { return nodes_[i]; }

    void setBoundaryType(int x, int y, int z, BoundaryType type) {
        int idx = getIndex(x, y, z);
        boundaryTypes_[idx] = type;
    }

    BoundaryType getBoundaryType(int x, int y, int z) const {
        int idx = getIndex(x, y, z);
        return boundaryTypes_[idx];
    }
    
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }

    int getNumOfPoints() const;

    int getIndex(int x, int y, int z) const;

    void getCoords(int index, int& x, int& y, int& z) const;
    void getCoords(int index, int& x, int& y) const;

};
