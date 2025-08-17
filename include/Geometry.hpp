#pragma once

#include <vector>
#include <cstdint>

enum class NodeType : uint8_t { Solid, Fluid };

class Geometry {
    
private:
    int nx_, ny_, nz_;
    std::vector<NodeType> nodes_;

public:

    Geometry(int nx, int ny, int nz)
        : nx_(nx), ny_(ny), nz_(nz),
          nodes_(nx * ny * nz, NodeType::Fluid
        ) {}

    void setSolid(int x, int y, int z)   { nodes_[getIndex(x, y, z)] = NodeType::Solid; }
    void setFluid(int x, int y, int z)   { nodes_[getIndex(x, y, z)] = NodeType::Fluid; }

    NodeType getNode(int x, int y, int z) const { return nodes_[getIndex(x, y, z)]; }
    NodeType getNode(int i) const { return nodes_[i]; }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }

    int getNumOfPoints() const;

    int getIndex(int x, int y, int z) const;

    void getCoords(int index, int& x, int& y, int& z) const;
    void getCoords(int index, int& x, int& y) const;

};
