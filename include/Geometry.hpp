#pragma once

#include "boundary/BoundaryModel.hpp"

#include <vector>
#include <cstdint>


/**
 * @brief Enum representing the type of a node in the domain.
 */
enum class NodeType : uint8_t { Solid, Fluid };


/**
 * @brief Class that defines the problem geometry.
 *
 * This class stores the dimensions of the computational domain, the type of each node
 * (solid or fluid), and the boundary type for each node. 
 */
class Geometry {
    
private:
    int nx_; ///< Number of nodes in x-direction
    int ny_; ///< Number of nodes in y-direction
    int nz_; ///< Number of nodes in z-direction
    std::vector<NodeType> nodes_; ///< Vector storing node types (Solid/Fluid)
    std::vector<BoundaryType> boundaryTypes_; ///< Vector storing boundary types for each node

public:
    /**
     * @brief Constructs a Geometry object with given dimensions.
     *
     * Initially, all nodes are set to Fluid and all boundary types are None.
     *
     * @param nx Number of nodes in X-direction
     * @param ny Number of nodes in Y-direction
     * @param nz Number of nodes in Z-direction
     */
    Geometry(int nx, int ny, int nz) : 
        nx_(nx), ny_(ny), nz_(nz),
        nodes_(nx * ny * nz, NodeType::Fluid),
        boundaryTypes_(nx * ny * nz, BoundaryType::None) {}

    /**
     * @brief Marks a node as Solid.
     * @param x X-coordinate of the node
     * @param y Y-coordinate of the node
     * @param z Z-coordinate of the node
     */
    void setSolid(int x, int y, int z)   { nodes_[getIndex(x, y, z)] = NodeType::Solid; }

    /**
     * @brief Marks a node as Fluid.
     * @param x X-coordinate of the node
     * @param y Y-coordinate of the node
     * @param z Z-coordinate of the node
     */
    void setFluid(int x, int y, int z)   { nodes_[getIndex(x, y, z)] = NodeType::Fluid; }

    /**
     * @brief Returns the type of a node given its 3D coordinates.
     * @param x X-coordinate
     * @param y Y-coordinate
     * @param z Z-coordinate
     * @return NodeType of the specified node
     */
    NodeType getNode(int x, int y, int z) const { return nodes_[getIndex(x, y, z)]; }

    /**
     * @brief Returns the type of a node given its linear index.
     * @param i Linear index of the node
     * @return NodeType of the specified node
     */
    NodeType getNode(int i) const { return nodes_[i]; }

    /**
     * @brief Sets the boundary type of a node.
     * @param x X-coordinate
     * @param y Y-coordinate
     * @param z Z-coordinate
     * @param type BoundaryType to assign
     */
    void setBoundaryType(int x, int y, int z, BoundaryType type) {
        int idx = getIndex(x, y, z);
        boundaryTypes_[idx] = type;
    }

    /**
     * @brief Returns the boundary type of a node.
     * @param x X-coordinate
     * @param y Y-coordinate
     * @param z Z-coordinate
     * @return BoundaryType of the node
     */
    BoundaryType getBoundaryType(int x, int y, int z) const {
        int id = getIndex(x, y, z);
        return boundaryTypes_[id];
    }

    /**
     * @brief Returns the boundary type of a node.
     * @param id Index of the node
     * @return BoundaryType of the node
     */
    BoundaryType getBoundaryType(int id) const {
        return boundaryTypes_[id];
    }

    /// Returns the number of nodes in X-direction
    int nx() const { return nx_; }
    
    /// Returns the number of nodes in Y-direction    
    int ny() const { return ny_; }

    /// Returns the number of nodes in Z-direction
    int nz() const { return nz_; }

    /**
     * @brief Returns the total number of points in the domain.
     * @return Total number of nodes (nx * ny * nz)
     */
    int getNumOfPoints() const;

    /**
     * @brief Converts 3D coordinates into a linear index.
     * @param x X-coordinate
     * @param y Y-coordinate
     * @param z Z-coordinate
     * @return Linear index corresponding to the coordinates
     */
    int getIndex(int x, int y, int z) const;

    /**
     * @brief Converts a linear index into 3D coordinates.
     * @param index Linear index
     * @param x Reference to store X-coordinate
     * @param y Reference to store Y-coordinate
     * @param z Reference to store Z-coordinate
     */
    void getCoords(int index, int& x, int& y, int& z) const;

    /**
     * @brief Converts a linear index into 2D coordinates (for 2D problems).
     * @param index Linear index
     * @param x Reference to store X-coordinate
     * @param y Reference to store Y-coordinate
     */
    void getCoords(int index, int& x, int& y) const;

    /**
     * @brief Returns the velocity index for the direction specified.
     * @param index Linear node index
     * @param dir Desired direction (x, y and z)
     * 
     *   dir = 0 -> u_x
     * 
     *   dir = 1 -> u_y
     * 
     *   dir = 2 -> u_z
     *
     * @return The index in the velocity vector
     */
    int getVelocityIndex(int index, int dir) const;

};
