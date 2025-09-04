#pragma once

#include "LatticeModel.hpp"

/**
 * @brief Concrete implementation of the D2Q9 lattice model.
 *
 * Inherits all method documentation from LatticeModel.
 */
class D2Q9 : public LatticeModel {

public:

    D2Q9() : LatticeModel(

    /* NUM_OF_DIM , NUM_OF_VEL */    2, 9,                                    
    /* Speed of Sound*/              1/sqrt(3),                              
    /* cx */                         {0,  1, -1,  0,  0,  1, -1, -1,  1},   
    /* cy */                         {0,  0,  0,  1, -1,  1, -1,  1, -1},  
    /* cz */                         std::vector<int>(9,0),
    /* Weights */                    {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0} 
    
    ) {}

    /// \copydoc LatticeModel::computeFields
    void computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const override;

    /// \copydoc LatticeModel::computeEquilibrium
    void computeEquilibrium(double* feq, const double rho, const double vx, const double vy, double vz) const override;

    /// \copydoc LatticeModel::computeMoments
    void computeMoments(const double* f, double* m) const override;

    /// \copydoc LatticeModel::computeEquilibriumMoments
    void computeEquilibriumMoments(double* meq, const double* m) const override;

    /// \copydoc LatticeModel::reconstructDistribution
    void reconstructDistribution(double* f, const double* m) const override;

    /// \copydoc LatticeModel::relaxationMatrix
    std::vector<double> relaxationMatrix(const double tau) const override;
};

