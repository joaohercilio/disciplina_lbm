#pragma once

#include "LatticeModel.hpp"

class D2Q5 : public LatticeModel {

public:
    D2Q5() : LatticeModel(

    /* NUM_OF_DIM , NUM_OF_VEL */    2, 5,                                    
    /* Speed of Sound*/              1/sqrt(3),                              
    /* cx */                         {0,  1, -1,  0,  0},   
    /* cy */                         {0,  0,  0,  1, -1},    
    /* cz */                         std::vector<int>(5,0),
    /* Weights */                    {2.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0} 

    ) {}

    void computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const override;
    void computeEquilibrium(double* feq, const double rho, const double vx, const double vy, double vz) const override;
};
