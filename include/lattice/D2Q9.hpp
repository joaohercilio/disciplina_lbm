#pragma once

#include "LatticeModel.hpp"
#include <iostream>

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

    void computeFields(const double* f, double& rho, double& vx, double& vy, double& /*vz*/) const override;
    void computeEquilibrium(double* feq, const double rho, const double vx, const double vy, double /*vz*/) const override;
};
