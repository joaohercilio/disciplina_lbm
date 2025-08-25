#pragma once

#include "LatticeModel.hpp"

class D3Q19 : public LatticeModel {

public:

    D3Q19() : LatticeModel (

    /* NUM_OF_DIM , NUM_OF_VEL */    3, 19,                                    
    /* Speed of Sound*/              1/sqrt(3),                              
    /* cx */                         { 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 },   
    /* cy */                         { 0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0, -1,  1, -1,  1 },  
    /* cz */                         { 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1 },
    /* Weights */                    {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 
                                      1.0/36.0, 1.0/36.0, 1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0 ,1.0/36.0} 
    ) {}

    void computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const override;
    void computeEquilibrium(double* feq, const double rho, const double vx, const double vy, double vz) const override;

};