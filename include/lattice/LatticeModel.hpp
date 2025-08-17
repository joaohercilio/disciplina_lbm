#pragma once

#include <vector>
#include <cmath>

class LatticeModel {

protected:
    int numOfDim_;
    int numOfVel_;
    double cs_;
    std::vector<int> cx_;
    std::vector<int> cy_;
    std::vector<int> cz_;
    std::vector<double> w_;

public:
    LatticeModel(
        int numOfDim, int numOfVel, double cs,
        std::vector<int> cx, std::vector<int> cy, std::vector<int> cz,
        std::vector<double> w
    ):    
        numOfDim_(numOfDim),
        numOfVel_(numOfVel),
        cs_(cs),
        cx_(cx),
        cy_(cy),
        cz_(cz),
        w_(w)
    {}
 
    virtual ~LatticeModel() = default;

    int getNumOfVel() const { return numOfVel_; }
    int getNumOfDim() const { return numOfDim_; }
    const std::vector<int>& getCx() const { return cx_; }
    const std::vector<int>& getCy() const { return cy_; }
    const std::vector<int>& getCz() const { return cz_; }    
    const std::vector<double>& getWeights() const { return w_; }
    double getCs() const { return cs_; }

    virtual void computeFields(const double* f, double& rho, double& vx, double& vy, double& vz) const = 0;
    virtual void computeEquilibrium(double* feq, const double rho, const double vx, const double vy, const double vz) const = 0;
};
