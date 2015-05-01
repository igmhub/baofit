// Created 24-Apr-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "baofit/MetalCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include <cmath>

namespace local = baofit;

local::MetalCorrelationModel::MetalCorrelationModel(AbsCorrelationModel *base)
: AbsCorrelationModel("Metal Correlation Model"), _base(base ? *base:*this) {
    // Define parameters
    _indexBase = _base.defineParameter("metal ampl0",1,0.1);
    _base.defineParameter("metal ampl1",1,0.1);
    _base.defineParameter("metal ampl2",1,0.1);
    _base.defineParameter("metal ampl3",1,0.1);
    _base.defineParameter("metal width0",6,0.6);
}

local::MetalCorrelationModel::~MetalCorrelationModel() { }

double local::MetalCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    double rpar = std::fabs(r*mu);
    double rperp = r*std::sqrt(1-mu*mu);
    double exppar, expperp, corr0(0), corr1(0), corr2(0), corr3(0);
    // Template 0
    double rpar0 = 59.533;
    double sigpar0 = _base.getParameterValue(_indexBase+4);
    double rperp0 = 2.0;
    double sigperp0 = 5.55309;
    double ampl0 = _base.getParameterValue(_indexBase);
    exppar = (rpar-rpar0)/sigpar0;
    expperp = (rperp-rperp0)/sigperp0;
    if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr0 = 1e-4*ampl0*std::exp(-exppar*exppar-expperp);
    // Template 1
    double rpar1 = 111.344;
    double sigpar1 = 4.88626;
    double rperp1 = 2.0;
    double sigperp1 = 4.13032;
    double ampl1 = _base.getParameterValue(_indexBase+1);
    exppar = (rpar-rpar1)/sigpar1;
    expperp = (rperp-rperp1)/sigperp1;
    if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr1 = 1e-4*ampl1*std::exp(-exppar*exppar-expperp);
    // Template 2
    double rpar2 = 134.998;
    double sigpar2 = 5.489;
    double rperp2 = 2.0;
    double sigperp2 = 6.30484;
    double ampl2 = _base.getParameterValue(_indexBase+2);
    exppar = (rpar-rpar2)/sigpar2;
    expperp = (rperp-rperp2)/sigperp2;
    if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr2 = 1e-4*ampl2*std::exp(-exppar*exppar-expperp);
    // Template 3
    double rpar3 = 174.525;
    double sigpar3 = 8.48942;
    double rperp3 = 2.0;
    double sigperp3 = 7.85636;
    double ampl3 = _base.getParameterValue(_indexBase+3);
    exppar = (rpar-rpar3)/sigpar3;
    expperp = (rperp-rperp3)/sigperp3;
    if(std::fabs(exppar)<5 && std::fabs(expperp)<10) corr3 = 1e-4*ampl3*std::exp(-exppar*exppar-expperp);
    // Add the contributions
    return corr0 + corr1 + corr2 + corr3;
}

void  local::MetalCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using metal correlation templates" << std::endl;
}