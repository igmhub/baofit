// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/XiCorrelationModel.h"

#include "likely/AbsBinning.h"

#include "boost/format.hpp"

namespace local = baofit;

local::XiCorrelationModel::XiCorrelationModel(likely::AbsBinningCPtr rbins)
: AbsCorrelationModel("Xi"), _rbins(rbins)
{
    // Create parameters at the center of each radial bin.
    boost::format pname("xi%d-%d");
    for(int ell = 0; ell <= 4; ell += 2) {
        double perror = 1;
        if(2 == ell) perror = 0.1;
        else if(4 == ell) perror = 0.01;
        for(int index = 0; index < _rbins->getNBins(); ++index) {
            double rval(_rbins->getBinCenter(index));
            defineParameter(boost::str(pname % ell % index),0,perror);
        }
    }
}

local::XiCorrelationModel::~XiCorrelationModel() { }

double local::XiCorrelationModel::evaluate(double r, double mu, double z,
std::vector<double> const &params) const {
    
    return 0;
}

double local::XiCorrelationModel::evaluate(double r, cosmo::Multipole multipole, double z,
std::vector<double> const &params) const {
    int index = _rbins->getBinIndex(r);
    if(multipole == cosmo::Quadrupole) index += _rbins->getNBins();
    else if(multipole == cosmo::Hexadecapole) index += 2*_rbins->getNBins();
    return 1e-6*params[index];
}
