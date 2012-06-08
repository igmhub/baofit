// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/XiCorrelationModel.h"

#include "likely/AbsBinning.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::XiCorrelationModel::XiCorrelationModel(likely::AbsBinningCPtr rbins, double zref)
: AbsCorrelationModel("Xi Correlation Model"), _rbins(rbins), _zref(zref)
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
    // Define linear bias parameters.
    defineParameter("alpha",3.8,0.3);
    defineParameter("beta",1.0,0.1);
    defineParameter("(1+beta)*bias",-0.34,0.03);
}

local::XiCorrelationModel::~XiCorrelationModel() { }

double local::XiCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    // Fetch linear bias parameters.
    double alpha = getParameterValue("alpha");
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // Calculate the Legendre polynomials we need.
    double mu2(mu*mu);
    double P2(1.5*mu2-0.5), P4(4.375*mu2*mu2-3.75*mu2+0.375);
    // Calculate the beta functions for each multipole.
    double C0(1 + beta*((2./3.) + beta/5.)), C2(4*beta*((1./3.) + beta/7.)), C4((8./35.)*beta*beta);
    // Find which radial bin we are in.
    int index(_rbins->getBinIndex(r));
    int nbins(_rbins->getNBins());
    // Combine the multipoles.
    return 1e-6*bias*bias*zfactor*
        (C0*getParameterValue(index) + C2*P2*getParameterValue(index+nbins) +
        C4*P4*getParameterValue(index+2*nbins));
}

double local::XiCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    // Fetch linear bias parameters.
    double alpha = getParameterValue("alpha");
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // Find which radial bin we are in and calculate the appropriate normalization factor.
    int index = _rbins->getBinIndex(r);
    double norm = 1e-6*bias*bias*zfactor;
    if(multipole == cosmo::Quadrupole) {
        index += _rbins->getNBins();
        norm *= (8./35.)*beta*beta;
    }
    else if(multipole == cosmo::Hexadecapole) {
        index += 2*_rbins->getNBins();
        norm *= 4*beta*((1./3.) + beta/7.);
    }
    else {
        norm *= 1 + beta*((2./3.) + beta/5.);
    }
    return norm*getParameterValue(index);
}

void  local::XiCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
}
