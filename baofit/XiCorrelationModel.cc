// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/XiCorrelationModel.h"

#include "likely/AbsBinning.h"
#include "likely/Interpolator.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::XiCorrelationModel::XiCorrelationModel(likely::AbsBinningCPtr rbins, double zref,
double rVetoMin, double rVetoMax, std::string const &method)
: AbsCorrelationModel("Xi Correlation Model"), _method(method), _zref(zref)
{
    // Create parameters at the center of each radial bin.
    boost::format pname("xi%d-%d");
    for(int ell = 0; ell <= 4; ell += 2) {
        double perror = 1;
        if(2 == ell) perror = 0.1;
        else if(4 == ell) perror = 0.01;
        for(int index = 0; index < rbins->getNBins(); ++index) {
            double rval(rbins->getBinCenter(index));
            if(rval > rVetoMin && rval < rVetoMax) continue;
            defineParameter(boost::str(pname % ell % index),0,perror);
            if(0 == ell) _rValues.push_back(rval);
        }
    }
    // Define linear bias parameters.
    defineParameter("alpha-bias",3.8,0.3);
    defineParameter("beta",1.0,0.1);
    defineParameter("(1+beta)*bias",-0.34,0.03);
    // Pick a normalization scale that gives parameter values of order one.
    _normScale = 1e-2;
}

local::XiCorrelationModel::~XiCorrelationModel() { }

void local::XiCorrelationModel::_initializeInterpolators() const {
    int index, npoints(_rValues.size());
    // Do we need to (re)initialize our xi0 interpolator?
    for(index = 0; index < npoints; ++index) {
        if(isParameterValueChanged(index)) break;
    }
    if(index < npoints) {
        _xiValues.resize(0);
        for(index = 0; index < npoints; ++index) {
            _xiValues.push_back(getParameterValue(index));
        }
        _xi0.reset(new likely::Interpolator(_rValues,_xiValues,_method));
    }
    // Do we need to (re)initialize our xi2 interpolator?
    for(index = npoints; index < 2*npoints; ++index) {
        if(isParameterValueChanged(index)) break;
    }
    if(index < 2*npoints) {
        _xiValues.resize(0);
        for(index = npoints; index < 2*npoints; ++index) {
            _xiValues.push_back(getParameterValue(index));
        }
        _xi2.reset(new likely::Interpolator(_rValues,_xiValues,_method));
    }
    // Do we need to (re)initialize our xi4 interpolator?
    for(index = 2*npoints; index < 3*npoints; ++index) {
        if(isParameterValueChanged(index)) break;
    }
    if(index < 3*npoints) {
        _xiValues.resize(0);
        for(index = 2*npoints; index < 3*npoints; ++index) {
            _xiValues.push_back(getParameterValue(index));
        }
        _xi4.reset(new likely::Interpolator(_rValues,_xiValues,_method));
    }
}

double local::XiCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    // Fetch linear bias parameters.
    double alpha = getParameterValue("alpha-bias");
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
    // Rebuild our interpolators, if necessary.
    if(anyChanged) _initializeInterpolators();
    // Combine the multipoles.
    return _normScale*bias*bias*zfactor*(C0*(*_xi0)(r) + C2*P2*(*_xi2)(r) + C4*P4*(*_xi4)(r))/(r*r);
}

double local::XiCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    // Fetch linear bias parameters.
    double alpha = getParameterValue("alpha-bias");
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // Rebuild our interpolators, if necessary.
    if(anyChanged) _initializeInterpolators();
    // Return the appropriately normalization multipole.
    double norm = _normScale*bias*bias*zfactor;
    if(multipole == cosmo::Hexadecapole) {
        return norm*(8./35.)*beta*beta*(*_xi4)(r)/(r*r);
    }
    else if(multipole == cosmo::Quadrupole) {
        return norm*4*beta*((1./3.) + beta/7.)*(*_xi2)(r)/(r*r);
    }
    else {
        return norm*(1 + beta*((2./3.) + beta/5.))*(*_xi0)(r)/(r*r);
    }
}

void  local::XiCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
    out << "Interpolating with " << _rValues.size() << " points covering " << _rValues[0] << " to "
        << _rValues[_rValues.size()-1] << " Mpc/h" << std::endl;
}
