// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"

#include "cosmo/RsdCorrelationFunction.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "likely/Interpolator.h"
#include "likely/function.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &broadbandName, double zref)
: AbsCorrelationModel(), _zref(zref)
{
    boost::format fileName("%s.%d.dat"),bbandName("%s%c.%d.dat");
    std::string method("cspline");
    cosmo::CorrelationFunctionPtr
        fid0 = likely::createFunctionPtr(likely::createInterpolator(boost::str(fileName % fiducialName % 0),method)),
        fid2 = likely::createFunctionPtr(likely::createInterpolator(boost::str(fileName % fiducialName % 2),method)),
        fid4 = likely::createFunctionPtr(likely::createInterpolator(boost::str(fileName % fiducialName % 4),method)),
        nw0 = likely::createFunctionPtr(likely::createInterpolator(boost::str(fileName % nowigglesName % 0),method)),
        nw2 = likely::createFunctionPtr(likely::createInterpolator(boost::str(fileName % nowigglesName % 2),method)),
        nw4 = likely::createFunctionPtr(likely::createInterpolator(boost::str(fileName % nowigglesName % 4),method)),
        bbc0 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % 'c' % 0),method)),
        bbc2 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % 'c' % 2),method)),
        bbc4 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % 'c' % 4),method)),
        bb10 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % '1' % 0),method)),
        bb12 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % '1' % 2),method)),
        bb14 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % '1' % 4),method)),
        bb20 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % '2' % 0),method)),
        bb22 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % '2' % 2),method)),
        bb24 = likely::createFunctionPtr(likely::createInterpolator(boost::str(bbandName % broadbandName % '2' % 4),method));
    _fid.reset(new cosmo::RsdCorrelationFunction(fid0,fid2,fid4));
    _nw.reset(new cosmo::RsdCorrelationFunction(nw0,nw2,nw4));
    _bbc.reset(new cosmo::RsdCorrelationFunction(bbc0,bbc2,bbc4));
    _bb1.reset(new cosmo::RsdCorrelationFunction(bb10,bb12,bb14));
    _bb2.reset(new cosmo::RsdCorrelationFunction(bb20,bb22,bb24));
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

double local::BaoCorrelationModel::evaluate(double r, double mu, double z,
std::vector<double> const &params) const {
    double alpha(params[0]), bb(params[1]), beta(params[2]), ampl(params[3]), scale(params[4]);
    double bias = bb/(1+beta);
    double xio(params[5]), a0(params[6]), a1(params[7]), a2(params[8]);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // Apply redshift-space distortion to each model component.
    _fid->setDistortion(beta);
    _nw->setDistortion(beta);
    _bbc->setDistortion(beta);
    _bb1->setDistortion(beta);
    _bb2->setDistortion(beta);
    // Calculate the peak contribution with scaled radius.
    double fid((*_fid)(r*scale,mu)), nw((*_nw)(r*scale,mu)); // scale cancels in mu
    double peak = ampl*(fid-nw);
    // Calculate the additional broadband contribution with no radius scaling.
    double bbc((*_bbc)(r,mu)), bb0((*_nw)(r,mu)), bb1((*_bb1)(r,mu)), bb2((*_bb2)(r,mu));
    double broadband = xio*bbc + (1+a0)*bb0 + a1*bb1 + a2*bb2;
    // Combine the peak and broadband components, with bias and redshift evolution.
    return bias*bias*zfactor*(peak + broadband);
}

double local::BaoCorrelationModel::evaluate(double r, double z,
std::vector<double> const &params) const {
    return 0;
}

std::vector<double> local::BaoCorrelationModel::evaluateMultipoles(double r,
std::vector<double> const &params) const {
    std::vector<double> pcopy(params);
    pcopy[0] = 0; // alpha = 0 to fix z = zref
    //pcopy[1] = 1; // fix bias = 1

    pcopy[2] = 0; // mu=0, beta = 0 gives xi = xi0
    double xia = evaluate(r,0,0,pcopy);
    
    pcopy[2] = +1; // mu=0, beta = +1 gives xi = (28/15)xi0 - (20/21)xi2 +(3/35)xi4
    double xib = evaluate(r,0,0,pcopy);

    pcopy[2] = -1; // mu=0, beta = -1 gives xi = (8/15)xi0 + (8/21)xi2 + (3/35)xi4
    double xic = evaluate(r,0,0,pcopy);
    
    // Solve for xi0,xi2,xi4
    std::vector<double> xi(3);
    xi[0] = xia;
    xi[1] = xia - (3./4.)*(xib - xic);
    xi[2] = (-32*xia + 10*xib + 25*xic)/3;
    return xi;
}
