// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "cosmo/RsdCorrelationFunction.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "likely/Interpolator.h"
#include "likely/function.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &broadbandName, double zref, double initialAmp, double initialScale,
    bool fixAlpha, bool fixLinear, bool fixBao, bool fixScale, bool noBBand)
: AbsCorrelationModel(), _zref(zref)
{
    // Define our parameters. The order here determines the order of elements in our
    // parameter vector for our evaluate(...) methods.
    defineParameter("alpha",3.8,0.3, fixAlpha || fixLinear);
    defineParameter("beta",1.0,0.1, fixLinear || (!fixBao && !noBBand));
    defineParameter("(1+beta)*bias",-0.34,0.03, fixLinear || (!fixBao && !noBBand));
    defineParameter("BAO amplitude", initialAmp,0.15,fixBao);
    defineParameter("BAO scale", initialScale,0.02,fixBao || fixScale);
    defineParameter("BBand xio",0,0.001, noBBand);
    defineParameter("BBand a0",0,0.2, noBBand);
    defineParameter("BBand a1",0,2, noBBand);
    defineParameter("BBand a2",0,2, noBBand);
    // Load the interpolation data we will use for each multipole of each model.
    std::string root(modelrootName);
    if(0 < root.size() && root[root.size()-1] != '/') root += '/';
    boost::format fileName("%s%s.%d.dat"),bbandName("%s%s%c.%d.dat");
    std::string method("cspline");
    cosmo::CorrelationFunctionPtr
        fid0 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % fiducialName % 0),method)),
        fid2 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % fiducialName % 2),method)),
        fid4 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % fiducialName % 4),method)),
        nw0 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 0),method)),
        nw2 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 2),method)),
        nw4 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 4),method)),
        bbc0 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % 'c' % 0),method)),
        bbc2 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % 'c' % 2),method)),
        bbc4 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % 'c' % 4),method)),
        bb10 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % '1' % 0),method)),
        bb12 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % '1' % 2),method)),
        bb14 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % '1' % 4),method)),
        bb20 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % '2' % 0),method)),
        bb22 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % '2' % 2),method)),
        bb24 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(bbandName % root % broadbandName % '2' % 4),method));
    // Create redshift-space distorted correlation function models from the multipole interpolators.
    _fid.reset(new cosmo::RsdCorrelationFunction(fid0,fid2,fid4));
    _nw.reset(new cosmo::RsdCorrelationFunction(nw0,nw2,nw4));
    _bbc.reset(new cosmo::RsdCorrelationFunction(bbc0,bbc2,bbc4));
    _bb1.reset(new cosmo::RsdCorrelationFunction(bb10,bb12,bb14));
    _bb2.reset(new cosmo::RsdCorrelationFunction(bb20,bb22,bb24));
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

double local::BaoCorrelationModel::evaluate(double r, double mu, double z,
std::vector<double> const &params) const {
    double alpha(params[0]), beta(params[1]), bb(params[2]), ampl(params[3]), scale(params[4]);
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
    double peak(0);
    if(ampl != 0) {
        double fid((*_fid)(r*scale,mu)), nw((*_nw)(r*scale,mu)); // scale cancels in mu
        peak = ampl*(fid-nw);
    }
    // Calculate the additional broadband contribution with no radius scaling.
    double broadband(0);
    if(xio != 0) broadband += xio*(*_bbc)(r,mu);
    if(1+a0 != 0) broadband += (1+a0)*(*_nw)(r,mu);
    if(a1 != 0) broadband += a1*(*_bb1)(r,mu);
    if(a2 != 0) broadband += a2*(*_bb2)(r,mu);
    // Combine the peak and broadband components, with bias and redshift evolution.
    return bias*bias*zfactor*(peak + broadband);
}

double local::BaoCorrelationModel::evaluate(double r, cosmo::Multipole multipole, double z,
std::vector<double> const &params) const {
    double alpha(params[0]), beta(params[1]), bb(params[2]), ampl(params[3]), scale(params[4]);
    double bias = bb/(1+beta);
    double xio(params[5]), a0(params[6]), a1(params[7]), a2(params[8]);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // No need to apply redshift-space distortion to each model component since we are
    // working in undistorted multipoles here.

    // Calculate the peak contribution with scaled radius.
    double peak(0);
    if(ampl != 0) {
        double fid((*_fid)(r*scale,multipole)), nw((*_nw)(r*scale,multipole));
        peak = ampl*(fid-nw);
    }
    // Calculate the additional broadband contribution with no radius scaling.
    double broadband(0);
    if(xio != 0) broadband += xio*(*_bbc)(r,multipole);
    if(1+a0 != 0) broadband += (1+a0)*(*_nw)(r,multipole);
    if(a1 != 0) broadband += a1*(*_bb1)(r,multipole);
    if(a2 != 0) broadband += a2*(*_bb2)(r,multipole);
    // Combine the peak and broadband components, with bias and redshift evolution.
    return bias*bias*zfactor*(peak + broadband);
}

/*
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
*/