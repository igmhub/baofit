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
    std::string const &broadbandName, double zref)
: AbsCorrelationModel("BAO Correlation Model"), _zref(zref)
{
    // Define our parameters. The order here determines the order of elements in our
    // parameter vector for our evaluate(...) methods.
    defineParameter("alpha",3.8,0.3);
    defineParameter("beta",1.0,0.1);
    defineParameter("(1+beta)*bias",-0.34,0.03);
    defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO scale",1,0.02);
    defineParameter("BBand xio",0,0.001);
    defineParameter("BBand a0",0,0.2);
    defineParameter("BBand a1",0,2);
    defineParameter("BBand a2",0,2);
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

double local::BaoCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    double alpha = getParameterValue("alpha");
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO scale");
    double xio = getParameterValue("BBand xio");
    double a0 = getParameterValue("BBand a0");
    double a1 = getParameterValue("BBand a1");
    double a2 = getParameterValue("BBand a2");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
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

double local::BaoCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    double alpha = getParameterValue("alpha");
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO scale");
    double xio = getParameterValue("BBand xio");
    double a0 = getParameterValue("BBand a0");
    double a1 = getParameterValue("BBand a1");
    double a2 = getParameterValue("BBand a2");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // Calculate the redshift-space distortion scale factor for this multipole.
    double rsdScale;
    if(multipole == cosmo::Hexadecapole) {
        rsdScale = (8./35.)*beta*beta;
    }
    else if(multipole == cosmo::Quadrupole) {
        rsdScale = 4*beta*((1./3.) + beta/7.);
    }
    else {
        rsdScale = 1 + beta*((2./3.) + beta/5.);
    }
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
    return bias*bias*zfactor*rsdScale*(peak + broadband);
}

void  local::BaoCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
}
