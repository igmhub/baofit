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
    if(zref < 0) {
        throw RuntimeError("BaoCorrelationModel: expected zref >= 0.");
    }
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
    defineParameter("BBand mono const",0,1e-4);
    defineParameter("BBand quad const",0,1e-4);
    defineParameter("BBand hexa const",0,1e-4);
    defineParameter("BBand mono 1/r",0,0.01);
    defineParameter("BBand quad 1/r",0,0.02);
    defineParameter("BBand hexa 1/r",0,0.04);
    defineParameter("BBand mono 1/(r*r)",0,0.6);
    defineParameter("BBand quad 1/(r*r)",0,1.2);
    defineParameter("BBand hexa 1/(r*r)",0,2.4);
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
    // Hardcode our scale parameter prior for a first test.
    _scalePriorMin = 0.85;
    _scalePriorMax = 1.15;
    _scalePriorNorm = 2*0.01*0.01; // 2*sigma^2
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

// Evaluates -log(prior(scale)) where the prior is Gaussian for scale < min or scale > max,
// and equal to one for min < scale < max.
double local::BaoCorrelationModel::_evaluatePrior(bool anyChanged) const {
    double scale = getParameterValue("BAO scale");
    if(scale < _scalePriorMin) {
        double diff(scale - _scalePriorMin);
        return diff*diff/_scalePriorNorm;
    }
    if(scale > _scalePriorMax) {
        double diff(scale - _scalePriorMax);
        return diff*diff/_scalePriorNorm;
    }
    return 0;
}

namespace baofit {
    // Define a function object class that simply returns a constant. This could also be done
    // with boost::lambda using (_1 = value), but I don't know how to create a lambda functor
    // on the heap so it can be used with the likely::createFunctionPtr machinery.
    class BaoCorrelationModel::Offset {
    public:
        Offset(double value) : _value(value) { }
        double operator()(double r) { return _value; }
    private:
        double _value;
    };
}

#include "likely/function_impl.h"

template cosmo::CorrelationFunctionPtr likely::createFunctionPtr<local::BaoCorrelationModel::Offset>
    (local::BaoCorrelationModel::OffsetPtr pimpl);

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
    double c0 = getParameterValue("BBand c0");
    double c2 = getParameterValue("BBand c2");
    double c4 = getParameterValue("BBand c4");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    // Calculate redshift evolution factor.
    double zfactor = std::pow((1+z)/(1+_zref),alpha);
    // Build a model with xi(ell=0,2,4) = c(ell).
    cosmo::RsdCorrelationFunction offsetsModel(
        likely::createFunctionPtr(OffsetPtr(new Offset(c0))),
        likely::createFunctionPtr(OffsetPtr(new Offset(c2))),
        likely::createFunctionPtr(OffsetPtr(new Offset(c4))));
    // Apply redshift-space distortion to each model component.
    _fid->setDistortion(beta);
    _nw->setDistortion(beta);
    _bbc->setDistortion(beta);
    _bb1->setDistortion(beta);
    _bb2->setDistortion(beta);
    offsetsModel.setDistortion(beta);
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
    // Calculate any additional DC offset contributions.
    double offsets = offsetsModel(r,mu);
    std::cout << "offsets " << c0 << ' ' << c2 << ' ' << c4 << " => " << offsets << std::endl;
    // Combine the peak and broadband components, with bias and redshift evolution.
    return bias*bias*zfactor*(peak + broadband + offsets);
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
    double rsdScale, offset, rsq(r*r);
    if(multipole == cosmo::Hexadecapole) {
        rsdScale = (8./35.)*beta*beta;
        offset = getParameterValue("BBand hexa const") + getParameterValue("BBand hexa 1/r")/r
            + getParameterValue("BBand hexa 1/(r*r)")/(r*r);
    }
    else if(multipole == cosmo::Quadrupole) {
        rsdScale = 4*beta*((1./3.) + beta/7.);
        offset = getParameterValue("BBand quad const") + getParameterValue("BBand quad 1/r")/r
            + getParameterValue("BBand quad 1/(r*r)")/(r*r);
    }
    else {
        rsdScale = 1 + beta*((2./3.) + beta/5.);
        offset = getParameterValue("BBand mono const") + getParameterValue("BBand mono 1/r")/r
            + getParameterValue("BBand mono 1/(r*r)")/(r*r);
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
    return bias*bias*zfactor*rsdScale*(peak + broadband + offset);
}

void  local::BaoCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
}
