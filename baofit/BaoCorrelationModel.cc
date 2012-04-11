// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"

#include "cosmo/TransferFunctionPowerSpectrum.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &broadbandName, double zref)
: AbsCorrelationModel(), _zref(zref)
{
    boost::format fileName("%s.%d.dat"),bbandName("%s%c.%d.dat");
    cosmo::CorrelationFunctionPtr
        fid0 = cosmo::createFunctionPtr(load(boost::str(fileName % fiducialName % 0))),
        fid2 = cosmo::createFunctionPtr(load(boost::str(fileName % fiducialName % 2))),
        fid4 = cosmo::createFunctionPtr(load(boost::str(fileName % fiducialName % 4))),
        nw0 = cosmo::createFunctionPtr(load(boost::str(fileName % nowigglesName % 0))),
        nw2 = cosmo::createFunctionPtr(load(boost::str(fileName % nowigglesName % 2))),
        nw4 = cosmo::createFunctionPtr(load(boost::str(fileName % nowigglesName % 4))),
        bbc0 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % 'c' % 0))),
        bbc2 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % 'c' % 2))),
        bbc4 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % 'c' % 4))),
        bb10 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '1' % 0))),
        bb12 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '1' % 2))),
        bb14 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '1' % 4))),
        bb20 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '2' % 0))),
        bb22 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '2' % 2))),
        bb24 = cosmo::createFunctionPtr(load(boost::str(bbandName % broadbandName % '2' % 4)));
    _fid.reset(new cosmo::RsdCorrelationFunction(fid0,fid2,fid4));
    _nw.reset(new cosmo::RsdCorrelationFunction(nw0,nw2,nw4));
    _bbc.reset(new cosmo::RsdCorrelationFunction(bbc0,bbc2,bbc4));
    _bb1.reset(new cosmo::RsdCorrelationFunction(bb10,bb12,bb14));
    _bb2.reset(new cosmo::RsdCorrelationFunction(bb20,bb22,bb24));
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

double local::BaoCorrelationModel::evaluate(double r, double mu, double z, std::vector<double> params) const {
    double alpha(params[0]), bias(params[1]), beta(params[2]), ampl(params[3]), scale(params[4]);
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

double local::BaoCorrelationModel::evaluate(double r, double z, std::vector<double> params) const {
    return 0;
}
