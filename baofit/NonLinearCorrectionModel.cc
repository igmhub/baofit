// Created 14-Jan-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "baofit/NonLinearCorrectionModel.h"
#include "baofit/AbsCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/Interpolator.h"

#include <cmath>

namespace local = baofit;

local::NonLinearCorrectionModel::NonLinearCorrectionModel(double zref, double sigma8, bool nlCorrection, bool nlCorrectionAlt)
: _zref(zref), _sigma8(sigma8), _nlCorrection(nlCorrection), _nlCorrectionAlt(nlCorrectionAlt)
{
    if(nlCorrection) _initialize();
}

void local::NonLinearCorrectionModel::_initialize() {
    // Interpolation data from Table 7 of http://arxiv.org/abs/1506.04519
    double zArray[]      = { 2.2000, 2.4000, 2.6000, 2.8000, 3.0000 };
    double qnlArray[]    = { 0.8670, 0.8510, 0.7810, 0.7730, 0.7920 };
    double kvArray[]     = { 1.1200, 1.1122, 1.2570, 1.2765, 1.2928 };
    double avArray[]     = { 0.5140, 0.5480, 0.6110, 0.6080, 0.5780 };
    double bvArray[]     = { 1.6000, 1.6100, 1.6400, 1.6500, 1.6300 };
    double kpArray[]     = { 19.400, 19.500, 21.100, 19.200, 17.100 };
    std::vector<double> z, qnl, kv, av, bv, kp;
    z.assign (zArray,zArray+5);
    qnl.assign (qnlArray,qnlArray+5);
    kv.assign (kvArray,kvArray+5);
    av.assign (avArray,avArray+5);
    bv.assign (bvArray,bvArray+5);
    kp.assign (kpArray,kpArray+5);
    _qnlInterpolator.reset(new likely::Interpolator(z,qnl,"linear"));
    _kvInterpolator.reset(new likely::Interpolator(z,kv,"linear"));
    _avInterpolator.reset(new likely::Interpolator(z,av,"linear"));
    _bvInterpolator.reset(new likely::Interpolator(z,bv,"linear"));
    _kpInterpolator.reset(new likely::Interpolator(z,kp,"linear"));
}

local::NonLinearCorrectionModel::~NonLinearCorrectionModel() { }

double local::NonLinearCorrectionModel::_evaluateNLCorrection(double k, double mu_k, double pk, double z) const {
    double growth, pecvelocity, pressure, nonlinearcorr;
    // Non-linear correction model of http://arxiv.org/abs/1506.04519
    if(_nlCorrection) {
        double qnl = (*_qnlInterpolator)(z);
        double kv = (*_kvInterpolator)(z);
        double av = (*_avInterpolator)(z);
        double bv = (*_bvInterpolator)(z);
        double kp = (*_kpInterpolator)(z);
        double sigma8Sim(0.8338);
        double pi(4*std::atan(1));
        pk = pk*(sigma8Sim/_sigma8)*(sigma8Sim/_sigma8)*redshiftEvolution(1.,-2.,z,_zref);
        double deltaSq = k*k*k*pk/(2*pi*pi);
        growth = qnl*deltaSq;
        pecvelocity = std::pow(k/kv,av)*std::pow(std::fabs(mu_k),bv);
        pressure = (k/kp)*(k/kp);
        nonlinearcorr = std::exp(growth*(1-pecvelocity)-pressure);
    }
    // Non-linear correction model of http://arxiv.org/abs/astro-ph/0108064
    // Parameter values at z=2.25 taken from Table 1.
    else if(_nlCorrectionAlt) {
        double knl(6.4), pnl(0.569), kpp(15.3), pp(2.01), kv0(1.22), pv(1.5), kvi(0.923), pvi(0.451);
        double kvel = kv0*std::pow(1+k/kvi,pvi);
        growth = std::pow(k/knl,pnl);
        pressure = std::pow(k/kpp,pp);
        pecvelocity = std::pow(std::fabs(k*mu_k)/kvel,pv);
        nonlinearcorr = std::exp(growth-pressure-pecvelocity);
    }
    else {
        nonlinearcorr = 1;
    }
    
    return nonlinearcorr;
}
