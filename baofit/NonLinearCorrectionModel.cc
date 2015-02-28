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
    // Interpolation data from Table 6 of http://arxiv.org/abs/xxxx.xxxxx
    double zArray[]      = { 2.2000, 2.4000, 2.6000, 2.8000, 3.0000 };
    double q1Array[]     = { 0.2000, 0.2300, 0.2600, 0.3700, 0.5000 };
    double q2Array[]     = { 0.2500, 0.2600, 0.2600, 0.2100, 0.1400 };
    double kvArray[]     = { 0.0034, 0.0131, 0.0356, 0.1310, 0.3120 };
    double avArray[]     = { 0.1300, 0.1660, 0.2050, 0.2770, 0.3480 };
    double bvArray[]     = { 1.5100, 1.5300, 1.5500, 1.5500, 1.5400 };
    double kpArray[]     = { 9.4200, 9.7900, 10.000, 10.400, 10.500 };
    std::vector<double> z, q1, q2, kv, av, bv, kp;
    z.assign (zArray,zArray+5);
    q1.assign (q1Array,q1Array+5);
    q2.assign (q2Array,q2Array+5);
    kv.assign (kvArray,kvArray+5);
    av.assign (avArray,avArray+5);
    bv.assign (bvArray,bvArray+5);
    kp.assign (kpArray,kpArray+5);
    _q1Interpolator.reset(new likely::Interpolator(z,q1,"linear"));
    _q2Interpolator.reset(new likely::Interpolator(z,q2,"linear"));
    _kvInterpolator.reset(new likely::Interpolator(z,kv,"linear"));
    _avInterpolator.reset(new likely::Interpolator(z,av,"linear"));
    _bvInterpolator.reset(new likely::Interpolator(z,bv,"linear"));
    _kpInterpolator.reset(new likely::Interpolator(z,kp,"linear"));
}

local::NonLinearCorrectionModel::~NonLinearCorrectionModel() { }

double local::NonLinearCorrectionModel::_evaluateNLCorrection(double k, double mu_k, double pk, double z) const {
    double growth, pecvelocity, pressure, nonlinearcorr;
    // Non-linear correction model of http://arxiv.org/abs/xxxx.xxxxx
    if(_nlCorrection) {
        double q1 = (*_q1Interpolator)(z);
        double q2 = (*_q2Interpolator)(z);
        double kv = (*_kvInterpolator)(z);
        double av = (*_avInterpolator)(z);
        double bv = (*_bvInterpolator)(z);
        double kp = (*_kpInterpolator)(z);
        double sigma8Sim(0.878);
        double pi(4*std::atan(1));
        pk = pk*(sigma8Sim/_sigma8)*(sigma8Sim/_sigma8)*redshiftEvolution(1.,-2.,z,_zref);
        double dk = k*k*k*pk/(2*pi*pi);
        growth = q1*dk + q2*dk*dk;
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
