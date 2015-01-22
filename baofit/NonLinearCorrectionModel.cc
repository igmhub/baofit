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
    double zArray[]      = { 2.200, 2.400, 2.600, 2.800, 3.000 };
    double qnlArray[]    = { 0.036, 0.035, 0.040, 0.036, 0.040 };
    double kvArray[]     = { 0.756, 0.787, 0.830, 0.891, 0.980 };
    double avArray[]     = { 0.454, 0.478, 0.490, 0.502, 0.500 };
    double bvArray[]     = { 1.520, 1.540, 1.560, 1.570, 1.560 };
    double kpArray[]     = { 12.50, 12.40, 12.00, 11.40, 10.70 };
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
    // Non-linear correction model of http://arxiv.org/abs/xxxx.xxxxx
    if(_nlCorrection) {
        //double qnl = (*_qnlInterpolator)(z);
        //double kv = (*_kvInterpolator)(z);
        //double av = (*_avInterpolator)(z);
        //double bv = (*_bvInterpolator)(z);
        //double kp = (*_kpInterpolator)(z);
        double qnl(0.036), kv(0.756), av(0.454), bv(1.52), kp(12.5);
        pk = pk*(0.88/_sigma8)*(0.88/_sigma8)*redshiftEvolution(1.,-2.,z,_zref);
        growth = k*k*k*pk*qnl;
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
