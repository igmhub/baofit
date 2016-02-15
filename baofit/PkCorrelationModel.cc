// Created 31-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/PkCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

#include "gsl/gsl_sf_expint.h"

#include <fstream>
#include <cmath>

namespace local = baofit;

local::PkCorrelationModel::PkCorrelationModel(std::string const &modelrootName, std::string const &nowigglesName,
double klo, double khi, int nk, int splineOrder, bool independentMultipoles, double zref,
bool crossCorrelation)
: AbsCorrelationModel("P(ell,k) Correlation Model"), _klo(klo), _nk(nk), _splineOrder(splineOrder),
_independentMultipoles(independentMultipoles)
{
    // Check inputs.
    if(klo >= khi) throw RuntimeError("PkCorrelationModel: expected khi > klo.");
    if(nk - splineOrder < 2) throw RuntimeError("PkCorrelationModel: expected nk - splineOrder >= 2.");
    if(zref < 0) throw RuntimeError("PkCorrelationModel: expected zref >= 0.");
    if(splineOrder != 0 && splineOrder != 1 && splineOrder != 3) {
        throw RuntimeError("PkCorrelationModel: only splineOrder = 0,1,3 are implemented so far.");
    }
    // Precompute useful quantities
    _dk = (khi - klo)/(nk - 1);
    _dk2 = _dk*_dk;
    _dk3 = _dk2*_dk;
    _dk4 = _dk2*_dk2;
    _rsave = -1;
    double pi(4*std::atan(1));
    _twopisq = 2*pi*pi;
    _sinInt.resize(nk);
    _sin.resize(nk);
    _cos.resize(nk);
    // Linear bias parameters
    _indexBase = 1 + _defineLinearBiasParameters(zref,crossCorrelation);
    // B-spline coefficients for each multipole.
    boost::format name("Pk b-%d-%d");
    for(int ell = 0; ell <= 4; ell += 2) {
        if(!_independentMultipoles && ell > 0) break;
        for(int j = 0; j <= nk - splineOrder - 2; ++j) {
            defineParameter(boost::str(name % ell % j),0,0.1);
        }
    }
    // Load the interpolation data for the specified no-wiggles model.
    std::string root(modelrootName);
    if(0 < root.size() && root[root.size()-1] != '/') root += '/';
    boost::format powerName("%s%s_matterpower.dat"),fileName("%s%s.%d.dat");
    std::string method("cspline");
    try {
        _nwPower = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(powerName % root % nowigglesName),method));
        _nw0 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 0),method));
        _nw2 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 2),method));
        _nw4 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 4),method));
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("PkCorrelationModel: error while reading model interpolation data.");
    }
    {
        /**
        // Debug output to compare with Mathematica results
        int j = 2;
        double kj = _klo + _dk*j, r = 100;
        _fillCache(r);
        std::cout << "EE0" << _splineOrder << '[' << kj << ',' << _dk << ",100] = "
            << _getE(j,r,cosmo::Monopole) << std::endl;
        std::cout << "EE2" << _splineOrder << '[' << kj << ',' << _dk << ",100] = "
            << _getE(j,r,cosmo::Quadrupole) << std::endl;
        std::cout << "EE4" << _splineOrder << '[' << kj << ',' << _dk << ",100] = "
            << _getE(j,r,cosmo::Hexadecapole) << std::endl;
        **/
    }
}

local::PkCorrelationModel::~PkCorrelationModel() { }

void local::PkCorrelationModel::_fillCache(double r) const {
    if(r == _rsave) return;
    for(int j = 0; j < _nk; ++j) {
        double kj = _klo + _dk*j;
        _sin[j] = std::sin(kj*r);
        _cos[j] = std::cos(kj*r);
        _sinInt[j] = gsl_sf_Si(kj*r);
    }
    _rsave = r;
}

double local::PkCorrelationModel::_getB(int j, double k) const {
    double kj = _klo + _dk*j;
    if(k <= kj) return 0;
    double t = (k - kj)/((_splineOrder+1)*_dk);
    if(t > 1) return 0;
    // Use symmetry to put 0 < t <= 1/2
    if(t > 0.5) t = 1-t;
    double result(0);
    if(_splineOrder == 3) {
        if(t < 0.25) {
            result = (32./3.)*t*t*t;
        }
        else { // t <= 1/2
            result = (2./3.) - 8*t*(1 - 4*t*(1 - t));
        }
    }
    else if(_splineOrder == 1) {
        result = 2*t;
    }
    else { // _splineOrder == 0
        return 1;
    }
    return result;
}

double local::PkCorrelationModel::_getE(int j, double r, cosmo::Multipole multipole) const {
    double kj = _klo + _dk*j;
    double r2(r*r),r3(r2*r),r5(r3*r2),kj2(kj*kj),kj3(kj2*kj),kj4(kj2*kj2),kj5(kj3*kj2),kj6(kj3*kj3);
    if(_splineOrder == 3) {
        if(multipole == cosmo::Monopole) {
            return (_sin[j] - 4*_sin[j+1] + 6*_sin[j+2] - 4*_sin[j+3] + _sin[j+4])/(_dk3*r5);
        }
        else if(multipole == cosmo::Quadrupole) {
            return -(3*kj*r*_cos[j] - 12*_dk*r*_cos[j+1] - 12*kj*r*_cos[j+1] + 36*_dk*r*_cos[j+2] + 18*kj*r*_cos[j+2] - 
                  36*_dk*r*_cos[j+3] - 12*kj*r*_cos[j+3] + 12*_dk*r*_cos[j+4] + 3*kj*r*_cos[j+4] + 5*_sin[j] - 
                  20*_sin[j+1] + 30*_sin[j+2] - 20*_sin[j+3] + 5*_sin[j+4] +
                  3*kj2*r2*_sinInt[j] - 
                  12*(_dk2+2*_dk*kj+kj2)*r2*_sinInt[j+1] +
                  72*_dk2*r2*_sinInt[j+2] + 72*_dk*kj*r2*_sinInt[j+2] + 18*kj2*r2*_sinInt[j+2] - 
                  108*_dk2*r2*_sinInt[j+3] - 72*_dk*kj*r2*_sinInt[j+3] - 12*kj2*r2*_sinInt[j+3] +
                  48*_dk2*r2*_sinInt[j+4] + 24*_dk*kj*r2*_sinInt[j+4] + 3*kj2*r2*_sinInt[j+4]
                  )/(2.*_dk3*r5);
        }
        else { // Hexadecapole
            return -(15*kj*r*_cos[j] - 60*_dk*r*_cos[j+1] - 60*kj*r*_cos[j+1] + 180*_dk*r*_cos[j+2] + 90*kj*r*_cos[j+2] - 
                  180*_dk*r*_cos[j+3] - 60*kj*r*_cos[j+3] + 60*_dk*r*_cos[j+4] + 15*kj*r*_cos[j+4] + 11*_sin[j] - 
                  44*_sin[j+1] + 66*_sin[j+2] - 44*_sin[j+3] + 11*_sin[j+4] + 
                  5*(14 + 3*kj2*r2)*_sinInt[j] -
                  20*(14 + 3*_dk2*r2 + 6*_dk*kj*r2 + 3*kj2*r2)*_sinInt[j+1] +
                  420*_sinInt[j+2] + 360*_dk2*r2*_sinInt[j+2] + 360*_dk*kj*r2*_sinInt[j+2] + 90*kj2*r2*_sinInt[j+2] -
                  280*_sinInt[j+3] - 540*_dk2*r2*_sinInt[j+3] - 360*_dk*kj*r2*_sinInt[j+3] - 60*kj2*r2*_sinInt[j+3] +
                  70*_sinInt[j+4] + 240*_dk2*r2*_sinInt[j+4] + 120*_dk*kj*r2*_sinInt[j+4] + 15*kj2*r2*_sinInt[j+4]
                  )/(4.*_dk3*r5);
        }
    }
    else if(_splineOrder == 1) {
        if(multipole == cosmo::Monopole) {
            return -((_sin[j] - 2*_sin[j+1] + _sin[j+2])/(_dk*r3));
        }
        else if(multipole == cosmo::Quadrupole) {
            return (_sin[j] - 2*_sin[j+1] + _sin[j+2] - 3*_sinInt[j] + 6*_sinInt[j+1] - 3*_sinInt[j+2])/(_dk*r3);
        }
        else { // Hexadecapole
            double tmp = 2*_dk2 + 3*_dk*kj + kj2;
            double kjj = kj + _dk, kjjj = kjj + _dk;
            return -(140*_dk4*kj*r*_cos[j] + 420*_dk3*kj2*r*_cos[j] + 455*_dk2*kj3*r*_cos[j] + 210*_dk*kj4*r*_cos[j] + 
                  35*kj5*r*_cos[j] - 280*_dk3*kj2*r*_cos[j+1] - 560*_dk2*kj3*r*_cos[j+1] - 
                  350*_dk*kj4*r*_cos[j+1] - 70*kj5*r*_cos[j+1] + 70*_dk3*kj2*r*_cos[j+2] + 
                  175*_dk2*kj3*r*_cos[j+2] + 140*_dk*kj4*r*_cos[j+2] + 35*kj5*r*_cos[j+2] - 
                  140*_dk4*_sin[j] - 420*_dk3*kj*_sin[j] - 455*_dk2*kj2*_sin[j] - 210*_dk*kj3*_sin[j] - 
                  35*kj4*_sin[j] + 8*_dk4*kj2*r2*_sin[j] + 24*_dk3*kj3*r2*_sin[j] + 
                  26*_dk2*kj4*r2*_sin[j] + 12*_dk*kj5*r2*_sin[j] + 2*kj6*r2*_sin[j] + 
                  280*_dk2*kj2*_sin[j+1] + 280*_dk*kj3*_sin[j+1] + 70*kj4*_sin[j+1] - 
                  16*_dk4*kj2*r2*_sin[j+1] - 48*_dk3*kj3*r2*_sin[j+1] - 
                  52*_dk2*kj4*r2*_sin[j+1] - 24*_dk*kj5*r2*_sin[j+1] - 
                  4*kj6*r2*_sin[j+1] - 35*_dk2*kj2*_sin[j+2] - 70*_dk*kj3*_sin[j+2] - 
                  35*kj4*_sin[j+2] + 8*_dk4*kj2*r2*_sin[j+2] + 
                  24*_dk3*kj3*r2*_sin[j+2] + 26*_dk2*kj4*r2*_sin[j+2] + 
                  12*_dk*kj5*r2*_sin[j+2] + 2*kj6*r2*_sin[j+2] + 
                  15*kj2*tmp*tmp*r2*_sinInt[j] - 
                  30*kj2*tmp*tmp*r2*_sinInt[j+1] + 
                  60*_dk4*kj2*r2*_sinInt[j+2] + 180*_dk3*kj3*r2*_sinInt[j+2] + 
                  195*_dk2*kj4*r2*_sinInt[j+2] + 90*_dk*kj5*r2*_sinInt[j+2] + 
                  15*kj6*r2*_sinInt[j+2])/(2.*_dk*kj2*kjj*kjj*kjjj*kjjj*r5);
        }
    }
    else { // _splineOrder == 0
        if(multipole == cosmo::Monopole) {
            return (_cos[j] - _cos[j+1])/r2;
        }
        else if(multipole == cosmo::Quadrupole) {
            return (-(kj*(_dk + kj)*r*_cos[j]) + kj*(_dk + kj)*r*_cos[j+1] + 3*((_dk + kj)*_sin[j] -
                kj*_sin[j+1]))/(kj*(_dk + kj)*r3);
        }
        else { // Hexadecapole
            double kjj = kj + _dk, kjj3 = kjj*kjj*kjj;
            return (kj*kjj3*r*(-35 + kj2*r2)*_cos[j] - kj3*(_dk + kj)*r*(-35 + _dk2*r2 + 2*_dk*kj*r2 + kj2*r2)*_cos[j+1] - 
                 5*(kjj3*(-7 + 2*kj2*r2)*_sin[j] - kj3*(-7 + 2*_dk2*r2 + 4*_dk*kj*r2 + 2*kj2*r2)*_sin[j+1]))/(kj3*kjj3*r5);
        }
        
    }
}

double local::PkCorrelationModel::_xi(double r, cosmo::Multipole multipole) const {
    // Evaluate the smooth baseline model.
    double xi(0), sign(1);
    int nj = _nk-_splineOrder-1, offset = _indexBase;
    switch(multipole) {
    case cosmo::Monopole:
        xi = (*_nw0)(r);
        break;
    case cosmo::Quadrupole:
        xi = (*_nw2)(r);
        if(_independentMultipoles) offset += nj;
        sign = -1;
        break;
    case cosmo::Hexadecapole:
        xi = (*_nw4)(r);
        if(_independentMultipoles) offset += 2*nj;
        break;
    }
    // Add the splined interpolation.
    for(int j = 0; j < nj; ++j) {
        xi += sign/_twopisq*getParameterValue(offset+j)*_getE(j,r,multipole);
    }
    return xi;
}

double local::PkCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged, int index) const {
    // Cache expensive sine integrals.
    _fillCache(r);
    // Calculate the Legendre weights.
    double muSq(mu*mu);
    double L0(1), L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;
    // Put the pieces together.
    return
        _getNormFactor(cosmo::Monopole,z)*L0*_xi(r,cosmo::Monopole) +
        _getNormFactor(cosmo::Quadrupole,z)*L2*_xi(r,cosmo::Quadrupole) +
        _getNormFactor(cosmo::Hexadecapole,z)*L4*_xi(r,cosmo::Hexadecapole);
}

double local::PkCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged, int index) const {
    // Cache expensive sine integrals.
    _fillCache(r);
    return _getNormFactor(multipole,z)*_xi(r,multipole);
}

double local::PkCorrelationModel::_evaluateKSpace(double k, double mu_k, double pk, double z) const { }

void  local::PkCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
}

void local::PkCorrelationModel::dump(std::string const &dumpName, double kmin, double kmax, int nk,
likely::Parameters const &params, double zref) {
    if(kmax <= kmin || nk <= 2) throw RuntimeError("PkCorrelationModel::dump: bad inputs (kmin,kmax,nk).");
    // Open the requested file.
    std::ofstream out(dumpName.c_str());
    // Load the requested parameters.
    updateParameterValues(params);
    // Loop over k values.
    double dk = (kmax-kmin)/(nk-1);
    int nj = _nk-_splineOrder-1;
    for(int ik = 0; ik < nk; ++ik) {
        double k = kmin + dk*ik;
        // Calculate the smooth model's k*P(ell,k) values.
        double kPk = k*(*_nwPower)(k);
        double kPk0 = _getNormFactor(cosmo::Monopole,zref)*kPk;
        double kPk2 = _getNormFactor(cosmo::Quadrupole,zref)*kPk;
        double kPk4 = _getNormFactor(cosmo::Hexadecapole,zref)*kPk;
        // Calculate the additive k*dP(ell,k) values.
        double kdPk0(0), kdPk2(0), kdPk4(0);
        for(int j = 0; j < nj; ++j) {
            double kdPk = _getB(j,k);
            double b0 = getParameterValue(_indexBase + j);
            double b2(b0),b4(b0);
            if(_independentMultipoles) {
                b2 = getParameterValue(_indexBase + nj + j);
                b4 = getParameterValue(_indexBase + 2*nj + j);
            }
            kdPk0 += _getNormFactor(cosmo::Monopole,zref)*b0*kdPk;
            kdPk2 += _getNormFactor(cosmo::Quadrupole,zref)*b2*kdPk;
            kdPk4 += _getNormFactor(cosmo::Hexadecapole,zref)*b4*kdPk;
        }
        out << k << ' ' << kPk0 << ' ' << kPk2 << ' ' << kPk4 << ' '
            << kdPk0 << ' ' << kdPk2 << ' ' << kdPk4 << std::endl;
    }
    out.close();
}
