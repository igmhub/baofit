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
double klo, double khi, int nk, int splineOrder, bool independentMultipoles, double zref)
: AbsCorrelationModel("P(ell,k) Correlation Model"), _klo(klo), _nk(nk), _splineOrder(splineOrder),
_independentMultipoles(independentMultipoles), _zref(zref)
{
    // Check inputs.
    if(klo >= khi) throw RuntimeError("PkCorrelationModel: expected khi > klo.");
    if(nk - splineOrder < 2) throw RuntimeError("PkCorrelationModel: expected nk - splineOrder >= 2.");
    if(zref < 0) throw RuntimeError("PkCorrelationModel: expected zref >= 0.");
    if(splineOrder != 1 && splineOrder != 3) {
        throw RuntimeError("PkCorrelationModel: only splineOrder = 1,3 is implemented so far.");
    }
    // Precompute useful quantities
    _dk = (khi - klo)/(nk - 1);
    _dk2 = _dk*_dk;
    _dk3 = _dk2*_dk;
    _dk4 = _dk3*_dk;
    _rsave = -1;
    double pi(4*std::atan(1));
    _twopisq = 2*pi*pi;
    _sinInt.resize(nk);
    _indexBase = -1;
    // Linear bias parameters
    defineParameter("beta",1.4,0.1);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    // Redshift evolution parameters
    defineParameter("gamma-bias",3.8,0.3);
    defineParameter("gamma-beta",0,0.1);
    // Multiplicative broadband distortion factors.
    defineParameter("Pk s-0",1,0.1);
    defineParameter("Pk s-2",1,0.1);
    defineParameter("Pk s-4",1,0.1);
    // B-spline coefficients for each multipole.
    boost::format name("Pk b-%d-%d");
    for(int ell = 0; ell <= 4; ell += 2) {
        if(!_independentMultipoles && ell > 0) break;
        for(int j = 0; j <= nk - splineOrder - 2; ++j) {
            int index = defineParameter(boost::str(name % ell % j),0,0.1);
            if(_indexBase == -1) _indexBase = index;
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
}

local::PkCorrelationModel::~PkCorrelationModel() { }

void local::PkCorrelationModel::_calculateNorm(double z) const {
    // Lookup the linear bias parameters.
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    double biasSq = bias*bias;
    // Calculate redshift evolution of bias and beta.
    double ratio = (1+z)/(1+_zref);
    biasSq *= std::pow(ratio,getParameterValue("gamma-bias"));
    beta *= std::pow(ratio,getParameterValue("gamma-beta"));
    // Calculate the linear bias normalization factors.
    _norm0 = biasSq*(1 + beta*(2./3. + (1./5.)*beta));
    _norm2 = biasSq*beta*(4./3. + (4./7.)*beta);
    _norm4 = biasSq*beta*beta*(8./35.);
    // Apply multiplicative distortion.
    _norm0 *= getParameterValue("Pk s-0");
    _norm2 *= getParameterValue("Pk s-2");
    _norm4 *= getParameterValue("Pk s-4");
}

void local::PkCorrelationModel::_fillSinIntegralCache(double r) const {
    if(r == _rsave) return;
    for(int j = 0; j < _nk; ++j) {
        double kj = _klo + _dk*j;
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
    else { // _splineOrder == 1
        result = 2*t;
    }
    return result;
}

double local::PkCorrelationModel::_getE(int j, double r, cosmo::Multipole multipole) const {
    double kj = _klo + _dk*j;
    double r2(r*r),r3(r2*r),r5(r3*r2),kj2(kj*kj),kj3(kj2*kj),kj4(kj2*kj2),kj5(kj3*kj2),kj6(kj3*kj3);
    if(_splineOrder == 3) {
        if(multipole == cosmo::Monopole) {
            return (std::sin(kj*r) + (6 - 8*std::cos(_dk*r))*std::sin((2*_dk + kj)*r) + std::sin((4*_dk + kj)*r))/(_dk3*r5);
        }
        else if(multipole == cosmo::Quadrupole) {
            return -(3*kj*r*std::cos(kj*r) - 12*_dk*r*std::cos((_dk + kj)*r) - 12*kj*r*std::cos((_dk + kj)*r) + 36*_dk*r*std::cos((2*_dk + kj)*r) + 18*kj*r*std::cos((2*_dk + kj)*r) - 
                  36*_dk*r*std::cos((3*_dk + kj)*r) - 12*kj*r*std::cos((3*_dk + kj)*r) + 12*_dk*r*std::cos((4*_dk + kj)*r) + 3*kj*r*std::cos((4*_dk + kj)*r) + 5*std::sin(kj*r) - 
                  20*std::sin((_dk + kj)*r) + 30*std::sin((2*_dk + kj)*r) - 20*std::sin((3*_dk + kj)*r) + 5*std::sin((4*_dk + kj)*r) +
                  3*kj2*r2*_sinInt[j] - 
                  12*(_dk2+2*_dk*kj+kj2)*r2*_sinInt[j+1] +
                  72*_dk2*r2*_sinInt[j+2] + 72*_dk*kj*r2*_sinInt[j+2] + 18*kj2*r2*_sinInt[j+2] - 
                  108*_dk2*r2*_sinInt[j+3] - 72*_dk*kj*r2*_sinInt[j+3] - 12*kj2*r2*_sinInt[j+3] +
                  48*_dk2*r2*_sinInt[j+4] + 24*_dk*kj*r2*_sinInt[j+4] + 3*kj2*r2*_sinInt[j+4]
                  )/(2.*_dk3*r5);
        }
        else { // Hexadecapole
            return -(15*kj*r*std::cos(kj*r) - 60*_dk*r*std::cos((_dk + kj)*r) - 60*kj*r*std::cos((_dk + kj)*r) + 180*_dk*r*std::cos((2*_dk + kj)*r) + 90*kj*r*std::cos((2*_dk + kj)*r) - 
                  180*_dk*r*std::cos((3*_dk + kj)*r) - 60*kj*r*std::cos((3*_dk + kj)*r) + 60*_dk*r*std::cos((4*_dk + kj)*r) + 15*kj*r*std::cos((4*_dk + kj)*r) + 11*std::sin(kj*r) - 
                  44*std::sin((_dk + kj)*r) + 66*std::sin((2*_dk + kj)*r) - 44*std::sin((3*_dk + kj)*r) + 11*std::sin((4*_dk + kj)*r) + 
                  5*(14 + 3*kj2*r2)*_sinInt[j] -
                  20*(14 + 3*_dk2*r2 + 6*_dk*kj*r2 + 3*kj2*r2)*_sinInt[j+1] +
                  420*_sinInt[j+2] + 360*_dk2*r2*_sinInt[j+2] + 360*_dk*kj*r2*_sinInt[j+2] + 90*kj2*r2*_sinInt[j+2] -
                  280*_sinInt[j+3] - 540*_dk2*r2*_sinInt[j+3] - 360*_dk*kj*r2*_sinInt[j+3] - 60*kj2*r2*_sinInt[j+3] +
                  70*_sinInt[j+4] + 240*_dk2*r2*_sinInt[j+4] + 120*_dk*kj*r2*_sinInt[j+4] + 15*kj2*r2*_sinInt[j+4]
                  )/(4.*_dk3*r5);
        }
    }
    else { // _splineOrder == 1
        if(multipole == cosmo::Monopole) {
            return -((std::sin(kj*r) - 2*std::sin((_dk + kj)*r) + std::sin((2*_dk + kj)*r))/(_dk*r3));
        }
        else if(multipole == cosmo::Quadrupole) {
            return (std::sin(kj*r) - 2*std::sin((_dk + kj)*r) + std::sin((2*_dk + kj)*r) - 3*_sinInt[j] + 6*_sinInt[j+1] - 3*_sinInt[j+2])/(_dk*r3);
        }
        else { // Hexadecapole
            double tmp = 2*_dk2 + 3*_dk*kj + kj2;
            double kjj = kj + _dk, kjjj = kjj + _dk;
            return -(140*_dk4*kj*r*std::cos(kj*r) + 420*_dk3*kj2*r*std::cos(kj*r) + 455*_dk2*kj3*r*std::cos(kj*r) + 210*_dk*kj4*r*std::cos(kj*r) + 
                  35*kj5*r*std::cos(kj*r) - 280*_dk3*kj2*r*std::cos((_dk + kj)*r) - 560*_dk2*kj3*r*std::cos((_dk + kj)*r) - 
                  350*_dk*kj4*r*std::cos((_dk + kj)*r) - 70*kj5*r*std::cos((_dk + kj)*r) + 70*_dk3*kj2*r*std::cos((2*_dk + kj)*r) + 
                  175*_dk2*kj3*r*std::cos((2*_dk + kj)*r) + 140*_dk*kj4*r*std::cos((2*_dk + kj)*r) + 35*kj5*r*std::cos((2*_dk + kj)*r) - 
                  140*_dk4*std::sin(kj*r) - 420*_dk3*kj*std::sin(kj*r) - 455*_dk2*kj2*std::sin(kj*r) - 210*_dk*kj3*std::sin(kj*r) - 
                  35*kj4*std::sin(kj*r) + 8*_dk4*kj2*r2*std::sin(kj*r) + 24*_dk3*kj3*r2*std::sin(kj*r) + 
                  26*_dk2*kj4*r2*std::sin(kj*r) + 12*_dk*kj5*r2*std::sin(kj*r) + 2*kj6*r2*std::sin(kj*r) + 
                  280*_dk2*kj2*std::sin((_dk + kj)*r) + 280*_dk*kj3*std::sin((_dk + kj)*r) + 70*kj4*std::sin((_dk + kj)*r) - 
                  16*_dk4*kj2*r2*std::sin((_dk + kj)*r) - 48*_dk3*kj3*r2*std::sin((_dk + kj)*r) - 
                  52*_dk2*kj4*r2*std::sin((_dk + kj)*r) - 24*_dk*kj5*r2*std::sin((_dk + kj)*r) - 
                  4*kj6*r2*std::sin((_dk + kj)*r) - 35*_dk2*kj2*std::sin((2*_dk + kj)*r) - 70*_dk*kj3*std::sin((2*_dk + kj)*r) - 
                  35*kj4*std::sin((2*_dk + kj)*r) + 8*_dk4*kj2*r2*std::sin((2*_dk + kj)*r) + 
                  24*_dk3*kj3*r2*std::sin((2*_dk + kj)*r) + 26*_dk2*kj4*r2*std::sin((2*_dk + kj)*r) + 
                  12*_dk*kj5*r2*std::sin((2*_dk + kj)*r) + 2*kj6*r2*std::sin((2*_dk + kj)*r) + 
                  15*kj2*tmp*tmp*r2*_sinInt[j] - 
                  30*kj2*tmp*tmp*r2*_sinInt[j+1] + 
                  60*_dk4*kj2*r2*_sinInt[j+2] + 180*_dk3*kj3*r2*_sinInt[j+2] + 
                  195*_dk2*kj4*r2*_sinInt[j+2] + 90*_dk*kj5*r2*_sinInt[j+2] + 
                  15*kj6*r2*_sinInt[j+2])/(2.*_dk*kj2*kjj*kjj*kjjj*kjj*r5);
        }
    }
}

double local::PkCorrelationModel::_xi(double r, cosmo::Multipole multipole) const {
    // Evaluate the smooth baseline model.
    double xi(0), sign(1), twopisq();
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

double local::PkCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    // Calculate normalization factors.
    _calculateNorm(z);
    // Cache expensive sine integrals.
    _fillSinIntegralCache(r);
    // Calculate the Legendre weights.
    double muSq(mu*mu);
    double L0(1), L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;
    // Put the pieces together.
    return _norm0*L0*_xi(r,cosmo::Monopole) + _norm2*L2*_xi(r,cosmo::Quadrupole) + _norm4*L4*_xi(r,cosmo::Hexadecapole);
}

double local::PkCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    // Calculate normalization factors.
    _calculateNorm(z);
    // Cache expensive sine integrals.
    _fillSinIntegralCache(r);
    if(multipole == cosmo::Hexadecapole) {
        return _norm4*_xi(r,cosmo::Hexadecapole);
    }
    else if(multipole == cosmo::Quadrupole) {
        return _norm2*_xi(r,cosmo::Quadrupole);
    }
    else { // Monopole
        return _norm0*_xi(r,cosmo::Monopole);
    }
}

void  local::PkCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
}

void local::PkCorrelationModel::dump(std::string const &dumpName, double kmin, double kmax, int nk,
likely::Parameters const &params, double zref) {
    if(kmax <= kmin || nk <= 2) throw RuntimeError("PkCorrelationModel::dump: bad inputs (kmin,kmax,nk).");
    // Open the requested file.
    std::ofstream out(dumpName.c_str());
    // Load the requested parameters.
    updateParameterValues(params);
    // Calculate normalization factors.
    _calculateNorm(zref);
    // Loop over k values.
    double dk = (kmax-kmin)/(nk-1);
    int nj = _nk-_splineOrder-1;
    for(int ik = 0; ik < nk; ++ik) {
        double k = kmin + dk*ik;
        // Calculate the smooth model's k*P(ell,k) values.
        double kPk = k*(*_nwPower)(k);
        double kPk0 = _norm0*kPk;
        double kPk2 = _norm2*kPk;
        double kPk4 = _norm4*kPk;
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
            kdPk0 += _norm0*b0*kdPk;
            kdPk2 += _norm2*b2*kdPk;
            kdPk4 += _norm4*b4*kdPk;
        }
        out << k << ' ' << kPk0 << ' ' << kPk2 << ' ' << kPk4 << ' '
            << kdPk0 << ' ' << kdPk2 << ' ' << kdPk4 << std::endl;
    }
    out.close();
}
