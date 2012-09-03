// Created 31-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/PkCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::PkCorrelationModel::PkCorrelationModel(std::string const &modelrootName, std::string const &nowigglesName,
double klo, double khi, int nk, int splineOrder, double zref)
: AbsCorrelationModel("P(ell,k) Correlation Model"), _klo(klo), _nk(nk), _splineOrder(splineOrder), _zref(zref)
{
    // Check inputs.
    if(klo >= khi) throw RuntimeError("PkCorrelationModel: expected khi > klo.");
    if(nk - splineOrder < 2) throw RuntimeError("PkCorrelationModel: expected nk - splineOrder >= 2.");
    if(zref < 0) throw RuntimeError("PkCorrelationModel: expected zref >= 0.");
    // Precompute useful quantities
    _dk = (khi - klo)/(nk - 1);
    double pi(4*std::atan(1));
    _twopisq = 2*pi*pi;
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

double local::PkCorrelationModel::_xi(double r, cosmo::Multipole multipole) const {
    // Evaluate the smooth baseline model.
    double xi(0), sign(1), twopisq();
    boost::format name;
    switch(multipole) {
    case cosmo::Monopole:
        xi = (*_nw0)(r);
        name = boost::format("Pk b-0-%d");
        break;
    case cosmo::Quadrupole:
        xi = (*_nw2)(r);
        name = boost::format("Pk b-2-%d");
        sign = -1;
        break;
    case cosmo::Hexadecapole:
        xi = (*_nw4)(r);
        name = boost::format("Pk b-4-%d");
        break;
    }
    // Add the splined interpolation.
    for(int j = 0; j <= _nk - _splineOrder - 2; ++j) {
        double kj = _klo + _dk*j;
        xi += sign/_twopisq*getParameterValue(boost::str(name % j))*_getE(kj,r,multipole);
    }
    return xi;
}

double local::PkCorrelationModel::_getE(double kj, double r, cosmo::Multipole multipole) const {
    return 0;
}

double local::PkCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    // Calculate normalization factors.
    _calculateNorm(z);
    // Calculate the Legendre weights.
    double muSq(mu*mu);
    double L0(1), L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;
    // Put the pieces together.
    return _norm0*L0*_xi(r,cosmo::Monopole) + _norm2*L2*_xi(r,cosmo::Quadrupole) + _norm4*L4*_xi(r,cosmo::Hexadecapole);
}

double local::PkCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    _calculateNorm(z);
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
