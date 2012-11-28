// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"
#include "baofit/RuntimeError.h"
#include "baofit/BroadbandModel.h"

#include "cosmo/RsdCorrelationFunction.h"
#include "cosmo/TransferFunctionPowerSpectrum.h"
#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &bbandAdd, std::string const &bbandMul, double bbandR0,
    double zref, bool anisotropic, bool decoupled)
: AbsCorrelationModel("BAO Correlation Model"), _anisotropic(anisotropic), _decoupled(decoupled)
{
    // Linear bias parameters
    _defineLinearBiasParameters(zref);
    // BAO peak parameters
    defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO alpha-iso",1,0.02);
    defineParameter("BAO alpha-parallel",1,0.1);
    defineParameter("BAO alpha-perp",1,0.1);
    defineParameter("gamma-scale",0,0.5);
    // Load the interpolation data we will use for each multipole of each model.
    std::string root(modelrootName);
    if(0 < root.size() && root[root.size()-1] != '/') root += '/';
    boost::format fileName("%s%s.%d.dat");
    std::string method("cspline");
    try {
        _fid0 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % fiducialName % 0),method));
        _fid2 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % fiducialName % 2),method));
        _fid4 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % fiducialName % 4),method));
        _nw0 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 0),method));
        _nw2 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 2),method));
        _nw4 = likely::createFunctionPtr(likely::createInterpolator(
            boost::str(fileName % root % nowigglesName % 4),method));
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("BaoCorrelationModel: error while reading model interpolation data.");
    }
    // Define our broadband distortion models, if any.
    if(bbandAdd.length() > 0) {
        _distortAdd.reset(new baofit::BroadbandModel("Additive broadband distortion",
            "dist add",bbandAdd,bbandR0,zref,this));
    }
    if(bbandMul.length() > 0) {
        _distortMul.reset(new baofit::BroadbandModel("Multiplicative broadband distortion",
            "dist mul",bbandMul,bbandR0,zref,this));
    }
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

double local::BaoCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {

    // Lookup parameter values by name.
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO alpha-iso");
    double scale_parallel = getParameterValue("BAO alpha-parallel");
    double scale_perp = getParameterValue("BAO alpha-perp");
    double gamma_scale = getParameterValue("gamma-scale");

    // Calculate redshift evolution of the scale parameters.
    scale = _redshiftEvolution(scale,gamma_scale,z);
    scale_parallel = _redshiftEvolution(scale_parallel,gamma_scale,z);
    scale_perp = _redshiftEvolution(scale_perp,gamma_scale,z);

    // Transform (r,mu) to (rBAO,muBAO) using the scale parameters.
    double rBAO, muBAO;
    if(_anisotropic) {
        double ap1(scale_parallel);
        double bp1(scale_perp);
        double musq(mu*mu);
        // Exact (r,mu) transformation
        double rscale = std::sqrt(ap1*ap1*musq + (1-musq)*bp1*bp1);
        rBAO = r*rscale;
        muBAO = mu*ap1/rscale;
        // Linear approximation, equivalent to multipole model below
        /*
        rBAO = r*(1 + (ap1-1)*musq + (bp1-1)*(1-musq));
        muBAO = mu*(1 + (ap1-bp1)*(1-musq));
        */
    }
    else {
        rBAO = r*scale;
        muBAO = mu;
    }

    // Calculate the cosmological prediction.
    double norm0 = _getNormFactor(cosmo::Monopole,z), norm2 = _getNormFactor(cosmo::Quadrupole,z),
        norm4 = _getNormFactor(cosmo::Hexadecapole,z);
    double muSq(muBAO*muBAO);
    double L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;
    double fid = norm0*(*_fid0)(rBAO) + norm2*L2*(*_fid2)(rBAO) + norm4*L4*(*_fid4)(rBAO);
    double nw = norm0*(*_nw0)(rBAO) + norm2*L2*(*_nw2)(rBAO) + norm4*L4*(*_nw4)(rBAO);
    double peak = ampl*(fid-nw);
    double smooth = nw;
    if(_decoupled) {
        // Recalculate the smooth cosmological prediction using (r,mu) instead of (rBAO,muBAO)
        double muSq(mu*mu);
        double L2 = (3*muSq - 1)/2., L4 = (35*muSq*muSq - 30*muSq + 3)/8.;    
        smooth = norm0*(*_nw0)(r) + norm2*L2*(*_nw2)(r) + norm4*L4*(*_nw4)(r);
    }
    double xi = peak + smooth;
    
    // Add broadband distortions, if any.
    if(_distortMul) xi *= 1 + _distortMul->_evaluate(r,mu,z,anyChanged);
    if(_distortAdd) {
        double distortion = _distortAdd->_evaluate(r,mu,z,anyChanged);
        // The additive distortion is multiplied by ((1+z)/(1+z0))^gamma_bias
        double gamma_bias = getParameterValue("gamma-bias");
        xi += _redshiftEvolution(distortion,gamma_bias,z);
    }

    return xi;
}

double local::BaoCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    /**
    double beta = getParameterValue("beta");
    double bb = getParameterValue("(1+beta)*bias");
    double gamma_bias = getParameterValue("gamma-bias");
    double gamma_beta = getParameterValue("gamma-beta");
    double ampl = getParameterValue("BAO amplitude");
    double scale = getParameterValue("BAO alpha-iso");
    double scale_parallel = getParameterValue("BAO alpha-parallel");
    double scale_perp = getParameterValue("BAO alpha-perp");
    double gamma_scale = getParameterValue("gamma-scale");
    // Calculate bias(zref) from beta(zref) and bb(zref).
    double bias = bb/(1+beta);
    // Calculate redshift evolution.
    double zratio((1+z)/(1+_zref));
    double zfactor = std::pow(zratio,gamma_bias);
    double scaleFactor = std::pow(zratio,gamma_scale);
    scale *= scaleFactor;
    scale_parallel *= scaleFactor;
    scale_perp *= scaleFactor;
    beta *= std::pow(zratio,gamma_beta);
    // Calculate the redshift-space distortion scale factor for this multipole.
    double rsdScale, bband2, rsq(r*r);
    if(multipole == cosmo::Hexadecapole) {
        rsdScale = (8./35.)*beta*beta;
        bband2 = getParameterValue("BBand2 hexa const") + getParameterValue("BBand2 hexa 1/r")/r
            + getParameterValue("BBand2 hexa 1/(r*r)")/(r*r);
    }
    else if(multipole == cosmo::Quadrupole) {
        rsdScale = 4*beta*((1./3.) + beta/7.);
        bband2 = getParameterValue("BBand2 quad const") + getParameterValue("BBand2 quad 1/r")/r
            + getParameterValue("BBand2 quad 1/(r*r)")/(r*r);
    }
    else {
        rsdScale = 1 + beta*((2./3.) + beta/5.);
        bband2 = getParameterValue("BBand2 mono const") + getParameterValue("BBand2 mono 1/r")/r
            + getParameterValue("BBand2 mono 1/(r*r)")/(r*r);
    }
    // Calculate the peak contribution with scaled radius.
    double cosmo(0);
    if(ampl != 0) {
        double fid, nw;
        if(_anisotropic) {
            double fid0 = (*_fid)(r,cosmo::Monopole), fid2 = (*_fid)(r,cosmo::Quadrupole), fid4 = (*_fid)(r,cosmo::Hexadecapole);
            double nw0 = (*_nw)(r,cosmo::Monopole), nw2 = (*_nw)(r,cosmo::Quadrupole), nw4 = (*_nw)(r,cosmo::Hexadecapole);

            double dr = 1;
            double fid0p = ((*_fid)(r+dr,cosmo::Monopole) - (*_fid)(r-dr,cosmo::Monopole))/(2*dr);
            double fid2p = ((*_fid)(r+dr,cosmo::Quadrupole) - (*_fid)(r-dr,cosmo::Quadrupole))/(2*dr);
            double fid4p = ((*_fid)(r+dr,cosmo::Hexadecapole) - (*_fid)(r-dr,cosmo::Hexadecapole))/(2*dr);
            double nw0p = ((*_nw)(r+dr,cosmo::Monopole) - (*_nw)(r-dr,cosmo::Monopole))/(2*dr);
            double nw2p = ((*_nw)(r+dr,cosmo::Quadrupole) - (*_nw)(r-dr,cosmo::Quadrupole))/(2*dr);
            double nw4p = ((*_nw)(r+dr,cosmo::Hexadecapole) - (*_nw)(r-dr,cosmo::Hexadecapole))/(2*dr);
        
            double a(scale_parallel-1);
            double b(scale_perp-1);

            // !! TODO: add hexadecapole terms below
            switch(multipole) {
            case cosmo::Monopole:
                fid = fid0 + (2./5.)*(a-b)*fid2 + (a+2*b)/3*r*fid0p + (2./15.)*(a-b)*r*fid2p;
                nw = nw0 + (2./5.)*(a-b)*nw2 + (a+2*b)/3*r*nw0p + (2./15.)*(a-b)*r*nw2p;
                break;
            case cosmo::Quadrupole:
                fid = fid2*(1 + (2./7.)*(a-b)) + (2./3.)*(a-b)*r*fid0p + (11*a+10*b)/21*r*fid2p;
                nw = nw2*(1 + (2./7.)*(a-b)) + (2./3.)*(a-b)*r*nw0p + (11*a+10*b)/21*r*nw2p;
                break;
            case cosmo::Hexadecapole:
                //throw RuntimeError("BaoCorrelationModel: anisotropic hexadecapole not implemented yet.");
                fid = nw = 0;
                break;
            }
        }
        else {
            fid = (*_fid)(r*scale,multipole);
            nw = (*_nw)(r*scale,multipole);
        }
        peak = ampl*(fid-nw);
    }
    // Combine the peak and broadband components, with bias and redshift evolution.
    return cosmo + bias*bias*zfactor*rsdScale*bband2;
    **/
    throw RuntimeError("BaoCorrelationModel: multipole evaluation temporarily disabled.");
    return 0;
}

void  local::BaoCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using " << (_anisotropic ? "anisotropic":"isotropic") << " BAO scales." << std::endl;
    out << "Scales apply to BAO peak " << (_decoupled ? "only." : "and cosmological broadband.") << std::endl;
    if(_distortAdd) _distortAdd->printToStream(out,formatSpec);
    if(_distortMul) _distortMul->printToStream(out,formatSpec);
}
