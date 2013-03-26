// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoCorrelationModel.h"
#include "baofit/RuntimeError.h"
#include "baofit/BroadbandModel.h"

#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

#include <cmath>

namespace local = baofit;

local::BaoCorrelationModel::BaoCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName,
    std::string const &distAdd, std::string const &distMul, double distR0,
    double zref, bool anisotropic, bool decoupled, bool crossCorrelation)
: AbsCorrelationModel("BAO Correlation Model"), _anisotropic(anisotropic), _decoupled(decoupled)
{
    // Linear bias parameters
    _indexBase = _defineLinearBiasParameters(zref,crossCorrelation);
    // BAO peak parameters (values can be retrieved efficiently as offsets from _indexBase)
    defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO alpha-iso",1,0.02);
    defineParameter("BAO alpha-parallel",1,0.1);
    defineParameter("BAO alpha-perp",1,0.1);
    defineParameter("gamma-scale",0,0.5);
    // quasar radiation parameters 
    defineParameter("Rad strength",0.,0.1); 
    defineParameter("Rad anisotropy",0.,0.1);
    defineParameter("Rad mean free path",200.,10.); // in Mpc/h
    defineParameter("Rad quasar lifetime",10.,0.1); // in Myr
    // by default, the radiation parameters are fixed
    configureFitParameters("fix[Rad*]=0");

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
    if(distAdd.length() > 0) {
        _distortAdd.reset(new baofit::BroadbandModel("Additive broadband distortion",
            "dist add",distAdd,distR0,zref,this));
    }
    if(distMul.length() > 0) {
        _distortMul.reset(new baofit::BroadbandModel("Multiplicative broadband distortion",
            "dist mul",distMul,distR0,zref,this));
    }
}

local::BaoCorrelationModel::~BaoCorrelationModel() { }

double local::BaoCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {

    // Lookup parameter values by name.
    double ampl = getParameterValue(_indexBase + 1); //("BAO amplitude");
    double scale = getParameterValue(_indexBase + 2); //"BAO alpha-iso");
    double scale_parallel = getParameterValue(_indexBase + 3); //("BAO alpha-parallel");
    double scale_perp = getParameterValue(_indexBase + 4); //("BAO alpha-perp");
    double gamma_scale = getParameterValue(_indexBase + 5); //("gamma-scale");

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
    double musq(muBAO*muBAO);
    double L2 = (-1+3*musq)/2., L4 = (3+musq*(-30+35*musq))/8.;
    double fid = norm0*(*_fid0)(rBAO) + norm2*L2*(*_fid2)(rBAO) + norm4*L4*(*_fid4)(rBAO);
    double nw = norm0*(*_nw0)(rBAO) + norm2*L2*(*_nw2)(rBAO) + norm4*L4*(*_nw4)(rBAO);
    double peak = ampl*(fid-nw);
    double smooth = nw;
    if(_decoupled) {
        // Recalculate the smooth cosmological prediction using (r,mu) instead of (rBAO,muBAO)
        double musq(mu*mu);
        double L2 = (-1+3*musq)/2., L4 = (3+musq*(-30+35*musq))/8.;
        smooth = norm0*(*_nw0)(r) + norm2*L2*(*_nw2)(r) + norm4*L4*(*_nw4)(r);
    }
    double xi = peak + smooth;
    
    // Add broadband distortions, if any.
    if(_distortMul) xi *= 1 + _distortMul->_evaluate(r,mu,z,anyChanged);
    if(_distortAdd) {
        double distortion = _distortAdd->_evaluate(r,mu,z,anyChanged);
        // The additive distortion is multiplied by ((1+z)/(1+z0))^gamma_bias
        double gamma_bias = getParameterValue(_indexBase - 1); //("gamma-bias");
        xi += _redshiftEvolution(distortion,gamma_bias,z);
    }

    // Lookup radiation parameters, also value by name.
    double rad_strength = getParameterValue(_indexBase + 6);
    double rad_aniso = getParameterValue(_indexBase + 7);
    double mean_free_path = getParameterValue(_indexBase + 8);
    double quasar_lifetime = getParameterValue(_indexBase + 9);

    // add quasar radiation effects (for cross-correlations only)
    // allways works with decoupled
    if(rad_strength>0 && r>0.){ 
        // isotropical radiation
        double rad = rad_strength/(r*r);
        // attenuation
        rad *= std::exp(-r/mean_free_path);
        // anisotropy
        rad *= (1 - rad_aniso*(1-mu*mu));
        // time effects 
        double ctd = r*(1-mu)/(1+z);
        rad *= std::exp(-ctd/quasar_lifetime);
        xi += rad;
    }

    return xi;
}

void  local::BaoCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using " << (_anisotropic ? "anisotropic":"isotropic") << " BAO scales." << std::endl;
    out << "Scales apply to BAO peak " << (_decoupled ? "only." : "and cosmological broadband.") << std::endl;
}
