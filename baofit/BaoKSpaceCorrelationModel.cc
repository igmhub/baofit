// Created 27-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BaoKSpaceCorrelationModel.h"
#include "baofit/RuntimeError.h"
#include "baofit/BroadbandModel.h"

#include "likely/Interpolator.h"
#include "likely/function_impl.h"

#include "cosmo/RuntimeError.h"
#include "cosmo/TabulatedPower.h"
#include "cosmo/DistortedPowerCorrelation.h"

#include "boost/format.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/bind.hpp"

#include <cmath>
#include <iostream>

namespace local = baofit;

local::BaoKSpaceCorrelationModel::BaoKSpaceCorrelationModel(std::string const &modelrootName,
    std::string const &fiducialName, std::string const &nowigglesName, double zref,
    double rmin, double rmax, double dilmin, double dilmax, double relerr, double abserr, int ellMax,
    std::string const &distAdd, std::string const &distMul, double distR0,
    bool anisotropic, bool decoupled,  bool nlBroadband, bool crossCorrelation, bool verbose)
: AbsCorrelationModel("BAO k-Space Correlation Model"), _dilmin(dilmin), _dilmax(dilmax),
_anisotropic(anisotropic), _decoupled(decoupled), _nlBroadband(nlBroadband),
_crossCorrelation(crossCorrelation), _verbose(verbose), _nWarnings(0), _maxWarnings(10)
{
    _setZRef(zref);
    // Linear bias parameters
    defineParameter("beta",1.4,0.1);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    defineParameter("gamma-bias",3.8,0.3);
    defineParameter("gamma-beta",0,0.1);
    // Non-linear broadening parameters
    _nlBase = defineParameter("SigmaNL-perp",3.26,0.3);
    defineParameter("1+f",2,0.1);
    // BAO peak parameters
    _baoBase = defineParameter("BAO amplitude",1,0.15);
    defineParameter("BAO alpha-iso",1,0.02);
    defineParameter("BAO alpha-parallel",1,0.1);
    defineParameter("BAO alpha-perp",1,0.1);
    defineParameter("gamma-scale",0,0.5);

    // Load the P(k) interpolation data we will use for each multipole of each model.
    std::string root(modelrootName);
    if(0 < root.size() && root[root.size()-1] != '/') root += '/';
    boost::format fileName("%s%s_matterpower.dat");
    bool extrapolateBelow(true),extrapolateAbove(true);
    double maxRelError(1e-3);
    cosmo::TabulatedPowerCPtr Pfid,Pnw;
    try {
        Pfid = cosmo::createTabulatedPower(boost::str(fileName % root % fiducialName),
            extrapolateBelow,extrapolateAbove,maxRelError);
        Pnw = cosmo::createTabulatedPower(boost::str(fileName % root % nowigglesName),
            extrapolateBelow,extrapolateAbove,maxRelError);
    }
    catch(cosmo::RuntimeError const &e) {
        throw RuntimeError("BaoKSpaceCorrelationModel: error while reading model interpolation data.");
    }
    // Internally, we use the peak = (fid - nw) and smooth = nw components to evaluate our model.
    cosmo::TabulatedPowerCPtr Ppk = Pfid->createDelta(Pnw);

    // Create smart pointers to our power spectra Ppk(k) and Pnw(fid)
    likely::GenericFunctionPtr PpkPtr =
        likely::createFunctionPtr<const cosmo::TabulatedPower>(Ppk);
    likely::GenericFunctionPtr PnwPtr =
        likely::createFunctionPtr<const cosmo::TabulatedPower>(Pnw);

    // Create a smart pointer to our k-space distortion model D(k,mu_k)
    cosmo::RMuFunctionCPtr distortionModelPtr(new cosmo::RMuFunction(boost::bind(
        &BaoKSpaceCorrelationModel::_evaluateKSpaceDistortion,this,_1,_2)));

    // Create our fiducial and no-wiggles models. We don't initialize our models
    // yet, and instead wait until we are first evaluated and have values for
    // our distortion parameters.
    if(rmin >= rmax) throw RuntimeError("BaoKSpaceCorrelationModel: expected rmin < rmax.");
    if(rmin <= 0) throw RuntimeError("BaoKSpaceCorrelationModel: expected rmin > 0.");
    if(dilmin > dilmax) throw RuntimeError("BaoKSpaceCorrelationModel: expected dilmin <= dilmax.");
    if(dilmin <= 0) throw RuntimeError("BaoKSpaceCorrelationModel: expected dilmin > 0.");
    // Expand the radial ranges needed for transforms to allow for the min/max dilation.
    rmin *= dilmin;
    rmax *= dilmax;
    // Space interpolation points at ~1 Mpc/h.
    int nr = (int)std::ceil(rmax-rmin); 
    double abspow(0);
    bool symmetric(true);
    // Xipk(r,mu) ~ D(k,mu_k)*Ppk(k)
    _Xipk.reset(new cosmo::DistortedPowerCorrelation(PpkPtr,distortionModelPtr,
        rmin,rmax,nr,ellMax,symmetric,relerr,abserr,abspow));
    // Xinw(r,mu) ~ D(k,mu_k)*Pnw(k)
    _Xinw.reset(new cosmo::DistortedPowerCorrelation(PnwPtr,distortionModelPtr,
        rmin,rmax,nr,ellMax,symmetric,relerr,abserr,abspow));

    // Define our r-space broadband distortion models, if any.
    if(distAdd.length() > 0) {
        _distortAdd.reset(new baofit::BroadbandModel("Additive broadband distortion",
            "dist add",distAdd,distR0,zref,this));
    }
    if(distMul.length() > 0) {
        _distortMul.reset(new baofit::BroadbandModel("Multiplicative broadband distortion",
            "dist mul",distMul,distR0,zref,this));
    }
}

local::BaoKSpaceCorrelationModel::~BaoKSpaceCorrelationModel() { }

double local::BaoKSpaceCorrelationModel::_evaluateKSpaceDistortion(double k, double mu_k) const {
    double mu2(mu_k*mu_k);
    // Calculate linear bias model
    double tmp = 1 + _betaz*mu2;
    double linear = tmp*tmp;
    // Calculate non-linear broadening
    double snl2 = _snlPar2*mu2 + _snlPerp2*(1-mu2);
    double nonlinear = std::exp(-0.5*snl2*k*k);
    // Put the pieces together
    return nonlinear*linear;
}

double local::BaoKSpaceCorrelationModel::_evaluate(double r, double mu, double z,
bool anyChanged) const {

    // Lookup linear bias parameters.
    double beta = getParameterValue(0);
    double bb = getParameterValue(1);
    // Calculate bias^2 from beta and bb.
    double bias = bb/(1+beta);
    double biasSq = bias*bias;

    // Lookup linear bias redshift evolution parameters.
    double gammaBias = getParameterValue(2);
    double gammaBeta = getParameterValue(3);
    // Apply redshift evolution
    biasSq = _redshiftEvolution(biasSq,gammaBias,z);
    _betaz = _redshiftEvolution(beta,gammaBeta,z);

    // Lookup non-linear broadening parameters.
    double snlPerp = getParameterValue(_nlBase);
    double snlPar = snlPerp*getParameterValue(_nlBase+1);
    _snlPerp2 = snlPerp*snlPerp;
    _snlPar2 = snlPar*snlPar;

    // Redo the transforms from (k,mu_k) to (r,mu) if necessary
    if(anyChanged) {
        bool nlChanged = isParameterValueChanged(_nlBase) || isParameterValueChanged(_nlBase+1);
        bool otherChanged = isParameterValueChanged(0);
        int nmu(20),minSamplesPerDecade(40);
        double margin(4), vepsMax(1e-1), vepsMin(1e-6);
        bool optimize(false),bypass(false),converged(true);
        if(!_Xipk->isInitialized()) {
            // Initialize the first time. This is when the automatic calculation of numerical
            // precision parameters takes place.
            _Xipk->initialize(nmu,minSamplesPerDecade,margin,vepsMax,vepsMin,optimize);
            if(_verbose) {
                std::cout << "-- Initialized peak k-space model:" << std::endl;
                _Xipk->printToStream(std::cout);
            }
        }
        else if(nlChanged || otherChanged) {
            // We are already initialized, so just redo the transforms.
            converged &= _Xipk->transform(bypass);
        }
        // Are we only applying non-linear broadening to the peak?
        if(!_nlBroadband) {
            _snlPerp2 = _snlPar2 = 0;
            nlChanged = false;
        }
        if(!_Xinw->isInitialized()) {
            // Initialize the first time. This is when the automatic calculation of numerical
            // precision parameters takes place.
            _Xinw->initialize(nmu,minSamplesPerDecade,margin,vepsMax,vepsMin,optimize);
            if(_verbose) {
                std::cout << "-- Initialized no-wiggles k-space model:" << std::endl;
                _Xinw->printToStream(std::cout);
            }
        }
        else if(nlChanged || otherChanged) {
            // We are already initialized, so just redo the transforms.
            converged &= _Xinw->transform(bypass);
        }
        if(!converged) {
            if(++_nWarnings <= _maxWarnings) {
                std::cout << "WARNING: transforms not converged with:" << std::endl;
                printCurrentValues(std::cout);
                if(_nWarnings == _maxWarnings) {
                    std::cout << "(will not print any more warnings like this)" << std::endl;
                }
            }
        }
    }

    // Lookup BAO peak parameter values.
    double ampl = getParameterValue(_baoBase);
    double scale = getParameterValue(_baoBase + 1);
    double scale_parallel = getParameterValue(_baoBase + 2);
    double scale_perp = getParameterValue(_baoBase + 3);
    double gamma_scale = getParameterValue(_baoBase + 4);

    // Transform (r,mu) to (rBAO,muBAO) using the scale parameters.
    double rBAO, muBAO;
    if(_anisotropic) {
        double apar = _redshiftEvolution(scale_parallel,gamma_scale,z);
        double aperp = _redshiftEvolution(scale_perp,gamma_scale,z);
        double musq(mu*mu);
        scale = std::sqrt(apar*apar*musq + aperp*aperp*(1-musq));
        rBAO = r*scale;
        muBAO = apar*mu/scale;
    }
    else {
        scale = _redshiftEvolution(scale,gamma_scale,z);
        rBAO = r*scale;
        muBAO = mu;
    }

    // Check dilation limits
    if(scale < _dilmin) {
        throw RuntimeError("BaoKSpaceCorrelationModel: hit min dilation limit.");
    }
    else if(scale > _dilmax) {
        throw RuntimeError("BaoKSpaceCorrelationModel: hit max dilation limit.");
    }

    // Calculate the cosmological predictions...
    // the peak model is always evaluated at (rBAO,muBAO)
    double peak = _Xipk->getCorrelation(rBAO,muBAO);
    // the decoupled option determines where we evaluate the smooth model
    double smooth = (_decoupled) ? _Xinw->getCorrelation(r,mu) : _Xinw->getCorrelation(rBAO,muBAO);
    // Combine the pieces with the appropriate normalization factors
    double xi = biasSq*(ampl*peak + smooth);
    
    // Add r-space broadband distortions, if any.
    if(_distortMul) xi *= 1 + _distortMul->_evaluate(r,mu,z,anyChanged);
    if(_distortAdd) {
        double distortion = _distortAdd->_evaluate(r,mu,z,anyChanged);
        // The additive distortion is multiplied by ((1+z)/(1+z0))^gammaBias
        xi += _redshiftEvolution(distortion,gammaBias,z);
    }

    return xi;
}

void  local::BaoKSpaceCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
    out << "Using " << (_anisotropic ? "anisotropic":"isotropic") << " BAO scales." << std::endl;
    out << "Scales apply to BAO peak " << (_decoupled ? "only." : "and cosmological broadband.") << std::endl;
}
