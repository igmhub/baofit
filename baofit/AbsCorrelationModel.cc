// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "cosmo/TransferFunctionPowerSpectrum.h" // for getMultipole(...)

#include "boost/bind.hpp"

#include <cmath>

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel(std::string const &name)
: FitModel(name), _indexBase(-1), _crossCorrelation(false), _combinedBias(false), _dvIndex(-1), _betabiasIndex(-1),
_nbins(0)
{ }

local::AbsCorrelationModel::~AbsCorrelationModel() { }

double local::AbsCorrelationModel::evaluate(double r, double mu, double z,
likely::Parameters const &params, int index) {
    bool anyChanged = updateParameterValues(params);
    _updateInternalParameters();
    if(_dvIndex >= 0) _applyVelocityShift(r,mu,z);
    if(_dvIndex >= 0 && _nbins > 0 && anyChanged) {
        double rbin, mubin, zbin;
        for(int i = 0; i < _nbins; ++i) {
            rbin = _rbin[i];
            mubin = _mubin[i];
            zbin = _zbin[i];
            _applyVelocityShift(rbin,mubin,zbin);
            _rbinShift[i] = rbin;
            _mubinShift[i] = mubin;
            _zbinShift[i] = zbin;
        }
    }
    double result = _evaluate(r,mu,z,anyChanged,index);
    resetParameterValuesChanged();
    return result;
}

double local::AbsCorrelationModel::evaluate(double r, cosmo::Multipole multipole, double z,
likely::Parameters const &params, int index) {
    bool anyChanged = updateParameterValues(params);
    double result = _evaluate(r,multipole,z,anyChanged,index);
    resetParameterValuesChanged();
    return result;
}

void local::AbsCorrelationModel::setCoordinates(std::vector<double> rbin, std::vector<double> mubin,
std::vector<double> zbin) {
    _rbin = rbin;
    _mubin = mubin;
    _zbin = zbin;
    _nbins = _rbin.size();
    if(_nbins != _mubin.size() || _nbins != _zbin.size()) {
        throw RuntimeError("AbsCorrelationModel::setCoordinates: coordinate vectors not the same size.");
    }
    _rbinShift = rbin;
    _mubinShift = mubin;
    _zbinShift = zbin;
}

double local::AbsCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged, int index) const {
    // Get a pointer to our (r,mu,z) evaluator. We need a typedef here to disambiguate the two
    // overloaded _evaluate methods.
    typedef double (AbsCorrelationModel::*fOfRMuZ)(double, double, double, bool, int) const;
    fOfRMuZ fptr(&AbsCorrelationModel::_evaluate);
    // Call our (r,mu,z) evaluator once with mu=0 and the input value of anyChanged so it can
    // do any necessary one-time calculations. Subsequent calls will use anyChanged = false.
    (this->*fptr)(r,0,z,anyChanged,index);
    // Create a smart pointer to a function object of mu with the other args (r,z,anyChanged,index) bound.
    likely::GenericFunctionPtr fOfMuPtr(
        new likely::GenericFunction(boost::bind(fptr,this,r,_1,z,false,index)));
    // Finally we have something we can pass to the generic multipole projection integrator.
    return cosmo::getMultipole(fOfMuPtr,(int)multipole);
}

void local::AbsCorrelationModel::_setZRef(double zref) {
    if(zref < 0) throw RuntimeError("AbsCorrelationModel: expected zref >= 0.");
    _zref = zref;
}

int local::AbsCorrelationModel::_defineLinearBiasParameters(double zref, bool crossCorrelation, bool combinedBias) {
    if(_indexBase >= 0) throw RuntimeError("AbsCorrelationModel: linear bias parameters already defined.");
    _setZRef(zref);
    // Linear bias parameters
    _indexBase = defineParameter("beta",1.4,0.1);
    _setBetaIndex(_indexBase + BETA);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    _setBbIndex(_indexBase + BB);
    // Redshift evolution parameters
    defineParameter("gamma-bias",3.8,0.3);
    _setGammaBiasIndex(_indexBase + GAMMA_BIAS);
    int last = defineParameter("gamma-beta",0,0.1);
    _setGammaBetaIndex(_indexBase + GAMMA_BETA);
    if(crossCorrelation) {
        _crossCorrelation = true;
        // Amount to shift each separation's line of sight velocity in km/s
        defineParameter("delta-v",0,10);
        _setDVIndex(_indexBase + DELTA_V);
        // We don't use beta2 and (1+beta2)*bias2 here since for galaxies or quasars
        // the combination beta2*bias2 = f = dln(G)/dln(a) is well constrained.
        defineParameter("bias2",3.6,0.1);
        _setBias2Index(_indexBase + BIAS2);
        last = defineParameter("beta2*bias2",1,0.05);
        _setBeta2Bias2Index(_indexBase + BB2);
    }
    else {
        // not really necessary since the ctor already does this and you cannot call this method
        // more than once
        _crossCorrelation = false;
    }
    if(combinedBias) {
        _combinedBias = true;
        _combBiasBase = defineParameter("beta*bias",-0.196,0.02);
        last = _combBiasBase;
        _setBetaBiasIndex(_combBiasBase);
        // Fix the parameter that is not being used
        configureFitParameters("fix[(1+beta)*bias]=0");
    }
    return last;
}

void local::AbsCorrelationModel::_updateInternalParameters() {
    _beta = getParameterValue(_betaIndex);
    _bias = getParameterValue(_bbIndex)/(1+_beta);
    _gammaBias = getParameterValue(_gammabiasIndex);
    _gammaBeta = getParameterValue(_gammabetaIndex);
    if(_dvIndex >= 0) {
        _bias2 = getParameterValue(_bias2Index);
        _beta2 = getParameterValue(_beta2bias2Index)/_bias2;
    }
    if(_betabiasIndex >= 0) {
        _bias = getParameterValue(_betabiasIndex)/_beta;
    }
}

void local::AbsCorrelationModel::_applyVelocityShift(double &r, double &mu, double z) {
    // Lookup value of delta_v
    double dv = getParameterValue(_dvIndex);
    // Convert dv in km/s to dpi in Mpc/h using a flat matter+lambda cosmology with OmegaLambda = 0.73
    double zp1 = 1+z;
    double dpi = (dv/100.)*(1+z)/std::sqrt(0.73+0.27*zp1*zp1*zp1);
    // Calculate the effect of changing pi by dpi in the separation
    double rnew = std::sqrt(r*r + 2*r*mu*dpi + dpi*dpi);
    double munew = (r*mu+dpi)/rnew;
    r = rnew;
    mu = munew;
}

double local::redshiftEvolution(double p0, double gamma, double z, double zref) {
    return p0*std::pow((1+z)/(1+zref),gamma);
}

double local::AbsCorrelationModel::_getNormFactor(cosmo::Multipole multipole, double z) const {
    if(_indexBase < 0) throw RuntimeError("AbsCorrelationModel: no linear bias parameters defined.");
    // Lookup the linear bias parameters at the reference redshift.
    double beta = getParameterValue(_indexBase + BETA);
    double bb = getParameterValue(_indexBase + BB);
    // Calculate bias from beta and bb.
    double bias = bb/(1+beta);
    if(_combinedBias) {
        double betabias = getParameterValue(_combBiasBase);
        bias = betabias/beta;
    }
    // For cross correlations, the linear and quadratic beta terms are independent and
    // the overall bias could be negative.
    double betaAvg,betaProd,biasSq;
    if(_crossCorrelation) {
        double bias2 = getParameterValue(_indexBase + BIAS2);
        double bb2 = getParameterValue(_indexBase + BB2);
        double beta2 = bb2/bias2;
        betaAvg = (beta + beta2)/2;
        betaProd = beta*beta2;
        biasSq = bias*bias2;
    }
    else {
        betaAvg = beta;
        betaProd = beta*beta;
        biasSq = bias*bias;
    }
    // Calculate redshift evolution of biasSq, betaAvg and betaProd.
    double gammaBias = getParameterValue(_indexBase + GAMMA_BIAS);
    double gammaBeta = getParameterValue(_indexBase + GAMMA_BETA);
    biasSq = redshiftEvolution(biasSq,gammaBias,z,_zref);
    betaAvg = redshiftEvolution(betaAvg,gammaBeta,z,_zref);
    betaProd = redshiftEvolution(betaProd,2*gammaBeta,z,_zref);
    // Return the requested normalization factor.
    switch(multipole) {
    case cosmo::Hexadecapole:
        return biasSq*betaProd*(8./35.);
    case cosmo::Quadrupole:
        return biasSq*((4./3.)*betaAvg + (4./7.)*betaProd);
    default:
        return biasSq*(1 + (2./3.)*betaAvg + (1./5.)*betaProd);
    }
}

double local::AbsCorrelationModel::_getRBin(int index) const {
    if(index < 0 || index >= _nbins) {
        throw RuntimeError("AbsCorrelationModel::getRBin: invalid index.");
    }
    return _rbinShift[index];
}

double local::AbsCorrelationModel::_getMuBin(int index) const {
    if(index < 0 || index >= _nbins) {
        throw RuntimeError("AbsCorrelationModel::getMuBin: invalid index.");
    }
    return _mubinShift[index];
}

double local::AbsCorrelationModel::_getZBin(int index) const {
    if(index < 0 || index >= _nbins) {
        throw RuntimeError("AbsCorrelationModel::getZBin: invalid index.");
    }
    return _zbinShift[index];
}

void  local::AbsCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    FitModel::printToStream(out,formatSpec);
    out << std::endl << "Reference redshift = " << _zref << std::endl;
}
