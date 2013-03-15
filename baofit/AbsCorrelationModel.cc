// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include <cmath>

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel(std::string const &name)
: FitModel(name), _indexBase(-1), _crossCorrelation(false)
{ }

local::AbsCorrelationModel::~AbsCorrelationModel() { }

double local::AbsCorrelationModel::evaluate(double r, double mu, double z,
likely::Parameters const &params) {
    bool anyChanged = updateParameterValues(params);
    _applyVelocityShift(r,mu,z);
    double result = _evaluate(r,mu,z,anyChanged);
    resetParameterValuesChanged();
    return result;
}

double local::AbsCorrelationModel::evaluate(double r, cosmo::Multipole multipole, double z,
likely::Parameters const &params) {
    bool anyChanged = updateParameterValues(params);
    double result = _evaluate(r,multipole,z,anyChanged);
    resetParameterValuesChanged();
    return result;
}

int local::AbsCorrelationModel::_defineLinearBiasParameters(double zref, bool crossCorrelation) {
    if(_indexBase >= 0) throw RuntimeError("AbsCorrelationModel: linear bias parameters already defined.");
    if(zref < 0) throw RuntimeError("AbsCorrelationModel: expected zref >= 0.");
    _zref = zref;
    // Linear bias parameters
    _indexBase = defineParameter("beta",1.4,0.1);
    defineParameter("(1+beta)*bias",-0.336,0.03);
    // Redshift evolution parameters
    defineParameter("gamma-bias",3.8,0.3);
    defineParameter("gamma-beta",0,0.1);    
    // Amount to shift each separation's line of sight velocity in km/s
    int last = defineParameter("delta-v",0,10);
    if(crossCorrelation) {
        _crossCorrelation = true;
        // We use don't use beta2 and (1+beta2)*bias2 here since for galaxies or quasars
        // so that this parameter corresponds directly to f = dln(G)/dln(a), which is what
        // we usually want to constrain when the second component is galaxies or quasars.
        defineParameter("bias2",1.4,0.1);
        last = defineParameter("beta2*bias2",-0.336,0.03);
    }
    else {
        // not really necessary since the ctor already does this and you cannot call this method
        // more than once
        _crossCorrelation = false;
    }
    return last;
}

void local::AbsCorrelationModel::_applyVelocityShift(double &r, double &mu, double z) {
    // Lookup value of delta_v
    double dv = getParameterValue(_indexBase + DELTA_V);
    // Convert dv in km/s to dpi in Mpc/h using a flat matter+lambda cosmology with OmegaLambda = 0.73
    double zp1 = 1+z;
    double dpi = (dv/100.)*(1+z)/std::sqrt(0.73+0.27*zp1*zp1*zp1);
    // Calculate the effect of changing pi by dpi in the separation
    double rnew = std::sqrt(r*r + 2*r*mu*dpi + dpi*dpi);
    double munew = (r*mu+dpi)/rnew;
    r = rnew;
    mu = munew;    
}

double local::AbsCorrelationModel::_redshiftEvolution(double p0, double gamma, double z) const {
    return p0*std::pow((1+z)/(1+_zref),gamma);
}

double local::AbsCorrelationModel::_getNormFactor(cosmo::Multipole multipole, double z) const {
    if(_indexBase < 0) throw RuntimeError("AbsCorrelationModel: no linear bias parameters defined.");
    // Lookup the linear bias parameters.
    double beta0 = getParameterValue(_indexBase + BETA);
    double bb0 = getParameterValue(_indexBase + BB);
    // Calculate bias from beta and bb.
    double bias0 = bb0/(1+beta0);
    double bias0Sq = bias0*bias0;
    // Calculate redshift evolution of bias and beta.
    double biasSq = _redshiftEvolution(bias0Sq,getParameterValue(_indexBase + GAMMA_BIAS),z);
    double beta = _redshiftEvolution(beta0,getParameterValue(_indexBase + GAMMA_BETA),z);
    // Build the combinations needed below.
    double betaAvg,betaProd;
    if(_crossCorrelation) {
        double bias2 = getParameterValue(_indexBase + BIAS2);
        double bb2 = getParameterValue(_indexBase + BB2);
        double beta2 = bb2/bias2;
        betaAvg = (beta + beta2)/2;
        betaProd = beta*beta2;
        biasSq = std::sqrt(biasSq)*bias2;
        // TODO: what about redshift evolution of extra RSD parameters?
    }
    else {
        betaAvg = beta;
        betaProd = beta*beta;
    }
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

void  local::AbsCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    FitModel::printToStream(out,formatSpec);
    if(_indexBase >= 0) out << std::endl << "Reference redshift = " << _zref << std::endl;        
}
