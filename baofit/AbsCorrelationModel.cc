// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include <cmath>

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel(std::string const &name)
: FitModel(name), _indexBase(-1)
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

int local::AbsCorrelationModel::_defineLinearBiasParameters(double zref) {
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
    // Return the requested normalization factor.
    switch(multipole) {
    case cosmo::Hexadecapole:
        return biasSq*beta*beta*(8./35.);
    case cosmo::Quadrupole:
        return biasSq*beta*(4./3. + (4./7.)*beta);
    default:
        return biasSq*(1 + beta*(2./3. + (1./5.)*beta));
    }
}

void  local::AbsCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    FitModel::printToStream(out,formatSpec);
    if(_indexBase >= 0) out << std::endl << "Reference redshift = " << _zref << std::endl;        
}
