// Created 31-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/PkCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include "likely/Interpolator.h"
#include "likely/function.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

namespace local = baofit;

local::PkCorrelationModel::PkCorrelationModel(std::string const &modelrootName, std::string const &nowigglesName,
double kmin, double kmax, int nk)
: AbsCorrelationModel("P(ell,k) Correlation Model")
{
    // Define the overall normalization parameters for each multipole.
    // Calculate initial values assuming the following linear bias parameter values:
    double bias = -0.14, beta = 1.4, gamma_bias = 3.8, gamma_beta = 0;
    defineParameter("Pk norm ell 0",bias*bias*(1 + (2./3.)*beta + (1/5.)*beta*beta),0.002);
    defineParameter("Pk norm ell 2",bias*bias*((4./3.)*beta + (4./7.)*beta*beta),0.002);
    defineParameter("Pk norm ell 4",bias*bias*((8./35.)*beta*beta),0.001);
    defineParameter("gamma Pk norm ell 0", gamma_bias,0.2);
    defineParameter("gamma Pk norm ell 2", gamma_bias + gamma_beta,0.2);
    defineParameter("gamma Pk norm ell 4", gamma_bias + 2*gamma_beta,0.2);
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
        throw RuntimeError("BaoCorrelationModel: error while reading model interpolation data.");
    }
}

local::PkCorrelationModel::~PkCorrelationModel() { }

double local::PkCorrelationModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    return 0;
}

double local::PkCorrelationModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    return 0;
}

void  local::PkCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    AbsCorrelationModel::printToStream(out,formatSpec);
}
