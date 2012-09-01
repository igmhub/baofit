// Created 31-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/PkCorrelationModel.h"

namespace local = baofit;

local::PkCorrelationModel::PkCorrelationModel(double kmin, double kmax, int nk)
: AbsCorrelationModel("P(ell,k) Correlation Model")
{ }

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
