// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel(std::string const &name)
: FitModel(name)
{ }

local::AbsCorrelationModel::~AbsCorrelationModel() { }

double local::AbsCorrelationModel::evaluate(double r, double mu, double z,
likely::Parameters const &params) {
    bool anyChanged = updateParameterValues(params);
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
