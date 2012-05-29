// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel() { }

local::AbsCorrelationModel::~AbsCorrelationModel() { }

void local::AbsCorrelationModel::defineParameter(std::string const &name,
double value, double error, bool fixed) {
    _parameters.push_back(likely::FitParameter(name,value,fixed ? 0 : error));
}

likely::FitParameters const &local::AbsCorrelationModel::getParameters() const {
    return _parameters;
}
