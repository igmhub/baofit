// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"

#include "likely/FitParameter.h"

#include <iostream>

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel(std::string const &name)
: _name(name)
{ }

local::AbsCorrelationModel::~AbsCorrelationModel() { }

void local::AbsCorrelationModel::defineParameter(std::string const &name, double value, double error) {
    _parameters.push_back(likely::FitParameter(name,value,error));
}

likely::FitParameters const &local::AbsCorrelationModel::getParameters() const {
    return _parameters;
}

int local::AbsCorrelationModel::getNParameters(bool onlyFloating) const {
    return onlyFloating ? likely::countFloatingFitParameters(_parameters) : _parameters.size();
}

void  local::AbsCorrelationModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    out << "Correlation Model \"" << _name << "\" has initial parameters:" << std::endl;
    printFitParametersToStream(_parameters,out,formatSpec);
}

void local::AbsCorrelationModel::configure(std::string const &script) {
    likely::modifyFitParameters(_parameters,script);
}
