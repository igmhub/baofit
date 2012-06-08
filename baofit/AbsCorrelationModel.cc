// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationModel.h"

namespace local = baofit;

local::AbsCorrelationModel::AbsCorrelationModel(std::string const &name)
: FitModel(name)
{ }

local::AbsCorrelationModel::~AbsCorrelationModel() { }
