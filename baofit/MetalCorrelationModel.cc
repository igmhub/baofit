// Created 24-Apr-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "baofit/MetalCorrelationModel.h"
#include "baofit/AbsCorrelationModel.h"
#include "baofit/RuntimeError.h"

#include <cmath>

namespace local = baofit;

local::MetalCorrelationModel::MetalCorrelationModel(int indexBase)
: _indexBase(indexBase) { }

local::MetalCorrelationModel::~MetalCorrelationModel() { }

double local::MetalCorrelationModel::_evaluateMetalCorrelation(double r, double mu) const { }