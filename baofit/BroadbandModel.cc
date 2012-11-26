// Created 26-Nov-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/BroadbandModel.h"

namespace local = baofit;

local::BroadbandModel::BroadbandModel(std::string const &name)
: AbsCorrelationModel(name)
{ }

local::BroadbandModel::~BroadbandModel() { }

double local::BroadbandModel::_evaluate(double r, double mu, double z, bool anyChanged) const {
    return 0;
}

double local::BroadbandModel::_evaluate(double r, cosmo::Multipole multipole, double z,
bool anyChanged) const {
    return 0;
}
