// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"

namespace local = baofit;

local::AbsCorrelationData::AbsCorrelationData(
likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3,
TransverseBinningType type)
: likely::BinnedData(axis1,axis2,axis3), _type(type)
{
}

local::AbsCorrelationData::AbsCorrelationData(std::vector<likely::AbsBinningCPtr> axes,
TransverseBinningType type)
: likely::BinnedData(axes), _type(type)
{
}

local::AbsCorrelationData::~AbsCorrelationData() { }

double local::AbsCorrelationData::getCosAngle(int index) const { return 0; }

cosmo::Multipole local::AbsCorrelationData::getMultipole(int index) const { return cosmo::Monopole; }