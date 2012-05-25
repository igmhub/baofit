// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"

namespace local = baofit;

local::AbsCorrelationData::AbsCorrelationData(
likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3)
: likely::BinnedData(axis1,axis2,axis3)
{
}

local::AbsCorrelationData::~AbsCorrelationData() { }
