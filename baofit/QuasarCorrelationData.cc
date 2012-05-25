// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/QuasarCorrelationData.h"

namespace local = baofit;

local::QuasarCorrelationData::QuasarCorrelationData(likely::BinnedDataCPtr data)
: AbsCorrelationData(data)
{
}

local::QuasarCorrelationData::~QuasarCorrelationData() { }

double local::QuasarCorrelationData::getRadius(int offset) const {
    return 0;
}

double local::QuasarCorrelationData::getCosAngle(int offset) const {
    return 0;
}

double local::QuasarCorrelationData::getRedshift(int offset) const {
    return 0;
}
