// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/BinnedData.h"

namespace local = baofit;

local::AbsCorrelationData::AbsCorrelationData(likely::BinnedDataCPtr data)
: _data(data)
{
}

local::AbsCorrelationData::~AbsCorrelationData() { }

void local::AbsCorrelationData::setData(likely::BinnedDataCPtr data) {
    if(!_data->isCongruent(*data)) {
        throw RuntimeError("AbsCorrelationData::setData: new and old data are not congruent.");
    }
    _data = data;
}

int local::AbsCorrelationData::getSize() const {
    return _data->getNBinsWithData();
}

double local::AbsCorrelationData::getData(int offset) const {
    return _data->getData(_data->getIndexAtOffset(offset));
}
