// Created 21-Oct-2014 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "baofit/MultiCorrelationData.h"
#include "baofit/RuntimeError.h"

namespace local = baofit;

local::MultiCorrelationData::MultiCorrelationData() { }

local::MultiCorrelationData::~MultiCorrelationData() { }

void local::MultiCorrelationData::addDataSet(AbsCorrelationDataCPtr data) {
	_datasets.push_back(data);
}

double local::MultiCorrelationData::getRadius(int index) const {
}

double local::MultiCorrelationData::getCosAngle(int index) const {
}

cosmo::Multipole local::MultiCorrelationData::getMultipole(int index) const {
}

double local::MultiCorrelationData::getRedshift(int index) const {
}
