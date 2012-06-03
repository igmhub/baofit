// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/MultipoleCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/AbsBinning.h"

#include <cmath>
#include <iostream>

namespace local = baofit;

local::MultipoleCorrelationData::MultipoleCorrelationData(likely::AbsBinningCPtr axis1,
likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3)
: AbsCorrelationData(axis1,axis2,axis3,Multipole)
{ }

local::MultipoleCorrelationData::MultipoleCorrelationData(std::vector<likely::AbsBinningCPtr> axes)
: AbsCorrelationData(axes,Multipole)
{
    if(axes.size() != 3) {
        throw RuntimeError("MultipoleCorrelationData: expected 3 axes.");
    }
}

local::MultipoleCorrelationData::~MultipoleCorrelationData() { }

local::MultipoleCorrelationData *local::MultipoleCorrelationData::clone(bool binningOnly) const {
    return binningOnly ?
        new MultipoleCorrelationData(getAxisBinning()) :
        new MultipoleCorrelationData(*this);
}

double local::MultipoleCorrelationData::getRadius(int index) const {
    _setIndex(index);
    return _rLast;
}

cosmo::Multipole local::MultipoleCorrelationData::getMultipole(int index) const {
    _setIndex(index);
    return _ellLast;
}

double local::MultipoleCorrelationData::getRedshift(int index) const {
    _setIndex(index);
    return _zLast;
}

void local::MultipoleCorrelationData::_setIndex(int index) const {
    if(index == _lastIndex) return;
    getBinCenters(index,_binCenter);
    _rLast = _binCenter[0];
    switch((int)std::floor(_binCenter[1]+0.5)) {
    case 0:
        _ellLast = cosmo::Monopole;
        break;
    case 2:
        _ellLast = cosmo::Quadrupole;
        break;
    case 4:
        _ellLast = cosmo::Hexadecapole;
        break;
    default:
        throw RuntimeError("MultipoleCorrelationData: invalid multipole value.");
    }
    _zLast = _binCenter[2];
    _lastIndex = index;
}

void local::MultipoleCorrelationData::dump(
std::ostream &out, cosmo::Multipole multipole, int zIndex) const {
    std::vector<likely::AbsBinningCPtr> binning = getAxisBinning();
    int nRadialBins(binning[0]->getNBins());
    int ellIndex(binning[1]->getBinIndex((double)multipole));
    std::vector<int> bin(3);
    for(int rIndex = 0; rIndex < nRadialBins; ++rIndex) {
        double rval(binning[0]->getBinCenter(rIndex));
        bin[0] = rIndex;
        int index = getIndex(bin);
        double cov = hasCovariance() ? getCovariance(index,index) : 0;
        out << rval << ' ' << getData(index) << ' ' << std::sqrt(cov) << std::endl;
    }
}
