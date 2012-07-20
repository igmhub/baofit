// Created 11-Jul-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/ComovingCorrelationData.h"
#include "baofit/RuntimeError.h"

namespace local = baofit;

local::ComovingCorrelationData::ComovingCorrelationData(likely::AbsBinningCPtr rBins,
likely::AbsBinningCPtr muBins, likely::AbsBinningCPtr zBins,
double rmin, double rmax, double muMin, double muMax, double rVetoMin, double rVetoMax) 
: AbsCorrelationData(rBins,muBins,zBins,Coordinate)
{
    _initialize(rmin,rmax,muMin,muMax,rVetoMin,rVetoMax);
}

local::ComovingCorrelationData::ComovingCorrelationData(std::vector<likely::AbsBinningCPtr> axes,
double rmin, double rmax, double muMin, double muMax, double rVetoMin, double rVetoMax) 
: AbsCorrelationData(axes,Coordinate)
{
    _initialize(rmin,rmax,muMin,muMax,rVetoMin,rVetoMax);
}

local::ComovingCorrelationData::~ComovingCorrelationData() { }

local::ComovingCorrelationData *local::ComovingCorrelationData::clone(bool binningOnly) const {
    return binningOnly ?
        new ComovingCorrelationData(getAxisBinning(),_rmin,_rmax,_muMin,_muMax,_rVetoMin,_rVetoMax) :
        new ComovingCorrelationData(*this);
}

void local::ComovingCorrelationData::_initialize(double rmin, double rmax, double muMin, double muMax,
double rVetoMin, double rVetoMax) {
    if(rmin >= rmax) {
        throw RuntimeError("ComovingCorrelationData: expected rmin < rmax.");
    }
    if(muMin >= muMax) {
        throw RuntimeError("MultipoleCorrelationData: expected mu-min < mu-max.");
    }
    _rmin = rmin;
    _rmax = rmax;
    _muMin = muMin;
    _muMax = muMax;
    if(rVetoMin > rVetoMax) {
        throw RuntimeError("ComovingCorrelationData: expected rVetoMin <= rVetoMax.");
    }
    _rVetoMin = rVetoMin;
    _rVetoMax = rVetoMax;
    _lastIndex = -1;
}

void local::ComovingCorrelationData::finalize() {
    std::set<int> keep;
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of r,mu,z at the center of this bin.
        int index(*iter);
        double r(getRadius(index)), mu(getCosAngle(index)), z(getRedshift(index));
        // Keep this bin in our pruned dataset?
        if(r >= _rmin && r < _rmax) {
            if(r <= _rVetoMin || r >= _rVetoMax) keep.insert(index);
        }
    }
    prune(keep);
    AbsCorrelationData::finalize();
}

void local::ComovingCorrelationData::_setIndex(int index) const {
    if(index == _lastIndex) return;
    getBinCenters(index,_binCenter);
    _lastIndex = index;
}

double local::ComovingCorrelationData::getRadius(int index) const {
    _setIndex(index);
    return _binCenter[0];
}

double local::ComovingCorrelationData::getCosAngle(int index) const {
    _setIndex(index);
    return _binCenter[1];
}

double local::ComovingCorrelationData::getRedshift(int index) const {
    _setIndex(index);
    return _binCenter[2];
}
