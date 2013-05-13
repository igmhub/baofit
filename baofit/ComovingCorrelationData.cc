// Created 11-Jul-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/ComovingCorrelationData.h"
#include "baofit/RuntimeError.h"

#include <cmath>

namespace local = baofit;

local::ComovingCorrelationData::ComovingCorrelationData(likely::BinnedGrid grid,
CoordinateSystem coordinateSystem) 
: AbsCorrelationData(grid,Coordinate), _coordinateSystem(coordinateSystem), _lastIndex(-1)
{
}

local::ComovingCorrelationData::~ComovingCorrelationData() { }

local::ComovingCorrelationData *local::ComovingCorrelationData::clone(bool binningOnly) const {
    ComovingCorrelationData *data = binningOnly ?
        new ComovingCorrelationData(getGrid(),_coordinateSystem) :
        new ComovingCorrelationData(*this);
    _cloneFinalCuts(*data);
    return data;
}

void local::ComovingCorrelationData::finalize() {
    std::set<int> keep;
    _applyFinalCuts(keep);
    prune(keep);
    AbsCorrelationData::finalize();
}

void local::ComovingCorrelationData::_setIndex(int index) const {
    if(index == _lastIndex) return;
    getGrid().getBinCenters(index,_binCenter);
    _lastIndex = index;
}

double local::ComovingCorrelationData::getRadius(int index) const {
    _setIndex(index);
    if(_coordinateSystem == PolarCoordinates) {
        return _binCenter[0];
    }
    else {
        double rpar = _binCenter[0], rperp = _binCenter[1];
        return std::sqrt(rpar*rpar+rperp*rperp);
    }
}

double local::ComovingCorrelationData::getCosAngle(int index) const {
    _setIndex(index);
    if(_coordinateSystem == PolarCoordinates) {
        return _binCenter[1];
    }
    else {
        double rpar = _binCenter[0], rperp = _binCenter[1];
        return rpar/std::sqrt(rpar*rpar+rperp*rperp);
    }
}

double local::ComovingCorrelationData::getRedshift(int index) const {
    _setIndex(index);
    return _binCenter[2];
}
