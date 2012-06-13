// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/MultipoleCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/AbsBinning.h"

#include <cmath>
#include <iostream>

namespace local = baofit;

local::MultipoleCorrelationData::MultipoleCorrelationData(likely::AbsBinningCPtr axis1,
likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3, double rmin, double rmax)
: AbsCorrelationData(axis1,axis2,axis3,Multipole)
{
    _initialize(rmin,rmax);
}

local::MultipoleCorrelationData::MultipoleCorrelationData(std::vector<likely::AbsBinningCPtr> axes,
double rmin, double rmax)
: AbsCorrelationData(axes,Multipole)
{
    if(axes.size() != 3) {
        throw RuntimeError("MultipoleCorrelationData: expected 3 axes.");
    }
    _initialize(rmin,rmax);
}

void local::MultipoleCorrelationData::_initialize(double rmin, double rmax) {
    if(rmin >= rmax) {
        throw RuntimeError("MultipoleCorrelationData: expected rmin < rmax.");
    }
    _rmin = rmin;
    _rmax = rmax;
    _lastIndex = -1;
}

local::MultipoleCorrelationData::~MultipoleCorrelationData() { }

local::MultipoleCorrelationData *local::MultipoleCorrelationData::clone(bool binningOnly) const {
    return binningOnly ?
        new MultipoleCorrelationData(getAxisBinning(),_rmin,_rmax) :
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

void local::MultipoleCorrelationData::finalize() {
    std::set<int> keep;
    std::vector<double> binCenter;
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(*iter);
        double r(getRadius(index));
        // Keep this bin in our pruned dataset?
        if(r >= _rmin && r < _rmax) {
            keep.insert(index);
        }
    }
    prune(keep);
    AbsCorrelationData::finalize();
}

void local::MultipoleCorrelationData::dump(std::ostream &out, int zIndex) const {
    std::vector<likely::AbsBinningCPtr> binning = getAxisBinning();
    int nRadialBins(binning[0]->getNBins()), nEllBins(binning[1]->getNBins());
    std::vector<int> bin(3);
    bin[2] = 0;
    for(int rIndex = 0; rIndex < nRadialBins; ++rIndex) {
        double rval(binning[0]->getBinCenter(rIndex));
        if(rval < _rmin) continue;
        if(rval >= _rmax) break;
        out << rval;
        bin[0] = rIndex;
        for(int ellIndex = 0; ellIndex < nEllBins; ++ellIndex) {
            bin[1] = ellIndex;
            int index = getIndex(bin);
            if(hasData(index)) {
                out << ' ' << getData(index);
                if(hasCovariance()) {
                    out << ' ' << std::sqrt(getCovariance(index,index));
                }
                else {
                    out << " 0";
                }
            }
            else {
                out << " 0 -1";
            }
        }
        out << std::endl;
    }
}
