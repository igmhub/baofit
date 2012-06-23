// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/MultipoleCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/AbsBinning.h"

#include <cmath>
#include <iostream>

namespace local = baofit;

local::MultipoleCorrelationData::MultipoleCorrelationData(likely::AbsBinningCPtr axis1,
likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3, double rmin, double rmax,
cosmo::Multipole ellmin, cosmo::Multipole ellmax)
: AbsCorrelationData(axis1,axis2,axis3,Multipole)
{
    _initialize(rmin,rmax,ellmin,ellmax);
}

local::MultipoleCorrelationData::MultipoleCorrelationData(std::vector<likely::AbsBinningCPtr> axes,
double rmin, double rmax, cosmo::Multipole ellmin, cosmo::Multipole ellmax)
: AbsCorrelationData(axes,Multipole)
{
    if(axes.size() != 3) {
        throw RuntimeError("MultipoleCorrelationData: expected 3 axes.");
    }
    _initialize(rmin,rmax,ellmin,ellmax);
}

void local::MultipoleCorrelationData::_initialize(double rmin, double rmax,
cosmo::Multipole ellmin, cosmo::Multipole ellmax) {
    if(rmin >= rmax) {
        throw RuntimeError("MultipoleCorrelationData: expected rmin < rmax.");
    }
    if(ellmin > ellmax) {
        throw RuntimeError("MultipoleCorrelationData: expected ellmin <= ellmax.");
    }
    _rmin = rmin;
    _rmax = rmax;
    _ellmin = ellmin;
    _ellmax = ellmax;
    _lastIndex = -1;
}

local::MultipoleCorrelationData::~MultipoleCorrelationData() { }

local::MultipoleCorrelationData *local::MultipoleCorrelationData::clone(bool binningOnly) const {
    return binningOnly ?
        new MultipoleCorrelationData(getAxisBinning(),_rmin,_rmax,_ellmin,_ellmax) :
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
    // We don't check for a valid enum type here on purpose, so that additional modes can be
    // included in the dataset for correct weighted, and then pruned out in the finalize step.
    _ellLast = static_cast<cosmo::Multipole>(std::floor(_binCenter[1]+0.5));
    _zLast = _binCenter[2];
    _lastIndex = index;
}

void local::MultipoleCorrelationData::finalize() {
    std::set<int> keep;
    std::vector<double> binCenter;
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(*iter), ell(getMultipole(index));
        double r(getRadius(index));
        // Keep this bin in our pruned dataset?
        if(r >= _rmin && r < _rmax && ell >= _ellmin && ell <= _ellmax) {
            keep.insert(index);
        }
    }
    prune(keep);
    AbsCorrelationData::finalize();
}

void local::MultipoleCorrelationData::dump(std::ostream &out, std::vector<double> const &weights) const {
    if(weights.size() > 0 && weights.size() != getNBinsWithData()) {
        throw RuntimeError("MultipoleCorrelationData::dump: unexpected size of weights.");
    }
    std::vector<likely::AbsBinningCPtr> binning = getAxisBinning();
    int nRadialBins(binning[0]->getNBins()), nEllBins(binning[1]->getNBins()), nZBins(binning[2]->getNBins());
    std::vector<int> bin(3);
    for(int rIndex = 0; rIndex < nRadialBins; ++rIndex) {
        double rval(binning[0]->getBinCenter(rIndex));
        if(rval < _rmin) continue;
        if(rval >= _rmax) break;
        out << rval;
        bin[0] = rIndex;
        for(int zIndex = 0; zIndex < nZBins; ++zIndex) {
            bin[2] = zIndex;
            for(int ellIndex = 0; ellIndex < nEllBins; ++ellIndex) {
                bin[1] = ellIndex;
                int index = getIndex(bin);
                double data(0),error(-1);
                if(hasData(index)) {
                    data = getData(index);
                    if(0 < weights.size()) {
                        error = weights[getOffsetForIndex(index)];
                        error = (error > 0) ? 1/std::sqrt(error) : -1;
                    }
                    else if(hasCovariance()) {
                        error = std::sqrt(getCovariance(index,index));
                    }
                }
                out << ' ' << data << ' ' << error;
            }
        }
        out << std::endl;
    }
}
