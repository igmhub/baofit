// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/MultipoleCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/AbsBinning.h"

#include <cmath>
#include <iostream>

namespace local = baofit;

local::MultipoleCorrelationData::MultipoleCorrelationData(likely::AbsBinningCPtr axis1,
likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3)
: AbsCorrelationData(axis1,axis2,axis3,Multipole), _lastIndex(-1)
{
}

local::MultipoleCorrelationData::MultipoleCorrelationData(std::vector<likely::AbsBinningCPtr> axes)
: AbsCorrelationData(axes,Multipole), _lastIndex(-1)
{
    if(axes.size() != 3) {
        throw RuntimeError("MultipoleCorrelationData: expected 3 axes.");
    }
}

local::MultipoleCorrelationData::~MultipoleCorrelationData() { }

local::MultipoleCorrelationData *local::MultipoleCorrelationData::clone(bool binningOnly) const {
    MultipoleCorrelationData *data = binningOnly ?
        new MultipoleCorrelationData(getAxisBinning()) : new MultipoleCorrelationData(*this);
    _cloneFinalCuts(*data);
    return data;
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
    // included in the dataset for correct weighting, and then pruned out in the finalize step.
    _ellLast = static_cast<cosmo::Multipole>(std::floor(_binCenter[1]+0.5));
    _zLast = _binCenter[2];
    _lastIndex = index;
}

void local::MultipoleCorrelationData::finalize() {
    std::set<int> keep;
    _applyFinalCuts(keep);
    prune(keep);
    AbsCorrelationData::finalize();
}

void local::MultipoleCorrelationData::dump(std::ostream &out, double rmin, double rmax,
std::vector<double> const &weights) const {
    if(weights.size() > 0 && weights.size() != getNBinsWithData()) {
        throw RuntimeError("MultipoleCorrelationData::dump: unexpected size of weights.");
    }
    std::vector<likely::AbsBinningCPtr> binning = getAxisBinning();
    int nRadialBins(binning[0]->getNBins()), nEllBins(binning[1]->getNBins()), nZBins(binning[2]->getNBins());
    std::vector<int> bin(3);
    for(int rIndex = 0; rIndex < nRadialBins; ++rIndex) {
        double rval(binning[0]->getBinCenter(rIndex));
        if(rval < rmin) continue;
        if(rval > rmax) break;
        out << rval;
        bin[0] = rIndex;
        for(int zIndex = 0; zIndex < nZBins; ++zIndex) {
            bin[2] = zIndex;
            for(int ellIndex = 0; ellIndex < nEllBins; ++ellIndex) {
                bin[1] = ellIndex;
                int index = getIndex(bin);
                double data(0),error(0);
                if(hasData(index)) {
                    data = getData(index);
                    if(0 < weights.size()) {
                        double weight = weights[getOffsetForIndex(index)];
                        if(weight != 0) error = 1/std::sqrt(std::fabs(weight));
                        if(weight < 0) error = -error;
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
