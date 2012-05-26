// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/QuasarCorrelationData.h"

#include "likely/BinnedData.h"
#include "cosmo/AbsHomogeneousUniverse.h"

#include <cmath>

namespace local = baofit;

local::QuasarCorrelationData::QuasarCorrelationData(
likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3,
cosmo::AbsHomogeneousUniversePtr cosmology)
: AbsCorrelationData(axis1,axis2,axis3), _cosmology(cosmology),
_finalized(false), _lastIndex(-1)
{
    _arcminToRad = 4*std::atan(1)/(60.*180.);
}

local::QuasarCorrelationData::~QuasarCorrelationData() { }

void local::QuasarCorrelationData::finalize(double rmin, double rmax, double llmin) {
    std::set<int> keep;
    std::vector<double> binCenter,binWidth;
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(*iter);
        double r(getRadius(index)), mu(getCosAngle(index)), z(getRedshift(index));
        double ll(_binCenter[0]);
        // Keep this bin in our pruned dataset?
        if(r >= rmin && r < rmax && ll >= llmin) {
            keep.insert(index);            
            // Remember these values.
            _rLookup.push_back(r);
            _muLookup.push_back(mu);
            _zLookup.push_back(z);
        }
    }
    prune(keep);
    _finalized = true;
}

void local::QuasarCorrelationData::transform(double ll, double sep, double dsep, double z,
double &r, double &mu) const {
    double ratio(std::exp(0.5*ll)),zp1(z+1);
    double z1(zp1/ratio-1), z2(zp1*ratio-1);
    double drLos = _cosmology->getLineOfSightComovingDistance(z2) -
        _cosmology->getLineOfSightComovingDistance(z1);
    // Calculate the geometrically weighted mean separation of this bin as
    // Integral[s^2,{s,smin,smax}]/Integral[s,{s,smin,smax}] = s + dsep^2/(12*s)
    double swgt = sep + (dsep*dsep/12)/sep;
    double drPerp = _cosmology->getTransverseComovingScale(z)*(swgt*_arcminToRad);
    double rsq = drLos*drLos + drPerp*drPerp;
    r = std::sqrt(rsq);
    mu = std::abs(drLos)/r;
}

void local::QuasarCorrelationData::_setIndex(int index) const {
    if(index == _lastIndex) return;
    getBinCenters(index,_binCenter);
    getBinWidths(index,_binWidth);
    _zLast = _binCenter[2];
    transform(_binCenter[0],_binCenter[1],_binWidth[1],_zLast,_rLast,_muLast);
    _lastIndex = index;
}

double local::QuasarCorrelationData::getRadius(int index) const {
    if(_finalized) return _rLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _rLast;
}

double local::QuasarCorrelationData::getCosAngle(int index) const {
    if(_finalized) return _muLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _muLast;
}

double local::QuasarCorrelationData::getRedshift(int index) const {
    if(_finalized) return _zLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _zLast;
}
