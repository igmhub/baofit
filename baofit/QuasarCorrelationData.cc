// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/QuasarCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "cosmo/AbsHomogeneousUniverse.h"

#include <cmath>

namespace local = baofit;

local::QuasarCorrelationData::QuasarCorrelationData(
likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3,
double rmin, double rmax, double muMin, double muMax, double llmin, bool fixCov, double rVetoMin, double rVetoMax,
cosmo::AbsHomogeneousUniversePtr cosmology)
: AbsCorrelationData(axis1,axis2,axis3,Coordinate)
{
  _initialize(rmin,rmax,muMin,muMax,llmin,fixCov,rVetoMin,rVetoMax,cosmology);
}

local::QuasarCorrelationData::QuasarCorrelationData(
std::vector<likely::AbsBinningCPtr> axes, double rmin, double rmax, double muMin, double muMax, double llmin,
bool fixCov, double rVetoMin, double rVetoMax, cosmo::AbsHomogeneousUniversePtr cosmology)
: AbsCorrelationData(axes,Coordinate)
{
    if(axes.size() != 3) {
        throw RuntimeError("QuasarCorrelationData: expected 3 axes.");
    }
    _initialize(rmin,rmax,muMin,muMax,llmin, fixCov, rVetoMin,rVetoMax,cosmology);
}

void local::QuasarCorrelationData::_initialize(double rmin, double rmax, double muMin, double muMax, double llmin,
bool fixCov, double rVetoMin, double rVetoMax, cosmo::AbsHomogeneousUniversePtr cosmology) {
    if(rmin >= rmax) {
        throw RuntimeError("QuasarCorrelationData: expected rmin < rmax.");
    }
    if(muMin >= muMax) {
        throw RuntimeError("MultipoleCorrelationData: expected mu-min < mu-max.");
    }
    _rmin = rmin;
    _rmax = rmax;
    _muMin = muMin;
    _muMax = muMax;
    _llmin = llmin;
    _fixCov = fixCov; 
    if(rVetoMin > rVetoMax) {
        throw RuntimeError("QuasarCorrelationData: expected rVetoMin <= rVetoMax.");
    }
    _rVetoMin = rVetoMin;
    _rVetoMax = rVetoMax;
    _cosmology = cosmology;
    _lastIndex = -1;
    _arcminToRad = 4*std::atan(1)/(60.*180.);    
}

local::QuasarCorrelationData::~QuasarCorrelationData() { }

local::QuasarCorrelationData *local::QuasarCorrelationData::clone(bool binningOnly) const {
    return binningOnly ?
        new QuasarCorrelationData(getAxisBinning(),_rmin,_rmax,_muMin,_muMax,_llmin,
				  _fixCov,_rVetoMin,_rVetoMax,_cosmology) :
        new QuasarCorrelationData(*this);
}

void local::QuasarCorrelationData::fixCovariance(double ll0, double c0, double c1, double c2) {

    if (!isCovarianceModifiable()) {
        throw RuntimeError("QuasarCorrelationData::fixCovariance: not modifiable.");
    }
    // Make sure that our our data vector is un-weighted.
    getData(*begin());

    // Save values in the outer loop, for re-use in the inner loop.
    std::vector<double> dll;
    dll.reserve(getNBinsWithData());
    std::vector<int> bin(3);

    // Loop over all bins.
    for(IndexIterator iter1 = begin(); iter1 != end(); ++iter1) {
        int i1(*iter1);
        // Calculate and save the value of ll - ll0 at the center of this bin.
        getBinCenters(i1,_binCenter);
        dll.push_back(_binCenter[0] - ll0);
        // Remember the indices of this 3D bin along our sep,z axes
        getBinIndices(i1,bin);
        int sepIndex(bin[1]), zIndex(bin[2]);
        // Loop over bins with index i2 <= i1
        for(IndexIterator iter2 = begin(); iter2 != end(); ++iter2) {
            int i2(*iter2);
            // we want to to iter1 end() but not beyond!!
            if (i2>i1) continue;
            // Check that this bin has the same sep,z indices
            getBinIndices(i2,bin);
            if(bin[1] != sepIndex || bin[2] != zIndex) continue;
            // Calculate (ll1 - ll0)*(ll2 - ll0) using cached values.
            double d = dll[i1]*dll[i2];
            // Update the covariance for (i1,i2)
            // magic constants are set by the requirement that for
            // a certain cov, you should add something that is "large"
            // but at the same time does not make numerical errors unbearable
            double C(getCovariance(i1,i2));
            C += c0 + c1*d + c2*d*d;
            setCovariance(i1,i2,C);
        }
    }
}

void local::QuasarCorrelationData::finalize() {

    // First fix Covariance
    if (_fixCov) fixCovariance();

    // Next do pruning
    std::set<int> keep;
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(*iter);
        double r(getRadius(index)), mu(getCosAngle(index)), z(getRedshift(index));
        double ll(_binCenter[0]);
        // Keep this bin in our pruned dataset?
        if(r >= _rmin && r < _rmax && mu >= _muMin && mu <= _muMax && ll >= _llmin) {
            if(r <= _rVetoMin || r >= _rVetoMax) {
                keep.insert(index);
                // Remember these values.
                _rLookup.push_back(r);
                _muLookup.push_back(mu);
                _zLookup.push_back(z);
            }
        }
    }
    prune(keep);
    AbsCorrelationData::finalize();
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
    _llLast = _binCenter[0];
    _sepLast = _binCenter[1];
    _zLast = _binCenter[2];
    
    transform(_llLast, _sepLast, _binWidth[1],_zLast,_rLast,_muLast);
    _lastIndex = index;
}

double local::QuasarCorrelationData::getRadius(int index) const {
    if(isFinalized()) return _rLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _rLast;
}

double local::QuasarCorrelationData::getCosAngle(int index) const {
    if(isFinalized()) return _muLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _muLast;
}

double local::QuasarCorrelationData::getLogLambda(int index) const {
  // not yet caching loglambda
  //    if(isFinalized()) return _llLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _llLast;
}

 double local::QuasarCorrelationData::getSeparation(int index) const {
  // not yet caching separation
   // if(isFinalized()) return _setLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _sepLast;
}

double local::QuasarCorrelationData::getRedshift(int index) const {
    if(isFinalized()) return _zLookup[getOffsetForIndex(index)];
    _setIndex(index);
    return _zLast;
}
