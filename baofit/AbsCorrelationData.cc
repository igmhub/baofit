// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/CovarianceMatrix.h"

#include "boost/lexical_cast.hpp"

#include <iostream>

namespace local = baofit;

local::AbsCorrelationData::AbsCorrelationData(
likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3,
TransverseBinningType type)
: likely::BinnedData(axis1,axis2,axis3), _type(type), _haveFinalCuts(false)
{
}

local::AbsCorrelationData::AbsCorrelationData(std::vector<likely::AbsBinningCPtr> axes,
TransverseBinningType type)
: likely::BinnedData(axes), _type(type), _haveFinalCuts(false)
{
}

local::AbsCorrelationData::~AbsCorrelationData() { }

double local::AbsCorrelationData::getCosAngle(int index) const { return 0; }

cosmo::Multipole local::AbsCorrelationData::getMultipole(int index) const { return cosmo::Monopole; }

void local::AbsCorrelationData::setFinalCuts(double rMin, double rMax, double rVetoMin, double rVetoMax,
double muMin, double muMax, cosmo::Multipole lMin, cosmo::Multipole lMax,
double zMin, double zMax) {
    if(rMin > rMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected r-min <= r-max.");
    if(rVetoMin > rVetoMax) {
        throw RuntimeError("AbsCorrelationData::setFinalCuts: expected rveto-min <= rveto-max.");
    }
    if(muMin > muMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected mu-min <= mu-max.");
    if(lMin > lMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected lmin <= lmax.");
    if(zMin > zMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected z-min <= z-max.");
    _rMin = rMin; _rMax = rMax;
    _rVetoMin = rVetoMin; _rVetoMax = rVetoMax;
    _muMin = muMin; _muMax = muMax;
    _lMin = lMin; _lMax = lMax;
    _zMin = zMin; _zMax = zMax;
    _haveFinalCuts = true;
}

void local::AbsCorrelationData::_cloneFinalCuts(AbsCorrelationData &other) const {
    other._rMin = _rMin; other._rMax = _rMax;
    other._rVetoMin = _rVetoMin; other._rVetoMax = _rVetoMax;
    other._muMin = _muMin; other._muMax = _muMax;
    other._lMin = _lMin; other._lMax = _lMax;
    other._zMin = _zMin; other._zMax = _zMax;
    other._haveFinalCuts = _haveFinalCuts;
}

void local::AbsCorrelationData::_applyFinalCuts(std::set<int> &keep) const {
    if(!_haveFinalCuts) throw RuntimeError("AbsCorrelationData: no final cuts specified yet.");
    if(!keep.empty()) throw RuntimeError("AbsCorrelationData: expected empty set.");
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(*iter);
        double r(getRadius(index)), z(getRedshift(index));
        // Keep this bin in our pruned dataset?
        if(r < _rMin || r > _rMax) continue;
        if(r > _rVetoMin && r < _rVetoMax) continue;
        if(z < _zMin || z > _zMax) continue;
        if(_type == Coordinate) {
            double mu(getCosAngle(index));
            if(mu < _muMin || mu > _muMax) continue;
        }
        else { // _type == Multipole
            int ell(getMultipole(index));
            if(ell < _lMin || ell > _lMax) continue;
        }
        // This bin passes all cuts so we keep it.
        keep.insert(index);
    }
}
