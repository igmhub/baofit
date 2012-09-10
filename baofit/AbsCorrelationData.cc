// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"
#include "baofit/RuntimeError.h"

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
