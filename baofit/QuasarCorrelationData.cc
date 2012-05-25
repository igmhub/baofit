// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/QuasarCorrelationData.h"

#include "likely/BinnedData.h"
#include "cosmo/AbsHomogeneousUniverse.h"

#include <cmath>

namespace local = baofit;

local::QuasarCorrelationData::QuasarCorrelationData(
likely::BinnedDataCPtr data, cosmo::AbsHomogeneousUniversePtr cosmology)
: AbsCorrelationData(data)
{
    _rLookup.reserve(getSize());
    _muLookup.reserve(getSize());
    _zLookup.reserve(getSize());
    // Loop over bins.
    double arcminToRad = 4*std::atan(1)/(60.*180.);
    std::vector<double> binCenter,binWidth;
    for(int offset = 0; offset < getSize(); ++offset) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(_getData().getIndexAtOffset(offset));
        _getData().getBinCenters(index,binCenter);
        double ll(binCenter[0]), sep(binCenter[1]), z(binCenter[2]);
        _getData().getBinWidths(index,binWidth);
        double ds(binWidth[1]);
        // Calculate the corresponding values of r,mu,z.
        double ratio(std::exp(0.5*ll)),zp1(z+1);
        double z1(zp1/ratio-1), z2(zp1*ratio-1);
        double drLos = cosmology->getLineOfSightComovingDistance(z2) -
            cosmology->getLineOfSightComovingDistance(z1);
        // Calculate the geometrically weighted mean separation of this bin as
        // Integral[s^2,{s,smin,smax}]/Integral[s,{s,smin,smax}] = s + ds^2/(12*s)
        double swgt = sep + (ds*ds/12)/sep;
        double drPerp = cosmology->getTransverseComovingScale(z)*(swgt*arcminToRad);
        double rsq = drLos*drLos + drPerp*drPerp;
        double r3d = std::sqrt(rsq);
        // Remember these values.
        _rLookup.push_back(r3d);
        _muLookup.push_back(std::abs(drLos)/r3d);
        _zLookup.push_back(z);
    }
}

local::QuasarCorrelationData::~QuasarCorrelationData() { }

double local::QuasarCorrelationData::getRadius(int offset) const {
    _checkOffset(offset);
    return _rLookup[offset];
}

double local::QuasarCorrelationData::getCosAngle(int offset) const {
    _checkOffset(offset);
    return _muLookup[offset];
}

double local::QuasarCorrelationData::getRedshift(int offset) const {
    _checkOffset(offset);
    return _zLookup[offset];
}
