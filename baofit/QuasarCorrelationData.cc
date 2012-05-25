// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/QuasarCorrelationData.h"

#include "likely/BinnedData.h"
#include "cosmo/AbsHomogeneousUniverse.h"

#include <cmath>

namespace local = baofit;

local::QuasarCorrelationData::QuasarCorrelationData(
likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3,
cosmo::AbsHomogeneousUniversePtr cosmology)
: AbsCorrelationData(axis1,axis2,axis3), _cosmology(cosmology)
{
}

void local::QuasarCorrelationData::finalize(double rmin, double rmax, double llmin) {
    std::set<int> keep;
    double arcminToRad = 4*std::atan(1)/(60.*180.);
    std::vector<double> binCenter,binWidth;
    // Loop over bins with data.
    for(IndexIterator iter = begin(); iter != end(); ++iter) {
        // Lookup the value of ll,sep,z at the center of this bin.
        int index(*iter);
        getBinCenters(index,binCenter);
        double ll(binCenter[0]), sep(binCenter[1]), z(binCenter[2]);
        getBinWidths(index,binWidth);
        double ds(binWidth[1]);
        // Calculate the corresponding values of r,mu,z.
        double ratio(std::exp(0.5*ll)),zp1(z+1);
        double z1(zp1/ratio-1), z2(zp1*ratio-1);
        double drLos = _cosmology->getLineOfSightComovingDistance(z2) -
            _cosmology->getLineOfSightComovingDistance(z1);
        // Calculate the geometrically weighted mean separation of this bin as
        // Integral[s^2,{s,smin,smax}]/Integral[s,{s,smin,smax}] = s + ds^2/(12*s)
        double swgt = sep + (ds*ds/12)/sep;
        double drPerp = _cosmology->getTransverseComovingScale(z)*(swgt*arcminToRad);
        double rsq = drLos*drLos + drPerp*drPerp;
        double r3d = std::sqrt(rsq);
        double mu = std::abs(drLos)/r3d;
        // Keep this bin in our pruned dataset?
        if(r3d >= rmin && r3d < rmax && ll >= llmin) {
            keep.insert(index);            
            // Remember these values.
            _rLookup.push_back(r3d);
            _muLookup.push_back(mu);
            _zLookup.push_back(z);
        }
    }
    prune(keep);
}

local::QuasarCorrelationData::~QuasarCorrelationData() { }

double local::QuasarCorrelationData::getRadius(int index) const {
    return _rLookup[getOffsetForIndex(index)];
}

double local::QuasarCorrelationData::getCosAngle(int index) const {
    return _muLookup[getOffsetForIndex(index)];
}

double local::QuasarCorrelationData::getRedshift(int index) const {
    return _zLookup[getOffsetForIndex(index)];
}
