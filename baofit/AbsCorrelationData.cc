// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/CovarianceMatrix.h"
#include "likely/AbsBinning.h"
#include "likely/RuntimeError.h"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"

#include <iostream>

namespace local = baofit;

local::AbsCorrelationData::AbsCorrelationData(likely::BinnedGrid grid, TransverseBinningType type)
: likely::BinnedData(grid), _type(type), _haveFinalCuts(false)
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


likely::BinnedGrid local::createCorrelationGrid(std::string const &axis1Bins, std::string const &axis2Bins,
    std::string const &axis3Bins, std::string const &axisLabels, bool verbose) {
    // Extract the comma-separated 3 axis labels
    std::vector<std::string> labels;
    boost::split(labels,axisLabels,boost::is_any_of(","));
    if(labels.size() != 3) {
        throw RuntimeError("createCorrelationGrid: expected 3 axis labels.");
    }
    // Parse the binning from the strings provided for each axis
    likely::AbsBinningCPtr axis1ptr,axis2ptr,axis3ptr;
    try {
        axis1ptr = likely::createBinning(axis1Bins);
        if(verbose) {
            std::cout << labels[0] << " bin centers:";
            int nbins = axis1ptr->getNBins();
            for(int bin = 0; bin < nbins; ++bin) {
                std::cout << (bin ? ',':' ') << axis1ptr->getBinCenter(bin);
            }
            std::cout << " (n = " << nbins << ")" << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("createCorrelationGrid: error in axis 1 (" + labels[0] + ") binning.");
    }
    try {
        axis2ptr = likely::createBinning(axis2Bins);
        if(verbose) {
            std::cout << labels[1] << " bin centers:";
            int nbins = axis2ptr->getNBins();
            for(int bin = 0; bin < nbins; ++bin) {
                std::cout << (bin ? ',':' ') << axis2ptr->getBinCenter(bin);
            }
            std::cout << " (n = " << nbins << ")" << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("createCorrelationGrid: error in axis 2 (" + labels[1] + ") binning.");
    }
    try {
        axis3ptr = likely::createBinning(axis3Bins);
        if(verbose) {
            std::cout << labels[2] << " bin centers:";
            int nbins = axis3ptr->getNBins();
            for(int bin = 0; bin < nbins; ++bin) {
                std::cout << (bin ? ',':' ') << axis3ptr->getBinCenter(bin);
            }
            std::cout << " (n = " << nbins << ")" << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("createCorrelationGrid: error in axis 3 (" + labels[2] + ") binning.");
    }
    // Return a BinnedGrid object for these 3 axes.
    return likely::BinnedGrid(axis1ptr,axis2ptr,axis3ptr);    
}
