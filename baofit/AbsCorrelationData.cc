// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/AbsCorrelationData.h"
#include "baofit/RuntimeError.h"

#include "likely/CovarianceMatrix.h"
#include "likely/AbsBinning.h"
#include "likely/RuntimeError.h"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"
#include "boost/smart_ptr.hpp"

#include <fstream>
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
double muMin, double muMax, double rperpMin, double rperpMax, double rparMin, double rparMax,
cosmo::Multipole lMin, cosmo::Multipole lMax, double zMin, double zMax) {
    if(rMin > rMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected r-min <= r-max.");
    if(rVetoMin > rVetoMax) {
        throw RuntimeError("AbsCorrelationData::setFinalCuts: expected rveto-min <= rveto-max.");
    }
    if(muMin > muMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected mu-min <= mu-max.");
    if(rperpMin > rperpMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected rperp-min <= rperp-max.");
    if(rparMin > rparMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected rpar-min <= rpar-max.");
    if(lMin > lMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected lmin <= lmax.");
    if(zMin > zMax) throw RuntimeError("AbsCorrelationData::setFinalCuts: expected z-min <= z-max.");
    _rMin = rMin; _rMax = rMax;
    _rVetoMin = rVetoMin; _rVetoMax = rVetoMax;
    _muMin = muMin; _muMax = muMax;
    _rperpMin = rperpMin; _rperpMax = rperpMax;
    _rparMin = rparMin; _rparMax = rparMax;
    _lMin = lMin; _lMax = lMax;
    _zMin = zMin; _zMax = zMax;
    _haveFinalCuts = true;
}

void local::AbsCorrelationData::_cloneFinalCuts(AbsCorrelationData &other) const {
    other._rMin = _rMin; other._rMax = _rMax;
    other._rVetoMin = _rVetoMin; other._rVetoMax = _rVetoMax;
    other._muMin = _muMin; other._muMax = _muMax;
    other._rperpMin = _rperpMin; other._rperpMax = _rperpMax;
    other._rparMin = _rparMin; other._rparMax = _rparMax;
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
            double rpar = r*mu;
            if(rpar < _rparMin || rpar > _rparMax) continue;
            double rperp = r*std::sqrt(1-mu*mu);
            if(rperp < _rperpMin || rperp > _rperpMax) continue;
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

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

baofit::AbsCorrelationDataPtr local::loadCorrelationData(std::string const &dataName,
baofit::AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted) {

    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData(dynamic_cast<AbsCorrelationData*>(prototype->clone(true)));

    // General stuff we will need for reading both files.
    std::string line;
    int lines;
    
    // import boost spirit parser symbols
    using qi::double_;
    using qi::int_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;

    // Loop over lines in the parameter file.
    std::string paramsName = dataName + (weighted ? ".wdata" : ".data");
    std::ifstream paramsIn(paramsName.c_str());
    if(!paramsIn.good()) throw RuntimeError("loadCorrelationData: Unable to open " + paramsName);
    lines = 0;
    int index;
    double data;
    std::vector<double> bin(3);
    while(std::getline(paramsIn,line)) {
        lines++;
        bin.resize(0);
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(index) = _1] >> double_[ref(data) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadCorrelationData: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        binnedData->setData(index,data,weighted);
    }
    paramsIn.close();
    int ndata = binnedData->getNBinsWithData();
    int nbins = binnedData->getGrid().getNBinsTotal();
    if(verbose) {
        std::cout << "Read " << ndata << " of " << nbins << " data values from "
            << paramsName << std::endl;
    }

    // Loop over lines in the (inverse) covariance file.
    std::string covName = dataName + (icov ? ".icov" : ".cov");
    std::ifstream covIn(covName.c_str());
    if(!covIn.good()) throw RuntimeError("loadCorrelationData: Unable to open " + covName);
    lines = 0;
    double value;
    int index1,index2;
    while(std::getline(covIn,line)) {
        lines++;
        bin.resize(0);
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(index1) = _1] >> int_[ref(index2) = _1] >> double_[ref(value) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadCorrelationData: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + covName);
        }
        // Check for invalid offsets.
        if(index1 < 0 || index2 < 0 || index1 >= nbins || index2 >= nbins ||
        !binnedData->hasData(index1) || !binnedData->hasData(index2)) {
            throw RuntimeError("loadCorrelationData: invalid covariance indices on line " +
                boost::lexical_cast<std::string>(lines) + " of " + covName);
        }
        // Add this covariance to our dataset.
        if(icov) {
            binnedData->setInverseCovariance(index1,index2,value);
        }
        else {
            binnedData->setCovariance(index1,index2,value);            
        }
    }
    covIn.close();
    if(verbose) {
        int ncov = (ndata*(ndata+1))/2;
        std::cout << "Read " << lines << " of " << ncov
            << " covariance values from " << covName << std::endl;
    }

    return binnedData;
}
