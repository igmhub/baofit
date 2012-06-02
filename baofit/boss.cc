// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/boss.h"
#include "baofit/types.h"
#include "baofit/RuntimeError.h"
#include "baofit/MultipoleCorrelationData.h"
#include "baofit/QuasarCorrelationData.h"

#include "cosmo/AbsHomogeneousUniverse.h"

#include "likely/UniformBinning.h"
#include "likely/UniformSampling.h"

#include "boost/lexical_cast.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

#include <fstream>
#include <iostream>
#include <cassert>

namespace local = baofit::boss;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

std::vector<double> local::twoStepSampling(
int nBins, double breakpoint,double dlog, double dlin, double eps) {
    assert(breakpoint > 0 && dlog > 0 && dlin > 0 && eps > 0);
    std::vector<double> samplePoints;
    // first sample is at zero.
    samplePoints.push_back(0);
    // next samples are uniformly spaced up to the breakpoint.
    int nUniform = std::floor(breakpoint/dlin);
    for(int k = 1; k <= nUniform; ++k) {
        samplePoints.push_back((k-0.5)*dlin);
    }
    // remaining samples are logarithmically spaced, with log-weighted bin centers.
    double ratio = std::log((breakpoint+dlog)/breakpoint);
    for(int k = 1; k < nBins-nUniform; ++k) {
        samplePoints.push_back(breakpoint*std::exp(ratio*(k-0.5)));
    }
    return samplePoints;
}

// Loads a binned correlation function in French format and returns a shared pointer to
// a MultipoleCorrelationData.
baofit::AbsCorrelationDataPtr
local::loadFrench(std::string dataName, double zref, bool verbose) {

    // Create the new BinnedData that we will fill.
    likely::AbsBinningCPtr
        rBins(new likely::UniformBinning(0,200,50)),
        ellBins(new likely::UniformSampling(0,0,1)), // only monopole for now
        zBins(new likely::UniformSampling(zref,zref,1));
    baofit::AbsCorrelationDataPtr
        binnedData(new baofit::MultipoleCorrelationData(rBins,ellBins,zBins));

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
    std::string paramsName(dataName + ".txt");
    std::ifstream paramsIn(paramsName.c_str());
    if(!paramsIn.good()) throw RuntimeError("loadFrench: Unable to open " + paramsName);
    lines = 0;
    double rval,mono,quad;
    std::vector<double> bin(3);
    while(std::getline(paramsIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                double_[ref(rval) = _1] >> double_[ref(mono) = _1] >> double_[ref(quad) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadFrench: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        bin[0] = rval;
        bin[2] = zref;
        bin[1] = 0;
        int monoIndex = binnedData->getIndex(bin);
        binnedData->setData(monoIndex,mono);
    }
    paramsIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " data values from " << paramsName << std::endl;
    }
    
    // Loop over lines in the covariance file.
    std::string covName = paramsName;
    int pos = covName.rfind('/',-1);
    covName.insert(pos+1,"cov_");
    std::ifstream covIn(covName.c_str());
    if(!covIn.good()) throw RuntimeError("Unable to open " + covName);
    lines = 0;
    int index1,index2;
    double cov;
    while(std::getline(covIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(index1) = _1] >> int_[ref(index2) = _1] >> double_[ref(cov) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadFrench: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + covName);
        }
        if(index1 <= index2 && index2 < 50) binnedData->setCovariance(index1,index2,cov);
    }
    covIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " covariance values from " << covName << std::endl;
    }
    return binnedData;
}

// Loads a binned correlation function in cosmolib format and returns a BinnedData object.
// The fast option disables regexp checks for valid numeric inputs.
baofit::AbsCorrelationDataPtr local::loadCosmolib(std::string dataName,
likely::AbsBinningCPtr llBins, likely::AbsBinningCPtr sepBins, likely::AbsBinningCPtr zBins,
double rmin, double rmax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology,
bool verbose, bool icov, bool fast) {

    // Create the new BinnedData that we will fill.
    baofit::AbsCorrelationDataPtr
        binnedData(new baofit::QuasarCorrelationData(llBins,sepBins,zBins,rmin,rmax,llmin,cosmology));

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
    std::string paramsName(dataName + ".params");
    std::ifstream paramsIn(paramsName.c_str());
    if(!paramsIn.good()) throw RuntimeError("loadCosmolib: Unable to open " + paramsName);
    lines = 0;
    double xi;
    std::vector<double> bin(3);
    while(std::getline(paramsIn,line)) {
        lines++;
        bin.resize(0);
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                double_[ref(xi) = _1] >> double_ >> "| Lya covariance 3D (" >>
                double_[push_back(ref(bin),_1)] >> ',' >> double_[push_back(ref(bin),_1)] >>
                ',' >> double_[push_back(ref(bin),_1)] >> ')'
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadCosmolib: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        int index = binnedData->getIndex(bin);
        binnedData->setData(index,xi);        
    }
    paramsIn.close();
    if(verbose) {
        std::cout << "Read " << binnedData->getNBinsWithData() << " of "
            << binnedData->getNBinsTotal() << " data values from " << paramsName << std::endl;
    }
    
    // Loop over lines in the covariance file.
    std::string covName(dataName + (icov ? ".icov" : ".cov"));
    std::ifstream covIn(covName.c_str());
    if(!covIn.good()) throw RuntimeError("Unable to open " + covName);
    lines = 0;
    double value;
    int offset1,offset2;
    while(std::getline(covIn,line)) {
        lines++;
        bin.resize(0);
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(offset1) = _1] >> int_[ref(offset2) = _1] >> double_[ref(value) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadCosmolib: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        // Add this covariance to our dataset.
        if(icov) value = -value; // !?! see line #388 of Observed2Point.cpp
        int index1 = *(binnedData->begin()+offset1), index2 = *(binnedData->begin()+offset2);
        if(icov) {
            binnedData->setInverseCovariance(index1,index2,value);
        }
        else {
            binnedData->setCovariance(index1,index2,value);
        }
    }
    covIn.close();
    if(verbose) {
        int ndata = binnedData->getNBinsWithData();
        int ncov = (ndata*(ndata+1))/2;
        std::cout << "Read " << lines << " of " << ncov
            << " covariance values from " << covName << std::endl;
    }

    // Check for zero values on the diagonal
    for(likely::BinnedData::IndexIterator iter = binnedData->begin();
    iter != binnedData->end(); ++iter) {
        int index = *iter;
        if(icov) {
            if(0 == binnedData->getInverseCovariance(index,index)) {
                binnedData->setInverseCovariance(index,index,1e-30);
            }
        }
        else {
            if(0 == binnedData->getCovariance(index,index)) {
                binnedData->setCovariance(index,index,1e40);
            }                
        }
    }
    // Compress our binned data to take advantage of a potentially sparse covariance matrix.
    binnedData->compress();
    return binnedData;
}
