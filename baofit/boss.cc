// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/boss.h"
#include "baofit/types.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/MultipoleCorrelationData.h"
#include "baofit/QuasarCorrelationData.h"
#include "baofit/ComovingCorrelationData.h"

#include "cosmo/AbsHomogeneousUniverse.h"

#include "likely/AbsBinning.h"
#include "likely/UniformBinning.h"
#include "likely/UniformSampling.h"
#include "likely/NonUniformSampling.h"
#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"

#include "boost/regex.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"
#include "boost/smart_ptr.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

namespace local = baofit::boss;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

// Reproduce the hybrid linear-log binning of cosmolib's ForestCovariance3DTheory_Xi::BinToWavelength_3D
std::vector<double> local::twoStepSampling(double breakpoint,double llmax,double dlog,double dlin) {
    if(!(breakpoint > 0 && dlog > 0 && dlin > 0 && llmax > breakpoint)) {
        throw RuntimeError("twoStepSampling: invalid parameters.");
    }
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
    int k = 1;
    while(1) {
        double llval = breakpoint*std::exp(ratio*(k-0.5));
        // Stop when we pass llmax
        if(llval > llmax) break;
        samplePoints.push_back(llval);
        k++;
    }
    return samplePoints;
}

baofit::AbsCorrelationDataPtr local::createComovingPrototype(ComovingCorrelationData::CoordinateSystem coords,
bool verbose, std::string const &axis1Bins, std::string const &axis2Bins, std::string const &axis3Bins) {
    // Parse the binning from the strings provided.
    likely::AbsBinningCPtr axis1ptr,axis2ptr,axis3ptr;
    try {
        axis1ptr = likely::createBinning(axis1Bins);
        if(verbose) {
            int nbins = axis1ptr->getNBins();
            if(coords == ComovingCorrelationData::CartesianCoordinates) {
                std::cout << "r_par bin centers:";
            }
            else {
                std::cout << "r bins centers:";
            }
            for(int bin = 0; bin < nbins; ++bin) {
                std::cout << (bin ? ',':' ') << axis1ptr->getBinCenter(bin);
            }
            std::cout << " (n = " << axis1ptr->getNBins() << ")" << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("createComovingPrototype: error in axis 1 binning.");
    }
    try {
        axis2ptr = likely::createBinning(axis2Bins);
        if(verbose) {
            int nbins = axis2ptr->getNBins();
            if(coords == ComovingCorrelationData::PolarCoordinates) {
                std::cout << "mu bin centers:";
            }
            else if(coords == ComovingCorrelationData::CartesianCoordinates) {
                std::cout << "r_perp bin centers:";
            }
            else {
                std::cout << "multipoles:";
            }
            for(int bin = 0; bin < nbins; ++bin) {
                std::cout << (bin ? ',':' ') << axis2ptr->getBinCenter(bin);
            }
            std::cout << " (n = " << axis2ptr->getNBins() << ")" << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("createComovingPrototype: error in axis 2 binning.");
    }
    try {
        axis3ptr = likely::createBinning(axis3Bins);
        if(verbose) {
            int nbins = axis3ptr->getNBins();
            std::cout << "z bin centers:";
            for(int bin = 0; bin < nbins; ++bin) {
                std::cout << (bin ? ',':' ') << axis3ptr->getBinCenter(bin);
            }
            std::cout << " (n = " << axis3ptr->getNBins() << ")" << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError("createComovingPrototype: error in axis 3 binning.");
    }

    likely::BinnedGrid grid(axis1ptr,axis2ptr,axis3ptr);
    baofit::AbsCorrelationDataPtr prototype(new baofit::ComovingCorrelationData(grid,coords));
    return prototype;    
}

// Creates a prototype QuasarCorrelationData with the specified binning and cosmology.
baofit::AbsCorrelationDataPtr local::createCosmolibPrototype(
double minsep, double dsep, int nsep, double minz, double dz, int nz,
double minll, double maxll, double dll, double dll2,
double llMin, double llMax, double sepMin, double sepMax,
bool fixCov, cosmo::AbsHomogeneousUniversePtr cosmology) {

    // Initialize the (logLambda,separation,redshift) binning from command-line params.
    likely::AbsBinningCPtr llBins,
        sepBins(new likely::UniformBinning(minsep,minsep+nsep*dsep,nsep)),
        zBins(new likely::UniformSampling(minz+0.5*dz,minz+(nz-0.5)*dz,nz));
    if(0 == dll2) {
        // If dll does not divide evenly into [minll,maxll], stick with minll,dll and adjust maxll.
        int nll = std::floor((maxll-minll)/dll+0.5);
        maxll = minll + dll*nll;
        llBins.reset(new likely::UniformBinning(minll,maxll,nll));
    }
    else {
        llBins.reset(new likely::NonUniformSampling(twoStepSampling(minll,maxll,dll,dll2)));
    }

    // Create the new BinnedData that we will fill.
    likely::BinnedGrid grid(llBins,sepBins,zBins);
    baofit::AbsCorrelationDataPtr
        prototype(new baofit::QuasarCorrelationData(grid,llMin,llMax,sepMin,sepMax,fixCov,cosmology));
    return prototype;
}

baofit::AbsCorrelationDataPtr local::loadSaved(std::string const &dataName,
baofit::AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted) {
    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData((baofit::QuasarCorrelationData *)(prototype->clone(true)));

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
    if(!paramsIn.good()) throw RuntimeError("loadSaved: Unable to open " + paramsName);
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
            throw RuntimeError("loadSaved: error reading line " +
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
    if(!covIn.good()) throw RuntimeError("loadSaved: Unable to open " + covName);
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
            throw RuntimeError("loadSaved: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + covName);
        }
        // Check for invalid offsets.
        if(index1 < 0 || index2 < 0 || index1 >= nbins || index2 >= nbins ||
        !binnedData->hasData(index1) || !binnedData->hasData(index2)) {
            throw RuntimeError("loadSaved: invalid covariance indices on line " +
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

// Loads a binned correlation function in cosmolib format and returns a BinnedData object.
baofit::AbsCorrelationDataPtr local::loadCosmolib(std::string const &dataName,
baofit::AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted,
int &reuseCovIndex, int reuseCov) {

    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData((baofit::QuasarCorrelationData *)(prototype->clone(true)));
    likely::BinnedGrid grid(binnedData->getGrid());

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
    double data,cinvData;
    std::vector<double> bin(3);
    while(std::getline(paramsIn,line)) {
        lines++;
        bin.resize(0);
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                double_[ref(data) = _1] >> double_[ref(cinvData) = _1] >> "| Lya covariance 3D (" >>
                double_[push_back(ref(bin),_1)] >> ',' >> double_[push_back(ref(bin),_1)] >>
                ',' >> double_[push_back(ref(bin),_1)] >> ')'
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadCosmolib: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        int index = grid.getIndex(bin);
        binnedData->setData(index, weighted ? cinvData : data, weighted);        
    }
    paramsIn.close();
    int ndata = binnedData->getNBinsWithData();
    if(verbose) {
        std::cout << "Read " << ndata << " of " << binnedData->getGrid().getNBinsTotal()
            << " data values from " << paramsName << std::endl;
    }

    // Do we need to reuse the covariance estimated for the first realization of this plate?
    std::string covName;
    if(reuseCov >= 0) {
        // Parse the data name.
        boost::regex namePattern("([a-zA-Z0-9/_\\.]+/)?([0-9]+_)([0-9]+)\\.cat\\.([0-9]+)");
        boost::match_results<std::string::const_iterator> what;
        if(!boost::regex_match(dataName,what,namePattern)) {
            throw RuntimeError("loadCosmolib: cannot parse name \"" + dataName + "\"");
	    }
        covName = what[1]+what[2]+boost::lexical_cast<std::string>(reuseCov)+".cat."+what[4];
    }
    else {
        covName = dataName;
    }
    covName += (icov ? ".icov" : ".cov");
    
    // Can we reuse a previously loaded covariance matrix?
    // Initialize a dictionary of dataset indices and covariance filenames.
    typedef std::map<std::string,int> CovarianceCache;
    static CovarianceCache covarianceCache;
    typedef CovarianceCache::value_type CovarianceCacheEntry;
    static int nextIndex(0);
    if(reuseCov) {
        CovarianceCache::const_iterator found = covarianceCache.find(covName);
        if(found == covarianceCache.end()) {
            covarianceCache.insert(CovarianceCacheEntry(covName,nextIndex));
        }
        else {
            reuseCovIndex = found->second;
            if(verbose) {
                std::cout << "Reusing covariance matrix for observation ["
                    << reuseCovIndex << "] from " << covName << std::endl;
            }
        }
    }
    nextIndex++;

    if(reuseCovIndex < 0) {
        // Loop over lines in the covariance file.
        std::ifstream covIn(covName.c_str());
        if(!covIn.good()) throw RuntimeError("loadCosmolib: Unable to open " + covName);
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
            // Check for invalid offsets.
            if(offset1 < 0 || offset2 < 0 || offset1 >= ndata || offset2 >= ndata) {
                throw RuntimeError("loadCosmolib: invalid covariance indices on line " +
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
    }
    
    return binnedData;
}
