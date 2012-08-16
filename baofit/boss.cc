// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/boss.h"
#include "baofit/types.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/MultipoleCorrelationData.h"
#include "baofit/QuasarCorrelationData.h"
#include "baofit/ComovingCorrelationData.h"

#include "cosmo/AbsHomogeneousUniverse.h"

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
#include <sstream>

namespace local = baofit::boss;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

// Creates a prototype MultipoleCorrelationData with the specified binning.
baofit::AbsCorrelationDataCPtr local::createDR9LRGPrototype(double zref, double rmin, double rmax,
double rVetoMin, double rVetoMax, std::string const &covName, bool verbose) {
    // Create the new BinnedData that we will fill.
    int nbins(50);
    likely::AbsBinningCPtr
        rBins(new likely::UniformBinning(2,202,nbins)),
        ellBins(new likely::UniformSampling(0,0,1)), // only monopole for now
        zBins(new likely::UniformSampling(zref,zref,1));
    baofit::AbsCorrelationDataPtr
        prototype(new baofit::MultipoleCorrelationData(rBins,ellBins,zBins,rmin,rmax,
            rVetoMin,rVetoMax,cosmo::Monopole,cosmo::Monopole));
        
    // Pre-fill each bin with zero values.
    for(int index = 0; index < nbins; ++index) prototype->setData(index,0);
    
    // Load the specified covariance matrix...
    std::string line;
    int lines;
    
    // import boost spirit parser symbols
    using qi::double_;
    using qi::int_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;

    std::ifstream covIn(covName.c_str());
    if(!covIn.good()) throw RuntimeError("createDR9LRG: unable to open " + covName);
    lines = 0;
    std::vector<double> covRow;
    while(std::getline(covIn,line)) {
        lines++;
        covRow.resize(0);
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                +double_[push_back(ref(covRow),_1)]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("createDR9LRG: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + covName);
        }
        if(covRow.size() != nbins) {
            throw RuntimeError("createDR9LRG: got unexpected number of covariances.");
        }
        int row(lines-1);
        for(int col = 0; col <= row; ++col) prototype->setCovariance(row,col,covRow[col]);
    }
    covIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " covariance values from " << covName << std::endl;
    }    
    
    return prototype;
}

baofit::AbsCorrelationDataPtr
local::loadDR9LRG(std::string const &dataName, baofit::AbsCorrelationDataCPtr prototype, bool verbose) {

    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData((baofit::MultipoleCorrelationData *)(prototype->clone(false)));
    
    // Lookup our reference redshift.
    double zref = prototype->getAxisBinning()[2]->getBinCenter(0);

    std::string line;
    int lines;
    
    // import boost spirit parser symbols
    using qi::double_;
    using qi::int_;
    using qi::_1;
    using phoenix::ref;
    using phoenix::push_back;

    // Loop over lines in the parameter file.
    std::ifstream paramsIn(dataName.c_str());
    if(!paramsIn.good()) throw RuntimeError("loadDR9LRG: Unable to open " + dataName);
    lines = 0;
    double rval,mono;
    std::vector<double> bin(3);
    while(std::getline(paramsIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                double_[ref(rval) = _1] >> double_[ref(mono) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadDR9LRG: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + dataName);
        }
        bin[0] = rval;
        bin[2] = zref;
        bin[1] = 0;
        try {
            int monoIndex = binnedData->getIndex(bin);
            binnedData->setData(monoIndex,mono);
        }
        catch(likely::RuntimeError const &e) {
            // The correlation function has radial bins that go beyond the coverage of
            // our covariance matrix, so stop reading when go beyond that coverage.
            lines--;
            break;
        }
    }
    paramsIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " data values from " << dataName << std::endl;
    }

    return binnedData;
}

// Creates a prototype MultipoleCorrelationData with the specified binning.
baofit::AbsCorrelationDataCPtr local::createFrenchPrototype(double zref, double rmin, double rmax,
double rVetoMin, double rVetoMax, cosmo::Multipole ellmin, cosmo::Multipole ellmax) {
    if(ellmin < cosmo::Monopole || ellmax > cosmo::Quadrupole || ellmin > ellmax) {
        throw RuntimeError("createFrenchPrototype: invalid ell range.");
    }
    // Create the new BinnedData that we will fill.
    likely::AbsBinningCPtr
        rBins(new likely::UniformBinning(0,200,50)),
        zBins(new likely::UniformSampling(zref,zref,1)),
        ellBins(new likely::UniformSampling(cosmo::Monopole,cosmo::Quadrupole,2));
    baofit::AbsCorrelationDataPtr
        prototype(new baofit::MultipoleCorrelationData(rBins,ellBins,zBins,rmin,rmax,
            rVetoMin,rVetoMax,ellmin,ellmax));
    return prototype;
}

// Loads a binned correlation function in French format and returns a shared pointer to
// a MultipoleCorrelationData.
baofit::AbsCorrelationDataPtr
local::loadFrench(std::string const &dataName, baofit::AbsCorrelationDataCPtr prototype,
bool verbose, bool unweighted, bool expanded, bool checkPosDef) {

    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData((baofit::MultipoleCorrelationData *)(prototype->clone(true)));
    
    // Lookup the number of radial bins.
    int nrbins = prototype->getAxisBinning()[0]->getNBins();
    
    // Lookup our reference redshift.
    double zref = prototype->getAxisBinning()[2]->getBinCenter(0);
    
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
        bool ok(false);
        if(expanded) {
            ok = qi::phrase_parse(line.begin(),line.end(),
                (
                    // data uses an expanded format: r mono dmono quad dquad ...
                    // https://trac.sdss3.org/wiki/BOSS/LyaForestsurvey/FPGAnalysis/Unblinding_details#DATAfiles
                    double_[ref(rval) = _1] >> double_[ref(mono) = _1] >> double_ >> double_[ref(quad) = _1]
                ),
                ascii::space);            
        }
        else {
            ok = qi::phrase_parse(line.begin(),line.end(),
                (
                    // mocks use the format: r mono quad ...
                    double_[ref(rval) = _1] >> double_[ref(mono) = _1] >> double_[ref(quad) = _1]
                ),
                ascii::space);
        }
        if(!ok) {
            throw RuntimeError("loadFrench: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        bin[0] = rval;
        bin[2] = zref;
        bin[1] = cosmo::Monopole;
        int monoIndex = binnedData->getIndex(bin);
        binnedData->setData(monoIndex,mono);
        bin[1] = cosmo::Quadrupole;
        int quadIndex = binnedData->getIndex(bin);
        binnedData->setData(quadIndex,quad);
    }
    paramsIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " data values from " << paramsName << std::endl;
    }
    
    if(!unweighted) {
        // Loop over lines in the covariance file.
        std::string covName = paramsName;
        int pos = covName.rfind('/',-1);
        covName.insert(pos+1,"cov_");
        std::ifstream covIn(covName.c_str());
        if(!covIn.good()) throw RuntimeError("Unable to open " + covName);
        lines = 0;
        int index1,index2;
        double cov;
        std::vector<int> bin1(3),bin2(3);
        bin1[2] = bin2[2] = 0;
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
            // Ignore entries above the diagonal since they are duplicates by symmetry.
            if(index1 > index2) continue;
            // Remap file indexing to a BinnedData global index.
            bin1[0] = index1 % nrbins;
            bin2[0] = index2 % nrbins;
            bin1[1] = index1/nrbins;
            bin2[1] = index2/nrbins;
            binnedData->setCovariance(binnedData->getIndex(bin1), binnedData->getIndex(bin2), cov);
        }
        covIn.close();
        if(verbose) {
            std::cout << "Read " << lines << " and stored "
                << binnedData->getCovarianceMatrix()->getNElements()
                << " covariance values from " << covName << std::endl;
        }
        if(checkPosDef) {
            // Check that the covariance is positive definite by triggering an inversion.
            try {
                binnedData->getInverseCovariance(0,0);
            }
            catch(likely::RuntimeError const &e) {
                std::cerr << "### Inverse covariance not positive-definite: " << covName << std::endl;
            }
        }
    }

    return binnedData;
}

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

baofit::AbsCorrelationDataCPtr local::createSectorsPrototype(double zref, double rmin, double rmax,
double muMin, double muMax, double rVetoMin, double rVetoMax) {
    // Initialize the fixed (r,mu,z) binning for this format.
    likely::AbsBinningCPtr
        rBins(new likely::UniformBinning(0,200,50)),
        muBins(new likely::UniformBinning(0,1,50)),
        zBins(new likely::UniformSampling(zref,zref,1));

    baofit::AbsCorrelationDataPtr
        prototype(new baofit::ComovingCorrelationData(rBins,muBins,zBins,rmin,rmax,
            muMin,muMax,rVetoMin,rVetoMax));

    return prototype;    
}

baofit::AbsCorrelationDataPtr local::loadSectors(std::string const &dataName,
baofit::AbsCorrelationDataCPtr prototype, bool verbose) {

    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData((baofit::ComovingCorrelationData *)(prototype->clone(true)));

    int length(0);
    boost::scoped_array<char> buffer;
    {
        // Get the length of the input file.
        std::string paramsName(dataName + ".dat");
        std::ifstream paramsIn(paramsName.c_str(),std::ios::binary);
        paramsIn.seekg(0, std::ios::end);
        length = paramsIn.tellg();
        if(length % 4 != 0) throw RuntimeError("loadSectors: unexpected file length.");

        // Read the input file into memory.
        buffer.reset(new char[length]);
        paramsIn.seekg(0, std::ios::beg);
        paramsIn.read(buffer.get(),length);
        paramsIn.close();
    }
    
    // Parse the input file in memory.
    int *iter = (int*)buffer.get(), *end = iter + length/4;
    while(iter != end) {
        int cmd = *iter++;
        int count = *iter++;
        if(0 == cmd) {
            // Load Cinv.data bin contents.
            bool weighted(true);
            for(int i = 0; i < count; ++i) {
                int index = *iter++;
                double cinvData = *((double*)iter);
                iter += 2;
                binnedData->setData(index, cinvData, weighted);
            }
        }
        else if(1 == cmd) {
            // Load inverse covariances.
            for(int i = 0; i < count; ++i) {
                int index1 = *iter++, index2 = *iter++;
                double cinv = *((double*)iter);
                iter += 2;
                binnedData->setInverseCovariance(index1, index2, cinv);
            }
        }
        else {
            throw RuntimeError("loadSectors: unexpected command byte in input file.");
        }
        if(iter > end) throw RuntimeError("loadSectors: internal parsing error.");
    }

    // Compress our binned data to take advantage of a very sparse (diagonal) covariance matrix.
    binnedData->compress();
    return binnedData;   
}

// Creates a prototype QuasarCorrelationData with the specified binning and cosmology.
baofit::AbsCorrelationDataCPtr local::createCosmolibPrototype(
double minsep, double dsep, int nsep, double minz, double dz, int nz,
double minll, double maxll, double dll, double dll2,
double rmin, double rmax, double muMin, double muMax, double rVetoMin, double rVetoMax, double llmin,
cosmo::AbsHomogeneousUniversePtr cosmology) {

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
    baofit::AbsCorrelationDataPtr
        prototype(new baofit::QuasarCorrelationData(llBins,sepBins,zBins,rmin,rmax,muMin,muMax,llmin,
            rVetoMin,rVetoMax,cosmology));

    return prototype;
}

// Loads a binned correlation function in cosmolib format and returns a BinnedData object.
baofit::AbsCorrelationDataPtr local::loadCosmolib(std::string const &dataName,
baofit::AbsCorrelationDataCPtr prototype, bool verbose, bool icov, bool weighted, int reuseCov,
bool checkPosDef) {

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
        int index = binnedData->getIndex(bin);
        binnedData->setData(index, weighted ? cinvData : data, weighted);        
    }
    paramsIn.close();
    if(false) {
        std::cout << "Read " << binnedData->getNBinsWithData() << " of "
            << binnedData->getNBinsTotal() << " data values from " << paramsName << std::endl;
    }
    

    // Do we need to reuse the covariance estimated for the first realization of this plate?
    std::string covName;
    if(reuseCov>=0) {
        // Parse the data name.
        boost::regex namePattern("([a-zA-Z0-9/_\\.]+/)?([0-9]+_)([0-9]+)\\.cat\\.([0-9]+)");
        boost::match_results<std::string::const_iterator> what;
        if(!boost::regex_match(dataName,what,namePattern)) {
            throw RuntimeError("loadCosmolib: cannot parse name \"" + dataName + "\"");
	}
	std::stringstream covnum;
	covnum << reuseCov;
        covName = what[1]+what[2]+covnum.str()+".cat."+what[4];
    }
    else {
        covName = dataName;
    }
    covName += (icov ? ".icov" : ".cov");

    // Loop over lines in the covariance file.
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
    if(false) {
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

    if(checkPosDef) {
        // Check that the covariance is positive definite by triggering an inversion.
        try {
            binnedData->getInverseCovariance(0,0);
            binnedData->getCovariance(0,0);
        }
        catch(likely::RuntimeError const &e) {
            std::cerr << "### Covariance not positive-definite: " << covName << std::endl;
        }
    }
    // Compress our binned data to take advantage of a potentially sparse covariance matrix.
    binnedData->compress();
    return binnedData;
}

baofit::AbsCorrelationDataCPtr local::createCosmolibXiPrototype(double minz, double dz, int nz,
double minr, double maxr, double nr, bool hasHexadecapole, double rmin, double rmax,
double rVetoMin, double rVetoMax, cosmo::Multipole ellmin, cosmo::Multipole ellmax) {
    if(ellmin > ellmax) {
        throw RuntimeError("createCosmolibXiPrototype: expected ellmin <= ellmax.");
    }
    // Create the new BinnedData that we will fill.
    std::vector<double> ellValues;
    ellValues.push_back(-1);
    ellValues.push_back(cosmo::Monopole);
    ellValues.push_back(cosmo::Quadrupole);
    if(hasHexadecapole) ellValues.push_back(cosmo::Hexadecapole);
    likely::AbsBinningCPtr
        rBins(new likely::UniformSampling(minr,maxr,nr)),
        ellBins(new likely::NonUniformSampling(ellValues)),
        zBins(new likely::UniformSampling(minz+0.5*dz,minz+(nz-0.5)*dz,nz));
    baofit::AbsCorrelationDataPtr
        prototype(new baofit::MultipoleCorrelationData(rBins,ellBins,zBins,rmin,rmax,
            rVetoMin,rVetoMax,ellmin,ellmax));
    return prototype;
}

baofit::AbsCorrelationDataPtr local::loadCosmolibXi(std::string const &dataName,
AbsCorrelationDataCPtr prototype, bool verbose, bool weighted, int reuseCov, bool checkPosDef) {
    
    // Create the new AbsCorrelationData that we will fill.
    baofit::AbsCorrelationDataPtr binnedData((baofit::MultipoleCorrelationData *)(prototype->clone(true)));

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
                double_[ref(data) = _1] >> double_[ref(cinvData) = _1] >> "| XiLR (" >>
                double_[push_back(ref(bin),_1)] >> ',' >> double_[push_back(ref(bin),_1)] >>
                ',' >> double_[push_back(ref(bin),_1)] >> ')'
            ),
            ascii::space);
        if(!ok) {
            throw RuntimeError("loadCosmolibXi: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
        // File contains (z,ell,r) but we need (r,ell,z)
        std::swap(bin[0],bin[2]);
        // File uses 1e30 for "large" r. Map this to r=200.
        if(bin[0] > 200) bin[0] = 200;
        int index = binnedData->getIndex(bin);
        binnedData->setData(index, weighted ? cinvData : data, weighted);
    }
    paramsIn.close();
    if(false) {
        std::cout << "Read " << binnedData->getNBinsWithData() << " of "
            << binnedData->getNBinsTotal() << " data values from " << paramsName << std::endl;
    }
    

    // Do we need to reuse the covariance estimated for the first realization of this plate?
    std::string covName;
    if(reuseCov>=0) {
        // Parse the data name.
        boost::regex namePattern("([a-zA-Z0-9/_\\.]+/)?([0-9]+_)([0-9]+)\\.cat\\.([0-9]+)");
        boost::match_results<std::string::const_iterator> what;
        if(!boost::regex_match(dataName,what,namePattern)) {
            throw RuntimeError("loadCosmolib: cannot parse name \"" + dataName + "\"");
        }
	std::stringstream covnum;
	covnum << reuseCov;
        covName = what[1]+what[2]+covnum.str()+".cat."+what[4];
    }
    else {
        covName = dataName;
    }
    covName += ".icov";

    // Loop over lines in the covariance file.
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
        value = -value; // !?! see line #388 of Observed2Point.cpp
        int index1 = *(binnedData->begin()+offset1), index2 = *(binnedData->begin()+offset2);
        binnedData->setInverseCovariance(index1,index2,value);
    }
    covIn.close();
    if(false) {
        int ndata = binnedData->getNBinsWithData();
        int ncov = (ndata*(ndata+1))/2;
        std::cout << "Read " << lines << " of " << ncov
            << " covariance values from " << covName << std::endl;
    }

    // Check for zero values on the diagonal
    for(likely::BinnedData::IndexIterator iter = binnedData->begin();
    iter != binnedData->end(); ++iter) {
        int index = *iter;
        if(0 == binnedData->getInverseCovariance(index,index)) {
            binnedData->setInverseCovariance(index,index,1e-30);
        }
    }

    if(checkPosDef) {
        // Check that the covariance is positive definite by triggering an inversion.
        try {
            binnedData->getInverseCovariance(0,0);
            binnedData->getCovariance(0,0);
        }
        catch(likely::RuntimeError const &e) {
            std::cerr << "### Inverse covariance not positive-definite: " << covName << std::endl;
        }
    }

    // Compress our binned data to take advantage of a potentially sparse covariance matrix.
    binnedData->compress();
    return binnedData;    
}
