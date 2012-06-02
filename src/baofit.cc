// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/baofit.h"
#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"
#include "boost/pointer_cast.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

namespace lk = likely;
namespace po = boost::program_options;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

std::vector<double> twoStepSampling(
int nBins, double breakpoint,double dlog, double dlin, double eps = 1e-3) {
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
void loadFrench(std::string dataName, bool verbose = true) {

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
    if(!paramsIn.good()) throw baofit::RuntimeError("loadFrench: Unable to open " + paramsName);
    lines = 0;
    double rval,mono,quad;
    while(std::getline(paramsIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                double_[ref(rval) = _1] >> double_[ref(mono) = _1] >> double_[ref(quad) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw baofit::RuntimeError("loadFrench: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + paramsName);
        }
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
    if(!covIn.good()) throw cosmo::RuntimeError("Unable to open " + covName);
    lines = 0;
    int index1,index2;
    double cov;
    likely::CovarianceMatrix Cboth(100), Cmono(50), Cquad(50);
    while(std::getline(covIn,line)) {
        lines++;
        bool ok = qi::phrase_parse(line.begin(),line.end(),
            (
                int_[ref(index1) = _1] >> int_[ref(index2) = _1] >> double_[ref(cov) = _1]
            ),
            ascii::space);
        if(!ok) {
            throw baofit::RuntimeError("loadFrench: error reading line " +
                boost::lexical_cast<std::string>(lines) + " of " + covName);
        }
        if(index1 <= index2) Cboth.setCovariance(index1,index2,cov);
        if(index1 <= index2 && index2 < 50) Cmono.setCovariance(index1,index2,cov);
        if(index1 <= index2 && index1 >= 50) Cquad.setCovariance(index1-50,index2-50,cov);
            
    }
    covIn.close();
    if(verbose) {
        std::cout << "Read " << lines << " covariance values from " << covName << std::endl;
    }
    // Try to invert the matrices...
    try {
        Cmono.getInverseCovariance(0,0);
        std::cout << "mono ok" << std::endl;
        Cquad.getInverseCovariance(0,0);
        std::cout << "quad ok" << std::endl;
        Cboth.getInverseCovariance(0,0);
        std::cout << "both ok" << std::endl;
    }
    catch(likely::RuntimeError const &e) {
        std::cout << "At least one covariance matrix is not positive definite." << std::endl;
    }
}

// Loads a binned correlation function in cosmolib format and returns a BinnedData object.
// The fast option disables regexp checks for valid numeric inputs.
baofit::QuasarCorrelationDataPtr loadCosmolib(std::string dataName,
    likely::AbsBinningCPtr llBins, likely::AbsBinningCPtr sepBins, likely::AbsBinningCPtr zBins,
    double rmin, double rmax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology,
    bool verbose, bool icov = false, bool fast = false) {

    // Create the new BinnedData that we will fill.
    baofit::QuasarCorrelationDataPtr
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
    if(!paramsIn.good()) throw baofit::RuntimeError("loadCosmolib: Unable to open " + paramsName);
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
            throw baofit::RuntimeError("loadCosmolib: error reading line " +
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
    if(!covIn.good()) throw cosmo::RuntimeError("Unable to open " + covName);
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
            throw baofit::RuntimeError("loadCosmolib: error reading line " +
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

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("BAO fitting");
    double OmegaMatter,hubbleConstant,zref,minll,dll,dll2,minsep,dsep,minz,dz,rmin,rmax,llmin;
    int nll,nsep,nz,ncontour,modelBins,maxPlates,bootstrapTrials,bootstrapSize,randomSeed;
    std::string modelrootName,fiducialName,nowigglesName,broadbandName,dataName,dumpName;
    double initialAmp,initialScale;
    std::string platelistName,platerootName,bootstrapSaveName,bootstrapCurvesName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.27),
            "Present-day value of OmegaMatter.")
        ("hubble-constant", po::value<double>(&hubbleConstant)->default_value(0.7),
            "Present-day value of the Hubble parameter h = H0/(100 km/s/Mpc).")
        ("modelroot", po::value<std::string>(&modelrootName)->default_value(""),
                "Common path to prepend to all model filenames.")
        ("fiducial", po::value<std::string>(&fiducialName)->default_value(""),
            "Fiducial correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("nowiggles", po::value<std::string>(&nowigglesName)->default_value(""),
            "No-wiggles correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("broadband", po::value<std::string>(&broadbandName)->default_value(""),
            "Broadband models will be read from <name>bb<x>.<ell>.dat with x=c,1,2 and ell=0,2,4.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift.")
        ("rmin", po::value<double>(&rmin)->default_value(0),
            "Minimum 3D comoving separation (Mpc/h) to use in fit.")
        ("rmax", po::value<double>(&rmax)->default_value(200),
            "Maximum 3D comoving separation (Mpc/h) to use in fit.")
        ("llmin", po::value<double>(&llmin)->default_value(0),
            "Minimum value of log(lam2/lam1) to use in fit.")
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D covariance data will be read from <data>.params and <data>.cov")
        ("platelist", po::value<std::string>(&platelistName)->default_value(""),
            "3D covariance data will be read from individual plate datafiles listed in this file.")
        ("plateroot", po::value<std::string>(&platerootName)->default_value(""),
            "Common path to prepend to all plate datafiles listed in the platelist.")
        ("max-plates", po::value<int>(&maxPlates)->default_value(0),
            "Maximum number of plates to load (zero uses all available plates).")
        ("fast-load", "Bypasses numeric input validation when reading data.")
        ("bootstrap-trials", po::value<int>(&bootstrapTrials)->default_value(0),
            "Number of bootstrap trials to run if a platelist was provided.")
        ("bootstrap-size", po::value<int>(&bootstrapSize)->default_value(0),
            "Size of each bootstrap trial or zero to use the number of plates.")
        ("bootstrap-save", po::value<std::string>(&bootstrapSaveName)->default_value("bstrials.txt"),
            "Name of file to write with results of each bootstrap trial.")
        ("bootstrap-curves", po::value<std::string>(&bootstrapCurvesName)->default_value(""),
            "Name of file to write individual bootstrap fit multipole curves to.")
        ("naive-covariance", "Uses the naive covariance matrix for each bootstrap trial.")
        ("null-hypothesis", "Applies theory offsets to simulate the null hypothesis.")
        ("random-seed", po::value<int>(&randomSeed)->default_value(1966),
            "Random seed to use for generating bootstrap samples.")
        ("minll", po::value<double>(&minll)->default_value(0.0002),
            "Minimum log(lam2/lam1).")
        ("dll", po::value<double>(&dll)->default_value(0.004),
            "log(lam2/lam1) binsize.")
        ("dll2", po::value<double>(&dll2)->default_value(0),
            "log(lam2/lam1) second binsize parameter for two-step binning.")
        ("nll", po::value<int>(&nll)->default_value(14),
            "Maximum number of log(lam2/lam1) bins.")
        ("minsep", po::value<double>(&minsep)->default_value(0),
            "Minimum separation in arcmins.")
        ("dsep", po::value<double>(&dsep)->default_value(10),
            "Separation binsize in arcmins.")
        ("nsep", po::value<int>(&nsep)->default_value(14),
            "Maximum number of separation bins.")
        ("minz", po::value<double>(&minz)->default_value(1.7),
            "Minimum redshift.")
        ("dz", po::value<double>(&dz)->default_value(1.0),
            "Redshift binsize.")
        ("nz", po::value<int>(&nz)->default_value(2),
            "Maximum number of redshift bins.")
        ("dump", po::value<std::string>(&dumpName)->default_value(""),
            "Filename for dumping fit results.")
        ("ncontour",po::value<int>(&ncontour)->default_value(0),
            "Number of contour points to calculate in BAO parameters.")
        ("model-bins", po::value<int>(&modelBins)->default_value(200),
            "Number of high-resolution uniform bins to use for dumping best fit model.")
        ("minos", "Runs MINOS to improve error estimates.")
        ("fix-linear", "Fix linear bias parameters alpha, bias, beta.")
        ("fix-bao", "Fix BAO scale and amplitude parameters.")
        ("fix-scale", "Fix BAO scale parameter (amplitude floating).")
        ("no-bband", "Do not add any broadband contribution to the correlation function.")
        ("initial-amp", po::value<double>(&initialAmp)->default_value(1),
            "Initial value for the BAO amplitude parameter.")
        ("initial-scale", po::value<double>(&initialScale)->default_value(1),
            "Initial value for the BAO scale parameter.")
        ;

    // Do the command line parsing now.
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose")), minos(vm.count("minos")), fastLoad(vm.count("fast-load")),
        fixLinear(vm.count("fix-linear")), fixBao(vm.count("fix-bao")), fixScale(vm.count("fix-scale")),
        noBBand(vm.count("no-bband")), fixCovariance(0 == vm.count("naive-covariance")),
        nullHypothesis(vm.count("null-hypothesis"));

    // Check for the required filename parameters.
    if(0 == dataName.length() && 0 == platelistName.length()) {
        std::cerr << "Missing required parameter --data or --platelist." << std::endl;
        return -1;
    }
    if(0 == fiducialName.length()) {
        std::cerr << "Missing required parameter --fiducial." << std::endl;
        return -1;
    }
    if(0 == nowigglesName.length()) {
        std::cerr << "Missing required parameter --nowiggles." << std::endl;
        return -1;
    }
    if(0 == broadbandName.length()) {
        std::cerr << "Missing required parameter --broadband." << std::endl;
        return -1;
    }

    // Initialize our analyzer.
    baofit::CorrelationAnalyzer analyzer(randomSeed,verbose);

    // Initialize the cosmology calculations we will need.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    baofit::AbsCorrelationModelCPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        cosmology.reset(new cosmo::LambdaCdmRadiationUniverse(OmegaMatter,0,hubbleConstant));
        
         // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
         model.reset(new baofit::BaoCorrelationModel(
             modelrootName,fiducialName,nowigglesName,broadbandName,zref,
             initialAmp,initialScale,fixLinear,fixBao,fixScale,noBBand));

        if(verbose) std::cout << "Cosmology initialized." << std::endl;
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during cosmology initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    analyzer.setModel(model);
    
    loadFrench("APC/CF_for_tim/CF_M02_000_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_001_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_002_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_003_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_004_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_005_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_006_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_007_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_008_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_009_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_010_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_011_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_012_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_013_JK2D_fits");
    loadFrench("APC/CF_for_tim/CF_M02_014_JK2D_fits");
    
    // Load the data we will fit.
    try {
        // Initialize the (logLambda,separation,redshift) binning from command-line params.
        likely::AbsBinningCPtr llBins,
            sepBins(new likely::UniformBinning(minsep,minsep+nsep*dsep,nsep)),
            zBins(new likely::UniformSampling(minz+0.5*dz,minz+(nz-0.5)*dz,nz));
        if(0 == dll2) {
            llBins.reset(new likely::UniformBinning(minll,minll+nll*dll,nll));
        }
        else {
            llBins.reset(new likely::NonUniformSampling(twoStepSampling(nll,minll,dll,dll2)));
        }
        // Initialize the dataset we will fill.
        if(0 < dataName.length()) {
            // Load a single dataset.
            analyzer.addData(loadCosmolib(dataName,llBins,sepBins,zBins,
                rmin,rmax,llmin,cosmology,verbose,false,fastLoad));
        }
        else {
            // Load individual plate datasets.
            std::string plateName;
            boost::format platefile("%s%s");
            platelistName = platerootName + platelistName;
            std::ifstream platelist(platelistName.c_str());
            if(!platelist.good()) {
                std::cerr << "Unable to open platelist file " << platelistName << std::endl;
                return -1;
            }
            while(platelist.good() && !platelist.eof()) {
                platelist >> plateName;
                if(platelist.eof()) break;
                if(!platelist.good()) {
                    std::cerr << "Error while reading platelist from " << platelistName << std::endl;
                    return -1;
                }
                std::string filename(boost::str(platefile % platerootName % plateName));
                analyzer.addData(loadCosmolib(filename,llBins,sepBins,zBins,
                    rmin,rmax,llmin,cosmology,verbose,true,fastLoad));
                if(analyzer.getNData() == maxPlates) break;
            }
            platelist.close();
        }
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Do the requested analysis.
    try {
        lk::FunctionMinimumPtr fmin = analyzer.fitCombined("mn2::vmetric");
        if(bootstrapTrials > 0) {
            analyzer.doBootstrapAnalysis("mn2::vmetric",fmin,bootstrapTrials,bootstrapSize,fixCovariance);
        }
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }

    // All done: normal exit.
    return 0;
}