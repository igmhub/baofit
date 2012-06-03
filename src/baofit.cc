// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/baofit.h"
#include "baofit/boss.h"
#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace lk = likely;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("BAO fitting");
    double OmegaMatter,hubbleConstant,zref,minll,dll,dll2,minsep,dsep,minz,dz,rmin,rmax,llmin;
    int nll,nsep,nz,maxPlates,bootstrapTrials,bootstrapSize,randomSeed; //,ncontour,modelBins
    std::string modelrootName,fiducialName,nowigglesName,broadbandName,dataName; //,dumpName
    double initialAmp,initialScale,zfrench;
    std::string platelistName,platerootName; //,bootstrapSaveName,bootstrapCurvesName
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
            "Reference redshift used by model correlation functions.")
        ("rmin", po::value<double>(&rmin)->default_value(0),
            "Minimum 3D comoving separation (Mpc/h) to use in fit.")
        ("rmax", po::value<double>(&rmax)->default_value(200),
            "Maximum 3D comoving separation (Mpc/h) to use in fit.")
        ("llmin", po::value<double>(&llmin)->default_value(0),
            "Minimum value of log(lam2/lam1) to use in fit.")
        ("french", "3D correlation data files are in the French format (default is cosmolib).")
        ("zfrench", po::value<double>(&zfrench)->default_value(2.30),
            "Reference redshift used in French 3D correlation data")
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D correlation data will be read from the specified file.")
        ("platelist", po::value<std::string>(&platelistName)->default_value(""),
            "3D correlation data will be read from individual plate datafiles listed in this file.")
        ("plateroot", po::value<std::string>(&platerootName)->default_value(""),
            "Common path to prepend to all plate datafiles listed in the platelist.")
        ("max-plates", po::value<int>(&maxPlates)->default_value(0),
            "Maximum number of plates to load (zero uses all available plates).")
        ("bootstrap-trials", po::value<int>(&bootstrapTrials)->default_value(0),
            "Number of bootstrap trials to run if a platelist was provided.")
        ("bootstrap-size", po::value<int>(&bootstrapSize)->default_value(0),
            "Size of each bootstrap trial or zero to use the number of plates.")
        /**
        ("bootstrap-save", po::value<std::string>(&bootstrapSaveName)->default_value("bstrials.txt"),
            "Name of file to write with results of each bootstrap trial.")
        ("bootstrap-curves", po::value<std::string>(&bootstrapCurvesName)->default_value(""),
            "Name of file to write individual bootstrap fit multipole curves to.")
        ("naive-covariance", "Uses the naive covariance matrix for each bootstrap trial.")
        ("null-hypothesis", "Applies theory offsets to simulate the null hypothesis.")
        **/
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
        /**
        ("dump", po::value<std::string>(&dumpName)->default_value(""),
            "Filename for dumping fit results.")
        ("ncontour",po::value<int>(&ncontour)->default_value(0),
            "Number of contour points to calculate in BAO parameters.")
        ("model-bins", po::value<int>(&modelBins)->default_value(200),
            "Number of high-resolution uniform bins to use for dumping best fit model.")
        ("minos", "Runs MINOS to improve error estimates.")
        **/
        ("fix-alpha", "Fix linear bias parameter alpha.")
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
    bool verbose(vm.count("verbose")), french(vm.count("french")), fixAlpha(vm.count("fix-alpha")),
        fixLinear(vm.count("fix-linear")), fixBao(vm.count("fix-bao")), fixScale(vm.count("fix-scale")),
        noBBand(vm.count("no-bband")), fixCovariance(0 == vm.count("naive-covariance"));
    // minos(vm.count("minos")), nullHypothesis(vm.count("null-hypothesis"))

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

    // Initialize the cosmology models we will use.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    baofit::AbsCorrelationModelCPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        cosmology.reset(new cosmo::LambdaCdmRadiationUniverse(OmegaMatter,0,hubbleConstant));
        
         // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
         model.reset(new baofit::BaoCorrelationModel(
             modelrootName,fiducialName,nowigglesName,broadbandName,zref,
             initialAmp,initialScale,fixAlpha,fixLinear,fixBao,fixScale,noBBand));

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
    
    // Load the data we will fit.
    try {
        
        baofit::AbsCorrelationDataCPtr prototype;
        if(french) {
            prototype = baofit::boss::createFrenchPrototype(zfrench);
        }
        else {
            prototype = baofit::boss::createCosmolibPrototype(
                minsep,dsep,nsep,minz,dz,nz,minll,dll,dll2,nll,rmin,rmax,llmin,cosmology);
        }
        
        // Initialize the dataset we will fill.
        if(0 < dataName.length()) {
            if(french) {
                analyzer.addData(baofit::boss::loadFrench(dataName,prototype,verbose));
            }
            else {
                // Load a single cosmolib dataset, assumed to provide cov instead of icov.
                analyzer.addData(baofit::boss::loadCosmolib(dataName,prototype,verbose,false));
            }
        }
        else {
            // Load individual plate datasets, assumed to provided icov instead of cov.
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
                analyzer.addData(baofit::boss::loadCosmolib(filename,prototype,verbose,true));
                if(analyzer.getNData() == maxPlates) break;
            }
            platelist.close();
        }
    }
    catch(baofit::RuntimeError const &e) {
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