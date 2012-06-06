// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/baofit.h"
#include "baofit/boss.h"
#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/smart_ptr.hpp"

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
    double initialAmp,initialScale;
    std::string platelistName,platerootName,modelConfig; //,bootstrapSaveName,bootstrapCurvesName
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
        ("dr9lrg", "3D correlation data files are in the BOSS DR9 LRG galaxy format.")
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
        ("model-config", po::value<std::string>(&modelConfig)->default_value(""),
            "Model parameters configuration script.")
        ("fix-alpha", "Fix linear bias parameter alpha.")
        ("fix-beta", "Fix linear bias parameter beta.")
        ("fix-bias", "Fix linear bias parameter (1+beta)*bias.")
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
        fixBias(vm.count("fix-bias")), fixBao(vm.count("fix-bao")), fixScale(vm.count("fix-scale")),
        noBBand(vm.count("no-bband")), fixCovariance(0 == vm.count("naive-covariance")),
        dr9lrg(vm.count("dr9lrg")), fixBeta(vm.count("fix-beta"));
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

    // Initialize the models we will use.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    baofit::AbsCorrelationModelPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        cosmology.reset(new cosmo::LambdaCdmRadiationUniverse(OmegaMatter,0,hubbleConstant));
        
        // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
        model.reset(new baofit::BaoCorrelationModel(
            modelrootName,fiducialName,nowigglesName,broadbandName,zref,
            initialAmp,initialScale,fixAlpha,fixBeta,fixBias,fixBao,fixScale,noBBand));
             
        // Configure our fit model parameters, if requested.
         if(0 < modelConfig.length()) model->configure(modelConfig);

        if(verbose) std::cout << "Models initialized." << std::endl;
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR during model initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << "ERROR during model initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    if(verbose) model->printToStream(std::cout);
    analyzer.setModel(model);
    
    // Load the data we will fit.
    double zdata;
    try {
        
        // Create a prototype of the binned data we will be loading.
        baofit::AbsCorrelationDataCPtr prototype;
        if(french) {
            zdata = 2.30;
            prototype = baofit::boss::createFrenchPrototype(zdata,rmin,rmax);
        }
        else if(dr9lrg) {
            zdata = 0.57;
            prototype = baofit::boss::createDR9LRGPrototype(zdata,rmin,rmax,
                "LRG/Sample4_North.cov",verbose);
        }
        else {
            zdata = 2.25;
            prototype = baofit::boss::createCosmolibPrototype(
                minsep,dsep,nsep,minz,dz,nz,minll,dll,dll2,nll,rmin,rmax,llmin,cosmology);
        }
        
        // Build a list of the data files we will read.
        std::vector<std::string> filelist;
        if(0 < dataName.length()) {
            // Load a single named file specified by --data.
            filelist.push_back(dataName);
        }
        else {
            // Load individual plate files specified by --plateroot and --platelist.
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
                filelist.push_back(boost::str(platefile % platerootName % plateName));
                if(filelist.size() == maxPlates) break;
            }
            platelist.close();
            if(verbose) {
                std::cout << "Read " << filelist.size() << " entries from " << platelistName << std::endl;
            }
        }
        
        // Load each file into our analyzer.
        for(std::vector<std::string>::const_iterator filename = filelist.begin();
        filename != filelist.end(); ++filename) {
            if(french) {
                analyzer.addData(baofit::boss::loadFrench(*filename,prototype,verbose));
            }
            else if(dr9lrg) {
                analyzer.addData(baofit::boss::loadDR9LRG(*filename,prototype,verbose));
            }
            else {
                // Add a cosmolib dataset, assumed to provided icov instead of cov.
                analyzer.addData(baofit::boss::loadCosmolib(*filename,prototype,verbose,true));
            }            
        }
    }
    catch(baofit::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }

    if(french || dr9lrg) {
        std::ofstream out("monopole.dat");
        boost::shared_ptr<baofit::MultipoleCorrelationData> combined =
            boost::dynamic_pointer_cast<baofit::MultipoleCorrelationData>(analyzer.getCombined());
        combined->dump(out,cosmo::Monopole);
        out.close();
    }
    
    // Do the requested analysis.
    try {
        lk::FunctionMinimumPtr fmin = analyzer.fitCombined("mn2::vmetric");        
        if(bootstrapTrials > 0) {
            analyzer.doBootstrapAnalysis("mn2::vmetric",fmin,bootstrapTrials,bootstrapSize,fixCovariance);
        }
        {
            // Dump the best-fit monopole model.
            std::ofstream out("fitmono.dat");
            analyzer.dump(out,fmin,cosmo::Monopole,100,rmin,rmax,zdata);
            out.close();
        }
        {
            // Dump the best-fit monopole model without its peak contribution.
            fmin->setParameterValue("BAO amplitude",0);
            std::ofstream out("fitmono-smooth.dat");
            analyzer.dump(out,fmin,cosmo::Monopole,100,rmin,rmax,zdata);
            out.close();
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