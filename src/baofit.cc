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
    
    // Configure option processing
    po::options_description allOptions("Fits cosmological data to measure baryon acoustic oscillations"),
        genericOptions("Generic options"),modelOptions("Model options"), dataOptions("Data options"),
        frenchOptions("French data options"), cosmolibOptions("Cosmolib data options"),
        analysisOptions("Analysis options");

    double OmegaMatter,hubbleConstant,zref,minll,dll,dll2,minsep,dsep,minz,dz,rmin,rmax,llmin;
    int nll,nsep,nz,maxPlates,bootstrapTrials,bootstrapSize,randomSeed,numXi,ndump,jackknifeDrop;
    std::string modelrootName,fiducialName,nowigglesName,broadbandName,dataName,
        platelistName,platerootName,modelConfig,iniName,refitConfig,minMethod;

    // Default values in quotes below are to avoid roundoff errors leading to ugly --help
    // messages. See http://stackoverflow.com/questions/1734916/
    genericOptions.add_options()
        ("help,h", "Prints this info and exits.")
        ("quiet,q", "Runs silently unless there is a problem.")
        ("ini-file,i", po::value<std::string>(&iniName)->default_value(""),
            "Loads options from specified INI file (command line has priority).")
        ;
    modelOptions.add_options()
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0.27,"0.27"),
            "Present-day value of OmegaMatter.")
        ("hubble-constant", po::value<double>(&hubbleConstant)->default_value(0.7,"0.7"),
            "Present-day value of the Hubble parameter h = H0/(100 km/s/Mpc).")
        ("fiducial", po::value<std::string>(&fiducialName)->default_value(""),
            "Fiducial correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("nowiggles", po::value<std::string>(&nowigglesName)->default_value(""),
            "No-wiggles correlation functions will be read from <name>.<ell>.dat with ell=0,2,4.")
        ("broadband", po::value<std::string>(&broadbandName)->default_value(""),
            "Broadband models will be read from <name>bb<x>.<ell>.dat with x=c,1,2 and ell=0,2,4.")
        ("modelroot", po::value<std::string>(&modelrootName)->default_value(""),
            "Common path to prepend to all model filenames.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift used by model correlation functions.")
        ("xi-model", "Uses experimental binned correlation model.")
        ("num-xi", po::value<int>(&numXi)->default_value(9),
            "Number of points from rmin-rmax to use for interpolating xi(r)")
        ("model-config", po::value<std::string>(&modelConfig)->default_value(""),
            "Model parameters configuration script.")
        ;
    dataOptions.add_options()
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D correlation data will be read from the specified file.")
        ("platelist", po::value<std::string>(&platelistName)->default_value(""),
            "3D correlation data will be read from individual plate datafiles listed in this file.")
        ("plateroot", po::value<std::string>(&platerootName)->default_value(""),
            "Common path to prepend to all plate datafiles listed in the platelist.")
        ("french", "3D correlation data files are in the French format (default is cosmolib).")
        ("dr9lrg", "3D correlation data files are in the BOSS DR9 LRG galaxy format.")
        ("max-plates", po::value<int>(&maxPlates)->default_value(0),
            "Maximum number of plates to load (zero uses all available plates).")
        ("check-posdef", "Checks that each covariance is positive-definite (slow).")
        ;
    frenchOptions.add_options()
        ("unweighted", "Does not read covariance data.")
        ("use-quad", "Uses quadrupole correlations.")
        ;
    cosmolibOptions.add_options()
        ("weighted", "Data vectors are inverse-covariance weighted.")
        ("reuse-cov", "Reuse covariance estimated for first-realization of each plate.")
        ("minll", po::value<double>(&minll)->default_value(0.0002,"0.0002"),
            "Minimum log(lam2/lam1).")
        ("dll", po::value<double>(&dll)->default_value(0.004,"0.004"),
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
        ("minz", po::value<double>(&minz)->default_value(1.7,"1.7"),
            "Minimum redshift.")
        ("dz", po::value<double>(&dz)->default_value(1.0,"1.0"),
            "Redshift binsize.")
        ("nz", po::value<int>(&nz)->default_value(2),
            "Maximum number of redshift bins.")
        ;
    analysisOptions.add_options()
        ("rmin", po::value<double>(&rmin)->default_value(0),
            "Minimum 3D comoving separation (Mpc/h) to use in fit.")
        ("rmax", po::value<double>(&rmax)->default_value(200),
            "Maximum 3D comoving separation (Mpc/h) to use in fit.")
        ("llmin", po::value<double>(&llmin)->default_value(0),
            "Minimum value of log(lam2/lam1) to use in fit (cosmolib only).")
        ("ndump", po::value<int>(&ndump)->default_value(100),
            "Number of points spanning [rmin,rmax] to use for dumping models (zero for no dumps).")
        ("refit-config", po::value<std::string>(&refitConfig)->default_value(""),
            "Script to modify parameters for refits.")
        ("fit-each", "Fits each observation separately.")
        ("bootstrap-trials", po::value<int>(&bootstrapTrials)->default_value(0),
            "Number of bootstrap trials to run if a platelist was provided.")
        ("bootstrap-size", po::value<int>(&bootstrapSize)->default_value(0),
            "Size of each bootstrap trial or zero to use the number of plates.")
        ("jackknife-drop", po::value<int>(&jackknifeDrop)->default_value(0),
            "Number of observations to drop from each jackknife sample (zero for no jackknife analysis)")
        ("random-seed", po::value<int>(&randomSeed)->default_value(1966),
            "Random seed to use for generating bootstrap samples.")
        ("min-method", po::value<std::string>(&minMethod)->default_value("mn2::vmetric"),
            "Minimization method to use for fitting.")
        ;

    allOptions.add(genericOptions).add(modelOptions).add(dataOptions)
        .add(frenchOptions).add(cosmolibOptions).add(analysisOptions);
    po::variables_map vm;

    // Parse command line options first.
    try {
        po::store(po::parse_command_line(argc, argv, allOptions), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << allOptions << std::endl;
        return 1;
    }
    // If an INI file was specified, load it now.
    if(0 < iniName.length()) {
        try {
            std::ifstream iniFile(iniName.c_str());
            po::store(po::parse_config_file(iniFile, allOptions), vm);
            iniFile.close();
            po::notify(vm);
        }
        catch(std::exception const &e) {
            std::cerr << "Unable to parse INI file options: " << e.what() << std::endl;
            return -1;
        }        
    }
    
    // Extract boolean options.
    bool verbose(0 == vm.count("quiet")), french(vm.count("french")), weighted(vm.count("weighted")),
        checkPosDef(vm.count("check-posdef")), fixCovariance(0 == vm.count("naive-covariance")),
        xiModel(vm.count("xi-model")), dr9lrg(vm.count("dr9lrg")), unweighted(vm.count("unweighted")),
        useQuad(vm.count("use-quad")), fitEach(vm.count("fit-each")), reuseCov(vm.count("reuse-cov"));

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
    baofit::CorrelationAnalyzer analyzer(randomSeed,minMethod,rmin,rmax,verbose);

    // Initialize the models we will use.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    baofit::AbsCorrelationModelPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        cosmology.reset(new cosmo::LambdaCdmRadiationUniverse(OmegaMatter,0,hubbleConstant));
        
        if(xiModel) {
            likely::AbsBinningCPtr rbins(new likely::UniformSampling(rmin,rmax,numXi));
            model.reset(new baofit::XiCorrelationModel(rbins,zref,"cspline"));
        }
        else {
            // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
            model.reset(new baofit::BaoCorrelationModel(
                modelrootName,fiducialName,nowigglesName,broadbandName,zref));
        }
             
        // Configure our fit model parameters, if requested.
         if(0 < modelConfig.length()) model->configureFitParameters(modelConfig);

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
            prototype = baofit::boss::createFrenchPrototype(zdata,rmin,rmax,useQuad);
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
                analyzer.addData(baofit::boss::loadFrench(*filename,prototype,
                    verbose,unweighted,checkPosDef));
            }
            else if(dr9lrg) {
                analyzer.addData(baofit::boss::loadDR9LRG(*filename,prototype,verbose));
            }
            else {
                // Add a cosmolib dataset, assumed to provided icov instead of cov.
                analyzer.addData(baofit::boss::loadCosmolib(*filename,prototype,
                    verbose,true,weighted,reuseCov,checkPosDef));
            }
        }
    }
    catch(baofit::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }
    analyzer.setZData(zdata);

    if(french || dr9lrg) {
        std::ofstream out("combined.dat");
        boost::shared_ptr<baofit::MultipoleCorrelationData> combined =
            boost::dynamic_pointer_cast<baofit::MultipoleCorrelationData>(analyzer.getCombined());
        combined->dump(out);
        out.close();
    }

    // Do the requested analyses...
    try {
        // Always fit the combined sample.
        lk::FunctionMinimumPtr fmin = analyzer.fitCombined();
        if(ndump > 0) {
            // Dump the best-fit monopole model.
            std::ofstream out("fit.dat");
            analyzer.dumpModel(out,fmin,ndump);
            out.close();
        }
        if(!xiModel && ndump > 0) {
            // Dump the best-fit monopole model with its peak contribution forced to zero.
            std::ofstream out("fit-smooth.dat");
            analyzer.dumpModel(out,fmin,ndump,"value[BAO amplitude]=0");
            out.close();
        }
        {
            // Dump the best-fit residuals for each data bin.
            std::ofstream out("residuals.dat");
            analyzer.dumpResiduals(out,fmin);
            out.close();
        }
        // Refit the combined sample, if requested.
        lk::FunctionMinimumPtr fmin2;
        if(0 < refitConfig.size()) {
            if(verbose) {
                std::cout << std::endl << "Re-fitting combined with: " << refitConfig << std::endl;
            }
            fmin2 = analyzer.fitCombined(refitConfig);
            if(ndump > 0) {
                // Dump the best-fit monopole model.
                std::ofstream out("refit.dat");
                analyzer.dumpModel(out,fmin2,ndump);
                out.close();
            }
        }
        // Perform a bootstrap analysis, if requested.
        if(bootstrapTrials > 0) {
            analyzer.doBootstrapAnalysis(bootstrapTrials,bootstrapSize,fixCovariance,
                fmin,fmin2,refitConfig,"bs.dat",ndump);
        }
        // Perform a jackknife analysis, if requested.
        if(jackknifeDrop > 0) {
            analyzer.doJackknifeAnalysis(jackknifeDrop,fmin,fmin2,refitConfig,"jk.dat",ndump);
        }
        // Fit each observation separately, if requested.
        if(fitEach) {
            analyzer.fitEach(fmin,fmin2,refitConfig,"each.dat",ndump);
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