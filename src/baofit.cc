// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/baofit.h"
#include "baofit/boss.h"
#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/foreach.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure option processing
    po::options_description allOptions("Fits cosmological data to measure baryon acoustic oscillations"),
        genericOptions("Generic options"),modelOptions("Model options"), dataOptions("Data options"),
        frenchOptions("French data options"), cosmolibOptions("Cosmolib data options"),
        analysisOptions("Analysis options");

    double OmegaMatter,hubbleConstant,zref,minll,maxll,dll,dll2,minsep,dsep,minz,dz,rmin,rmax,
        rVetoWidth,rVetoCenter,xiRmin,xiRmax,muMin,muMax,kloSpline,khiSpline,toymcScale,saveICovScale,
        zMin,zMax,llMin,llMax,sepMin,sepMax,distR0;
    int nsep,nz,maxPlates,bootstrapTrials,bootstrapSize,randomSeed,ndump,jackknifeDrop,lmin,lmax,
        mcmcSave,mcmcInterval,toymcSamples,xiNr,reuseCov,nSpline,splineOrder,bootstrapCovTrials,
        projectModesNKeep;
    std::string modelrootName,fiducialName,nowigglesName,dataName,xiPoints,toymcConfig,
        platelistName,platerootName,iniName,refitConfig,minMethod,xiMethod,outputPrefix,altConfig,
        fixModeScales,distAdd,distMul;
    std::vector<std::string> modelConfig;

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
        ("modelroot", po::value<std::string>(&modelrootName)->default_value(""),
            "Common path to prepend to all model filenames.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift used by model correlation functions.")
        ("dist-add", po::value<std::string>(&distAdd)->default_value(""),
            "Parameterization to use for additive broadband distortion.")
        ("dist-mul", po::value<std::string>(&distMul)->default_value(""),
            "Parameterization to use for multiplicative broadband distortion.")
        ("dist-r0", po::value<double>(&distR0)->default_value(100),
            "Comoving separation in Mpc/h used to normalize broadband radial factors.")
        ("n-spline", po::value<int>(&nSpline)->default_value(0),
            "Number of spline knots to use spanning (klo,khi).")
        ("klo-spline", po::value<double>(&kloSpline)->default_value(0.02,"0.02"),
            "Minimum k in h/Mpc for k P(k) B-spline.")
        ("khi-spline", po::value<double>(&khiSpline)->default_value(0.2,"0.2"),
            "Maximum k in h/Mpc for k P(k) B-spline.")
        ("order-spline", po::value<int>(&splineOrder)->default_value(3),
            "Order of B-spline in k P(k).")
        ("multi-spline", "Fits independent parameters for each multipole.")
        ("xi-points", po::value<std::string>(&xiPoints)->default_value(""),
            "Comma-separated list of r values (Mpc/h) to use for interpolating r^2 xi(r)")
        ("xi-method", po::value<std::string>(&xiMethod)->default_value("cspline"),
            "Interpolation method to use in r^2 xi(r), use linear or cspline.")
        ("model-config", po::value<std::vector<std::string> >(&modelConfig)->composing(),
            "Model parameters configuration script (option can appear multiple times).")
        ("alt-config", po::value<std::string>(&altConfig)->default_value(""),
            "Parameter adjustments for dumping alternate best-fit model.")
        ("anisotropic", "Uses anisotropic scale parameters instead of an isotropic scale.")
        ("decoupled", "Only applies scale factors to BAO peak and not cosmological broadband.")
        ;
    dataOptions.add_options()
        ("data", po::value<std::string>(&dataName)->default_value(""),
            "3D correlation data will be read from the specified file.")
        ("platelist", po::value<std::string>(&platelistName)->default_value(""),
            "3D correlation data will be read from individual plate datafiles listed in this file.")
        ("plateroot", po::value<std::string>(&platerootName)->default_value(""),
            "Common path to prepend to all plate datafiles listed in the platelist.")
        ("dr9lrg", "3D correlation data files are in the BOSS DR9 LRG galaxy format.")
        ("max-plates", po::value<int>(&maxPlates)->default_value(0),
            "Maximum number of plates to load (zero uses all available plates).")
        ("check-posdef", "Checks that each covariance is positive-definite (slow).")
        ("save-data", "Saves the combined (unweighted) data after final cuts.")
        ("save-icov", "Saves the inverse covariance of the combined data after final cuts.")
        ("save-icov-scale", po::value<double>(&saveICovScale)->default_value(1),
            "Scale factor applied to inverse covariance elements when using save-icov.")
        ("fix-aln-cov", "Fixes covariance matrix of points in 'aln' parametrization")
        ("fix-mode-scales", po::value<std::string>(&fixModeScales)->default_value(""),
            "Fixes covariance matrix using mode scales from the specified file.")
        ("project-modes-keep", po::value<int>(&projectModesNKeep)->default_value(0),
            "Projects combined data onto the largest (nkeep>0) or smallest (nkeep<0) variance modes.")
        ;
    frenchOptions.add_options()
        ("french", "Correlation data files are in the French format (default is cosmolib).")
        ("unweighted", "Does not read covariance data.")
        ("expanded", "Data uses the expanded format.")
        ("use-quad", "Uses quadrupole correlations.")
        ("sectors", "Correlation data file in r-mu format by angular sector.")
        ;
    cosmolibOptions.add_options()
        ("weighted", "Data vectors are inverse-covariance weighted.")
        ("reuse-cov", po::value<int>(&reuseCov)->default_value(-1),
	        "Reuse covariance estimated for n-th realization of each plate (if >=0).")
        ("minll", po::value<double>(&minll)->default_value(0.0002,"0.0002"),
            "Minimum log(lam2/lam1).")
        ("maxll", po::value<double>(&maxll)->default_value(0.02,"0.02"),
            "Maximum log(lam2/lam1).")
        ("dll", po::value<double>(&dll)->default_value(0.004,"0.004"),
            "log(lam2/lam1) binsize.")
        ("dll2", po::value<double>(&dll2)->default_value(0),
            "log(lam2/lam1) second binsize parameter for two-step binning.")
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
        ("saved-format", "Cosmolib data in format written by our save options.")
        ("xi-format", "Cosmolib data in Xi format.")
        ("xi-rmin", po::value<double>(&xiRmin)->default_value(0),
            "Minimum separation in Mpc/h (Xi format only).")
        ("xi-rmax", po::value<double>(&xiRmax)->default_value(200),
            "Minimum separation in Mpc/h (Xi format only).")
        ("xi-nr", po::value<int>(&xiNr)->default_value(41),
            "Number of separation values equally spaced from minr-maxr (Xi format only).")
        ("xi-hexa", "Has hexadecapole (Xi format only).")
        ;
    analysisOptions.add_options()
        ("rmin", po::value<double>(&rmin)->default_value(0),
            "Final cut on minimum 3D comoving separation (Mpc/h) to use in fit.")
        ("rmax", po::value<double>(&rmax)->default_value(200),
            "Final cut on maximum 3D comoving separation (Mpc/h) to use in fit.")
        ("rveto-width", po::value<double>(&rVetoWidth)->default_value(0),
            "Full width (Mpc/h) of co-moving separation window to veto in fit (zero for no veto).")
        ("rveto-center", po::value<double>(&rVetoCenter)->default_value(114),
            "Center (Mpc/h) of co-moving separation window to veto in fit.")
        ("mu-min", po::value<double>(&muMin)->default_value(0),
            "Final cut on minimum value of mu = rL/r to use in the fit (coordinate data only).")
        ("mu-max", po::value<double>(&muMax)->default_value(1),
            "Final cut on maximum value of mu = rL/r to use in the fit (coordinate data only).")
        ("z-min", po::value<double>(&zMin)->default_value(0.),
            "Final cut on minimum value of redshift (coordinate data only).")
        ("z-max", po::value<double>(&zMax)->default_value(10.),
            "Final cut on maximum value of redshift (coordinate data only).")
        ("ll-min", po::value<double>(&llMin)->default_value(0),
            "Minimum value of log(lam2/lam1) to use in fit (Aln2 data only).")
        ("ll-max", po::value<double>(&llMax)->default_value(1),
            "Maximum value of log(lam2/lam1) to use in fit (Aln2 data only).")
        ("sep-min", po::value<double>(&sepMin)->default_value(0),
            "Minimum value of separation in arcmins to use in fit (Aln2 data only).")
        ("sep-max", po::value<double>(&sepMax)->default_value(200),
            "Maximum value of separation in arcmins to use in fit (Aln2 data only).")
        ("lmin", po::value<int>(&lmin)->default_value(0),
            "Final cut on minimum multipole ell (0,2,4) to use in fit (multipole data only).")
        ("lmax", po::value<int>(&lmax)->default_value(2),
            "Final cut on maximum multipole ell (0,2,4) to use in fit (multipole data only).")
        ("output-prefix", po::value<std::string>(&outputPrefix)->default_value(""),
            "Prefix to use for all analysis output files.")
        ("ndump", po::value<int>(&ndump)->default_value(100),
            "Number of points spanning [rmin,rmax] to use for dumping models (zero for no dumps).")
        ("decorrelated", "Combined data is saved with decorrelated errors.")
        ("no-initial-fit", "Skips initial fit to combined sample.")
        ("scalar-weights", "Combine plates using scalar weights instead of Cinv weights.")
        ("refit-config", po::value<std::string>(&refitConfig)->default_value(""),
            "Script to modify parameters for refits.")
        ("compare-each", "Compares each observation to the combined data, before final cuts.")
        ("compare-each-final", "Compares each observation to the combined data, after final cuts.")
        ("fit-each", "Fits each observation separately.")
        ("bootstrap-trials", po::value<int>(&bootstrapTrials)->default_value(0),
            "Number of bootstrap trials to run if a platelist was provided.")
        ("bootstrap-size", po::value<int>(&bootstrapSize)->default_value(0),
            "Size of each bootstrap trial or zero to use the number of plates.")
        ("bootstrap-cov-trials", po::value<int>(&bootstrapCovTrials)->default_value(0),
            "Number of bootstrap trials for estimating and saving combined covariance.")
        ("jackknife-drop", po::value<int>(&jackknifeDrop)->default_value(0),
            "Number of observations to drop from each jackknife sample (zero for no jackknife)")
        ("mcmc-save", po::value<int>(&mcmcSave)->default_value(0),
            "Number of Markov chain Monte Carlo samples to save (zero for no MCMC analysis)")
        ("mcmc-interval", po::value<int>(&mcmcInterval)->default_value(10),
            "Interval for saving MCMC trials (larger for less correlations and longer running time)")
        ("toymc-samples", po::value<int>(&toymcSamples)->default_value(0),
            "Number of toy MC samples to generate and fit.")
        ("toymc-config", po::value<std::string>(&toymcConfig)->default_value(""),
            "Fit parameter configuration to apply before generating toy MC samples.")
        ("toymc-save", "Saves first generated toy MC sample.")
        ("toymc-scale", po::value<double>(&toymcScale)->default_value(1),
            "Scales the covariance used for toy MC noise sampling (but not fitting).")
        ("random-seed", po::value<int>(&randomSeed)->default_value(1966),
            "Random seed to use for generating bootstrap samples.")
        ("min-method", po::value<std::string>(&minMethod)->default_value("mn2::vmetric"),
            "Minimization method to use for fitting.")
        ;

    allOptions.add(genericOptions).add(modelOptions).add(dataOptions)
        .add(frenchOptions).add(cosmolibOptions).add(analysisOptions);
    po::variables_map vm;

    // Parse command line options first so they override anything in an INI file
    // (except for --model-config)
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
    // Make a copy of any command-line model-config options.
    std::vector<std::string> modelConfigSave = modelConfig;
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
    // Shift any command-line model-config options after any INI file options.
    int nSave = modelConfigSave.size();
    std::copy(modelConfig.begin()+nSave,modelConfig.end(),modelConfig.begin());
    std::copy(modelConfigSave.begin(),modelConfigSave.end(),modelConfig.end()-nSave);
    
    // Extract boolean options.
    bool verbose(0 == vm.count("quiet")), french(vm.count("french")), weighted(vm.count("weighted")),
        checkPosDef(vm.count("check-posdef")), fixCovariance(0 == vm.count("naive-covariance")),
        dr9lrg(vm.count("dr9lrg")), unweighted(vm.count("unweighted")), anisotropic(vm.count("anisotropic")),
        fitEach(vm.count("fit-each")), xiHexa(vm.count("xi-hexa")), savedFormat(vm.count("saved-format")),
        xiFormat(vm.count("xi-format")), decorrelated(vm.count("decorrelated")),
        toymcSave(vm.count("toymc-save")), expanded(vm.count("expanded")), sectors(vm.count("sectors")),
        saveICov(vm.count("save-icov")), multiSpline(vm.count("multi-spline")),
        fixAlnCov(vm.count("fix-aln-cov")), saveData(vm.count("save-data")),
        scalarWeights(vm.count("scalar-weights")), noInitialFit(vm.count("no-initial-fit")),
        compareEach(vm.count("compare-each")), compareEachFinal(vm.count("compare-each-final")),
        decoupled(vm.count("decoupled"));

    // Check for the required filename parameters.
    if(0 == dataName.length() && 0 == platelistName.length()) {
        std::cerr << "Missing required parameter --data or --platelist." << std::endl;
        return -1;
    }

    // Check for valid multipole options.
    if(lmin != cosmo::Monopole && lmin != cosmo::Quadrupole && lmin != cosmo::Hexadecapole) {
        std::cerr << "Expected 0,2,4 for lmin but got " << lmin << std::endl;
        return -1;
    }
    cosmo::Multipole ellmin = static_cast<cosmo::Multipole>(lmin);
    if(lmax != cosmo::Monopole && lmax != cosmo::Quadrupole && lmax != cosmo::Hexadecapole) {
        std::cerr << "Expected 0,2,4 for lmax but got " << lmax << std::endl;
        return -1;
    }
    cosmo::Multipole ellmax = static_cast<cosmo::Multipole>(lmax);

    // Calculate veto window.
    double rVetoMin = rVetoCenter - 0.5*rVetoWidth, rVetoMax = rVetoCenter + 0.5*rVetoWidth;

    // Initialize our analyzer.
    likely::Random::instance()->setSeed(randomSeed);
    baofit::CorrelationAnalyzer analyzer(minMethod,rmin,rmax,verbose,scalarWeights);

    // Initialize the fit model we will use.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    baofit::AbsCorrelationModelPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        cosmology.reset(new cosmo::LambdaCdmRadiationUniverse(OmegaMatter,0,hubbleConstant));
        
        if(nSpline > 0) {
            model.reset(new baofit::PkCorrelationModel(modelrootName,nowigglesName,
                kloSpline,khiSpline,nSpline,splineOrder,multiSpline,zref));
        }
        else if(xiPoints.length() > 0) {
            model.reset(new baofit::XiCorrelationModel(xiPoints,zref,xiMethod));
        }
        else {
            // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
            model.reset(new baofit::BaoCorrelationModel(
                modelrootName,fiducialName,nowigglesName,distAdd,distMul,distR0,zref,anisotropic,decoupled));
        }
             
        // Configure our fit model parameters by applying all model-config options in turn,
        // starting with those in the INI file and ending with any command-line options.
        BOOST_FOREACH(std::string const &config, modelConfig) {
            model->configureFitParameters(config);
        }

        if(verbose) std::cout << "Model initialized." << std::endl;
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR during model initialization:\n  " << e.what() << std::endl;
        return -2;
    }
    if(verbose) model->printToStream(std::cout);
    analyzer.setModel(model);
    
    // Load the data we will fit.
    double zdata;
    baofit::AbsCorrelationDataCPtr combined;
    try {
        
        // Create a prototype of the binned data we will be loading.
        baofit::AbsCorrelationDataPtr prototype;
        if(french) {
            zdata = 2.30;
            prototype = baofit::boss::createFrenchPrototype(zdata);
        }
        else if(sectors) {
            zdata = 2.30;
            prototype = baofit::boss::createSectorsPrototype(zdata);
        }
        else if(dr9lrg) {
            zdata = 0.57;
            prototype = baofit::boss::createDR9LRGPrototype(zdata,"LRG/Sample4_North.cov",verbose);
        }
        else if(xiFormat) {
            zdata = 2.25;
            prototype = baofit::boss::createCosmolibXiPrototype(minz,dz,nz,xiRmin,xiRmax,xiNr,xiHexa);
        }
        else { // default is cosmolib (saved) format
            zdata = 2.25;
            prototype = baofit::boss::createCosmolibPrototype(
                minsep,dsep,nsep,minz,dz,nz,minll,maxll,dll,dll2,llMin,llMax,sepMin,sepMax,
                fixAlnCov,cosmology);
        }
        // Set the final cuts that have not already been specified in the prototype ctors above.
        prototype->setFinalCuts(rmin,rmax,rVetoMin,rVetoMax,muMin,muMax,ellmin,ellmax,zMin,zMax);
        
        // Load the vector of mode scale corrections to apply, if any
        std::vector<double> modeScales;
        if(fixModeScales.length() > 0) {
            std::ifstream in(fixModeScales.c_str());
            std::vector<std::vector<double> > vectors(1);
            int nread = likely::readVectors(in, vectors, false);
            in.close();
            modeScales = vectors[0];
            if(verbose) {
                std::cout << "Read " << modeScales.size()
                    << " mode scales from " << fixModeScales << std::endl;
            }
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
                return -3;
            }
            while(platelist.good() && !platelist.eof()) {
                platelist >> plateName;
                if(platelist.eof()) break;
                if(!platelist.good()) {
                    std::cerr << "Error while reading platelist from " << platelistName << std::endl;
                    return -3;
                }
                filelist.push_back(boost::str(platefile % platerootName % plateName));
                if(filelist.size() == maxPlates) break;
            }
            platelist.close();
            if(verbose) {
                std::cout << "Read " << filelist.size() << " entries from "
                    << platelistName << std::endl;
            }
        }
        
        // Load each file into our analyzer.
        for(std::vector<std::string>::const_iterator filename = filelist.begin();
        filename != filelist.end(); ++filename) {
            baofit::AbsCorrelationDataPtr data;
            int reuseCovIndex(-1);
            if(french) {
                data = baofit::boss::loadFrench(*filename,prototype,
                    verbose,unweighted,expanded);
            }
            else if(sectors) {
                data = baofit::boss::loadSectors(*filename,prototype,verbose);
            }
            else if(dr9lrg) {
                data = baofit::boss::loadDR9LRG(*filename,prototype,verbose);
            }
            else if(xiFormat) {
                data = baofit::boss::loadCosmolibXi(*filename,prototype,
                    verbose,weighted,reuseCov);
            }
            else {
                // Add a cosmolib dataset, assumed to provided icov instead of cov.
                if(savedFormat) {
                    data = baofit::boss::loadCosmolibSaved(*filename,prototype,verbose); 
                }
                else {
                    data = baofit::boss::loadCosmolib(*filename,prototype,
                        verbose,true,weighted,reuseCovIndex,reuseCov);
                }
            }
            if(checkPosDef && !data->getCovarianceMatrix()->isPositiveDefinite()) {
                std::cerr << "!!! Covariance matrix not positive-definite for "
                    << *filename << std::endl;
            }
            if(reuseCovIndex < 0 && modeScales.size() > 0) {
                if(verbose) std::cout << "Correcting mode scales..." << std::endl;
                data->rescaleEigenvalues(modeScales);
                {
                    std::string outName = outputPrefix + "fixed.data";
                    std::ofstream out(outName.c_str());
                    data->saveData(out);
                    out.close();
                }
                {
                    std::string outName = outputPrefix + "fixed.icov";
                    std::ofstream out(outName.c_str());
                    data->saveInverseCovariance(out);
                    out.close();
                }
            }
            analyzer.addData(data,reuseCovIndex);
        }
        // Specify the nominal redshift associated with the data.
        analyzer.setZData(zdata);
        // Initialize combined as a read-only pointer to the finalized data to fit...
        if(projectModesNKeep != 0) {
            // Project onto eigenmodes before finalizing.
            if(verbose) std::cout << "Projecting onto modes with nkeep = " << projectModesNKeep << std::endl;
            baofit::AbsCorrelationDataPtr beforeCuts = analyzer.getCombined(false,false);
            beforeCuts->projectOntoModes(projectModesNKeep);
            beforeCuts->finalize();
            combined = beforeCuts;
        }
        else {
            // Fetch the combined data after final cuts.
            combined = analyzer.getCombined(verbose);
        }
        // Check that the combined covariance is positive definite.
        if(!combined->getCovarianceMatrix()->isPositiveDefinite()) {
            std::cerr << "Combined covariance matrix is not positive definite." << std::endl;
            return -3;
        }
        // Save the combined (unweighted) data, if requested.
        if(saveData) {
            std::string outName = outputPrefix + "save.data";
            std::ofstream out(outName.c_str());
            combined->saveData(out);
            out.close();
        }
        // Save the combined inverse covariance, if requested.
        if(saveICov) {
            std::string outName = outputPrefix + "save.icov";
            std::ofstream out(outName.c_str());
            combined->saveInverseCovariance(out);
            out.close();
        }
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -3;
    }

    // Do the requested analyses...
    try {
        if(compareEach && analyzer.getNData() > 1) {
            // Compare each observation to the combined average before final cuts.
            std::cout << "Comparing each observation with combined before final cuts:" << std::endl;
            analyzer.compareEach(outputPrefix + "compare.dat",false);
        }
        if(compareEachFinal && analyzer.getNData() > 1) {
            // Compare each observation to the combined average after final cuts.
            std::cout << "Comparing each observation with combined after final cuts:" << std::endl;
            analyzer.compareEach(outputPrefix + "final_compare.dat",true);
        }
        // Fit the combined sample or use the initial model-config.
        likely::FunctionMinimumPtr fmin;
        if(noInitialFit) {
            baofit::CorrelationFitter fitter(combined,model);
            fmin = fitter.guess();
        }
        else {
            fmin = analyzer.fitSample(combined);
        }
        // Dump the fit parameters in model-config format.
        {
            std::string outName = outputPrefix + "fit.config";
            std::ofstream out(outName.c_str());
            out << likely::fitParametersToScript(fmin->getFitParameters());
            out.close();
        }
        // Dump the fit parameter values and diagonal errors.
        {
            std::string outName = outputPrefix + "save.pars";
            std::ofstream out(outName.c_str());
            fmin->saveParameters(out);
            out.close();
        }
        // Dump the fit parameter covariance matrix.
        if(fmin->hasCovariance()) {
            std::string outName = outputPrefix + "save.pcov";
            std::ofstream out(outName.c_str());
            fmin->saveFloatingParameterCovariance(out);
            out.close();
        }
        // Print out some extra info if this fit has floating "BAO alpha-*"
        // and "gamma-alpha" parameters.
        std::cout << std::endl;
        analyzer.printScaleZEff(fmin,zref,"BAO alpha-iso");
        analyzer.printScaleZEff(fmin,zref,"BAO alpha-parallel");
        analyzer.printScaleZEff(fmin,zref,"BAO alpha-perp");
        // Dump the combined multipole data points with decorrelated errors, if possible.
        if(french || dr9lrg || xiFormat) {
            std::string outName = outputPrefix + "combined.dat";
            std::ofstream out(outName.c_str());
            boost::shared_ptr<const baofit::MultipoleCorrelationData> combinedMultipoles =
                boost::dynamic_pointer_cast<const baofit::MultipoleCorrelationData>(combined);
            std::vector<double> dweights;
            if(decorrelated) {
                // Calculate the decorrelated weights for the combined fit.
                analyzer.getDecorrelatedWeights(combined,fmin->getParameters(),dweights);
            }
            combinedMultipoles->dump(out,rmin,rmax,dweights);
            out.close();
        }
        if(ndump > 0) {
            // Dump the best-fit model.
            std::string outName = outputPrefix + "fit.dat";
            std::ofstream out(outName.c_str());
            analyzer.dumpModel(out,fmin->getFitParameters(),ndump);
            out.close();
        }
        if(ndump > 0 && altConfig.length() > 0) {
            // Dump an alternate best-fit model with some parameters modified (e.g., no BAO features)
            std::string outName = outputPrefix + "alt.dat";
            std::ofstream out(outName.c_str());
            analyzer.dumpModel(out,fmin->getFitParameters(),ndump,altConfig);
            out.close();
        }
        if(ndump > 0 && nSpline > 0) {
            // Dump the P(k) model corresponding to our best fit xi(r).
            std::string outName = outputPrefix + "pk.dat";
            boost::shared_ptr<baofit::PkCorrelationModel> pkModel =
                boost::dynamic_pointer_cast<baofit::PkCorrelationModel>(model);
            pkModel->dump(outName,0.001,0.35,ndump,fmin->getParameters(),zref);
        }
        {
            // Dump the best-fit residuals for each data bin.
            std::string outName = outputPrefix + "residuals.dat";
            std::ofstream out(outName.c_str());
            analyzer.dumpResiduals(out,fmin,combined);
            out.close();
        }
        if(altConfig.length() > 0) {
            // Dump the best-fit residuals for each data bin using an alternate model.
            std::string outName = outputPrefix + "altresiduals.dat";
            std::ofstream out(outName.c_str());
            analyzer.dumpResiduals(out,fmin,combined,altConfig);
            out.close();
        }
        // Calculate and save a bootstrap estimate of the (unfinalized) combined covariance
        // matrix, if requested.
        if(bootstrapCovTrials > 0) {
            if(verbose) std::cout << "Estimating combined covariance with bootstrap..." << std::endl;
            // Although we will only save icov, we still need a copy of the unfinalized combined data
            // in order to get the indexing right.
            bool verbose(false),finalized(false);
            baofit::AbsCorrelationDataPtr copy = analyzer.getCombined(verbose,finalized);
            copy->setCovarianceMatrix(
                analyzer.estimateCombinedCovariance(bootstrapCovTrials, outputPrefix + "bs_cov_work.dat"));
            // Save the inverse covariance if we have enough statistics for a positive definite estimate.
            if(copy->getCovarianceMatrix()->isPositiveDefinite()) {
                std::string outName = outputPrefix + "bs.icov";
                std::ofstream out(outName.c_str());
                copy->saveInverseCovariance(out);
                out.close();
            }
            else {
                std::cout << "Bootstrap covariance matrix estimate is not positive definite." << std::endl;
            }
        }
        // Generate a Markov-chain for marginalization, if requested.
        if(mcmcSave > 0) {
            std::string outName = outputPrefix + "mcmc.dat";
            analyzer.generateMarkovChain(mcmcSave,mcmcInterval,fmin,outName,ndump);
        }
        // Refit the combined sample, if requested.
        likely::FunctionMinimumPtr fmin2;
        if(0 < refitConfig.size()) {
            if(verbose) {
                std::cout << std::endl << "Re-fitting combined with: " << refitConfig << std::endl;
            }
            fmin2 = analyzer.fitSample(combined,refitConfig);
            if(ndump > 0) {
                // Dump the best-fit model.
                std::string outName = outputPrefix + "refit.dat";
                std::ofstream out(outName.c_str());
                analyzer.dumpModel(out,fmin2->getFitParameters(),ndump);
                out.close();
            }
            std::cout << "Delta ChiSquare = "
                << 2*(fmin2->getMinValue() - fmin->getMinValue()) << std::endl;
        }
        // Generate and fit MC samples, if requested.
        if(toymcSamples > 0) {
            std::string outName = outputPrefix + "toymc.dat";
            std::string toymcSaveName;
            if(toymcSave) toymcSaveName = outputPrefix + "toymcsave.data";
            analyzer.doToyMCSampling(toymcSamples,toymcConfig,toymcSaveName,toymcScale,
                fmin,fmin2,refitConfig,outName,ndump);
        }
        // Perform a bootstrap analysis, if requested.
        if(bootstrapTrials > 0) {
            std::string outName = outputPrefix + "bs.dat";
            analyzer.doBootstrapAnalysis(bootstrapTrials,bootstrapSize,fixCovariance,
                fmin,fmin2,refitConfig,outName,ndump);
        }
        // Perform a jackknife analysis, if requested.
        if(jackknifeDrop > 0) {
            std::string outName = outputPrefix + "jk.dat";
            analyzer.doJackknifeAnalysis(jackknifeDrop,fmin,fmin2,refitConfig,outName,ndump);
        }
        // Fit each observation separately, if requested.
        if(fitEach) {
            std::string outName = outputPrefix + "each.dat";
            analyzer.fitEach(fmin,fmin2,refitConfig,outName,ndump);
        }
    }
    catch(std::runtime_error const &e) {
        std::cerr << "ERROR during analysis:\n  " << e.what() << std::endl;
        return -4;
    }
    // All done: normal exit.
    return 0;
}
