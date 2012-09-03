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

    double OmegaMatter,hubbleConstant,zref,minll,maxll,dll,dll2,minsep,dsep,minz,dz,rmin,rmax,llmin,
        rVetoWidth,rVetoCenter,xiRmin,xiRmax,muMin,muMax,kloSpline,khiSpline;
    int nsep,nz,maxPlates,bootstrapTrials,bootstrapSize,randomSeed,ndump,jackknifeDrop,lmin,lmax,
      mcmcSave,mcmcInterval,mcSamples,xiNr,reuseCov,nSpline,splineOrder;
    std::string modelrootName,fiducialName,nowigglesName,broadbandName,dataName,xiPoints,mcConfig,
        platelistName,platerootName,iniName,refitConfig,minMethod,xiMethod,outputPrefix,altConfig;
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
        ("broadband", po::value<std::string>(&broadbandName)->default_value(""),
            "Broadband models will be read from <name>bb<x>.<ell>.dat with x=c,1,2 and ell=0,2,4.")
        ("modelroot", po::value<std::string>(&modelrootName)->default_value(""),
            "Common path to prepend to all model filenames.")
        ("zref", po::value<double>(&zref)->default_value(2.25),
            "Reference redshift used by model correlation functions.")
        ("n-spline", po::value<int>(&nSpline)->default_value(0),
            "Number of spline knots to use spanning (klo,khi).")
        ("klo-spline", po::value<double>(&kloSpline)->default_value(0.02,"0.02"),
            "Minimum k in h/Mpc for k P(k) B-spline.")
        ("khi-spline", po::value<double>(&khiSpline)->default_value(0.2,"0.2"),
            "Maximum k in h/Mpc for k P(k) B-spline.")
        ("order-spline", po::value<int>(&splineOrder)->default_value(3),
            "Order of B-spline in k P(k).")
        ("xi-points", po::value<std::string>(&xiPoints)->default_value(""),
            "Comma-separated list of r values (Mpc/h) to use for interpolating r^2 xi(r)")
        ("xi-method", po::value<std::string>(&xiMethod)->default_value("cspline"),
            "Interpolation method to use in r^2 xi(r), use linear or cspline.")
        ("model-config", po::value<std::vector<std::string> >(&modelConfig)->composing(),
            "Model parameters configuration script (option can appear multiple times).")
        ("alt-config", po::value<std::string>(&altConfig)->default_value(""),
            "Parameter adjustments for dumping alternate best-fit model.")
        ("anisotropic", "Uses anisotropic a,b parameters instead of isotropic scale.")
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
        ("save-icov", "Saves the inverse covariance of the combined data after final cuts.")
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
        ("demo-format", "Cosmolib data in demo format.")
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
        ("llmin", po::value<double>(&llmin)->default_value(0),
            "Minimum value of log(lam2/lam1) to use in fit (multipole data only).")
        ("lmin", po::value<int>(&lmin)->default_value(0),
            "Final cut on minimum multipole ell (0,2,4) to use in fit (multipole data only).")
        ("lmax", po::value<int>(&lmax)->default_value(2),
            "Final cut on maximum multipole ell (0,2,4) to use in fit (multipole data only).")
        ("output-prefix", po::value<std::string>(&outputPrefix)->default_value(""),
            "Prefix to use for all analysis output files.")
        ("ndump", po::value<int>(&ndump)->default_value(100),
            "Number of points spanning [rmin,rmax] to use for dumping models (zero for no dumps).")
        ("decorrelated", "Combined data is saved with decorrelated errors.")
        ("refit-config", po::value<std::string>(&refitConfig)->default_value(""),
            "Script to modify parameters for refits.")
        ("fit-each", "Fits each observation separately.")
        ("bootstrap-trials", po::value<int>(&bootstrapTrials)->default_value(0),
            "Number of bootstrap trials to run if a platelist was provided.")
        ("bootstrap-size", po::value<int>(&bootstrapSize)->default_value(0),
            "Size of each bootstrap trial or zero to use the number of plates.")
        ("jackknife-drop", po::value<int>(&jackknifeDrop)->default_value(0),
            "Number of observations to drop from each jackknife sample (zero for no jackknife analysis)")
        ("mcmc-save", po::value<int>(&mcmcSave)->default_value(0),
            "Number of Markov chain Monte Carlo samples to save (zero for no MCMC analysis)")
        ("mcmc-interval", po::value<int>(&mcmcInterval)->default_value(10),
            "Interval for saving MCMC trials (larger for less correlations and longer running time)")
        ("mc-samples", po::value<int>(&mcSamples)->default_value(0),
            "Number of MC samples to generate and fit.")
        ("mc-config", po::value<std::string>(&mcConfig)->default_value(""),
            "Fit parameter configuration to apply before generating samples.")
        ("mc-save", "Saves first generated MC sample.")
        ("random-seed", po::value<int>(&randomSeed)->default_value(1966),
            "Random seed to use for generating bootstrap samples.")
        ("min-method", po::value<std::string>(&minMethod)->default_value("mn2::vmetric"),
            "Minimization method to use for fitting.")
        ;

    allOptions.add(genericOptions).add(modelOptions).add(dataOptions)
        .add(frenchOptions).add(cosmolibOptions).add(analysisOptions);
    po::variables_map vm;

    // Parse command line options first so they override anything in an INI file (except for --model-config)
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
        fitEach(vm.count("fit-each")), xiHexa(vm.count("xi-hexa")), demoFormat(vm.count("demo-format")),
        xiFormat(vm.count("xi-format")), decorrelated(vm.count("decorrelated")), mcSave(vm.count("mc-save")),
        expanded(vm.count("expanded")), sectors(vm.count("sectors")), saveICov(vm.count("save-icov"));

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
    baofit::CorrelationAnalyzer analyzer(minMethod,rmin,rmax,verbose);

    // Initialize the fit model we will use.
    cosmo::AbsHomogeneousUniversePtr cosmology;
    baofit::AbsCorrelationModelPtr model;
    try {
        // Build the homogeneous cosmology we will use.
        cosmology.reset(new cosmo::LambdaCdmRadiationUniverse(OmegaMatter,0,hubbleConstant));
        
        if(nSpline > 0) {
            model.reset(new baofit::PkCorrelationModel(modelrootName,nowigglesName,
                kloSpline,khiSpline,nSpline,splineOrder,zref));
        }
        else if(xiPoints.length() > 0) {
            model.reset(new baofit::XiCorrelationModel(xiPoints,zref,xiMethod));
        }
        else {
            // Build our fit model from tabulated ell=0,2,4 correlation functions on disk.
            model.reset(new baofit::BaoCorrelationModel(
                modelrootName,fiducialName,nowigglesName,broadbandName,zref,anisotropic));
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
    try {
        
        // Create a prototype of the binned data we will be loading.
        baofit::AbsCorrelationDataCPtr prototype;
        if(french) {
            zdata = 2.30;
            prototype = baofit::boss::createFrenchPrototype(zdata,rmin,rmax,rVetoMin,rVetoMax,ellmin,ellmax);
        }
        else if(sectors) {
            zdata = 2.30;
            prototype = baofit::boss::createSectorsPrototype(zdata,rmin,rmax,muMin,muMax,rVetoMin,rVetoMax);
        }
        else if(dr9lrg) {
            zdata = 0.57;
            prototype = baofit::boss::createDR9LRGPrototype(zdata,rmin,rmax,rVetoMin,rVetoMax,
                "LRG/Sample4_North.cov",verbose);
        }
        else if(xiFormat) {
            zdata = 2.25;
            prototype = baofit::boss::createCosmolibXiPrototype(minz,dz,nz,xiRmin,xiRmax,xiNr,xiHexa,
                rmin,rmax,rVetoMin,rVetoMax,ellmin,ellmax);
        }
        else { // default is cosmolib (demo) format
            zdata = 2.25;
            prototype = baofit::boss::createCosmolibPrototype(
                minsep,dsep,nsep,minz,dz,nz,minll,maxll,dll,dll2,rmin,rmax,muMin,muMax,
                rVetoMin,rVetoMax,llmin,cosmology);
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
                    verbose,unweighted,expanded,checkPosDef));
            }
            else if(sectors) {
                analyzer.addData(baofit::boss::loadSectors(*filename,prototype,verbose));
            }
            else if(dr9lrg) {
                analyzer.addData(baofit::boss::loadDR9LRG(*filename,prototype,verbose));
            }
            else if(xiFormat) {
                analyzer.addData(baofit::boss::loadCosmolibXi(*filename,prototype,
                    verbose,weighted,reuseCov,checkPosDef));
            }
            else {
                // Add a cosmolib dataset, assumed to provided icov instead of cov.
                if(demoFormat) {
                    analyzer.addData(baofit::boss::loadCosmolibDemo(*filename,prototype,verbose));                    
                }
                else {
                    analyzer.addData(baofit::boss::loadCosmolib(*filename,prototype,
                        verbose,true,weighted,reuseCov,checkPosDef));
                }
            }
        }
    }
    catch(baofit::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }
    analyzer.setZData(zdata);

    // Do the requested analyses...
    try {
        // Always fit the combined sample.
        likely::FunctionMinimumPtr fmin = analyzer.fitCombined();
        // Print out some extra info if this fit has floating "BAO alpha-*" and "gamma-alpha" parameters.
        std::cout << std::endl;
        analyzer.printScaleZEff(fmin,zref,"BAO alpha-iso");
        analyzer.printScaleZEff(fmin,zref,"BAO alpha-parallel");
        analyzer.printScaleZEff(fmin,zref,"BAO alpha-perp");
        // Dump the combined multipole data points with decorrelated errors, if possible.
        if(french || dr9lrg || xiFormat) {
            std::string outName = outputPrefix + "combined.dat";
            std::ofstream out(outName.c_str());
            boost::shared_ptr<baofit::MultipoleCorrelationData> combined =
                boost::dynamic_pointer_cast<baofit::MultipoleCorrelationData>(analyzer.getCombined());
            std::vector<double> dweights;
            if(decorrelated) {
                // Calculate the decorrelated weights for the combined fit.
                analyzer.getDecorrelatedWeights(analyzer.getCombined(),fmin->getParameters(),dweights);
            }
            combined->dump(out,dweights);
            out.close();
        }
        // Save the combined inverse covariance, if requested.
        if(saveICov) {
            baofit::AbsCorrelationDataPtr combined = analyzer.getCombined();
            std::string outName = outputPrefix + "icov.dat";
            std::ofstream out(outName.c_str());
            for(likely::BinnedData::IndexIterator iter1 = combined->begin(); iter1 != combined->end(); ++iter1) {
                int index1(*iter1);
                // Save all diagonal elements.
                out << index1 << ' ' << index1 << ' ' << combined->getInverseCovariance(index1,index1) << std::endl;
                // Loop over pairs with index2 > index1
                for(likely::BinnedData::IndexIterator iter2 = iter1; ++iter2 != combined->end();) {
                    int index2(*iter2);
                    // Only save non-zero off-diagonal elements.
                    double Cinv(combined->getInverseCovariance(index1,index2));
                    if(Cinv != 0) out << index1 << ' ' << index2 << ' ' << Cinv << std::endl;
                }
            }
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
        {
            // Dump the best-fit residuals for each data bin.
            std::string outName = outputPrefix + "residuals.dat";
            std::ofstream out(outName.c_str());
            analyzer.dumpResiduals(out,fmin);
            out.close();
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
            fmin2 = analyzer.fitCombined(refitConfig);
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
        if(mcSamples > 0) {
            std::string outName = outputPrefix + "mc.dat";
            std::string mcSaveName;
            if(mcSave) mcSaveName = outputPrefix + "mcsave.dat";
            analyzer.doMCSampling(mcSamples,mcConfig,mcSaveName,fmin,fmin2,refitConfig,outName,ndump);
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
        std::cerr << "ERROR during fit:\n  " << e.what() << std::endl;
        return -2;
    }
    // All done: normal exit.
    return 0;
}
