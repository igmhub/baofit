// Created 31-Jan-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/baofit.h"
#include "cosmo/cosmo.h"
#include "likely/likely.h"
// the following are not part of the public API, so not included by likely.h
#include "likely/MinuitEngine.h"
#include "likely/EngineRegistry.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"

#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"
#include "boost/foreach.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/pointer_cast.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <set>
#include <algorithm>

// Declare bindings to BLAS,LAPACK routines we need
extern "C" {
    // http://www.netlib.org/lapack/double/dpptrf.f
    void dpptrf_(char const *uplo, int const *n, double *ap, int *info);
    // http://www.netlib.org/lapack/double/dpptri.f
    void dpptri_(char const *uplo, int const *n, double *ap, int *info);
    // http://netlib.org/blas/dspmv.f
    void dspmv_(char const *uplo, int const *n, double const *alpha, double const *ap,
        double const *x, int const *incx, double const *beta, double *y, int const *incy);
    // http://www.netlib.org/blas/dsymm.f
    void dsymm_(char const *side, char const *uplo, int const *m, int const *n,
        double const *alpha, double const *a, int const *lda, double const *b,
        int const *ldb, double const *beta, double *c, int const *ldc);
}

namespace lk = likely;
namespace po = boost::program_options;

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

void getDouble(std::string::const_iterator const &begin, std::string::const_iterator const &end,
    double &value) {
    // Use boost::spirit::parse instead of the easier boost::lexical_cast since this is
    // a bottleneck when reading many files. For details, see:
    // http://tinodidriksen.com/2011/05/28/cpp-convert-string-to-double-speed/
    std::string tokenString(begin,end);
    char const *tokenPtr = tokenString.c_str();
    boost::spirit::qi::parse(tokenPtr, &tokenPtr[tokenString.size()],
        boost::spirit::qi::double_, value);    
}

void getInt(std::string::const_iterator const &begin, std::string::const_iterator const &end,
    int &value) {
    std::string tokenString(begin,end);
    value = std::atoi(tokenString.c_str());        
}

// Loads a binned correlation function in cosmolib format and returns a BinnedData object.
// The fast option disables regexp checks for valid numeric inputs.
baofit::QuasarCorrelationDataPtr loadCosmolib(std::string dataName,
    likely::AbsBinningCPtr llBins, likely::AbsBinningCPtr sepBins, likely::AbsBinningCPtr zBins,
    double rmin, double rmax, double llmin, cosmo::AbsHomogeneousUniversePtr cosmology,
    bool verbose, bool icov = false, bool fast = false) {
    // Create the new BinnedData.
    baofit::QuasarCorrelationDataPtr
        binnedData(new baofit::QuasarCorrelationData(llBins,sepBins,zBins,rmin,rmax,llmin,cosmology));
    // General stuff we will need for reading both files.
    std::string line;
    int lineNumber(0);
    // Capturing regexps for positive integer and signed floating-point constants.
    std::string ipat("(0|(?:[1-9][0-9]*))"),fpat("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)");
    if(fast) {
        // Replace validation patterns with simple non-whitespace groups.
        ipat = "(\\S+)";
        fpat = "(\\S+)";
    }
    boost::match_results<std::string::const_iterator> what;
    // Loop over lines in the parameter file.
    std::string paramsName(dataName + ".params");
    std::ifstream paramsIn(paramsName.c_str());
    if(!paramsIn.good()) throw cosmo::RuntimeError("Unable to open " + paramsName);
    boost::regex paramPattern(
        boost::str(boost::format("\\s*%s\\s+%s\\s*\\| Lya covariance 3D \\(%s,%s,%s\\)\\s*")
        % fpat % fpat % fpat % fpat % fpat));
    std::vector<double> axisValues(3);
    while(paramsIn.good() && !paramsIn.eof()) {
        std::getline(paramsIn,line);
        if(paramsIn.eof()) break;
        if(!paramsIn.good()) {
            throw cosmo::RuntimeError("Unable to read line " +
                boost::lexical_cast<std::string>(lineNumber));
        }
        lineNumber++;
        // Parse this line with a regexp.
        if(!boost::regex_match(line,what,paramPattern)) {
            throw cosmo::RuntimeError("Badly formatted params line " +
                boost::lexical_cast<std::string>(lineNumber) + ": '" + line + "'");
        }
        // Expected tokens are [0] value [1] Cinv*d (ignored) [2] logLambda [3] separation [4] redshift.
        int nTokens(5);
        std::vector<double> token(nTokens);
        for(int tok = 0; tok < nTokens; ++tok) {
            getDouble(what[tok+1].first,what[tok+1].second,token[tok]);
        }
        // Add this bin to our dataset.
        //!!addData(token[0],token[2],token[3],token[4]);
        axisValues[0] = token[2];
        axisValues[1] = token[3];
        axisValues[2] = token[4];
        int index = binnedData->getIndex(axisValues);
        binnedData->setData(index,token[0]);        
    }
    //!!finalizeData();
    paramsIn.close();
    if(verbose) {
        std::cout << "Read " << binnedData->getNBinsWithData() << " of "
            << binnedData->getNBinsTotal() << " data values from " << paramsName << std::endl;
    }
    // Loop over lines in the covariance file.
    std::string covName(dataName + (icov ? ".icov" : ".cov"));
    std::ifstream covIn(covName.c_str());
    if(!covIn.good()) throw cosmo::RuntimeError("Unable to open " + covName);
    boost::regex covPattern(boost::str(boost::format("\\s*%s\\s+%s\\s+%s\\s*")
        % ipat % ipat % fpat));
    lineNumber = 0;
    double value;
    int offset1,offset2;
    while(covIn.good() && !covIn.eof()) {
        std::getline(covIn,line);
        if(covIn.eof()) break;
        if(!covIn.good()) {
            throw cosmo::RuntimeError("Unable to read line " +
                boost::lexical_cast<std::string>(lineNumber));
        }
        lineNumber++;
        // Parse this line with a regexp.
        if(!boost::regex_match(line,what,covPattern)) {
            throw cosmo::RuntimeError("Badly formatted cov line " +
                boost::lexical_cast<std::string>(lineNumber) + ": '" + line + "'");
        }
        getInt(what[1].first,what[1].second,offset1);
        getInt(what[2].first,what[2].second,offset2);
        getDouble(what[3].first,what[3].second,value);
        // Add this covariance to our dataset.
        if(icov) value = -value; // !?! see line #388 of Observed2Point.cpp
        //!!addCovariance(offset1,offset2,value,icov);
        int index1 = *(binnedData->begin()+offset1), index2 = *(binnedData->begin()+offset2);
        if(icov) {
            binnedData->setInverseCovariance(index1,index2,value);
        }
        else {
            binnedData->setCovariance(index1,index2,value);
        }
    }
    //!!finalizeCovariance(icov);
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
    covIn.close();
    if(verbose) {
        int ndata = binnedData->getNBinsWithData();
        int ncov = (ndata*(ndata+1))/2;
        std::cout << "Read " << lineNumber << " of " << ncov
            << " covariance values from " << covName << std::endl;
    }
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
    
    // Load the data we will fit.
    /*
    baofit::QuasarCorrelationDataPtr binnedData;
    likely::BinnedDataResampler resampler(randomSeed);
    */
    baofit::CorrelationAnalyzer analyzer(randomSeed,verbose);
    analyzer.setModel(model);
    
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
            analyzer.addData(
                loadCosmolib(dataName,llBins,sepBins,zBins,rmin,rmax,llmin,cosmology,verbose,false,fastLoad));
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
            //!!binnedData.reset(new baofit::QuasarCorrelationData(llBins,sepBins,zBins,cosmology));
            while(platelist.good() && !platelist.eof()) {
                platelist >> plateName;
                if(platelist.eof()) break;
                if(!platelist.good()) {
                    std::cerr << "Error while reading platelist from " << platelistName << std::endl;
                    return -1;
                }
                std::string filename(boost::str(platefile % platerootName % plateName));
                analyzer.addData(
                    loadCosmolib(filename,llBins,sepBins,zBins,rmin,rmax,llmin,cosmology,verbose,true,fastLoad));
                /**
                plateBinned->compress();
                resampler.addObservation(
                    boost::dynamic_pointer_cast<const likely::BinnedData>(plateBinned));
                **/
                if(analyzer.getNData() == maxPlates) break;
            }
            platelist.close();
            //!!binnedData = boost::dynamic_pointer_cast<baofit::QuasarCorrelationData>(resampler.combined());
        }
        
        //!!binnedData->finalize(rmin,rmax,llmin);

    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << "ERROR while reading data:\n  " << e.what() << std::endl;
        return -2;
    }
    
    // Minimize the -log(Likelihood) function.
    try {
        lk::FunctionMinimumPtr fmin = analyzer.fitCombined("mn2::vmetric");

        /**
        // Fit the combined dataset.
        baofit::CorrelationFitter fitEngine(binnedData,model);
        lk::FunctionMinimumPtr fitResult = fitEngine.fit("mn2::vmetric");
        fitResult->printToStream(std::cout);
        **/

        if(bootstrapTrials > 0) {
            analyzer.doBootstrapAnalysis("mn2::vmetric",fmin,bootstrapTrials,bootstrapSize,fixCovariance);
        }

        /**
        // Perform bootstrap trials, if requested
        int nPlates(resampler.getNObservations()), nInvalid(0);
        if(0 < bootstrapTrials && 0 < nPlates) {
            if(0 == bootstrapSize) bootstrapSize = nPlates;
            baofit::QuasarCorrelationDataPtr bsData;
            // Lookup the values of floating parameters from the combined fit.
            std::vector<double> baseline = fitResult->getParameters(true);
            // Initialize statistics accumulators.
            int nstats = baseline.size()+1;
            boost::scoped_array<lk::WeightedAccumulator> stats(new lk::WeightedAccumulator[nstats]);
            likely::CovarianceAccumulator accumulator(nstats);

            for(int trial = 0; trial < bootstrapTrials; ++trial) {
                bsData = boost::dynamic_pointer_cast<baofit::QuasarCorrelationData>(
                    resampler.bootstrap(bootstrapSize,fixCovariance));

            //int trial(0);
            //while(bsData = boost::dynamic_pointer_cast<baofit::QuasarCorrelationData>(
            //resampler.jackknife(1,trial++))) {

                bsData->finalize(rmin,rmax,llmin);
                baofit::CorrelationFitter bsFitEngine(bsData,model);
                lk::FunctionMinimumPtr bsMin = bsFitEngine.fit("mn2::vmetric");
                if(bsMin->getStatus() == likely::FunctionMinimum::OK) {
                    // Lookup the fitted values of floating parameters.
                    std::vector<double> pvalues = bsMin->getParameters(true);
                    for(int par = 0; par < pvalues.size(); ++par) {
                        // Accumulate statistics for this parameter.
                        stats[par].accumulate(pvalues[par]);
                        // Calculate differences from the baseline fit result (to minimize
                        // roundoff error when accumulating covariance statistics).
                        pvalues[par] -= baseline[par];
                    }
                    // Include the fit chiSquare (relative to the baseline value) in
                    // our statistics.
                    stats[nstats-1].accumulate(bsMin->getMinValue());
                    pvalues.push_back(bsMin->getMinValue() - fitResult->getMinValue());
                    accumulator.accumulate(pvalues);
                }
                else {
                    nInvalid++;
                }
                if(trial == bootstrapTrials-1 || (verbose && (0 == (trial+1)%10))) {
                    std::cout << "Completed " << (trial+1) << " bootstrap trials (" << nInvalid
                        << " invalid)" << std::endl;
                }
            }
            std::vector<std::string> labels(fitResult->getNames(true));
            labels.push_back("ChiSquare");
            boost::format resultFormat("%20s = %12.6f +/- %12.6f\n");
            std::cout << std::endl << "Bootstrap Results:" << std::endl;
            for(int stat = 0; stat < nstats; ++stat) {
                std::cout << resultFormat % labels[stat] % stats[stat].mean() % stats[stat].error();
            }
            std::cout << std::endl << "Bootstrap Errors & Correlations:" << std::endl;
            accumulator.getCovariance()->printToStream(std::cout,true,"%12.6f",labels);
        }
        **/

/**        
        lk::GradientCalculatorPtr gcptr;
        LyaBaoLikelihood nll(binnedData,model,rmin,rmax,fixLinear,fixBao,fixScale,noBBand,
            initialAmp,initialScale);
        lk::FunctionPtr fptr(new lk::Function(boost::ref(nll)));

        int npar(nll.getNPar());
        lk::AbsEnginePtr engine = lk::getEngine("mn2::vmetric",fptr,gcptr,model->getParameters());
        lk::MinuitEngine &minuit = dynamic_cast<lk::MinuitEngine&>(*engine);        
        lk::MinuitEngine::StatePtr initialState(new ROOT::Minuit2::MnUserParameterState());
        nll.initialize(initialState);
        std::cout << *initialState;
        
        ROOT::Minuit2::MnStrategy strategy(1); // lo(0),med(1),hi(2)
        ROOT::Minuit2::MnMigrad fitter((ROOT::Minuit2::FCNBase const&)minuit,*initialState,strategy); 

        int maxfcn = 100*npar*npar;
        double edmtol = 0.1;
        ROOT::Minuit2::FunctionMinimum fmin = fitter(maxfcn,edmtol);
        
        if(minos) {
            ROOT::Minuit2::MnMinos minosError((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            for(int ipar = 0; ipar < npar; ++ipar) {
                std::pair<double,double> error = minosError(ipar,maxfcn);
                std::cout << "MINOS error[" << ipar << "] = +" << error.second
                    << ' ' << error.first << std::endl;
            }
        }
        
        std::cout << fmin;
        std::cout << fmin.UserCovariance();
        std::cout << fmin.UserState().GlobalCC();
        
        // Remember the best-fit parameters and errors.
        std::vector<double> bestParams = fmin.UserParameters().Params(),
            bestErrors = fmin.UserParameters().Errors();
        double bestFval = fmin.Fval();
        
        std::vector<ContourPoints> contourData;
        if(ncontour > 0) {
            if(verbose) {
                std::cout << "Calculating contours with " << ncontour << " points..." << std::endl;
            }
            // 95% CL (see http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node33.html)
            // Calculate in mathematica using:
            // Solve[CDF[ChiSquareDistribution[2], x] == 0.95, x]
            nll.setErrorScale(5.99146);
            fmin = fitter(maxfcn,edmtol);
            ROOT::Minuit2::MnContours contours95((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            // Parameter indices: 1=bias, 2=beta, 3=BAO amp, 4=BAO scale, 5=bband a1/10, 6=bband a2/1000
            contourData.push_back(contours95(5,6,ncontour));
            contourData.push_back(contours95(4,6,ncontour));
            contourData.push_back(contours95(1,6,ncontour));
            contourData.push_back(contours95(5,3,ncontour));
            contourData.push_back(contours95(4,3,ncontour));
            contourData.push_back(contours95(1,3,ncontour));            
            contourData.push_back(contours95(5,2,ncontour));
            contourData.push_back(contours95(4,2,ncontour));
            contourData.push_back(contours95(1,2,ncontour));
            // 68% CL
            nll.setErrorScale(2.29575);
            fmin = fitter(maxfcn,edmtol);
            ROOT::Minuit2::MnContours contours68((ROOT::Minuit2::FCNBase const&)minuit,fmin,strategy);
            contourData.push_back(contours68(5,6,ncontour));
            contourData.push_back(contours68(4,6,ncontour));
            contourData.push_back(contours68(1,6,ncontour));
            contourData.push_back(contours68(5,3,ncontour));
            contourData.push_back(contours68(4,3,ncontour));
            contourData.push_back(contours68(1,3,ncontour));            
            contourData.push_back(contours68(5,2,ncontour));
            contourData.push_back(contours68(4,2,ncontour));
            contourData.push_back(contours68(1,2,ncontour));
            // reset
            nll.setErrorScale(1);
        }
        
        // Simulate the null hypothesis by applying theory offsets to each plate, if requested.
        if(nullHypothesis) {
            std::vector<double> nullParams(bestParams);
            nullParams[3] = 0; // BAO peak amplitude
            BOOST_FOREACH(LyaDataPtr plate, plateData) {
                plate->applyTheoryOffsets(model,bestParams,nullParams);
            }
        }
        
        int nplates(plateData.size()), nInvalid(0);
        if(0 < bootstrapTrials && 0 < nplates) {
            lk::Random &random(lk::Random::instance());
            random.setSeed(randomSeed);
            boost::scoped_array<lk::WeightedAccumulator>
                accumulators(new lk::WeightedAccumulator[npar+1]);
            if(0 == bootstrapSize) bootstrapSize = nplates;
            std::ofstream out(bootstrapSaveName.c_str());
            out << "trial nuniq alpha bias beta amp scale xio a0 a1 a2 chisq" << std::endl;
            boost::scoped_ptr<std::ofstream> curvesOut;
            if(0 < bootstrapCurvesName.length()) {
                curvesOut.reset(new std::ofstream(bootstrapCurvesName.c_str()));
            }
            for(int k = 0; k < bootstrapTrials; ++k) {
                // First, decide how many copies of each plate to use in this trial.
                std::vector<double> counter(nplates,0);
                for(int p = 0; p < bootstrapSize; ++p) {
                    // Pick a random plate to use in this trial.
                    int index = (int)std::floor(random.getUniform()*nplates);
                    // Keep track of how many times we use this plate.
                    counter[index]++;
                }
                // Next, build the dataset for this trial.
                data->reset();
                for(int index = 0; index < nplates; ++index) {
                    int repeat = counter[index];
                    if(0 < repeat) data->add(*plateData[index],repeat);
                }
                data->finalize(!naiveCovariance);
                // Count total number of different plates used.
                int nuniq = nplates - std::count(counter.begin(),counter.end(),0);
                // Reset parameters to their initial values.
                initialState.reset(new ROOT::Minuit2::MnUserParameterState());
                nll.initialize(initialState);
                // Do the fit.
                ROOT::Minuit2::MnMigrad
                    bsFitter((ROOT::Minuit2::FCNBase const&)minuit,*initialState,strategy);
                fmin = bsFitter(maxfcn,edmtol);
                if(fmin.IsValid()) {
                    // Save the fit results and accumulate bootstrap stats for each parameter.
                    out << k << ' ' << nuniq << ' ';
                    std::vector<double> params = fmin.UserParameters().Params();
                    for(int i = 0; i < npar; ++i) {
                        double value = params[i];
                        accumulators[i].accumulate(value);
                        out << value << ' ';
                    }
                    out << fmin.Fval() << std::endl;
                    accumulators[npar].accumulate(fmin.Fval());
                    // Output curves of the best-fit multipoles if requested.
                    if(curvesOut) {
                        // A generic correlation model does not know how to calculate multipoles, so
                        // we must dynamically cast to our more specialized BAO model type.
                        boost::shared_ptr<const baofit::BaoCorrelationModel>
                            baoModel(boost::dynamic_pointer_cast<const baofit::BaoCorrelationModel>(model));
                        boost::format fmt(" %.3e %.3e %.3e");
                        double dr(1); // Mpc/h
                        int nr = 1+std::floor((rmax-rmin)/dr);
                        for(int i = 0; i < nr; ++i) {
                            double r = rmin + i*dr;
                            std::vector<double> xi = baoModel->evaluateMultipoles(r,params);
                            *curvesOut << fmt % xi[0] % xi[1] % xi[2];
                        }
                        *curvesOut << std::endl;
                    }
                }
                else {
                    nInvalid++;
                }
                if(verbose && (0 == (k+1)%10)) {
                    std::cout << "Completed " << (k+1) << " bootstrap trials (" << nInvalid
                        << " invalid)" << std::endl;
                }
            }
            if(curvesOut) curvesOut->close();
            out.close();
            for(int i = 0; i < npar; ++i) {
                std::cout << i << ' ' << accumulators[i].mean() << " +/- "
                    << accumulators[i].error() << "\t\t[ " << bestParams[i] << " +/- "
                    << bestErrors[i] << " ]" << std::endl;
            }
            std::cout << "minChiSq = " << accumulators[npar].mean() << " +/- "
                << accumulators[npar].error() << "\t\t[ " << bestFval << " ]" << std::endl;
        }
        
        if(dumpName.length() > 0) {
            if(verbose) std::cout << "Dumping fit results to " << dumpName << std::endl;
            //nll.dump(dumpName,fmin.UserParameters().Params(),contourData,modelBins);
        }
**/
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