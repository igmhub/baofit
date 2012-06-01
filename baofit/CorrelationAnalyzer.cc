// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/CorrelationAnalyzer.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/CorrelationFitter.h"

#include "likely/CovarianceAccumulator.h"
#include "likely/WeightedAccumulator.h"
#include "likely/FunctionMinimum.h"
#include "likely/CovarianceMatrix.h"

#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"

namespace local = baofit;

local::CorrelationAnalyzer::CorrelationAnalyzer(int randomSeed, bool verbose)
: _resampler(randomSeed), _verbose(verbose)
{ }

local::CorrelationAnalyzer::~CorrelationAnalyzer() { }

void local::CorrelationAnalyzer::addData(AbsCorrelationDataCPtr data) {
    _resampler.addObservation(boost::dynamic_pointer_cast<const likely::BinnedData>(data));
}

likely::FunctionMinimumPtr local::CorrelationAnalyzer::fitCombined(std::string const &method) const {
    AbsCorrelationDataPtr combined =
        boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(_resampler.combined());
    combined->finalize();
    CorrelationFitter fitter(combined,_model);
    likely::FunctionMinimumPtr fmin = fitter.fit(method);
    if(_verbose) fmin->printToStream(std::cout);
    return fmin;
}

int local::CorrelationAnalyzer::doBootstrapAnalysis(std::string const &method,
likely::FunctionMinimumPtr fmin, int bootstrapTrials, int bootstrapSize, bool fixCovariance) const {
    if(getNData() <= 1) {
        throw RuntimeError("CorrelationAnalyzer::doBootstrapAnalysis: need > 1 observation.");
    }
    if(0 == bootstrapSize) bootstrapSize = getNData();
    baofit::AbsCorrelationDataPtr bsData;
    // Lookup the values of floating parameters from the combined fit.
    std::vector<double> baseline = fmin->getParameters(true);
    // Initialize statistics accumulators.
    int nstats = baseline.size()+1;
    boost::scoped_array<likely::WeightedAccumulator> stats(new likely::WeightedAccumulator[nstats]);
    likely::CovarianceAccumulator accumulator(nstats);
    int nInvalid(0);

    for(int trial = 0; trial < bootstrapTrials; ++trial) {
        // Generate the next boostrap sample.
        bsData = boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(
            _resampler.bootstrap(bootstrapSize,fixCovariance));
        bsData->finalize();
        // Fit the sample.
        baofit::CorrelationFitter bsFitEngine(bsData,_model);
        likely::FunctionMinimumPtr bsMin = bsFitEngine.fit(method);
        // Accumulate the fit results.
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
            pvalues.push_back(bsMin->getMinValue() - fmin->getMinValue());
            accumulator.accumulate(pvalues);
        }
        else {
            nInvalid++;
        }
        // Print periodic updates while the analysis is running.
        if(trial == bootstrapTrials-1 || (_verbose && (0 == (trial+1)%10))) {
            std::cout << "Completed " << (trial+1) << " bootstrap trials (" << nInvalid
                << " invalid)" << std::endl;
        }
    }
    // Print a summary of the analysis results.
    std::vector<std::string> labels(fmin->getNames(true));
    labels.push_back("ChiSquare");
    boost::format resultFormat("%20s = %12.6f +/- %12.6f\n");
    std::cout << std::endl << "Bootstrap Results:" << std::endl;
    for(int stat = 0; stat < nstats; ++stat) {
        std::cout << resultFormat % labels[stat] % stats[stat].mean() % stats[stat].error();
    }
    std::cout << std::endl << "Bootstrap Errors & Correlations:" << std::endl;
    accumulator.getCovariance()->printToStream(std::cout,true,"%12.6f",labels);
    
    return nInvalid;
}
