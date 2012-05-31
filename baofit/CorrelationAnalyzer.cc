// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/CorrelationAnalyzer.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/CorrelationFitter.h"

namespace local = baofit;

local::CorrelationAnalyzer::CorrelationAnalyzer(int randomSeed)
: _resampler(randomSeed)
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
    return fitter.fit(method);
}
