// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_ANALYZER
#define BAOFIT_CORRELATION_ANALYZER

#include "baofit/types.h"

#include "likely/BinnedDataResampler.h"

namespace baofit {
	class CorrelationAnalyzer {
	public:
		CorrelationAnalyzer(int randomSeed, bool verbose = true);
		virtual ~CorrelationAnalyzer();
		// Set the verbose level during analysis.
        void setVerbose(bool value);
		// Adds a new correlation data object to this analyzer.
        void addData(AbsCorrelationDataCPtr data);
        // Returns the number of data objects added to this analyzer.
        int getNData() const;
        // Sets the correlation model to use.
        void setModel(AbsCorrelationModelCPtr model);
        // Fits the combined correlation data aadded to this analyzer and returns
        // the estimated function minimum.
        likely::FunctionMinimumPtr fitCombined(std::string const &method) const;
        // Performs a bootstrap analysis...
        // Returns the number of fits to bootstrap samples that failed.
        int doBootstrapAnalysis(std::string const &method,likely::FunctionMinimumPtr fmin,
            int bootstrapTrials, int bootstrapSize = 0, bool fixCovariance = true) const;
	private:
        bool _verbose;
        likely::BinnedDataResampler _resampler;
        AbsCorrelationModelCPtr _model;
	}; // CorrelationAnalyzer
	
    inline void CorrelationAnalyzer::setVerbose(bool value) { _verbose = value; }
    inline int CorrelationAnalyzer::getNData() const { return _resampler.getNObservations(); }
    inline void CorrelationAnalyzer::setModel(AbsCorrelationModelCPtr model) { _model = model; }

} // baofit

#endif // BAOFIT_CORRELATION_ANALYZER
