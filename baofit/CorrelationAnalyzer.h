// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_ANALYZER
#define BAOFIT_CORRELATION_ANALYZER

#include "baofit/types.h"

#include "likely/BinnedDataResampler.h"

namespace baofit {
	class CorrelationAnalyzer {
	public:
		CorrelationAnalyzer(int randomSeed);
		virtual ~CorrelationAnalyzer();
		// Adds a new correlation data object to this analyzer.
        void addData(AbsCorrelationDataCPtr data);
        // Returns the number of data objects added to this analyzer.
        int getNData() const;
        // Sets the correlation model to use.
        void setModel(AbsCorrelationModelCPtr model);
        // Fits the combined correlation data aadded to this analyzer and returns
        // the estimated function minimum.
        likely::FunctionMinimumPtr fitCombined(std::string const &method) const;
	private:
        likely::BinnedDataResampler _resampler;
        AbsCorrelationModelCPtr _model;
	}; // CorrelationAnalyzer
	
    inline int CorrelationAnalyzer::getNData() const { return _resampler.getNObservations(); }
    inline void CorrelationAnalyzer::setModel(AbsCorrelationModelCPtr model) { _model = model; }

} // baofit

#endif // BAOFIT_CORRELATION_ANALYZER
