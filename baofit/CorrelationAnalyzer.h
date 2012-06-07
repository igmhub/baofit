// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_ANALYZER
#define BAOFIT_CORRELATION_ANALYZER

#include "baofit/types.h"

#include "cosmo/types.h"

#include "likely/BinnedDataResampler.h"

#include <iosfwd>

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
        // Returns a shared pointer to the combined correlation data added to this
        // analyzer, after it has been finalized.
        AbsCorrelationDataPtr getCombined() const;
        // Fits the combined correlation data aadded to this analyzer and returns
        // the estimated function minimum.
        likely::FunctionMinimumPtr fitCombined(std::string const &method) const;
        // Performs a bootstrap analysis and returns the number of fits to bootstrap
        // samples that failed.
        int doBootstrapAnalysis(std::string const &method,likely::FunctionMinimumPtr fmin,
            int bootstrapTrials, int bootstrapSize = 0, bool fixCovariance = true) const;
        // Dumps the data, prediction, and diagonal error for each bin of the combined
        // data set to the specified output stream. The fit result is assumed to correspond
        // to model that is currently associated with this analyzer. Use the optional script
        // to modify the parameters used in the model.
        void dumpResiduals(std::ostream &out, likely::FunctionMinimumPtr fmin,
            std::string const &script = "") const;
        // Dumps the model predictions for the specified fit result to the specified
        // output stream. The fit result is assumed to correspond to model that is
        // currently associated with this analyzer. Use the optional script to modify
        // the parameters used in the model.
        void dumpModel(std::ostream &out, likely::FunctionMinimumPtr fmin,
            cosmo::Multipole multipole, int nr, double rmin, double rmax, double zval,
            std::string const &script = "") const;
        
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
