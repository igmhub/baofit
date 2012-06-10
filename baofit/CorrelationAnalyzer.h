// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_ANALYZER
#define BAOFIT_CORRELATION_ANALYZER

#include "baofit/types.h"

#include "cosmo/types.h"

#include "likely/BinnedDataResampler.h"

#include <iosfwd>

namespace baofit {
    // Accumulates correlation data and manages its analysis.
	class CorrelationAnalyzer {
	public:
	    // Creates a new analyzer using the specified random seed and minimization method.
		CorrelationAnalyzer(int randomSeed, std::string const &method, bool verbose = true);
		virtual ~CorrelationAnalyzer();
		// Set the verbose level during analysis.
        void setVerbose(bool value);
		// Adds a new correlation data object to this analyzer.
        void addData(AbsCorrelationDataCPtr data);
        // Returns the number of data objects added to this analyzer.
        int getNData() const;
        // Sets the correlation model to use.
        void setModel(AbsCorrelationModelPtr model);
        // Returns a shared pointer to the combined correlation data added to this
        // analyzer, after it has been finalized. If verbose, prints out the number
        // of bins with data before and after finalizing the data.
        AbsCorrelationDataPtr getCombined(bool verbose = false) const;
        // Fits the combined correlation data aadded to this analyzer and returns
        // the estimated function minimum. Use the optional config script to modify
        // the initial parameter configuration used for the fit (any changes do not
        // propagate back to the model or modify subsequent fits).
        likely::FunctionMinimumPtr fitCombined(std::string const &config = "") const;
        // Performs a bootstrap analysis and returns the number of fits to bootstrap
        // samples that failed. Specify a non-zero bootstrapSize to generate trials with
        // a number of observations different than getNData(). Specify a refitConfig script
        // to fit each bootstrap sample twice: first with the default model config, then
        // with the refit config script applied. In this case, a trial is only considered
        // successful if both fits succeed. Setting fixCovariance to false means that
        // fits will use a covariance matrix that does not correctly account for double
        // counting. See likely::BinnedDataResampler::bootstrap for details.
        int doBootstrapAnalysis(likely::FunctionMinimumPtr fmin, int bootstrapTrials,
            int bootstrapSize = 0, std::string const &refitConfig = "", bool fixCovariance = true) const;
        // Dumps the data, prediction, and diagonal error for each bin of the combined
        // data set to the specified output stream. The fit result is assumed to correspond
        // to model that is currently associated with this analyzer. Use the optional script
        // to modify the parameters used in the model.
        void dumpResiduals(std::ostream &out, likely::FunctionMinimumPtr fmin,
            std::string const &script = "") const;
        // Dumps the model predictions for the specified fit result to the specified
        // output stream. The fit result is assumed to correspond to model that is
        // currently associated with this analyzer. Use the optional script to modify
        // the parameters used in the model. By default, values are output as
        // "rval mono quad hexa" on separate lines. With oneLine = true, values of
        // "mono quad hexa" are concatenated onto a single line.
        void dumpModel(std::ostream &out, likely::FunctionMinimumPtr fmin,
            int nr, double rmin, double rmax, double zval, std::string const &script = "",
            bool oneLine = false) const;
        
	private:
        std::string _method;
        bool _verbose;
        likely::BinnedDataResampler _resampler;
        AbsCorrelationModelPtr _model;
	}; // CorrelationAnalyzer
	
    inline void CorrelationAnalyzer::setVerbose(bool value) { _verbose = value; }
    inline int CorrelationAnalyzer::getNData() const { return _resampler.getNObservations(); }
    inline void CorrelationAnalyzer::setModel(AbsCorrelationModelPtr model) { _model = model; }

} // baofit

#endif // BAOFIT_CORRELATION_ANALYZER
