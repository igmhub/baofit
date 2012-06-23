// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_FITTER
#define BAOFIT_CORRELATION_FITTER

#include "baofit/AbsCorrelationData.h"
#include "baofit/types.h"
#include "likely/types.h"

#include <vector>

namespace baofit {
	class CorrelationFitter {
	// Manages a correlation function fit.
	public:
	    // Creates a new fitter for the specified data and model.
		CorrelationFitter(AbsCorrelationDataCPtr data, AbsCorrelationModelPtr model);
		virtual ~CorrelationFitter();
		// Changes the error scale definition. The default value of 1 corresponds to the
		// usual 1-sigma errors.
        void setErrorScale(double scale);
        // Returns chiSquare/2 for the specified model parameter values.
        double operator()(likely::Parameters const &params) const;
        // Performs the fit and returns an estimate of the function minimum. Use the optional
        // config parameter to provide a script that will modify the initial parameter values
        // and errors (including fixed/floating) for this fit only.
        likely::FunctionMinimumPtr fit(std::string const &methodName, std::string const &config = "") const;
        // Generates nchain*interval Markov chain MC samples and fills the vector provided with the parameters
        // every interval samples. Uses the input fmin to initialize the MCMC proposal function and determine
        // which parameters are floating. On return, fmin is updated based on accumulated statistics.
        void mcmc(likely::FunctionMinimumPtr fmin, int nchain, int interval, std::vector<double> &samples) const;
	private:
        AbsCorrelationData::TransverseBinningType _type;
        AbsCorrelationDataCPtr _data;
        AbsCorrelationModelPtr _model;
        double _errorScale;
	}; // CorrelationFitter
} // baofit

#endif // BAOFIT_CORRELATION_FITTER
