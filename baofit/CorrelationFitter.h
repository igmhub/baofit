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
        // Fills the vector provided with the model prediction for the specified parameter values.
        void getPrediction(likely::Parameters const &params, std::vector<double> &prediction) const;
        // Returns chiSquare/2 for the specified model parameter values.
        double operator()(likely::Parameters const &params) const;
        // Performs the fit and returns an estimate of the function minimum. Use the optional
        // config parameter to provide a script that will modify the initial parameter values
        // and errors (including fixed/floating) for this fit only.
        likely::FunctionMinimumPtr fit(std::string const &methodName, std::string const &config = "") const;
        // Generates nchain*interval Markov chain MC samples and fills the vector provided with the parameters
        // every interval samples. Uses the input fmin, if one is provided, to initialize the MCMC proposal
        // function and determine which parameters are floating. If !fmin, then uses FitModel::guessMinimum
        // instead, which uses our model's initial fit parameter configuration. The MCMC chain is generated
        // without any periodic updates to the proposal function's covariance estimate. Returns the function
        // minimum actually used for MCMC (after any updates).
        likely::FunctionMinimumCPtr mcmc(likely::FunctionMinimumCPtr fmin, int nchain, int interval,
            std::vector<double> &samples) const;
	private:
        AbsCorrelationData::TransverseBinningType _type;
        AbsCorrelationDataCPtr _data;
        AbsCorrelationModelPtr _model;
        double _errorScale;
	}; // CorrelationFitter
} // baofit

#endif // BAOFIT_CORRELATION_FITTER
