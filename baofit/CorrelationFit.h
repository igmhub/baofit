// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_FIT
#define BAOFIT_CORRELATION_FIT

#include "baofit/types.h"
#include "likely/types.h"

namespace baofit {
	class CorrelationFit {
	// Manages a correlation function fit.
	public:
	    // Creates a new fit manager for the specified data and model.
		CorrelationFit(AbsCorrelationDataCPtr data, AbsCorrelationModelCPtr model);
		virtual ~CorrelationFit();
		// Changes the error scale definition. The default value of 1 corresponds to the
		// usual 1-sigma errors.
        void setErrorScale(double scale);
        // Returns chiSquare/2 for the specified model parameter values.
        double operator()(likely::Parameters const &params) const;
        // Performs the fit and returns an estimate of the function minimum.
        likely::FunctionMinimumPtr fit(std::string const &methodName) const;
	private:
        AbsCorrelationDataCPtr _data;
        AbsCorrelationModelCPtr _model;
        double _errorScale;
	}; // CorrelationFit
} // baofit

#endif // BAOFIT_CORRELATION_FIT
