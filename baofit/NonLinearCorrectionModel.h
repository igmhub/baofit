// Created 14-Jan-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_BAO_NON_LINEAR_CORRECTION_MODEL
#define BAOFIT_BAO_NON_LINEAR_CORRECTION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "likely/Interpolator.h"

namespace baofit {
	// Represents a non-linear correction in the flux power spectrum model
	class NonLinearCorrectionModel : public AbsCorrelationModel {
	public:
	    // Creates a non-linear correction model. Bool parameters determine which, if any,
	    // of the two available models to use. Sigma8 parameter is used to rescale one of
	    // the model parameters.
	    NonLinearCorrectionModel(double sigma8, bool nlCorrection = false, bool nlCorrectionAlt = false);
	    virtual ~NonLinearCorrectionModel();
	protected:
	    // Returns the non-linear correction in k space at point (k,mu_k). One of the models
	    // depends also on the linear matter power spectrum P(k). The redshift z is used
	    // to interpolate in the fixed model parameters and to calculate the linear growth factor.
	    double _evaluateNLCorrection(double k, double mu_k, double pk, double z) const;
	private:
	    double _sigma8;
	    bool _nlCorrection, _nlCorrectionAlt;
	    void _initialize();
	    likely::Interpolator _qnlInterpolator, _kvInterpolator, _avInterpolator, _bvInterpolator, _kpInterpolator;
    }; // NonLinearCorrectionModel
} // baofit

#endif // BAOFIT_BAO_NON_LINEAR_CORRECTION_MODEL
