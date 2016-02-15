// Created 14-Jan-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_BAO_NON_LINEAR_CORRECTION_MODEL
#define BAOFIT_BAO_NON_LINEAR_CORRECTION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "likely/types.h"

namespace baofit {
	// Represents a non-linear correction in the flux power spectrum model
	class NonLinearCorrectionModel : public AbsCorrelationModel {
	public:
	    // Creates a non-linear correction model. Bool parameters determine which, if any,
	    // of the two available models to use. The parameter sigma8 is used to rescale one
	    // of the model parameters. The parameter zref is used for calculating the linear
	    // growth factor.
	    NonLinearCorrectionModel(double zref, double sigma8, bool nlCorrection = false,
	        bool nlCorrectionAlt = false);
	    virtual ~NonLinearCorrectionModel();
	    // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
	protected:
	    virtual double _evaluate(double r, double mu, double z, bool anyChanged, int index) const;
	    // Returns the non-linear correction in k space at point (k,mu_k). One of the models
	    // depends also on the linear matter power spectrum P(k). The redshift z is used to
	    // interpolate in the fixed model parameters and to calculate the linear growth factor.
	    virtual double _evaluateKSpace(double k, double mu_k, double pk, double z) const;
	private:
	    double _zref, _sigma8;
	    bool _nlCorrection, _nlCorrectionAlt;
	    void _initialize();
	    mutable likely::InterpolatorPtr _qnlInterpolator, _kvInterpolator, _avInterpolator, _bvInterpolator, _kpInterpolator;
    }; // NonLinearCorrectionModel
} // baofit

#endif // BAOFIT_BAO_NON_LINEAR_CORRECTION_MODEL
