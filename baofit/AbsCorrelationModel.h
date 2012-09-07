// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_MODEL
#define BAOFIT_ABS_CORRELATION_MODEL

#include "likely/FitModel.h"

#include "cosmo/types.h"

#include <string>
#include <vector>

namespace baofit {
	class AbsCorrelationModel : public likely::FitModel {
	// Represents an abstract parameterized model of a two-point correlation function.
	public:
	    // Creates a new model with the specified name.
		AbsCorrelationModel(std::string const &name);
		virtual ~AbsCorrelationModel();
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h. Updates our current parameter values.
        double evaluate(double r, double mu, double z, likely::Parameters const &params);
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z. Updates our current parameter values.
        double evaluate(double r, cosmo::Multipole multipole, double z, likely::Parameters const &params);
    protected:
        // The public methods above call these protected methods after making parameter values
        // and changes available via our base class' getParameterValue() and isParameterValueChanged()
        // methods. Any registered changes to parameter values are reset after calling any of these.
        virtual double _evaluate(double r, double mu, double z, bool changed) const = 0;
        virtual double _evaluate(double r, cosmo::Multipole multipole, double z, bool changed) const = 0;
        // Defines the standard set of linear bias parameters used by _getNormFactor below. Returns
        // the index of the last parameter defined.
        int _defineLinearBiasParameters(double zref);
        // Evaluates the redshift evolution p(z) of a parameter for which p(zref)=p0 according to
        // p(z) = ((1+z)/(1+zref))^gamma.
        double _redshiftEvolution(double p0, double gamma, double z) const;
        // Updates the multipole normalization factors b^2(z)*C_ell(beta(z)) returned by getNormFactor(ell).
        double _getNormFactor(cosmo::Multipole multipole, double z) const;
    private:
        int _indexBase;
        double _zref;
        mutable double _normFactor0, _normFactor2, _normFactor4;
	}; // AbsCorrelationModel
} // baofit

#endif // BAOFIT_ABS_CORRELATION_MODEL
