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
		// be provided in Mpc/h.
        virtual double evaluate(double r, double mu, double z,
            std::vector<double> const &params) const = 0;
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z.
        virtual double evaluate(double r, cosmo::Multipole multipole, double z,
            std::vector<double> const &params) const = 0;
	}; // AbsCorrelationModel
} // baofit

#endif // BAOFIT_ABS_CORRELATION_MODEL
