// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_MODEL
#define BAOFIT_ABS_CORRELATION_MODEL

#include "likely/FitParameter.h"

#include <vector>

namespace baofit {
	class AbsCorrelationModel {
	// Represents an abstract parameterized model of a two-point correlation function.
	public:
		AbsCorrelationModel();
		virtual ~AbsCorrelationModel();
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double evaluate(double r, double mu, double z, std::vector<double> const &params) const = 0;
        // Returns the azimuthally averaged monopole correlation function evaluated at
        // comoving separation r in Mpc/h, with an average redshift z.
        virtual double evaluate(double r, double z, std::vector<double> const &params) const = 0;
        // Returns a reference to our model's fit parameters.
        likely::FitParameters const &getParameters() const;
        // Returns the number of model parameters (including both floating and fixed).
        int getNParameters() const;
    protected:
        void defineParameter(std::string const &name, double value, double error, bool fixed);
	private:
        likely::FitParameters _parameters;
	}; // AbsCorrelationModel

} // baofit

#endif // BAOFIT_ABS_CORRELATION_MODEL
