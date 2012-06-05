// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_MODEL
#define BAOFIT_ABS_CORRELATION_MODEL

#include "likely/FitParameter.h"

#include "cosmo/types.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace baofit {
	class AbsCorrelationModel {
	// Represents an abstract parameterized model of a two-point correlation function.
	public:
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
        // Returns a reference to our model's fit parameters.
        likely::FitParameters const &getParameters() const;
        // Returns the number of model parameters (including both floating and fixed).
        int getNParameters() const;
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
    protected:
        void defineParameter(std::string const &name, double value, double error, bool fixed);
	private:
        std::string _name;
        likely::FitParameters _parameters;
	}; // AbsCorrelationModel

} // baofit

#endif // BAOFIT_ABS_CORRELATION_MODEL
