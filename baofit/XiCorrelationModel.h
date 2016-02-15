// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_XI_CORRELATION_MODEL
#define BAOFIT_XI_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "likely/types.h"
#include "likely/Integrator.h"

#include <string>
#include <vector>

namespace baofit {
	// Represents a two-point correlation model parameterized as an interpolation in each multipole.
	class XiCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a new interpolating correlation model. The input points should be a comma-separated
        // list of r values (in Mpc/h) where spline points will be created. The interpolation method
        // should either be "linear" or "cspline". Creates independent parameters for each multipole
        // if requested, or else relative normalizations are fixed by Kaiser theory.
		XiCorrelationModel(std::string const &points, std::string const &method, bool independentMultipoles,
            double zref,  bool crossCorrelation = false);
		virtual ~XiCorrelationModel();
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
        // Saves the multipole parameter values and covariance to files <prefix>multipole.data
        // and <prefix>multipole.cov in a format suitable for input to a multipole fit. Also
        // prints out the config parameters necessary for using these files. Method is non-const
        // since it calls updateParameterValues. This should do the right thing even if some
        // multipole params are fixed or some non-multipole params are floating.
        void saveMultipolesAsData(std::string const &prefix, likely::FunctionMinimumCPtr fmin);
	protected:
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double _evaluate(double r, double mu, double z, bool anyChanged, int index) const;
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z.
        virtual double _evaluate(double r, cosmo::Multipole multipole, double z, bool anyChanged, int index) const;
        virtual double _evaluateKSpace(double k, double mu_k, double pk, double z) const;
	private:
        std::string _method;
        bool _independentMultipoles;
        int _indexBase;
        std::vector<double> _rValues;
        mutable std::vector<double> _xiValues;
        mutable likely::InterpolatorPtr _xi0, _xi2, _xi4;
        void _initializeInterpolators() const;
        double _xi2Integrand(double r) const;
        double _xi4Integrand(double r) const;
        likely::Integrator::IntegrandPtr _xi2IntegrandPtr, _xi4IntegrandPtr;
	}; // XiCorrelationModel



} // baofit

#endif // BAOFIT_XI_CORRELATION_MODEL
