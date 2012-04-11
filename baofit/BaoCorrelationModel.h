// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BAO_CORRELATION_MODEL
#define BAOFIT_BAO_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "cosmo/types.h"

#include <string>
#include <vector>

namespace baofit {
	class BaoCorrelationModel : public AbsCorrelationModel {
	// Represents a two-point correlation model parameterized for measuring the scale
	// and significance of a BAO peak.
	public:
		BaoCorrelationModel(std::string const &fiducialName, std::string const &nowigglesName,
            std::string const &broadbandName, double zref);
		virtual ~BaoCorrelationModel();
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in h/Mpc.
        virtual double evaluate(double r, double mu, double z, std::vector<double> const &params) const;
        // Returns the azimuthally averaged monopole correlation function evaluated at
        // comoving separation r in h/Mpc, with an average redshift z.
        virtual double evaluate(double r, double z, std::vector<double> const &params) const;
        // Returns a vector of ell=0,2,4 multipoles for the specified co-moving distance r in Mpc/h
        // and fit parameters. In order to avoid duplicating the code in evaluate(), we call
        // evaluate() with three different values of beta and solve for the multipoles.
        std::vector<double> evaluateMultipoles(double r, std::vector<double> const &params) const;
	private:
        double _zref;
        cosmo::RsdCorrelationFunctionPtr _fid, _nw, _bbc, _bb1, _bb2;
	}; // BaoCorrelationModel
} // baofit

#endif // BAOFIT_BAO_CORRELATION_MODEL
