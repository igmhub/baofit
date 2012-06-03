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
	    // Creates a new model using the specified tabulated correlation functions.
		BaoCorrelationModel(std::string const &modelrootName,
		    std::string const &fiducialName, std::string const &nowigglesName,
            std::string const &broadbandName, double zref, double initialAmp, double initialScale,
            bool fixAlpha, bool fixLinear, bool fixBao, bool fixScale, bool noBBand);
		virtual ~BaoCorrelationModel();
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double evaluate(double r, double mu, double z,
            std::vector<double> const &params) const;
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z.
        virtual double evaluate(double r, cosmo::Multipole multipole, double z,
            std::vector<double> const &params) const;
	private:
        double _zref;
        cosmo::RsdCorrelationFunctionPtr _fid, _nw, _bbc, _bb1, _bb2;
	}; // BaoCorrelationModel
} // baofit

#endif // BAOFIT_BAO_CORRELATION_MODEL
