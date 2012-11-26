// Created 26-Nov-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BROADBAND_MODEL
#define BAOFIT_BROADBAND_MODEL

#include "baofit/AbsCorrelationModel.h"

namespace baofit {
    // Represents a smooth two-point correlation model parameterized in powers of comoving separation
    // and redshift, and multipoles of mu = r.z
	class BroadbandModel : public AbsCorrelationModel {
	public:
		BroadbandModel(std::string const &name);
		virtual ~BroadbandModel();
	protected:
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double _evaluate(double r, double mu, double z, bool anyChanged) const;
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z.
        virtual double _evaluate(double r, cosmo::Multipole multipole, double z, bool anyChanged) const;
	private:
	}; // BroadbandModel
} // baofit

#endif // BAOFIT_BROADBAND_MODEL
