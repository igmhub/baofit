// Created 26-Nov-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BROADBAND_MODEL
#define BAOFIT_BROADBAND_MODEL

#include "baofit/AbsCorrelationModel.h"

namespace baofit {
    // Represents a smooth two-point correlation model parameterized in powers of comoving separation
    // and redshift, and multipoles of mu = r.z
	class BroadbandModel : public AbsCorrelationModel {
	public:
		BroadbandModel(std::string const &name, std::string const &paramSpec);
		virtual ~BroadbandModel();
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
	protected:
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double _evaluate(double r, double mu, double z, bool anyChanged) const;
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z.
        virtual double _evaluate(double r, cosmo::Multipole multipole, double z, bool anyChanged) const;
	private:
        int _rIndexMin,_rIndexMax,_rIndexStep;
        int _muIndexMin,_muIndexMax,_muIndexStep;
        int _zIndexMin,_zIndexMax,_zIndexStep;
	}; // BroadbandModel
} // baofit

#endif // BAOFIT_BROADBAND_MODEL
