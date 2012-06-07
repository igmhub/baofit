// Created 06-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_XI_CORRELATION_MODEL
#define BAOFIT_XI_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "likely/types.h"

namespace baofit {
	class XiCorrelationModel : public AbsCorrelationModel {
	public:
		XiCorrelationModel(likely::AbsBinningCPtr rbins);
		virtual ~XiCorrelationModel();
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
        likely::AbsBinningCPtr _rbins;
	}; // XiCorrelationModel
} // baofit

#endif // BAOFIT_XI_CORRELATION_MODEL
