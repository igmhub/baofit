// Created 31-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_PK_CORRELATION_MODEL
#define BAOFIT_PK_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "cosmo/types.h"

#include <vector>

namespace baofit {
	class PkCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a new correlation model parameterized as the sum of the specified tabulated smooth model
	    // with a uniform B-spline in k*P(ell,k) added to each k-space multipole, with nk uniformly spaced
	    // knots spanning the range (klo,khi) in h/Mpc.
		PkCorrelationModel(std::string const &modelrootName, std::string const &nowigglesName,
		    double klo, double khi, int nk, int splineOrder, double zref);
		virtual ~PkCorrelationModel();
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
        void _calculateNorm(double z) const;
        double _xi(double r, cosmo::Multipole multipole) const;
        double _getE(int j, double r, cosmo::Multipole multipole) const;
        void _fillSinIntegralCache(double r) const;
	    mutable std::vector<double> _sinInt;
        mutable double _rsave;
        int _nk, _splineOrder;
        double _klo, _dk, _dk2, _dk3, _zref, _twopisq;
        cosmo::PowerSpectrumPtr _nwPower;
	    cosmo::CorrelationFunctionPtr _nw0,_nw2,_nw4;
        mutable double _norm0, _norm2, _norm4;
	}; // PkCorrelationModel
} // baofit

#endif // BAOFIT_PK_CORRELATION_MODEL
