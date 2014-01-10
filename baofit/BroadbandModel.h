// Created 26-Nov-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BROADBAND_MODEL
#define BAOFIT_BROADBAND_MODEL

#include "baofit/AbsCorrelationModel.h"

namespace baofit {
    // Represents a smooth two-point correlation model parameterized in powers of comoving separation r,
    // multipoles of mu = r.z, the parallel (rP) and transverse (rT) components of r, and redshift z.
	class BroadbandModel : public AbsCorrelationModel {
	public:
	    // Creates a broadband model with parameters named "<tag> z%d mu%d r%d rP%d rT%d" where each
	    // %d slot is filled with a range of indices specified in paramSpec using python array notation:
	    //   n => only index is n
	    //   n1:n2 => indices n1...n2
	    //   n1:n2:dn => indices n1...n2 in steps of dn
		// For the r-axis only, dn < 0 is allowed and implies fractional steps of 1/(-dn).
	    // The set of axes that will be used is determined by one of the following prefixes:
	    //   r,mu=
	    //   rP,rT=
	    //   r,mu,z=
	    //   rP,rT,z=
	    //   r,mu,rP,rT,z=
	    // For example, the string "r,mu=-2:0,0:4:2" specifies a polynomial in (r/r0)^n P(ell,mu) with
	    // n=-2,-1,0 and ell = 0,2,4. If the prefix is missing and 3 axes are specified, then
	    // the r,mu,z syntax is assumed by default, for backwards compatibility. The parameters r0,z0
	    // are used to normalize the polynomial expansion factors as follows:
	    //
	    //   (r/r0 - q)^n, (rP/r0 - q)^n, (rT/r0 - q)^n, ((1+z)/(1+z0))^n
	    //
	    // where q = 0 if n <= 0 or else q = +1, e.g. ((r-r0)/r0)^2.
		BroadbandModel(std::string const &name, std::string const &tag, std::string const &paramSpec,
		    double r0, double z0, AbsCorrelationModel *base = 0);
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
        int _indexBase;
        int _rIndexMin,_rIndexMax,_rIndexStep,_rIndexDenom;
        int _muIndexMin,_muIndexMax,_muIndexStep;
        int _rPIndexMin,_rPIndexMax,_rPIndexStep;
        int _rTIndexMin,_rTIndexMax,_rTIndexStep;
        int _zIndexMin,_zIndexMax,_zIndexStep;
        double _r0, _z0;
        AbsCorrelationModel &_base;
	}; // BroadbandModel
    double legendreP(int ell, double mu);
} // baofit

#endif // BAOFIT_BROADBAND_MODEL
