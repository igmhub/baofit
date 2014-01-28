// Created 27-Jan-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BAO_K_SPACE_CORRELATION_MODEL
#define BAOFIT_BAO_K_SPACE_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"
#include "baofit/types.h"

#include "cosmo/types.h"

#include <string>

namespace baofit {
	// Represents a two-point correlation model derived from tabulated power spectra (with and
	// without wiggles), a parameterized anisotropic k-space distortion model D(k,mu_k), and
	// anisotropic BAO scale parameters. Optional multiplicative and/or additive broadband
	// distortion can also be added in r space.
	class BaoKSpaceCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a new model using the specified tabulated power spectra at the specified
	    // reference redshift.
		BaoKSpaceCorrelationModel(std::string const &modelrootName,
		    std::string const &fiducialName, std::string const &nowigglesName,
            std::string const &distAdd, std::string const &distMul, double distR0,
            double zref, bool anisotropic = false, bool decoupled = false,
            bool crossCorrelation = false);
		virtual ~BaoKSpaceCorrelationModel();
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
        // Evaluates our k-space distortion model D(k,mu_k) using our current parameter values.
        double evaluateKSpaceDistortion(double k, double mu_k) const;
	protected:
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double _evaluate(double r, double mu, double z, bool anyChanged) const;
	private:
        AbsCorrelationModelPtr _distortAdd, _distortMul;
        bool _anisotropic, _decoupled;
        int _indexBase;
        cosmo::CorrelationFunctionPtr _fid0, _fid2, _fid4, _nw0, _nw2, _nw4;
        cosmo::TabulatedPowerCPtr _Pfid, _Pnw;
        cosmo::DistortedPowerCorrelationPtr _Xifid, _Xinw;
	}; // BaoKSpaceCorrelationModel
} // baofit

#endif // BAOFIT_BAO_K_SPACE_CORRELATION_MODEL
