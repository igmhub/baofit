// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BAO_CORRELATION_MODEL
#define BAOFIT_BAO_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"
#include "baofit/types.h"

#include "cosmo/types.h"

#include <string>

namespace baofit {
	// Represents a two-point correlation model parameterized in terms of the relative scale and amplitude
	// of a BAO peak, with possible broadband distortion.
	class BaoCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a new model using the specified tabulated correlation functions at the specified
	    // reference redshift.
		BaoCorrelationModel(std::string const &modelrootName,
		    std::string const &fiducialName, std::string const &nowigglesName,
            std::string const &metalModelName, std::string const &distAdd,
            std::string const &distMul, double distR0, double zref,
            double OmegaMatter, bool anisotropic = false, bool decoupled = false,
            bool metalModel = false, bool metalModelInterpolate = false,
            bool metalCIV = false, bool toyMetal = false, bool combinedBias = false,
            bool combinedScale = false, bool crossCorrelation = false);
		virtual ~BaoCorrelationModel();
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
	protected:
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
        virtual double _evaluate(double r, double mu, double z, bool anyChanged, int index) const;
        virtual double _evaluateKSpace(double k, double mu_k, double pk, double z) const;
        virtual int _getIndexBase() const;
	private:
        AbsCorrelationModelPtr _metalCorr, _distortAdd, _distortMul;
        bool _anisotropic, _decoupled, _combinedScale;
        int _indexBase, _combScaleBase;
        cosmo::CorrelationFunctionPtr _fid0, _fid2, _fid4, _nw0, _nw2, _nw4;
	}; // BaoCorrelationModel
} // baofit

#endif // BAOFIT_BAO_CORRELATION_MODEL
