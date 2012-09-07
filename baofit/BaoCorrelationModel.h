// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_BAO_CORRELATION_MODEL
#define BAOFIT_BAO_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "cosmo/types.h"

#include <string>
#include <vector>

namespace baofit {
	// Represents a two-point correlation model parameterized in terms of the relative scale and amplitude
	// of a BAO peak, with possible broadband distortion.
	class BaoCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a new model using the specified tabulated correlation functions at the specified
	    // reference redshift. Set scalePriorMin < scalePriorMax to enable an optional top-hat prior
	    // on the BAO scale.
		BaoCorrelationModel(std::string const &modelrootName,
		    std::string const &fiducialName, std::string const &nowigglesName,
            std::string const &broadbandName, double zref, bool anisotropic = false);
		virtual ~BaoCorrelationModel();
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
        double _zref;
        bool _anisotropic;
        cosmo::RsdCorrelationFunctionPtr _fid, _nw, _bbc, _bb1, _bb2;
        class BBand2;
        typedef boost::shared_ptr<BBand2> BBand2Ptr;
	}; // BaoCorrelationModel
} // baofit

#endif // BAOFIT_BAO_CORRELATION_MODEL
