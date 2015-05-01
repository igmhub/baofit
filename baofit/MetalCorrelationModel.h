// Created 24-Apr-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_METAL_CORRELATION_MODEL
#define BAOFIT_METAL_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

namespace baofit {
	// Represents an r-space model for metal line correlations
	class MetalCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a metal correlation model.
	    MetalCorrelationModel(AbsCorrelationModel *base = 0);
	    virtual ~MetalCorrelationModel();
	    // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
	protected:
	    // Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h. Uses templates derived from metals added to mock data.
	    virtual double _evaluate(double r, double mu, double z, bool anyChanged) const;
	private:
	int _indexBase;
	AbsCorrelationModel &_base;
    }; // MetalCorrelationModel
} // baofit

#endif // BAOFIT_METAL_CORRELATION_MODEL