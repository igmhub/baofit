// Created 24-Apr-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_METAL_CORRELATION_MODEL
#define BAOFIT_METAL_CORRELATION_MODEL

#include "likely/types.h"

namespace baofit {
	// Represents a model for metal line correlations
	class MetalCorrelationModel {
	public:
	    // Creates a metal correlation model.
	    MetalCorrelationModel(int indexBase);
	    virtual ~MetalCorrelationModel();
	    // Returns the metal correlation in r space at point (r,mu).
	    double _evaluateMetalCorrelation(double r, double mu) const;
	protected:
	private:
	int _indexBase;
    }; // MetalCorrelationModel
} // baofit

#endif // BAOFIT_METAL_CORRELATION_MODEL