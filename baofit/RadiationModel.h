// Created 5-May-2015 by Ignasi Perez-Rafols (University of Barcelona) <iprafols@icc.ub.edu>

#ifndef BAOFIT_RADIATION_MODEL
#define BAOFIT_RADIATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include <cmath>
#include <complex>

namespace baofit {
	// Represents a k-space model for quasar radiation effects
	class RadiationModel : public AbsCorrelationModel {
	public:
	    // Creates a quasar radiation model.
	    RadiationModel(AbsCorrelationModel *base = 0);
	    virtual ~RadiationModel();
	    // Prints a multi-line description of this object to the specified output stream.
            virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
            // Returns the correlation function evaluated in k-space where (k,mu_k) is
            // the pair separation and z is their average redshift. The separation k should
            // be provided in Mpc/h.            
            double _evaluateRadiation(double k, double mu_k, double z) const;
	protected:
	    // Returns the correlation function evaluated in redshift space where (r,mu) is
	    // the pair separation and z is their average redshift. The separation k should
		// be provided in Mpc/h.
	    virtual double _evaluate(double r, double mu, double z, bool anyChanged) const;
        
	private:
	int _indexBase;
	AbsCorrelationModel &_base;
    }; // RadiationModel
} // baofit

#endif // BAOFIT_RADIATION_MODEL
