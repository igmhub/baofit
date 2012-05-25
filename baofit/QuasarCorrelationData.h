// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_QUASAR_CORRELATION_DATA
#define BAOFIT_QUASAR_CORRELATION_DATA

#include "baofit/AbsCorrelationData.h"

#include "cosmo/types.h"

namespace baofit {
	class QuasarCorrelationData : public AbsCorrelationData {
	// Represents a quasar transmission-fraction (F) correlation function binned in observed
	// coordinates log(lambda2/lambda1), angular separation between lines of sight, and
	// average absorption redshift.
	public:
	    // Creates a new object using the specified cosmology to map the observed coordinates
	    // into co-moving coordinates, and pruning the data to rmin <= r < rmax (in Mpc/h) and
	    // log(lambda2/lambda1) > llmin.
		QuasarCorrelationData(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2,
		    likely::AbsBinningCPtr axis3, cosmo::AbsHomogeneousUniversePtr cosmology);
		virtual ~QuasarCorrelationData();
        // Returns the 3D radius in Mpc/h associated with 0 <= offset < getSize().
        virtual double getRadius(int offset) const;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with 0 <= offset < getSize().
        // Returns zero if the data is azimuthally averaged.
        virtual double getCosAngle(int offset) const;
        // Returns the redshift associated with 0 <= offset < getSize().
        virtual double getRedshift(int offset) const;
        // Finalize a quasar dataset by...
        void finalize(double rmin, double rmax, double llmin);
	private:
        cosmo::AbsHomogeneousUniversePtr _cosmology;
        std::vector<double> _rLookup, _muLookup, _zLookup;
	}; // QuasarCorrelationData
} // baofit

#endif // BAOFIT_QUASAR_CORRELATION_DATA
