// Created 02-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_MULTIPOLE_CORRELATION_DATA
#define BAOFIT_MULTIPOLE_CORRELATION_DATA

#include "baofit/AbsCorrelationData.h"

#include <iosfwd>

namespace baofit {
	class MultipoleCorrelationData : public AbsCorrelationData {
	// Represents a 3D correlation function with its transverse degree of freedom projected
	// into multipole space.
	public:
	    //Creates a new multipole correlation dataset with the specified binning and final
	    // cuts rmin <= r < rmax, ellmin <= ell <= ellmax.
		MultipoleCorrelationData(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2,
		    likely::AbsBinningCPtr axis3);
        MultipoleCorrelationData(std::vector<likely::AbsBinningCPtr> axes);
		virtual ~MultipoleCorrelationData();
		// Polymorphic shallow copy so this type of data can be used with likely::BinnedDataResampler.
        virtual MultipoleCorrelationData *clone(bool binningOnly = false) const;
        // Returns the 3D radius in Mpc/h associated with the specified global index.
        virtual double getRadius(int index) const;
        // Returns the multipole (0,2,4) associated with the specified global index.
        virtual cosmo::Multipole getMultipole(int index) const;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const;
        // Finalize a multipole dataset by pruning to the final cuts specified in our
        // constructor. No further changes to our "shape" are possible after finalizing.
        // See the documentation for BinnedData::finalize() for details.
        virtual void finalize();
        // Dumps all available multipoles to the specified output stream. Uses the weights
        // provided for error bars, if weights.size() > 0, or else uses the covariance matrix
        // diagonal elements. Output format is separation followed by "value error" pairs for each
        // redshift slice (outer loop) and multipole (inner loop), for a total 1+2*nz*nell columns.
        void dump(std::ostream &out, double rmin, double rmax,
            std::vector<double> const &weights) const;
	private:
        void _setIndex(int index) const;
        mutable int _lastIndex;
        mutable std::vector<double> _binCenter;
        mutable double _rLast, _zLast;
        mutable cosmo::Multipole _ellLast;
	}; // MultipoleCorrelationData
} // baofit

#endif // BAOFIT_MULTIPOLE_CORRELATION_DATA
