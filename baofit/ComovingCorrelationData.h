// Created 11-Jul-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_COMOVING_CORRELATION_DATA
#define BAOFIT_COMOVING_CORRELATION_DATA

#include "baofit/AbsCorrelationData.h"

#include <vector>

namespace baofit {
	class ComovingCorrelationData : public AbsCorrelationData {
	// Represents 3D correlation data in comoving coordinates (r,mu,z).
	public:
		ComovingCorrelationData(likely::AbsBinningCPtr rBins, likely::AbsBinningCPtr muBins,
		    likely::AbsBinningCPtr zBins, double rmin, double rmax, double rVetoMin, double rVetoMax);
		ComovingCorrelationData(std::vector<likely::AbsBinningCPtr> axes, double rmin, double rmax,
		    double rVetoMin, double rVetoMax);
		virtual ~ComovingCorrelationData();
		// Polymorphic shallow copy so this type of data can be used with likely::BinnedDataResampler.
        virtual ComovingCorrelationData *clone(bool binningOnly = false) const;
        // Returns the 3D radius in Mpc/h associated with the specified global index.
        virtual double getRadius(int index) const;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with the specified global index.
        virtual double getCosAngle(int index) const;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const;
        // Finalize a comoving dataset by pruning to the limits specified in our constructor.
        // No further changes to our "shape" are possible after finalizing. See the documentation
        // for BinnedData::finalize() for details.
        virtual void finalize();
	private:
        double _rmin, _rmax, _rVetoMin, _rVetoMax;
        void _initialize(double rmin, double rmax, double rVetoMin, double rVetoMax);
        void _setIndex(int index) const;
        mutable int _lastIndex;
        mutable std::vector<double> _binCenter;
	}; // ComovingCorrelationData
} // baofit

#endif // BAOFIT_COMOVING_CORRELATION_DATA
