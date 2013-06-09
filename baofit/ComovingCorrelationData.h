// Created 11-Jul-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_COMOVING_CORRELATION_DATA
#define BAOFIT_COMOVING_CORRELATION_DATA

#include "baofit/AbsCorrelationData.h"

#include <vector>

namespace baofit {
	class ComovingCorrelationData : public AbsCorrelationData {
	// Represents 3D correlation data in comoving coordinates (r,mu,z) or (rpar,rperp,z).
	public:
        enum CoordinateSystem { PolarCoordinates, CartesianCoordinates, MultipoleCoordinates };
		ComovingCorrelationData(likely::BinnedGrid grid,
		    CoordinateSystem coordinateSystem = PolarCoordinates);
		virtual ~ComovingCorrelationData();
		// Polymorphic shallow copy so this type of data can be used with likely::BinnedDataResampler.
        virtual ComovingCorrelationData *clone(bool binningOnly = false) const;
        // Returns the 3D radius in Mpc/h associated with the specified global index.
        virtual double getRadius(int index) const;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with the specified global index.
        virtual double getCosAngle(int index) const;
        // Returns the multipole (0,2,4) associated with the specified global index.
        virtual cosmo::Multipole getMultipole(int index) const;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const;
        // Finalize a comoving dataset by pruning to the limits specified in our constructor.
        // No further changes to our "shape" are possible after finalizing. See the documentation
        // for BinnedData::finalize() for details.
        virtual void finalize();
	private:
        CoordinateSystem _coordinateSystem;
        void _setIndex(int index) const;
        mutable int _lastIndex;
        mutable std::vector<double> _binCenter;
	}; // ComovingCorrelationData
} // baofit

#endif // BAOFIT_COMOVING_CORRELATION_DATA
