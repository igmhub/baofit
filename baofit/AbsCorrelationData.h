// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_DATA
#define BAOFIT_ABS_CORRELATION_DATA

#include "likely/BinnedData.h"
#include "likely/types.h"

namespace baofit {
	class AbsCorrelationData : public likely::BinnedData {
	// Represents data binned in variables that map to the (r,mu,z) coordinates
	// used by an AbsCorrelationModel.
	public:
        enum TransverseBinningType { Coordinate, Multipole };
		AbsCorrelationData(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2,
		    likely::AbsBinningCPtr axis3, TransverseBinningType type);
        AbsCorrelationData(std::vector<likely::AbsBinningCPtr> axes, TransverseBinningType type);
		virtual ~AbsCorrelationData();
		// Returns the type of transverse binning: Coordinate means that a mu value can
		// be associated with each bin, and Multipole means that an ell value can be
		// associated with each bin.
        TransverseBinningType getTransverseBinningType() const;
        // Returns the 3D radius in Mpc/h associated with the specified global index.
        virtual double getRadius(int index) const = 0;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with the specified global index.
        // Will only be called if getTransverseBinningType() returns Coordinate.
        virtual double getCosAngle(int index) const;
        // Returns the multipole associated with the specified global index.
        // Will only be called if getTransverseBinningType() returns Multipole.
        virtual int getMultipole(int index) const;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const = 0;
    private:
        TransverseBinningType _type;
	}; // AbsCorrelationData
	
	inline AbsCorrelationData::TransverseBinningType
	    AbsCorrelationData::getTransverseBinningType() const { return _type; }

} // baofit

#endif // BAOFIT_ABS_CORRELATION_DATA
