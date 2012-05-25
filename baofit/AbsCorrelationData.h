// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_DATA
#define BAOFIT_ABS_CORRELATION_DATA

#include "likely/types.h"

namespace baofit {
	class AbsCorrelationData {
	// Represents an abstract mapping from binned data to the (r,mu,z) coordinates
	// used by an AbsCorrelationModel.
	public:
		AbsCorrelationData(likely::BinnedDataCPtr data);
		virtual ~AbsCorrelationData();
		// Replaces our binned data. Throws a RuntimeError unless the new data is
		// congruent with the original data.
        void setData(likely::BinnedDataCPtr data);
        // Returns the number of bins with data.
        int getSize() const;
        // Returns the data associated with 0 <= offset < getSize().
        double getData(int offset) const;
        // Returns the 3D radius in Mpc/h associated with 0 <= offset < getSize().
        virtual double getRadius(int offset) const = 0;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with 0 <= offset < getSize().
        // Returns zero if the data is azimuthally averaged.
        virtual double getCosAngle(int offset) const = 0;
        // Returns the redshift associated with 0 <= offset < getSize().
        virtual double getRedshift(int offset) const = 0;
	private:
        likely::BinnedDataCPtr _data;
	}; // AbsCorrelationData
} // baofit

#endif // BAOFIT_ABS_CORRELATION_DATA
