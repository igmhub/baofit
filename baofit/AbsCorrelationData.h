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
		AbsCorrelationData(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2,
		    likely::AbsBinningCPtr axis3);
        AbsCorrelationData(std::vector<likely::AbsBinningCPtr> axes);
		virtual ~AbsCorrelationData();
        // Returns the 3D radius in Mpc/h associated with the specified global index.
        virtual double getRadius(int index) const = 0;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with the specified global index.
        // Returns zero if the data is azimuthally averaged.
        virtual double getCosAngle(int index) const = 0;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const = 0;
        // Finalizes a dataset before fitting to a model.
        virtual void finalize();
	}; // AbsCorrelationData	
} // baofit

#endif // BAOFIT_ABS_CORRELATION_DATA
