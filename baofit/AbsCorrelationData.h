// Created 24-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_DATA
#define BAOFIT_ABS_CORRELATION_DATA

#include "likely/BinnedData.h"
#include "likely/types.h"

#include "cosmo/types.h"

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
        virtual cosmo::Multipole getMultipole(int index) const;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const = 0;
        // Records the final cuts that should be applied when this dataset is finalized.
        // It is up to subclasses to actually implement these cuts using the protected
        // _applyFinalCuts method in their finalize() implementation. Throws a RuntimeError
        // if any cuts don't make sense. Note that (lMin,lMax) will be ignored if
        // getTransverseBinningType() returns Coordinate and, otherwise, (muMin,muMax)
        // will be ignored.
        void setFinalCuts(double rMin, double rMax, double rVetoMin, double rVetoMax,
            double muMin, double muMax, cosmo::Multipole lMin, cosmo::Multipole lMax,
            double zMin, double zMax);
    private:
        TransverseBinningType _type;
        double _rMin,_rMax,_rVetoMin,_rVetoMax,_muMin,_muMax,_zMin,_zMax;
        cosmo::Multipole _lMin,_lMax;
        bool _haveFinalCuts;
	}; // AbsCorrelationData
	
	inline AbsCorrelationData::TransverseBinningType
	    AbsCorrelationData::getTransverseBinningType() const { return _type; }

} // baofit

#endif // BAOFIT_ABS_CORRELATION_DATA
