// Created 21-Oct-2014 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_MULTI_CORRELATION_DATA
#define BAOFIT_MULTI_CORRELATION_DATA

#include "baofit/types.h"

#include "baofit/AbsCorrelationData.h"

#include <vector>

namespace baofit {
	class MultiCorrelationData {
	// Represents 3D correlation data.
	public:
		MultiCorrelationData();
		virtual ~MultiCorrelationData();
		// Adds data sets.
		void addDataSet(AbsCorrelationDataCPtr data);
		// Returns the 3D radius in Mpc/h associated with the specified global index.
        virtual double getRadius(int index) const;
        // Returns the cosine of the angle between the separation vector and
        // the line of sight (aka mu) associated with the specified global index.
        virtual double getCosAngle(int index) const;
        // Returns the multipole (0,2,4) associated with the specified global index.
        virtual cosmo::Multipole getMultipole(int index) const;
        // Returns the redshift associated with the specified global index.
        virtual double getRedshift(int index) const;
		// Returns the number of data sets.
		int getNDataSets() const;
	private:
		std::vector<AbsCorrelationDataCPtr> _datasets;
	}; // MultiCorrelationData
	
	inline int MultiCorrelationData::getNDataSets() const { return _datasets.size(); }
	
} // baofit

#endif // BAOFIT_MULTI_CORRELATION_DATA
