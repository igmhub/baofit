// Created 21-Oct-2014 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_MULTI_CORRELATION_DATA
#define BAOFIT_MULTI_CORRELATION_DATA

#include "baofit/AbsCorrelationData.h"

namespace baofit {
	class MultiCorrelationData : public AbsCorrelationData {
	// Represents 3D correlation data
	public:
		MultiCorrelationData();
		virtual ~MultiCorrelationData();
	}; // MultiCorrelationData
} // baofit

#endif // BAOFIT_MULTI_CORRELATION_DATA
