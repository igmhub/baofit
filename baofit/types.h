// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_TYPES
#define BAOFIT_TYPES

#include "boost/tr1/memory.hpp"

namespace baofit {
    
    class AbsCorrelationModel;
    typedef std::tr1::shared_ptr<const AbsCorrelationModel> AbsCorrelationModelCPtr;
    
} // baofit

#endif // BAOFIT_TYPES
