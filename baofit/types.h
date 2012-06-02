// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_TYPES
#define BAOFIT_TYPES

#include "boost/smart_ptr.hpp"

namespace baofit {
    
    class AbsCorrelationModel;
    typedef boost::shared_ptr<AbsCorrelationModel> AbsCorrelationModelPtr;
    typedef boost::shared_ptr<const AbsCorrelationModel> AbsCorrelationModelCPtr;
    
    class AbsCorrelationData;
    typedef boost::shared_ptr<AbsCorrelationData> AbsCorrelationDataPtr;    
    typedef boost::shared_ptr<const AbsCorrelationData> AbsCorrelationDataCPtr;    

} // baofit

#endif // BAOFIT_TYPES
