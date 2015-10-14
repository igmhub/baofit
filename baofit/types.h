// Created 06-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_TYPES
#define BAOFIT_TYPES

#include "boost/smart_ptr.hpp"

namespace baofit {
    
    class AbsCorrelationModel;
    typedef boost::shared_ptr<AbsCorrelationModel> AbsCorrelationModelPtr;
    typedef boost::shared_ptr<const AbsCorrelationModel> AbsCorrelationModelCPtr;
    
    class NonLinearCorrectionModel;
    typedef boost::shared_ptr<NonLinearCorrectionModel> NonLinearCorrectionModelPtr;
    typedef boost::shared_ptr<const NonLinearCorrectionModel> NonLinearCorrectionModelCPtr;
    
    class AbsCorrelationData;
    typedef boost::shared_ptr<AbsCorrelationData> AbsCorrelationDataPtr;    
    typedef boost::shared_ptr<const AbsCorrelationData> AbsCorrelationDataCPtr;
    
    class RadiationModel;
    typedef boost::shared_ptr<RadiationModel> RadiationModelPtr;
    typedef boost::shared_ptr<const RadiationModel> RadiationModelCPtr;
    

} // baofit

#endif // BAOFIT_TYPES
