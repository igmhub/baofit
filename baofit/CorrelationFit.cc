// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/CorrelationFit.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/AbsCorrelationModel.h"

#include "likely/AbsEngine.h"

namespace local = baofit;

local::CorrelationFit::CorrelationFit(AbsCorrelationDataCPtr data, AbsCorrelationModelCPtr model)
: _data(data), _model(model), _errorScale(1)
{
}

local::CorrelationFit::~CorrelationFit() { }

void local::CorrelationFit::setErrorScale(double scale) {
    if(scale <= 0) {
        throw RuntimeError("CorrelationFit::setErrorScale: expected scale > 0.");
    }
    _errorScale = scale;
}

double local::CorrelationFit::operator()(likely::Parameters const &params) const {
    // Check that we have the expected number of parameters.
    if(params.size() != _model->getNParameters()) {
        throw RuntimeError("CorrelationFit: got unexpected number of parameters.");
    }
    // Loop over the dataset bins.
    std::vector<double> pred;
    pred.reserve(_data->getNBinsWithData());
    static int offset(0);
    for(baofit::AbsCorrelationData::IndexIterator iter = _data->begin(); iter != _data->end(); ++iter) {
        int index(*iter);
        double r = _data->getRadius(index);
        double mu = _data->getCosAngle(index);
        double z = _data->getRedshift(index);
        double predicted = _model->evaluate(r,mu,z,params);
        pred.push_back(predicted);
        //!!DK
        if(offset++ < 5) {
            std::cout << "NEW rr,mu,z = " << r << ',' << mu << ',' << z
                << " obs=" << _data->getData(index) << ", pred=" << predicted << std::endl;
        }
        //!!DK
    }
    // UP=0.5 is already hardcoded so we need a factor of 2 here since we are
    // calculating a chi-square. Apply an additional factor of _errorScale to
    // allow different error contours to be calculated.
    return 0.5*_data->chiSquare(pred)/_errorScale;
}

likely::FunctionMinimumPtr local::CorrelationFit::fit(std::string const &methodName) const {
    likely::FunctionPtr fptr(new likely::Function(*this));
    return likely::findMinimum(fptr,_model->getParameters(),methodName);
}