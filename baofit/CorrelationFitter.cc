// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/CorrelationFitter.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationModel.h"

#include "likely/AbsEngine.h"
#include "likely/FitParameter.h"
#include "likely/FunctionMinimum.h"
#include "likely/MarkovChainEngine.h"

#include "boost/bind.hpp"
#include "boost/ref.hpp"

namespace local = baofit;

local::CorrelationFitter::CorrelationFitter(AbsCorrelationDataCPtr data, AbsCorrelationModelPtr model)
: _data(data), _model(model), _errorScale(1), _type(data->getTransverseBinningType())
{
    if(!data || 0 == data->getNBinsWithData()) {
        throw RuntimeError("CorrelationFitter: need some data to fit.");
    }
    if(!model) {
        throw RuntimeError("CorrelationFitter: need a model to fit.");
    }
}

local::CorrelationFitter::~CorrelationFitter() { }

void local::CorrelationFitter::setErrorScale(double scale) {
    if(scale <= 0) {
        throw RuntimeError("CorrelationFitter::setErrorScale: expected scale > 0.");
    }
    _errorScale = scale;
}

void local::CorrelationFitter::getPrediction(likely::Parameters const &params,
std::vector<double> &prediction) const {
    prediction.reserve(_data->getNBinsWithData());
    prediction.resize(0);
    for(baofit::AbsCorrelationData::IndexIterator iter = _data->begin(); iter != _data->end(); ++iter) {
        int index(*iter);
        double z = _data->getRedshift(index);
        double r = _data->getRadius(index);
        double predicted;
        if(_type == AbsCorrelationData::Coordinate) {
            double mu = _data->getCosAngle(index);
            predicted = _model->evaluate(r,mu,z,params);
        }
        else {
            cosmo::Multipole multipole = _data->getMultipole(index);
            predicted = _model->evaluate(r,multipole,z,params);
        }
        prediction.push_back(predicted);
    }    
}

double local::CorrelationFitter::operator()(likely::Parameters const &params) const {
    // Check that we have the expected number of parameters.
    if(params.size() != _model->getNParameters()) {
        throw RuntimeError("CorrelationFitter: got unexpected number of parameters.");
    }
    // Calculate the prediction vector for these parameter values.
    std::vector<double> pred;
    getPrediction(params,pred);
    // Scale chiSquare by 0.5 since the likely minimizer expects a -log(likelihood).
    // Add any model prior on the parameters. The additional factor of _errorScale
    // is to allow arbitrary error contours to be calculated a la MNCONTOUR.
    return (0.5*_data->chiSquare(pred) + _model->evaluatePrior(params))/_errorScale;
}

likely::FunctionMinimumPtr local::CorrelationFitter::fit(std::string const &methodName,
std::string const &config) const {
    likely::FunctionPtr fptr(new likely::Function(*this));
    return _model->findMinimum(fptr,methodName,config);
}

namespace baofit {
    // A simple MCMC callback that appends sample and fval to samples.
    void mcmcCallback(std::vector<double> &samples, std::vector<double> const &sample, double fval) {
        samples.insert(samples.end(),sample.begin(),sample.end());
        samples.push_back(fval);
    }
}

void local::CorrelationFitter::mcmc(likely::FunctionMinimumPtr fmin, int nchain, int interval,
std::vector<double> &samples) const {
    likely::FunctionPtr fptr(new likely::Function(*this));
    likely::FitParameters params(fmin->getFitParameters());
    int npar(params.size());
    samples.reserve(nchain*npar);
    samples.resize(0);
    likely::MarkovChainEngine engine(fptr,likely::GradientCalculatorPtr(),params,"saunter");
    int ntrial(nchain*interval);
    likely::MarkovChainEngine::Callback callback = boost::bind(mcmcCallback,boost::ref(samples),_1,_3);
    engine.generate(fmin,ntrial,ntrial,callback,interval);
}
