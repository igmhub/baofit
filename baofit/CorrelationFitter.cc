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

#include <iostream>
#include <iomanip>
#include <cmath>

namespace local = baofit;

local::CorrelationFitter::CorrelationFitter(AbsCorrelationDataCPtr data, AbsCorrelationModelPtr model,
int covSampleSize)
: _data(data), _model(model), _errorScale(1), _type(data->getTransverseBinningType())
{
    if(!data || 0 == data->getNBinsWithData()) {
        throw RuntimeError("CorrelationFitter: need some data to fit.");
    }
    if(!model) {
        throw RuntimeError("CorrelationFitter: need a model to fit.");
    }
    int n = _data->getNBinsWithData();
    if(0 < covSampleSize) {
        if(covSampleSize <= n+2) {
            throw RuntimeError("CorrelationFitter: cannot fit with covSampleSize <= nbins+2");
        }
        // Correct for the mean bias in the inverse of an unbiased covariance estimate.
        // See eqn (25) of http://arxiv.org/abs/1212.4359
        _icovScale = (covSampleSize - n - 2.)/(covSampleSize - 1.);
    }
    else {
        _icovScale = 1;
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
    // Add any model priors on the parameters. The additional factor of _errorScale
    // is to allow arbitrary error contours to be calculated a la MNCONTOUR.
    double chisq = (0.5*_icovScale*_data->chiSquare(pred) + _model->evaluatePriors())/_errorScale;
    std::cout << std::setprecision(10) << _model->getParameterValue("beta") << " " << _model->getParameterValue("(1+beta)*bias") << " " << _model->getParameterValue("cont-kc") << " " << _model->getParameterValue("BAO alpha-parallel") << " " << _model->getParameterValue("BAO alpha-perp") << " " << 2*chisq << std::endl;
    return chisq;
    //return (0.5*_icovScale*_data->chiSquare(pred) + _model->evaluatePriors())/_errorScale;
}

likely::FunctionMinimumPtr local::CorrelationFitter::fit(std::string const &methodName,
std::string const &config) const {
    likely::FunctionPtr fptr(new likely::Function(*this));
    return _model->findMinimum(fptr,methodName,config);
}

likely::FunctionMinimumPtr local::CorrelationFitter::guess() const {
    likely::FunctionPtr fptr(new likely::Function(*this));
    return _model->guessMinimum(fptr);
}

namespace baofit {
    // A simple MCMC callback that appends sample and fval to samples.
    void mcmcCallback(std::vector<double> &samples, std::vector<double> const &sample, double fval) {
        static int ncall(0);
        samples.insert(samples.end(),sample.begin(),sample.end());
        samples.push_back(fval);
        if(++ncall % 10 == 0) std::cout << "Saved " << ncall << " MCMC trials." << std::endl;
    }
}

void local::CorrelationFitter::mcmc(likely::FunctionMinimumCPtr fminStart, int nchain, int interval,
std::vector<double> &samples) const {
    likely::FunctionPtr fptr(new likely::Function(*this));
    // Use a non-const copy of the input function minimum since the generate method below wants to
    // update it (but we will ignore the updates).
    likely::FunctionMinimumPtr fmin(new likely::FunctionMinimum(*fminStart));
    likely::FitParameters params(fmin->getFitParameters());
    int npar(params.size());
    samples.reserve(nchain*npar);
    samples.resize(0);
    likely::MarkovChainEngine engine(fptr,likely::GradientCalculatorPtr(),params,"saunter");
    int ntrial(nchain*interval);
    likely::MarkovChainEngine::Callback callback = boost::bind(mcmcCallback,boost::ref(samples),_1,_3);
    engine.generate(fmin,ntrial,ntrial,callback,interval);
}
