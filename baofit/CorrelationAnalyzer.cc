// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/CorrelationAnalyzer.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/AbsCorrelationModel.h"
#include "baofit/CorrelationFitter.h"

#include "likely/FunctionMinimum.h"
#include "likely/FitParameterStatistics.h"

#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/foreach.hpp"
#include "boost/math/special_functions/gamma.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace local = baofit;

local::CorrelationAnalyzer::CorrelationAnalyzer(std::string const &method, double rmin, double rmax, bool verbose)
: _method(method), _rmin(rmin), _rmax(rmax), _verbose(verbose)
{
    if(rmin >= rmax) {
        throw RuntimeError("CorrelationAnalyzer: expected rmin < rmax.");
    }
}

local::CorrelationAnalyzer::~CorrelationAnalyzer() { }

void local::CorrelationAnalyzer::setZData(double zdata) {
    if(zdata < 0) {
        throw RuntimeError("CorrelationAnalyzer: expected zdata >= 0.");        
    }
    _zdata = zdata;
}

void local::CorrelationAnalyzer::addData(AbsCorrelationDataCPtr data) {
    _resampler.addObservation(boost::dynamic_pointer_cast<const likely::BinnedData>(data));
}

local::AbsCorrelationDataPtr local::CorrelationAnalyzer::getCombined(bool verbose) const {
    AbsCorrelationDataPtr combined =
        boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(_resampler.combined());
    int nbefore = combined->getNBinsWithData();
    combined->finalize();
    if(verbose) {
        int nafter = combined->getNBinsWithData();
        std::cout << "Combined data has " << nafter << " (" << nbefore
            << ") bins with data after (before) finalizing." << std::endl;
    }
    return combined;    
}

likely::FunctionMinimumPtr local::CorrelationAnalyzer::fitCombined(std::string const &config) const {
    AbsCorrelationDataCPtr combined = getCombined(true);
    CorrelationFitter fitter(combined,_model);
    likely::FunctionMinimumPtr fmin = fitter.fit(_method,config);
    if(_verbose) {
        double chisq = 2*fmin->getMinValue();
        int nbins = combined->getNBinsWithData();
        int npar = fmin->getNParameters(true);
        double prob = 1 - boost::math::gamma_p((nbins-npar)/2,chisq/2);
        std::cout << std::endl << "Results of combined fit: chiSquare / dof = " << chisq << " / ("
            << nbins << '-' << npar << "), prob = " << prob << std::endl << std::endl;
        fmin->printToStream(std::cout);
    }
    return fmin;
}

namespace baofit {
    class CorrelationAnalyzer::AbsSampler {
    public:
        virtual AbsCorrelationDataCPtr nextSample() = 0;
    };
    class CorrelationAnalyzer::JackknifeSampler : public CorrelationAnalyzer::AbsSampler {
    public:
        JackknifeSampler(int ndrop, likely::BinnedDataResampler const &resampler)
        : _ndrop(ndrop), _seqno(0), _resampler(resampler) { }
        virtual AbsCorrelationDataCPtr nextSample() {
            AbsCorrelationDataPtr sample = boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(
                _resampler.jackknife(_ndrop,_seqno++));
            if(sample) sample->finalize();
            return sample;
        }
    private:
        int _ndrop, _seqno;
        likely::BinnedDataResampler const &_resampler;
    };
    class CorrelationAnalyzer::BootstrapSampler : public CorrelationAnalyzer::AbsSampler {
    public:
        BootstrapSampler(int trials, int size, bool fix, likely::BinnedDataResampler const &resampler)
        : _trials(trials), _size(size), _fix(fix), _resampler(resampler), _next(0) { }
        virtual AbsCorrelationDataCPtr nextSample() {
            AbsCorrelationDataPtr sample;
            if(++_next <= _trials) {
                sample = boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(
                    _resampler.bootstrap(_size,_fix));
                sample->finalize();
            }
            return sample;
        }
    private:
        int _trials, _size, _next;
        bool _fix;
        likely::BinnedDataResampler const &_resampler;
    };
    class CorrelationAnalyzer::EachSampler : public CorrelationAnalyzer::AbsSampler {
    public:
        EachSampler(likely::BinnedDataResampler const &resampler)
        : _resampler(resampler), _next(0) { }
        virtual AbsCorrelationDataCPtr nextSample() {
            AbsCorrelationDataPtr sample;
            if(++_next <= _resampler.getNObservations()) {
                sample = boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(
                    _resampler.getObservationCopy(_next-1));
                sample->finalize();
            }
            return sample;
        }
    private:
        int _next;
        likely::BinnedDataResampler const &_resampler;
    };
}

int local::CorrelationAnalyzer::doJackknifeAnalysis(int jackknifeDrop, likely::FunctionMinimumPtr fmin,
likely::FunctionMinimumPtr fmin2, std::string const &refitConfig, std::string const &saveName,
int nsave) const {
    if(jackknifeDrop <= 0) {
        throw RuntimeError("CorrelationAnalyzer::doJackknifeAnalysis: expected jackknifeDrop > 0.");
    }
    CorrelationAnalyzer::JackknifeSampler sampler(jackknifeDrop,_resampler);
    return doSamplingAnalysis(sampler, "Jackknife", fmin, fmin2, refitConfig, saveName, nsave);
}

int local::CorrelationAnalyzer::doBootstrapAnalysis(int bootstrapTrials, int bootstrapSize,
bool fixCovariance, likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
std::string const &refitConfig, std::string const &saveName, int nsave) const {
    if(bootstrapTrials <= 0) {
        throw RuntimeError("CorrelationAnalyzer::doBootstrapAnalysis: expected bootstrapTrials > 0.");
    }
    if(bootstrapSize < 0) {
        throw RuntimeError("CorrelationAnalyzer::doBootstrapAnalysis: expected bootstrapSize >= 0.");
    }
    if(0 == bootstrapSize) bootstrapSize = getNData();
    CorrelationAnalyzer::BootstrapSampler sampler(bootstrapTrials,bootstrapSize,fixCovariance,_resampler);
    return doSamplingAnalysis(sampler, "Bootstrap", fmin, fmin2, refitConfig, saveName, nsave);
}

int local::CorrelationAnalyzer::fitEach(likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
std::string const &refitConfig, std::string const &saveName, int nsave) const {
    CorrelationAnalyzer::EachSampler sampler(_resampler);
    return doSamplingAnalysis(sampler, "Individual", fmin, fmin2, refitConfig, saveName, nsave);    
}

int local::CorrelationAnalyzer::doSamplingAnalysis(CorrelationAnalyzer::AbsSampler &sampler,
std::string const &method, likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
std::string const &refitConfig, std::string const &saveName, int nsave) const {
    if(getNData() <= 1) {
        throw RuntimeError("CorrelationAnalyzer::doSamplingAnalysis: need > 1 observation.");
    }
    if(nsave < 0) {
        throw RuntimeError("CorrelationAnalyzer::doSamplingAnalysis: expected nsave >= 0.");
    }
    if((!fmin2 && 0 < refitConfig.size()) || !!fmin2 && 0 == refitConfig.size()) {
        throw RuntimeError("CorrelationAnalyzer::doSamplingAnalysis: inconsistent refit parameters.");
    }
    // Try to open the specified save file.
    boost::scoped_ptr<std::ofstream> save;
    if(0 < saveName.length()) {
        save.reset(new std::ofstream(saveName.c_str()));
        // Print a header consisting of the number of parameters, the number of dump points,
        // and the number of fits (1 = no-refit, 2 = with refit)
        *save << fmin->getNParameters() << ' ' << nsave << ' ' << (fmin2 ? 2:1) << std::endl;
        // Print the errors in fmin,fmin2.
        BOOST_FOREACH(double pvalue, fmin->getErrors()) {
            *save << pvalue << ' ';
        }
        if(fmin2) {
            BOOST_FOREACH(double pvalue, fmin2->getErrors()) {
                *save << pvalue << ' ';
            }
        }
        *save << std::endl;
        // The first line encodes the inputs fmin,fmin2 just like each sample below, for reference.
        BOOST_FOREACH(double pvalue, fmin->getParameters()) {
            *save << pvalue << ' ';
        }
        *save << 2*fmin->getMinValue() << ' ';
        if(fmin2) {
            BOOST_FOREACH(double pvalue, fmin2->getParameters()) {
                *save << pvalue << ' ';
            }
            *save << 2*fmin2->getMinValue() << ' ';
        }
        if(nsave > 0) {
            dumpModel(*save,fmin,nsave,"",true);
            if(fmin2) dumpModel(*save,fmin2,nsave,"",true);
        }
        *save << std::endl;
    }
    baofit::AbsCorrelationDataCPtr sample;
    // Initialize the parameter value statistics accumulators we will need.
    likely::FitParameterStatisticsPtr refitStats,
        fitStats(new likely::FitParameterStatistics(fmin->getFitParameters()));
    if(fmin2) {
        refitStats.reset(new likely::FitParameterStatistics(fmin2->getFitParameters()));
    }
    int nInvalid(0);
    // Loop over samples.
    int nsamples(0);
    while(sample = sampler.nextSample()) {
        // Fit the sample.
        baofit::CorrelationFitter fitEngine(sample,_model);
        likely::FunctionMinimumPtr sampleMin = fitEngine.fit(_method);
        bool ok = (sampleMin->getStatus() == likely::FunctionMinimum::OK);
        // Refit the sample if requested and the first fit succeeded.
        likely::FunctionMinimumPtr sampleMinRefit;
        if(ok && fmin2) {
            sampleMinRefit = fitEngine.fit(_method,refitConfig);
            // Did this fit succeed also?
            if(sampleMinRefit->getStatus() != likely::FunctionMinimum::OK) ok = false;
        }
        if(ok) {
            // Accumulate the fit results.
            fitStats->update(sampleMin);
            if(refitStats) refitStats->update(sampleMinRefit);
            // Save the fit results, if requested.
            if(save) {
                // Save fit parameter values and chisq.
                BOOST_FOREACH(double pvalue, sampleMin->getParameters()) {
                    *save << pvalue << ' ';
                }
                *save << 2*sampleMin->getMinValue() << ' ';
                // Save any refit parameter values and chisq.
                if(refitStats) {
                    BOOST_FOREACH(double pvalue, sampleMinRefit->getParameters()) {
                        *save << pvalue << ' ';
                    }
                    *save << 2*sampleMinRefit->getMinValue() << ' ';                    
                }
                // Save best-fit model multipoles, if requested.
                if(nsave > 0) {
                    dumpModel(*save,sampleMin,nsave,"",true);
                    if(refitStats) dumpModel(*save,sampleMinRefit,nsave,"",true);
                }
                *save << std::endl;
            }
        }
        else {
            nInvalid++;
        }
        // Print periodic updates while the analysis is running.
        nsamples++;
        if(_verbose && (0 == nsamples%10)) {
            std::cout << "Analyzed " << nsamples << " samples (" << nInvalid << " invalid)" << std::endl;
        }
    }
    // Print a summary of the analysis results.
    std::cout << std::endl << "== " << method << " Fit Results:" << std::endl;
    fitStats->printToStream(std::cout);
    if(refitStats) {
        std::cout << std::endl << "== " << method << " Re-Fit Results:" << std::endl;
        refitStats->printToStream(std::cout);        
    }
    // Close our save file if necessary.
    if(save) save->close();
    return nInvalid;
}

void local::CorrelationAnalyzer::dumpResiduals(std::ostream &out, likely::FunctionMinimumPtr fmin,
std::string const &script) const {
    if(getNData() == 0) {
        throw RuntimeError("CorrelationAnalyzer::dumpResiduals: no observations have been added.");
    }
    AbsCorrelationDataCPtr combined = getCombined();
    AbsCorrelationData::TransverseBinningType type = combined->getTransverseBinningType();
    // Get a copy of the the parameters at this minimum.
    likely::FitParameters parameters(fmin->getFitParameters());
    // Should check that these parameters are "congruent" (have same names?) with model params.
    // assert(parameters.isCongruent(model...))
    // Modify the parameters using the specified script, if any.
    if(0 < script.length()) likely::modifyFitParameters(parameters, script);    
    // Get the parameter values (floating + fixed)
    likely::Parameters parameterValues;
    likely::getFitParameterValues(parameters,parameterValues);
    // Loop over 3D bins in the combined dataset.
    std::vector<double> centers;
    for(likely::BinnedData::IndexIterator iter = combined->begin(); iter != combined->end(); ++iter) {
        int index(*iter);
        out << index;
        combined->getBinCenters(index,centers);
        BOOST_FOREACH(double center, centers) {
            out << ' ' << center;
        }
        double data = combined->getData(index);
        double error = combined->hasCovariance() ? std::sqrt(combined->getCovariance(index,index)) : 0;
        double z = combined->getRedshift(index);
        double r = combined->getRadius(index);
        double predicted;
        if(type == AbsCorrelationData::Coordinate) {
            double mu = combined->getCosAngle(index);
            predicted = _model->evaluate(r,mu,z,parameterValues);
            out  << ' ' << r << ' ' << mu << ' ' << z;
        }
        else {
            cosmo::Multipole multipole = combined->getMultipole(index);
            predicted = _model->evaluate(r,multipole,z,parameterValues);
            out  << ' ' << r << ' ' << (int)multipole << ' ' << z;
        }
        out << ' ' << predicted << ' ' << data << ' ' << error << std::endl;
    }
}

void local::CorrelationAnalyzer::dumpModel(std::ostream &out, likely::FunctionMinimumPtr fmin,
int ndump, std::string const &script, bool oneLine) const {
    if(ndump <= 1) {
        throw RuntimeError("CorrelationAnalyzer::dump: expected ndump > 1.");
    }
    // Get a copy of the the parameters at this minimum.
    likely::FitParameters parameters(fmin->getFitParameters());
    // Should check that these parameters are "congruent" (have same names?) with model params.
    // assert(parameters.isCongruent(model...))
    // Modify the parameters using the specified script, if any.
    if(0 < script.length()) likely::modifyFitParameters(parameters, script);
    // Get the parameter values (floating + fixed)
    likely::Parameters parameterValues;
    likely::getFitParameterValues(parameters,parameterValues);
    // Loop over the specified radial grid.
    double dr((_rmax - _rmin)/(ndump-1));
    for(int rIndex = 0; rIndex < ndump; ++rIndex) {
        double rval(_rmin + dr*rIndex);
        double mono = _model->evaluate(rval,cosmo::Monopole,_zdata,parameterValues);
        double quad = _model->evaluate(rval,cosmo::Quadrupole,_zdata,parameterValues);
        double hexa = _model->evaluate(rval,cosmo::Hexadecapole,_zdata,parameterValues);
        // Output the model predictions for this radius in the requested format.
        if(!oneLine) out << rval;
        out << ' ' << mono << ' ' << quad << ' ' << hexa;
        if(!oneLine) out << std::endl;
    }
}
