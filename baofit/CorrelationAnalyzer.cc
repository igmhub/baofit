// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "baofit/CorrelationAnalyzer.h"
#include "baofit/RuntimeError.h"
#include "baofit/AbsCorrelationData.h"
#include "baofit/AbsCorrelationModel.h"
#include "baofit/CorrelationFitter.h"

#include "likely/FunctionMinimum.h"
#include "likely/FitParameter.h"
#include "likely/CovarianceMatrix.h"
#include "likely/CovarianceAccumulator.h"
#include "likely/FitParameterStatistics.h"

#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/foreach.hpp"
#include "boost/utility.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/math/special_functions/gamma.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iterator>

namespace local = baofit;

local::CorrelationAnalyzer::CorrelationAnalyzer(std::string const &method, double rmin, double rmax,
int covSampleSize, bool verbose, bool scalarWeights)
: _method(method), _rmin(rmin), _rmax(rmax), _covSampleSize(covSampleSize),
_verbose(verbose), _resampler(scalarWeights)
{
    if(rmin >= rmax) {
        throw RuntimeError("CorrelationAnalyzer: expected rmin < rmax.");
    }
}

local::CorrelationAnalyzer::~CorrelationAnalyzer() { }

int local::CorrelationAnalyzer::addData(AbsCorrelationDataCPtr data, int reuseCovIndex) {
    return _resampler.addObservation(
        boost::dynamic_pointer_cast<const likely::BinnedData>(data),reuseCovIndex);
}

local::AbsCorrelationDataPtr local::CorrelationAnalyzer::getCombined(bool verbose, bool finalized) const {
    AbsCorrelationDataPtr combined =
        boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(_resampler.combined());
    int nbefore = combined->getNBinsWithData();
    if(finalized) combined->finalize();
    if(verbose && finalized) {
        int nafter = combined->getNBinsWithData();
        std::cout << "Combined data has " << nafter << " (" << nbefore
            << ") bins with data after (before) finalizing." << std::endl;
    }
    return combined;    
}

void local::CorrelationAnalyzer::compareEach(std::string const &saveName, bool finalized) const {
    if(_resampler.usesScalarWeights()) {
        throw RuntimeError("CorrelationAnalyzer::compareEach: not supported with scalar weights.");
    }
    // Open the output file.
    std::ofstream out(saveName.c_str());
    // Get our combined data to use as a reference
    baofit::AbsCorrelationDataCPtr refData = getCombined(false,finalized);
    // Get the combined data's covariance
    likely::CovarianceMatrixCPtr Cref = refData->getCovarianceMatrix();
    int nbins = refData->getNBinsWithData();
    // Initialize storage.
    std::vector<double> eigenvalues,eigenvectors,chi2modes;
    // Initialize formats.
    boost::format summary("%4d %.6lf %8.1lf %10.4lf\n"), oneValue(" %.4lg");
    // Loop over observations.
    std::cout << "   N     Prob     Chi2   log|C|/n" << std::endl;
    for(int obsIndex = 0; obsIndex < _resampler.getNObservations(); ++obsIndex) {
        AbsCorrelationDataPtr observation =
            boost::dynamic_pointer_cast<baofit::AbsCorrelationData>(
            _resampler.getObservationCopy(obsIndex));
        if(finalized) observation->finalize();
        // Calculate log(|C|)/nbins
        likely::CovarianceMatrixPtr Csub(new likely::CovarianceMatrix(*observation->getCovarianceMatrix()));
        double logDet = Csub->getLogDeterminant()/nbins;
        // Subtract the reference data covariance.
        for(int row = 0; row < nbins; ++row) {
            for(int col = 0; col <= row; ++col) {
                double value(Csub->getCovariance(row,col)-Cref->getCovariance(row,col));
                Csub->setCovariance(row,col,value);
            }
        }
        // Subtract the reference data vector.
        std::vector<double> delta;
        for(likely::BinnedData::IndexIterator iter = refData->begin(); iter != refData->end(); ++iter) {
            delta.push_back(observation->getData(*iter) - refData->getData(*iter));
        }
        // Calculate the chi-square of this observation relative to the combined data
        double chi2 = Csub->chiSquareModes(delta,eigenvalues,eigenvectors,chi2modes);
        // Calculate the corresponding chi-square probability and print a summary for this observation.
        double prob = 1 - boost::math::gamma_p(nbins/2.,chi2/2);
        std::cout << summary % obsIndex % prob % chi2 % logDet;
        // Save the contributions of each mode to the output file.
        out << obsIndex << ' ' << logDet << ' ' << chi2;
        for(int i = 0; i < nbins; ++i) out << oneValue % chi2modes[i];
        out << std::endl;
    }
    out.close();
}

bool local::CorrelationAnalyzer::printScaleZEff(likely::FunctionMinimumCPtr fmin, double zref,
std::string const &scaleName) const {
    likely::FitParameters params = fmin->getFitParameters();
    std::vector<std::string> pnames;
    bool onlyFloating(true);
    likely::getFitParameterNames(params,pnames,onlyFloating);
    // Is "gamma-scale" a floating parameter of this fit?
    std::vector<std::string>::iterator found = std::find(pnames.begin(),pnames.end(),"gamma-scale");
    if(found == pnames.end()) return false;
    int gammaIndex = std::distance(pnames.begin(),found);
    // Is scaleName a floating parameter of this fit?
    found = std::find(pnames.begin(),pnames.end(),scaleName);
    if(found == pnames.end()) return false;
    int scaleIndex = std::distance(pnames.begin(),found);
    // Look up the fit results for gamma,scale
    likely::Parameters pvalues = fmin->getParameters(onlyFloating);
    likely::Parameters perrors = fmin->getErrors(onlyFloating);
    double scale(pvalues[scaleIndex]), dscale(perrors[scaleIndex]);
    double gamma(pvalues[gammaIndex]), dgamma(perrors[gammaIndex]);
    // Look up the (scale,gamma) covariance.
    likely::CovarianceMatrixCPtr cov = fmin->getCovariance();
    if(!cov) return false;
    double rho(cov->getCovariance(scaleIndex,gammaIndex)/(dscale*dgamma));
    // Calculate zeff
    double a = dscale/(scale*dgamma), b = 1/(2*gamma);
    double logz = -b - a*rho + std::sqrt(b*b - a*a*(1-rho*rho));
    double zEff = std::exp(logz)*(1 + zref) - 1;
    // Calculate the scale at zref.
    double ratio = (1+zEff)/(1+zref);
    double evol = std::pow(ratio,gamma);
    double scaleEff = scale*evol;
    // Calculate the scale error at zref.
    double logRatio = std::log(ratio);
    double tmp = scale*dgamma*logRatio;
    double dscaleEff = evol*std::sqrt(dscale*dscale + 2*rho*scale*dscale*dgamma*logRatio + tmp*tmp);
    // Print results.
    std::vector<double> errors(1,dscaleEff);
    std::cout << boost::format("%18s") % scaleName << "(zeff = " << boost::format("%.3f") % zEff << ") = "
        << likely::roundValueWithError(scaleEff,errors) << std::endl;
    return true;
}

likely::FunctionMinimumPtr local::CorrelationAnalyzer::fitSample(
AbsCorrelationDataCPtr sample, std::string const &config) const {
    CorrelationFitter fitter(sample,_model,_covSampleSize);
    likely::FunctionMinimumPtr fmin = fitter.fit(_method,config);
    if(_verbose) {
        double chisq = 2*fmin->getMinValue();
        int nbins = sample->getNBinsWithData();
        int npar = fmin->getNParameters(true);
        if(chisq > 0 && nbins > npar) {
            double prob = 1 - boost::math::gamma_p((nbins-npar)/2.,chisq/2);
            std::cout << std::endl << "Fit results: chiSquare / dof = " << chisq << " / ("
                << nbins << '-' << npar << "), prob = " << prob << ", log(det(Covariance)) = "
                << sample->getCovarianceMatrix()->getLogDeterminant() << std::endl << std::endl;
        }
        else {
            std::cout << std::endl << "Fit is ill-conditioned with nbins = " << nbins
                << ", npar = " << npar << ", chiSquare = " << chisq << std::endl << std::endl;
        }
        fmin->printToStream(std::cout);
    }
    return fmin;
}

void local::CorrelationAnalyzer::dumpChisquare(std::ostream &out, likely::FunctionMinimumPtr fmin,
AbsCorrelationDataCPtr combined) const {
    double chisq = 2*fmin->getMinValue();
    int nbins = combined->getNBinsWithData();
    int npar = fmin->getNParameters(true);
    if(chisq > 0 && nbins > npar) {
        double prob = 1 - boost::math::gamma_p((nbins-npar)/2.,chisq/2);
        out << boost::lexical_cast<std::string>(chisq) << ' ' << boost::lexical_cast<std::string>(nbins) << ' '
            << boost::lexical_cast<std::string>(npar) << ' ' << boost::lexical_cast<std::string>(prob);
    }
    else {
        out << boost::lexical_cast<std::string>(chisq) << ' ' << boost::lexical_cast<std::string>(nbins) << ' '
            << boost::lexical_cast<std::string>(npar) << ' ' << boost::lexical_cast<std::string>(-1);
    }
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
    class CorrelationAnalyzer::ToyMCSampler : public CorrelationAnalyzer::AbsSampler {
    public:
        ToyMCSampler(int ngen, AbsCorrelationDataPtr prototype, std::vector<double> truth,
        std::string const &filename)
        : _remaining(ngen), _prototype(prototype), _truth(truth), _first(true), _filename(filename) { }
        virtual AbsCorrelationDataCPtr nextSample() {
            AbsCorrelationDataPtr sample;
            if(_remaining-- > 0) {
                // Generate a noise vector sampling from the prototype's covariance.
                _prototype->getCovarianceMatrix()->sample(_noise);
                // Clone our prototype (which only copies the covariance smart pointer, not
                // the whole matrix)
                sample.reset((AbsCorrelationData*)_prototype->clone());
                // Overwrite the bin values with truth+noise
                std::vector<double>::const_iterator nextTruth(_truth.begin()), nextNoise(_noise.begin());
                for(likely::BinnedData::IndexIterator iter = _prototype->begin();
                iter != _prototype->end(); ++iter) {
                    sample->setData(*iter,(*nextTruth++) + (*nextNoise++));
                }
                // We don't finalize here because the prototype should already be finalized.
                // Save this file?
                if(_first && _filename.length() > 0) {
                    std::ofstream out(_filename.c_str());
                    sample->saveData(out);
                    out.close();
                }
                _first = false;
            }
            return sample;
        }
    private:
        int _remaining;
        bool _first;
        std::string _filename;
        AbsCorrelationDataPtr _prototype;
        std::vector<double> _truth, _noise;
    };
}

int local::CorrelationAnalyzer::doJackknifeAnalysis(int jackknifeDrop, likely::FunctionMinimumPtr fmin,
likely::FunctionMinimumPtr fmin2, std::string const &refitConfig, std::string const &saveName,
int nsave, double zsave) const {
    if(jackknifeDrop <= 0) {
        throw RuntimeError("CorrelationAnalyzer::doJackknifeAnalysis: expected jackknifeDrop > 0.");
    }
    if(getNData() <= 1) {
        throw RuntimeError("CorrelationAnalyzer::doJackknifeAnalysis: need > 1 observation.");
    }
    CorrelationAnalyzer::JackknifeSampler sampler(jackknifeDrop,_resampler);
    return doSamplingAnalysis(sampler, "Jackknife", fmin, fmin2, refitConfig, saveName, nsave, zsave);
}

int local::CorrelationAnalyzer::doBootstrapAnalysis(int bootstrapTrials, int bootstrapSize,
bool fixCovariance, likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
std::string const &refitConfig, std::string const &saveName, int nsave, double zsave) const {
    if(bootstrapTrials <= 0) {
        throw RuntimeError("CorrelationAnalyzer::doBootstrapAnalysis: expected bootstrapTrials > 0.");
    }
    if(bootstrapSize < 0) {
        throw RuntimeError("CorrelationAnalyzer::doBootstrapAnalysis: expected bootstrapSize >= 0.");
    }
    if(getNData() <= 1) {
        throw RuntimeError("CorrelationAnalyzer::doBootstrapAnalysis: need > 1 observation.");
    }
    if(0 == bootstrapSize) bootstrapSize = getNData();
    CorrelationAnalyzer::BootstrapSampler sampler(bootstrapTrials,bootstrapSize,fixCovariance,_resampler);
    return doSamplingAnalysis(sampler, "Bootstrap", fmin, fmin2, refitConfig, saveName, nsave, zsave);
}

int local::CorrelationAnalyzer::fitEach(likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
std::string const &refitConfig, std::string const &saveName, int nsave, double zsave) const {
    CorrelationAnalyzer::EachSampler sampler(_resampler);
    return doSamplingAnalysis(sampler, "Individual", fmin, fmin2, refitConfig, saveName, nsave, zsave);    
}

int local::CorrelationAnalyzer::doToyMCSampling(int ngen, std::string const &mcConfig,
std::string const &mcSaveFile, double varianceScale, likely::FunctionMinimumPtr fmin,
likely::FunctionMinimumPtr fmin2, std::string const &refitConfig,
std::string const &saveName, int nsave, double zsave) const {
    if(ngen <= 0) {
        throw RuntimeError("CorrelationAnalyzer::doMCSampling: expected ngen > 0.");
    }
    // Get a copy of our (unfinalized!) combined dataset to use as a prototype.
    AbsCorrelationDataPtr prototype = getCombined(false,false);
    if(!prototype->hasCovariance()) {
        throw RuntimeError("CorrelationAnalyzer::doMCSampling: no covariance available.");
    }
    // Scale the prototype covariance, if requested.
    if(varianceScale <= 0) {
        throw RuntimeError("CorrelationAnalyzer::doMCSampling: expected varianceScale > 0.");
    }
    if(varianceScale != 1) {
        likely::CovarianceMatrixPtr covariance(new likely::CovarianceMatrix(*prototype->getCovarianceMatrix()));
        covariance->applyScaleFactor(varianceScale);
        prototype->setCovarianceMatrix(covariance);
    }
    // Finalize now, after any covariance scaling.
    prototype->finalize();
    // Configure the fit parameters for generating the truth vector.
    likely::FitParameters parameters = fmin->getFitParameters();
    likely::modifyFitParameters(parameters,mcConfig);
    std::vector<double> pvalues;
    likely::getFitParameterValues(parameters,pvalues);
    // Build a fitter to calculate the truth vector.
    CorrelationFitter fitter(prototype, _model, _covSampleSize);
    // Calculate the truth vector.
    std::vector<double> truth;
    fitter.getPrediction(pvalues,truth);
    // Build the sampler for this analysis.
    CorrelationAnalyzer::ToyMCSampler sampler(ngen,prototype,truth,mcSaveFile);
    return doSamplingAnalysis(sampler, "MonteCarlo", fmin, fmin2, refitConfig, saveName, nsave, zsave);
}

namespace baofit {
    // An implementation class to save the results of a sampling analysis in a standard format.
    class SamplingOutput : public boost::noncopyable {
    public:
        SamplingOutput(likely::FunctionMinimumCPtr fmin, likely::FunctionMinimumCPtr fmin2,
        std::string const &saveName, int nsave, double zsave, CorrelationAnalyzer const &parent)
        : _nsave(nsave), _zsave(zsave), _parent(parent) {
            if(0 < saveName.length()) {
                _save.reset(new std::ofstream(saveName.c_str()));
                // Print a header consisting of the number of parameters, the number of dump points,
                // and the number of fits (1 = no-refit, 2 = with refit)
                *_save << fmin->getNParameters() << ' ' << _nsave << ' ' << (fmin2 ? 2:1) << std::endl;
                // Print the errors in fmin,fmin2.
                BOOST_FOREACH(double pvalue, fmin->getErrors()) {
                    *_save << pvalue << ' ';
                }
                if(fmin2) {
                    BOOST_FOREACH(double pvalue, fmin2->getErrors()) {
                        *_save << pvalue << ' ';
                    }
                }
                *_save << std::endl;
                // The first line encodes the inputs fmin,fmin2 just like each sample below, for reference.
                BOOST_FOREACH(double pvalue, fmin->getParameters()) {
                    *_save << pvalue << ' ';
                }
                *_save << 2*fmin->getMinValue() << ' ';
                if(fmin2) {
                    BOOST_FOREACH(double pvalue, fmin2->getParameters()) {
                        *_save << pvalue << ' ';
                    }
                    *_save << 2*fmin2->getMinValue() << ' ';
                }
                if(_nsave > 0) {
                    _parent.dumpModel(*_save,fmin->getFitParameters(),_nsave,_zsave,"",true);
                    if(fmin2) _parent.dumpModel(*_save,fmin2->getFitParameters(),_nsave,_zsave,"",true);
                }
                *_save << std::endl;
            }            
        }
        ~SamplingOutput() {
            if(_save) _save->close();
        }
        void saveSample(likely::FitParameters parameters, double fval,
        likely::FitParameters parameters2 = likely::FitParameters(), double fval2 = 0) {
            if(!_save) return;
            // Save fit parameter values and chisq.
            likely::Parameters pvalues;
            likely::getFitParameterValues(parameters,pvalues);
            BOOST_FOREACH(double pvalue, pvalues) {
                *_save << pvalue << ' ';
            }
            // Factor of 2 converts -logL to chiSquare.
            *_save << 2*fval << ' ';
            // Save alternate fit parameter values and chisq, if any.
            if(parameters2.size() > 0) {
                likely::getFitParameterValues(parameters2,pvalues);
                BOOST_FOREACH(double pvalue, pvalues) {
                    *_save << pvalue << ' ';
                }
                *_save << 2*fval2 << ' ';
            }
            // Save best-fit model multipoles, if requested.
            if(_nsave > 0) {
                _parent.dumpModel(*_save,parameters,_nsave,_zsave,"",true);
                if(parameters2.size() > 0) _parent.dumpModel(*_save,parameters2,_nsave,_zsave,"",true);
            }
            *_save << std::endl;            
       }
    private:
        int _nsave;
        double _zsave;
        CorrelationAnalyzer const &_parent;
        boost::scoped_ptr<std::ofstream> _save;
    };
}

int local::CorrelationAnalyzer::doSamplingAnalysis(CorrelationAnalyzer::AbsSampler &sampler,
std::string const &method, likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
std::string const &refitConfig, std::string const &saveName, int nsave, double zsave) const {
    if(nsave < 0) {
        throw RuntimeError("CorrelationAnalyzer::doSamplingAnalysis: expected nsave >= 0.");
    }
    if((!fmin2 && 0 < refitConfig.size()) || (!!fmin2 && 0 == refitConfig.size())) {
        throw RuntimeError("CorrelationAnalyzer::doSamplingAnalysis: inconsistent refit parameters.");
    }
    SamplingOutput output(fmin,fmin2,saveName,nsave,zsave,*this);
    baofit::AbsCorrelationDataCPtr sample;
    // Make a copy of the initial parameters.
    likely::FitParameters initParams = fmin->getFitParameters();
    likely::FitParameters initParams2;
    if(fmin2) initParams2 = fmin2->getFitParameters();
    // Initialize the parameter value statistics accumulators we will need.
    likely::FitParameterStatisticsPtr refitStats,
        fitStats(new likely::FitParameterStatistics(fmin->getFitParameters()));
    if(fmin2) {
        refitStats.reset(new likely::FitParameterStatistics(fmin2->getFitParameters()));
    }
    int nInvalid(0);
    // Loop over samples.
    int nsamples(0);
    // Use double parentheses below to tell clang that the '=' is not a typo.
    while((sample = sampler.nextSample())) {
        // Fit the sample.
        baofit::CorrelationFitter fitEngine(sample,_model,_covSampleSize);
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
            // Accumulate the fit results if the fit was successful.
            bool onlyFloating(true);
            fitStats->update(sampleMin->getParameters(onlyFloating),sampleMin->getMinValue());
            if(refitStats) refitStats->update(sampleMinRefit->getParameters(onlyFloating),sampleMinRefit->getMinValue());
            // Save the fit results, if requested.
            output.saveSample(sampleMin->getFitParameters(),sampleMin->getMinValue(),
                sampleMinRefit ? sampleMinRefit->getFitParameters() : likely::FitParameters(),
                sampleMinRefit ? sampleMinRefit->getMinValue() : 0);
        }
        else {
            nInvalid++;
            // Save the fit results, if requested.
            output.saveSample(initParams,0,sampleMinRefit ? initParams2 : likely::FitParameters(),0);
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
    return nInvalid;
}

void local::CorrelationAnalyzer::generateMarkovChain(int nchain, int interval, likely::FunctionMinimumCPtr fmin,
std::string const &saveName, int nsave, double zsave) const {
    if(nchain <= 0) {
        throw RuntimeError("CorrelationAnalyzer::generateMarkovChain: expected nchain > 0.");
    }
    if(interval < 0) {
        throw RuntimeError("CorrelationAnalyzer::generateMarkovChain: expected interval >= 0.");        
    }
    // Create a fitter to calculate the likelihood.
    AbsCorrelationDataCPtr combined = getCombined(true);
    CorrelationFitter fitter(combined,_model,_covSampleSize);
    // Generate the MCMC chains, saving the results in a vector.
    std::vector<double> samples;
    fitter.mcmc(fmin, nchain, interval, samples);
    // Output the results and accumulate statistics.
    SamplingOutput output(fmin,likely::FunctionMinimumCPtr(),saveName,nsave,zsave,*this);
    likely::FitParameters parameters(fmin->getFitParameters());
    likely::FitParameterStatistics paramStats(parameters);
    int npar = parameters.size();
    bool onlyFloating(true);
    std::vector<double>::const_iterator iter(samples.begin()), next;
    for(int i = 0; i < nchain; ++i) {
        // Copy the parameter values for this MCMC sample into pvalues.
        next = iter+npar;
        likely::Parameters pvalues(iter,next);
        // Get the chi-square corresponding to these parameter values.
        double fval = *next++;
        iter = next;
        // Output this sample.
        likely::setFitParameterValues(parameters,pvalues);
        output.saveSample(parameters,fval);
        // Accumulate stats on floating parameters.
        likely::Parameters pfloating;
        likely::getFitParameterValues(parameters,pfloating,onlyFloating);
        paramStats.update(pfloating,fval);
    }
    paramStats.printToStream(std::cout);
}

int local::CorrelationAnalyzer::parameterScan(likely::FunctionMinimumCPtr fmin,
AbsCorrelationDataCPtr sample, std::string const &saveName, int nsave, double zsave) const {
    int nfits(0);
    // Initialize the grid we will be sampling.
    likely::FitParameters params = fmin->getFitParameters();
    likely::BinnedGrid grid = likely::getFitParametersGrid(params);
    // Prepare to output samples at each grid point.
    SamplingOutput output(fmin,likely::FunctionMinimumCPtr(),saveName,nsave,zsave,*this);
    // Loop over grid points.
    for(likely::BinnedGrid::Iterator iter = grid.begin(); iter != grid.end(); ++iter) {
        nfits++;
        std::string config = likely::getFitParametersGridConfig(params,grid,iter);
        if(_verbose) {
            std::cout << std::endl << "-- Performing scan step " << nfits << " of " << grid.getNBinsTotal()
                << " using " << config << std::endl;
        }
        // Do the fit.
        try {
            likely::FunctionMinimumCPtr gridMin = fitSample(sample,config);
            // Add the results to our output file.
            output.saveSample(gridMin->getFitParameters(),gridMin->getMinValue());
        }
        catch(std::runtime_error const &e) {
            std::cerr << "ERROR while fitting:\n  " << e.what() << std::endl;
            // Use a copy of initial parameters with the config for this step applied.
            likely::FitParameters pcopy(params);
            likely::modifyFitParameters(pcopy,config);
            // Save results with fval = 0
            output.saveSample(pcopy,0.);
        }
    }
    return nfits;
}

void local::CorrelationAnalyzer::dumpResiduals(std::ostream &out, likely::FunctionMinimumPtr fmin,
AbsCorrelationDataCPtr combined, std::string const &script, bool dumpGradients) const {
    if(getNData() == 0) {
        throw RuntimeError("CorrelationAnalyzer::dumpResiduals: no observations have been added.");
    }
    likely::BinnedGrid grid = combined->getGrid();
    AbsCorrelationData::TransverseBinningType type = combined->getTransverseBinningType();
    // Get a copy of the the parameters at this minimum.
    likely::FitParameters parameters(fmin->getFitParameters());
    // Should check that these parameters are "congruent" (have same names?) with model params.
    // assert(parameters.isCongruent(model...))
    // Modify the parameters using the specified script, if any.
    if(0 < script.length()) likely::modifyFitParameters(parameters, script);
    // Get the parameter values (floating + fixed)
    likely::Parameters parameterValues,parameterErrors;
    likely::getFitParameterValues(parameters,parameterValues);
    int npar = fmin->getNParameters();
    if(dumpGradients) likely::getFitParameterErrors(parameters,parameterErrors);
    // Loop over 3D bins in the combined dataset.
    std::vector<double> centers;
    for(likely::BinnedData::IndexIterator iter = combined->begin(); iter != combined->end(); ++iter) {
        int index(*iter);
        out << index;
        if(combined->useCustomGrid()) {
            combined->getCustomBinCenters(index,centers);
        }
        else {
            grid.getBinCenters(index,centers);
        }
        BOOST_FOREACH(double center, centers) {
            out << ' ' << center;
        }
        double data = combined->getData(index);
        double error = combined->hasCovariance() ? std::sqrt(combined->getCovariance(index,index)) : 0;
        double z = combined->getRedshift(index);
        double r = combined->getRadius(index);
        double predicted,mu;
        cosmo::Multipole multipole;
        if(type == AbsCorrelationData::Coordinate) {
            mu = combined->getCosAngle(index);
            predicted = _model->evaluate(r,mu,z,parameterValues);
            out  << ' ' << r << ' ' << mu << ' ' << z;
        }
        else {
            multipole = combined->getMultipole(index);
            predicted = _model->evaluate(r,multipole,z,parameterValues);
            out  << ' ' << r << ' ' << (int)multipole << ' ' << z;
        }
        out << ' ' << predicted << ' ' << data << ' ' << error;
        if(dumpGradients) {
            for(int ipar = 0; ipar < npar; ++ipar) {
                double gradient(0), dpar(0.1*parameterErrors[ipar]);
                if(dpar > 0) {
                    double p0 = parameterValues[ipar];
                    parameterValues[ipar] = p0 + 0.5*dpar;
                    double predHi = (type == AbsCorrelationData::Coordinate) ?
                        _model->evaluate(r,mu,z,parameterValues) :
                        _model->evaluate(r,multipole,z,parameterValues);
                    parameterValues[ipar] = p0 - 0.5*dpar;
                    double predLo = (type == AbsCorrelationData::Coordinate) ?
                        _model->evaluate(r,mu,z,parameterValues) :
                        _model->evaluate(r,multipole,z,parameterValues);
                    gradient = (predHi - predLo)/dpar;
                    parameterValues[ipar] = p0;
                }
                out << ' ' << gradient;
            }
        }
        out << std::endl;
    }
}

void local::CorrelationAnalyzer::dumpModel(std::ostream &out, likely::FitParameters parameters,
int ndump, double zdump, std::string const &script, bool oneLine) const {
    if(ndump <= 1) {
        throw RuntimeError("CorrelationAnalyzer::dump: expected ndump > 1.");
    }
    if(zdump < 0) {
        // Use the redshift of the first bin of the first subsample, if none is specified
        AbsCorrelationDataCPtr subsample =
            boost::dynamic_pointer_cast<const baofit::AbsCorrelationData>(_resampler.getObservation(0));
        zdump = subsample->getRedshift(0);
    }
    if(_verbose && !oneLine) {
        std::cout << "Dumping fit model at " << ndump << " points using z = " << zdump << std::endl;
    }

    // Should check that parameters are "congruent" (have same names?) with model params.
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
        double mono = _model->evaluate(rval,cosmo::Monopole,zdump,parameterValues);
        double quad = _model->evaluate(rval,cosmo::Quadrupole,zdump,parameterValues);
        double hexa = _model->evaluate(rval,cosmo::Hexadecapole,zdump,parameterValues);
        // Output the model predictions for this radius in the requested format.
        if(!oneLine) out << rval;
        out << ' ' << boost::lexical_cast<std::string>(mono) << ' '
            << boost::lexical_cast<std::string>(quad) << ' '
            << boost::lexical_cast<std::string>(hexa);
        if(!oneLine) out << std::endl;
    }
}

void local::CorrelationAnalyzer::getDecorrelatedWeights(AbsCorrelationDataCPtr data,
likely::Parameters const &params, std::vector<double> &dweights) const {
    CorrelationFitter fitter(data, _model,_covSampleSize);
    std::vector<double> prediction;
    fitter.getPrediction(params,prediction);
    data->getDecorrelatedWeights(prediction,dweights);
}

namespace baofit {
    bool accumulationCallback(likely::CovarianceAccumulatorCPtr accumulator) {
        std::cout << "accumulated " << accumulator->count() << " samples." << std::endl;
        return true;
    }
}

likely::CovarianceMatrixPtr
local::CorrelationAnalyzer::estimateCombinedCovariance(int nSamples, std::string const &filename) const {
    // Print a message every 10 samples during accumulation.
    likely::CovarianceAccumulatorPtr accumulator = _resampler.estimateCombinedCovariance(nSamples,
        likely::BinnedDataResampler::AccumulationCallback(accumulationCallback),10);
    // Save the work in progress if a filename was specified.
    if(filename.length() > 0) {
        std::cout << "saving work in progress to " << filename << std::endl;
        std::ofstream out(filename.c_str());
        accumulator->dump(out);
        out.close();
    }
    // Return the estimate covariance (which might not be positive definite)
    return accumulator->getCovariance();
}
